"""End-to-end CLI, file schema, units, projections and ledger regression checks."""
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile


def rows(path):
    with path.open() as file:
        return list(csv.DictReader(file))


def run(binary, parent, population=3, duration=.25):
    result = subprocess.run([binary, "--seed", "0", "--threads", "8", "--duration", str(duration),
                             "--dt", ".1", "--population", str(population), "--output", str(parent)],
                            text=True, capture_output=True)
    assert result.returncode == 0, result.stdout+result.stderr
    directory = sorted(parent.glob("run-*"))[-1]
    metadata = json.loads((directory/"metadata.json").read_text())
    assert metadata["seed"] == 0 and metadata["workers"] == 8
    assert len(metadata["source_sha256"]) == 64
    assert metadata["units"]["mass"] == "fg"
    assert metadata["duration_s"] == duration
    agent_files = sorted(directory.glob("agents_*.csv"))
    assert len(agent_files) == (4 if duration == .25 else 2)
    assert len(list(directory.glob("field_*.csv"))) == len(agent_files)
    for file in agent_files:
        agents = rows(file)
        assert len(agents) == population
        assert len({a["agent_id"] for a in agents}) == population
        for agent in agents:
            for key, value in agent.items():
                assert math.isfinite(float(value)), (key, value)
            assert float(agent["biomass_fg"]) > 0
            assert float(agent["eps_fg"]) >= 0
    field_files = sorted(directory.glob("field_*.csv"))
    final = rows(field_files[-1])
    assert len(final) == math.prod(metadata["shape"])*len(metadata["species"])
    for record in final:
        assert float(record["time_s"]) == duration
        assert float(record["concentration"]) >= 0
        assert record["concentration_unit"] == record["amount_unit"]+"/um^3"
        assert math.isclose(float(record["amount"]), float(record["concentration"])*float(record["volume_um3"]), rel_tol=1e-14, abs_tol=1e-20)
    projection = rows(sorted(directory.glob("depth_average_*.csv"))[-1])
    sums = {}
    for record in final:
        key = (record["x_um"], record["y_um"], record["species"])
        sums[key] = sums.get(key, 0)+float(record["concentration"])
    for record in projection:
        key = (record["x_um"], record["y_um"], record["species"])
        assert math.isclose(float(record["depth_average_concentration"]), sums[key]/metadata["shape"][2], rel_tol=1e-14, abs_tol=1e-20)
    diagnostics = json.loads(sorted(directory.glob("diagnostics_*.json"))[-1].read_text())
    assert diagnostics["time_s"] == duration
    for species in diagnostics["species"]:
        scale = abs(species["external"])+abs(species["reservoir"])+abs(species["reaction_source"])
        assert abs(species["budget_residual"]) < 1e-16+1e-8*scale
    if metadata["scenario"] == "chemotaxis":
        for record in rows(field_files[0]):
            expected = (.01+.0001*(metadata["extent_um"][1]-float(record["y_um"])))*1e-3
            assert math.isclose(float(record["concentration"]), expected, rel_tol=1e-14)
    if population:
        parameters = rows(directory/"parameters.csv")
        assert parameters and all(p["unit"] and p["provenance"] and p["status"] for p in parameters)
    return directory


def main():
    reference = Path(__file__).resolve().parent/"reference"
    provenance = json.loads((reference/"provenance.json").read_text())
    assert hashlib.sha256((reference/"coupled.csv").read_bytes()).hexdigest() == provenance["csv_sha256"]
    assert hashlib.sha256((reference.parent/"generate_reference.py").read_bytes()).hexdigest() == provenance["generator_sha256"]
    with tempfile.TemporaryDirectory(prefix="bnsim-output-validation-") as temp:
        for binary in sys.argv[1:]:
            parent = Path(temp)/(Path(binary).name+"_"+"long%-directory-"*8)
            first = run(binary, parent)
            second = run(binary, parent, population=0, duration=.025)
            assert first != second
            for args in [["--dt", "0"], ["--dt", "nan"], ["--threads", "0"], ["--threads", "99999"],
                         ["--duration", "inf"], ["--population", "-1"], ["--seed", "-1"], ["--output"],
                         ["--unknown", "1"], ["--dt", "0.1", "--dt", "0.2"]]:
                result = subprocess.run([binary, *args], text=True, capture_output=True)
                assert result.returncode != 0, args
                assert "BNSim:" in result.stderr, result.stderr
    print("Both examples: checked arguments, empty/mixed output, units, timestamps, unique directories, gradients, projections and budgets passed")


if __name__ == "__main__":
    main()
