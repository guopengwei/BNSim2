"""Independent amount/stoichiometry implementation; never imports BNSim.

Regenerate with Python, NumPy and SciPy. Radau and DOP853 must agree before
fixtures are written. Normal CTest runs use the committed CSV and need no SciPy.
"""
import csv
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
OUT = HERE / "reference"
TIMES = np.array([0, 0.1, 1, 5, 20, 60])
NAMES = "A1 R1 RA1 C1 RsaL A2 R2 RA2 C2 QSI1 RQSI1 QSI2 QSI3".split()
# Independent stoichiometric matrix: Las binding/dimerization, Rhl binding/
# dimerization and competitive binding. Chemical concentrations below use uM.
STOICH = np.zeros((13, 5))
for reaction, entries in enumerate([
    {0: -1, 1: -1, 2: 1}, {2: -2, 3: 1},
    {5: -1, 6: -1, 7: 1}, {7: -2, 8: 1}, {1: -1, 9: -1, 10: 1},
]):
    for species, coefficient in entries.items():
        STOICH[species, reaction] = coefficient


def rhs(_time, y):
    external, biomass, _eps, amounts = y[:6], y[6], y[7], y[8:]
    volume = biomass / 150
    a, r, ra, c, rsal, a2, r2, ra2, c2, q1, rq1, q2, q3 = amounts / volume * 1e6
    flux = np.array([
        0.7*a*r-0.2*ra, 0.3*ra**2-0.1*c,
        0.6*a2*r2-0.25*ra2, 0.2*ra2**2-0.08*c2,
        0.5*q1*r-0.12*rq1,
    ])
    intracellular = STOICH @ flux
    synthesis = np.array([
        0.002+0.02*c/(0.2+c)*0.2/(0.2+rsal),
        (0.002+0.02*c/(0.2+c))/(1+q3/0.3),
        0.002+0.02*c/(0.2+c),
        0.002+0.02*c2/(0.2+c2),
        0.002+0.02*c2/(0.2+c2)*c/(0.2+c),
    ])
    intracellular[[0, 1, 4, 5, 6]] += synthesis
    intracellular[[0, 1, 4, 5, 6]] -= np.array([0.03*a, 0.04*r, 0.02*rsal, 0.02*a2, 0.03*r2])
    intracellular[0] -= 0.2*a*q2/(1+q2)
    intracellular *= volume*1e-6
    outside = np.zeros(6)
    for internal, external_species, permeability in zip([0, 5, 9, 11, 12], [1, 2, 3, 4, 5], [.2, .25, .3, .35, .4]):
        transfer = permeability*(amounts[internal]-external[external_species]/1000*volume)
        intracellular[internal] -= transfer
        outside[external_species] += transfer
    nutrient = external[0]/1000
    saturation = nutrient/(0.2+nutrient)
    allocation = 0.01*c/(0.1+c)
    dx = ((0.03-allocation)*saturation-0.002)*biomass
    deps = allocation*saturation*biomass
    outside[0] = -0.03*saturation*biomass/0.8
    return np.r_[outside, dx, deps, intracellular]


def main():
    y0 = np.r_[np.array([.4, .5e-6, .2e-6, .8e-6, .4e-6, .6e-6])*1000,
               150, 5, np.array([1, .7, .2, .1, .05, .8, .6, .1, .03, .4, .05, .2, .3])*1e-6]
    atol = np.full(len(y0), 1e-20)
    atol[[0, 6, 7]] = 1e-11
    solutions = [solve_ivp(rhs, [0, 60], y0, method=method, t_eval=TIMES,
                           rtol=2e-12, atol=atol) for method in ("Radau", "DOP853")]
    for sol in solutions:
        assert sol.success, sol.message
        assert np.all(sol.y >= 0)
    delta = np.abs(solutions[0].y-solutions[1].y)
    tolerance = np.maximum(atol[:, None]*100, np.abs(solutions[0].y)*1e-9)
    assert np.all(delta < tolerance), np.max(delta/tolerance)
    OUT.mkdir(exist_ok=True)
    with (OUT/"coupled.csv").open("w", newline="") as file:
        writer = csv.writer(file, lineterminator="\n")
        writer.writerow(["time_s", "substrate_fg", "external_A1_fmol", "external_A2_fmol",
                         "external_QSI1_fmol", "external_QSI2_fmol", "external_QSI3_fmol",
                         "biomass_fg", "eps_fg"]+[name+"_fmol" for name in NAMES])
        for t, values in zip(TIMES, solutions[0].y.T):
            writer.writerow([format(v, ".17g") for v in [t, *values]])
    metadata = {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                "solvers": ["Radau", "DOP853"], "rtol": 2e-12, "atol": atol.tolist(),
                "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "csv_sha256": hashlib.sha256((OUT/"coupled.csv").read_bytes()).hexdigest(),
                "max_scaled_solver_disagreement": float(np.max(delta/tolerance)),
                "status": "numerical reference for specified equations; not empirical validation"}
    (OUT/"provenance.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
