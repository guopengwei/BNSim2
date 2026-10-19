# Validation execution record

Workspace date: 2026-10-18. Baseline revision: `203f29c`; validation was performed
on the modified working tree (`203f29c-dirty`), not on the unchanged baseline.
Simulator/CMake source SHA-256, as recorded in the example manifests:

```
7aa49c7e1838a4f31604917acfbf6cae566b6522abd15a836dfc9387be28ccae
```

Toolchain: Linux x86_64, GNU C++ 15.2.0, C++17, CMake 4.4.4. CTest selected
Python 3.12.12 for CLI/output verification. CMake was installed into a temporary
validation tools directory because it was absent initially. Builds used
`-Wall -Wextra -Wpedantic -Werror`; the Release checks also ran with `NDEBUG`.

| Build/check | Observed result |
| --- | --- |
| Release | 7/7 CTest entries passed, including all 42 C++ regression cases |
| AddressSanitizer + UndefinedBehaviorSanitizer | 7/7 entries passed; `detect_leaks=1`, `halt_on_error=1`; no findings |
| ThreadSanitizer | 7/7 entries passed; `halt_on_error=1`; no findings or suppressions |
| Worker reproducibility | Complete mutable states identical after each step with 1, 2 and 8 workers |
| Independent biochemical reference | All six times/all 21 reference components within committed tolerances |
| Source/reference consistency | Generator and CSV SHA-256 checks passed |
| Diff whitespace check | `git diff --check` passed |

The C++ cases comprise 5 input/ownership checks, 8 transport checks, 15 biology
checks, 5 mechanics/layer checks, 7 motility checks and 2 reproducibility checks.
The output test additionally exercises both executables with 3 and 0 cells,
8 workers, short final steps, invalid CLI values, long percent-containing paths,
unique run directories, manifest fields, field units, physical gradients,
projections and conservation budgets.

The fixed-seed motility ensemble contained 256 cells per condition over 40 s.
Mean x displacements were -13.4818 um (standard error 8.27987 um) in a uniform
field and +184.321 um (standard error 9.27001 um) in a rising attractant gradient.
This verifies qualitative directional response in this numerical scenario;
these are simulated values, not measured biological performance.

Reference generation used Python 3.14.4, NumPy 2.5.4 and SciPy 1.18.1. Radau
and DOP853 agreed with maximum scaled discrepancy 0.01038 against the generator's
allowed threshold of 1. Full provenance is in
[`tests/reference/provenance.json`](../tests/reference/provenance.json).

## Saved example runs

Both runs used seed 0, 8 workers, 16 cells, requested dt 0.1 s and end time 2.35 s.
Each generated 25 snapshots including time zero. The final step was shortened.
Outputs are deliberately ignored by Git but remain in the current workspace:

- `runs/chemotaxis/run-00000000/`: 1000 field rows per snapshot; final material
  budget residual exactly zero.
- `runs/qs_biofilm/run-00000000/`: 2592 field rows per snapshot; largest final
  residual approximately 1.05e-9 fg in the substrate budget, at roundoff scale
  relative to its approximately 216000-fg inventory.

An additional readback checked every final field amount/concentration/volume,
projection, timestamp, agent radius and wall bound. All quantities were finite
and nonnegative where required. Both examples completed successfully with 16
agents and final time exactly equal to the requested double value 2.35.

## Limits of this record

Remote CI was configured, not executed from this workspace. In particular, the
macOS job has not been locally observed. Cross-platform bitwise equality is not
claimed. Sanitizer results apply to the exercised scenarios, not every possible
input. Parameter and biological calibration remain provisional; open empirical
requirements are V01–V04 in the [issue register](issue_register.md).
