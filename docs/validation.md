# Reproducible validation

## Standard build and tests

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
# Individual checks and numerical statistics:
./build/bnsim_tests
```

The six C++ groups are `validation`, `transport`, `biology`, `mechanics`,
`motility` and `reproducibility`. The seventh CTest entry, `examples_and_outputs`,
runs both programs, validates their data and exercises rejected CLI inputs.
Assertions remain active in Release builds. The complete inventory and defect
mapping is in the [issue register](issue_register.md).

Transport checks include exact zero-decay invariance, analytical decay,
closed-system conservation, material discontinuities, no-flux walls, reservoir
replenishment and noncubic geometry. A no-flux cosine eigenmode is compared
against the discrete analytical eigenvalue for temporal refinement, and the
continuum eigenvalue for spatial refinement. The required observed error ratios
are >1.9 under timestep halving and >3.8 under grid-spacing halving with matching
timestep refinement.

Biochemistry checks cover all binding invariants, analytical two-compartment
exchange, shared availability, starvation, growth dilution, EPS allocation,
division, death, shedding, feedback and inhibitors. Failure-path tests check
complete rollback and joined workers. Budget residual tolerances include a
roundoff allowance proportional to throughput, necessary when final substrate
is nearly depleted.

Motility checks use fixed seeds and statistical tolerances: isotropic heading
moments; tangent-Gaussian orientation correlation; truncated gamma moments
against independent quadrature; exponential dwell means, survival probability
and two-state transition probability; stationary tumbles; adaptation; and
population drift in an attractant gradient. Their finite-sample tolerances are
not experimental confidence bounds on the biological model.

Reproducibility compares fields, membership, agents, motion state, entire random
engines, time, IDs and diagnostics exactly after every step across 1, 2 and 8
workers in a mixed world containing growth, division, shedding, dead cells and
inert EPS. Model parameters are inherited identically; separate checks verify
daughter concentrations and parameter inheritance.

## Independent chemical trajectories

`tests/reference/coupled.csv` contains six samples from 0 to 60 seconds for a
finite extracellular compartment coupled to a growing cell, both QS circuits,
all inhibitors and all five exchange fluxes. The reference generator expresses
reaction coupling through an independently assembled stoichiometric matrix and
uses uM for intracellular evaluation before converting to amounts. It neither
imports the simulator nor reuses its solver.

```sh
# Optional, only to regenerate references; not a CTest dependency.
python3 -m venv .venv-reference
.venv-reference/bin/pip install numpy scipy
.venv-reference/bin/python tests/generate_reference.py
```

Both SciPy Radau (implicit) and DOP853 (explicit eighth order) must agree within
the generator's scaled threshold before it writes the fixture. Regeneration
records Python/NumPy/SciPy versions, solver tolerances, source/CSV SHA-256 hashes
and maximum scaled disagreement in `tests/reference/provenance.json`.
The C++ solver is compared at all six times and across every extracellular,
intracellular and biomass/EPS component. This is numerical verification of the
specified equations, not an independently measured biological time course.
[SciPy's solver documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
describes the two reference methods and tolerance definitions.

## Runtime sanitizers

```sh
cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug -DBNSIM_SANITIZER=address
cmake --build build-asan --parallel 4
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  ctest --test-dir build-asan --output-on-failure

cmake -S . -B build-tsan -DCMAKE_BUILD_TYPE=Debug -DBNSIM_SANITIZER=thread
cmake --build build-tsan --parallel 4
TSAN_OPTIONS=halt_on_error=1 ctest --test-dir build-tsan --output-on-failure
```

The address configuration enables ASan and UBSan together, with leak detection
on Linux. ThreadSanitizer is a separate build. Empty worlds, one/excessive worker
counts, mixed agents, growth, birth/removal and exception paths are in these
same suites. The CI workflow builds on Linux and macOS and runs Linux sanitizers.
Its disposable TSan runner sets mmap entropy to a shadow-memory-compatible
value; it does not disable the sanitizer or suppress simulator reports.

## Execution record

See [validation_results.md](validation_results.md) for the actual local run,
toolchain, counts and remaining verification limits. A configured CI workflow
is not evidence that its remote Linux/macOS jobs have already executed.
