# Model and numerical contract

This is a corrected implementation of the reduced mechanisms in repository
revision `203f29c`, with explicit closures for unfinished behavior. It is not a
reimplementation of every equation or experiment in the cited papers.

## Quantities and ownership

Internal length, time and material units are micrometers, seconds and
femtograms. Signaling species use femtomoles. One liter is 10^15 um^3, so
1 M = 1 fmol/um^3, 1 uM = 10^-6 fmol/um^3, 1 mM = 10^-3 fmol/um^3,
and 1 g/L = 1 fg/um^3. `Species::quantity` distinguishes mass and molar fields;
mapping a mass species to a signaling network is rejected. A substrate cannot
be converted from molarity to mass without specifying its molecular weight.

Species, cells, embedded networks and grids have exclusive value ownership in
each `Universe`. A grid contains non-owning stable IDs, not owning cell pointers.
IDs are monotonic and never reused, including after explicit removal. Arrays of
intracellular amounts are in fmol. Concentrations are derived as `n / Vcell` at
each reaction evaluation; EPS volume is excluded from the cell interior.

The field carrier volume is `dx*dy*dz` in every voxel. Cell occupation does not
subtract solvent volume. This is a continuum approximation; diffusion in the
biomass layer represents reduced transport, without an explicit porosity model.

## Step order and failures

For a step of duration h:

1. Transport and decay for h/2.
2. Coupled extracellular chemistry, intracellular reactions and growth for h.
3. Motility, deaths, conservative divisions and EPS shedding.
4. Contact relaxation, membership commit and layer/reservoir updates.
5. Transport and decay for h/2; validate; advance time and step index.

All work occurs on a candidate state, committed only on success. Exceptions
include stage/voxel or timestep context. Integration and mechanics have bounded
iteration counts; they never report success after silently discarding a failure.
`run_until` shortens the final step. Configuration, geometry, resource limits,
species mappings, quantities and finite values are checked before use.

## Finite-volume transport

Positions use floor-based indexing into a uniform rectangular grid; spacings
are floating point. The coordinate exactly at a domain upper face maps to the
last voxel. Initial particle surfaces must lie inside the domain.

For one interface between voxels a and b, the amount transferred from a to b is

```
F_ab = h * D_ab / dx_axis^2 * (N_a - N_b)
D_ab = 2 D_a D_b / (D_a + D_b)
```

The harmonic mean is zero if either diffusivity is zero and is evaluated with
an overflow-resistant expression. Every interface appears once. All fluxes use
an immutable field; each face writes a separate result. Their signed changes
are reduced in fixed face order. Each substep satisfies
`h * max_i(sum_j D_ij / dx_ij^2) <= 0.9`. The solver inserts substeps or fails
if the configured limit would be exceeded. It does not clamp negative amounts.

Within each transport call, exact first-order decay is split half before and
half after diffusion: `N <- N exp(-k h)`. Decay remains active when diffusion is
disabled. Zero decay is exactly the identity. Closed walls have no faces to
the outside, hence no chemical flux. The diffusion update is forward Euler:
overall temporal order is at most one despite symmetric outer splitting.
Spatial order is two in smooth homogeneous regions. Interfaces and moving
layer classifications can reduce observed convergence order.

## Fixed reservoirs and budgets

Material type chooses diffusivity/decay; it never implies a reservoir. An
explicit per-species reservoir is configured through `set_reservoir`. Layered
bulk voxels also provide reservoirs. Explicit targets override bulk targets.
Fixed concentrations are restored/accounted during transport and held constant
in the local ODE. Removal of a reservoir leaves its current amount in place.
Decay in a reservoir is recorded along with the replenishment it requires.

Diagnostics record signed external edits/removals, signed reservoir supply,
positive extracellular decay and positive reaction sources/sinks. For each
species, including its mapped intracellular signal/inhibitor moieties,

```
total = external + reservoir + reaction_source - reaction_sink - decay
```

`external` includes initial field quantities and material imported in created
cells. Division and shedding are internal transfers, not external sources.
Las signal inventory includes `A1 + RA1 + 2*C1`; Rhl uses `A2 + RA2 + 2*C2`;
competitive inhibitor inventory includes `Q1 + RQ1`. Intracellular molecules
in dead cells remain in the inventory. Receptors and RsaL are not extracellular
species; their binding invariants are checked separately.

These are floating-point identities up to summation/integration roundoff.
Residual tolerances should scale with total material throughput, not only the
possibly vanishing final amount. Diagnostics report raw residuals without
subtracting or clipping numerical discrepancies.

## QS circuits

Concentrations in this section are intracellular, in fmol/um^3. Each rate below
is multiplied by the instantaneous cell volume to obtain amount derivatives.
The array names and ordering are defined in `include/bnsim/model.hpp`.

Define reversible reaction fluxes:

```
J1 = k1*A1*R1 - k2*RA1       H1 = k4*RA1^2 - k5*C1
J2 = k8*A2*R2 - k9*RA2       H2 = k11*RA2^2 - k12*C2
JQ = kon*Q1*R1 - koff*RQ1
```

The two complex derivatives are `dRA1 = J1-2*H1`, `dC1=H1`, and
`dRA2=J2-2*H2`, `dC2=H2`. Signal and receptor derivatives subtract the
corresponding J. The competition reaction subtracts JQ from both R1 and Q1 and
adds it to RQ1, restoring the receptor on unbinding. Thus reversible binding
conserves ligand, receptor and inhibitor moieties, including dimer factors of 2.

With `sat(c,K)=c/(K+c)`, synthesis terms are:

```
P_A1   = basal[0] + induced[0] sat(C1,KL1) * KS1/(KS1+RsaL)
P_R1   = (basal[1] + induced[1] sat(C1,KR1)) / (1+Q3/Ksupp)
P_RsaL = basal[2] + induced[2] sat(C1,KS)
P_A2   = basal[3] + induced[3] sat(C2,KL2)
P_R2   = basal[4] + induced[4] sat(C2,KL3) sat(C1,KL4)
```

Disabling RsaL feedback replaces its multiplicative suppression by 1; RsaL is
still synthesized. Free A1, R1, RsaL, A2 and R2 decay at k0, k3, k6, k7 and k10,
respectively. Complex degradation is absent, as in this legacy mechanism.
`activated` is true when a live cell's C1 concentration reaches its configured
threshold. Both Las and Rhl are enabled for every configured QS network.

The legacy inhibitor identities are retained:

| Mapping | Mechanism |
| --- | --- |
| QSI1 | Reversible competition for free LasR, through JQ |
| QSI2 | Catalytic A1 degradation, `kh*A1*Q2/(1+g*Q2)` |
| QSI3 | Suppression of LasR production, through `1/(1+Q3/Ksupp)` |

QSI2 and QSI3 are not consumed by their actions. A missing inhibitor mapping
disables its mechanism and exchange; nonzero initial amounts then fail
validation. Both free signals and each enabled free inhibitor exchange with
the voxel through `F = permeability * (n_inside - C_outside*Vcell)`. The same
flux is subtracted inside and added outside, in the same coupled solve.

[Wei et al. (2016)](https://www.nature.com/articles/srep35136) supports the
Las feedback/binding structure, Monod/EPS connection and layered environment.
Its signal-production inhibition targets LasI; the requested legacy QSI3 here
targets receptor production. Its degradation target is extracellular AHL;
legacy QSI2 here acts intracellularly after membrane exchange. Restored Rhl and
RsaL circuits come from the repository's commented equations. These differences
preclude claiming a reproduction of that paper's calibrated treatment results.

## Growth, EPS, division and death

The corrected legacy closure is:

```
f = S/(Ks+S)
p = pmax*C1/(Kc+C1)
dX/dt = ((umax-p)*f-maintenance)*X
dEPS/dt = p*f*X
substrate uptake = umax*f*X/yield
maintenance loss = maintenance*X
yield loss = umax*f*X*(1/yield-1)
```

Here X and EPS are fg, S is fg/um^3, C1 is fmol/um^3, and rates are 1/s.
Require `0 <= pmax <= umax` and `0 < yield <= 1`. With no QS network p=0.
For a closed system, substrate + biomass + EPS + accumulated maintenance loss
+ accumulated yield loss is conserved. Maintenance consumes existing biomass
even at zero substrate. QS protein/signal synthesis does not consume a
separately modeled carbon/amino-acid pool: this is an explicit missing metabolic
cost and not a calibrated elemental balance.

Biomass volume is X/rhoX; total volume adds EPS/rhoEPS. Spherical radii use
`cbrt(3*V/(4*pi))`. Growth therefore changes radius and dilutes intracellular
concentrations without removing intracellular amounts. All cell reactions and
extracellular consumption in a voxel are integrated simultaneously, preventing
cells from independently consuming the same available substrate.

At the lifecycle stage, a cell at or above its division mass splits biomass,
EPS and every intracellular amount exactly in half. Daughters inherit densities,
parameters and concentrations. Oversized descendants divide again, within
resource limits. Positions separate deterministically; each new ID gets a new
random stream. EPS at its shedding threshold transfers entirely to an inert
particle whose radius comes from transferred mass/density. A cell at or below
the death mass becomes inactive but retains all material and mechanical volume.
Lifecycle thresholds are checked at stage boundaries, not localized as ODE
events; timestep refinement is necessary for event timing.

## Local integration

Each occupied voxel's ODE vector contains extracellular amounts, all resident
active cell amounts/biomass/EPS, and source/sink counters. The solver is adaptive
Dormand–Prince 5(4), with a maximum component error norm using separate molecular
and mass absolute tolerances and a relative tolerance. Negative or nonfinite
trial stages/endpoints are rejected, with bounded retries and a minimum step.
The solver is explicit and can be inefficient or fail for very stiff parameters;
it reports the failure instead of returning an unreliable trajectory.

Dormand–Prince has negative internal weights: a zero-initialized reaction
cascade can have a negative intermediate value at every positive step size.
When that happens, a step-doubled SSPRK3 startup uses only validated nonnegative
stages, with its own error estimate `|two_half_steps - full_step|/7`. After
startup, Dormand–Prince resumes. No accepted chemical amount is clipped or
arbitrarily seeded. `ode_startup` records these steps. The exact polynomial
chain `x'=1, y'=x^2, z'=y^2` verifies this path, and independently computed
Radau/DOP853 trajectories verify the full coupled reaction system.

## Reduced chemotaxis and motility

The retained receptor model is `A = 1/(1+exp(F))`, where

```
F = N * [alpha*(m0-m) - log((1+L/KA)/(1+L/KI))]
dm/dt = kr*(1-A) - kb*A
Y = Ya*A
k_run_to_tumble = w*exp(-g/4 + g*Y/(2*(Kd+Y)))
k_tumble_to_run = w*exp(g/4 - g*Y/(2*(Kd+Y)))
```

The logistic and logarithms are evaluated stably. Adaptation uses RK4 and
applies the timestep once. The ligand is sampled from the cell's committed
voxel after the first half transport stage and held for the movement stage.
Subintervals no larger than `motility.max_step` freeze motor rates while
integrating adaptation. A residual unit-exponential hazard persists across
subintervals and user steps; switches occur when accumulated hazard reaches
that sample. Constant-rate dwell times are exact exponential draws. For varying
activity the hazard is approximated by these subintervals, requiring refinement.
Multiple switches per user step are possible.

Only run intervals contribute swimming displacement. During a run, two
independent tangent-plane Gaussian increments have variance `2*Dr*h` each.
The exponential map preserves a unit heading; angular mean-square displacement
approaches `4*Dr*h`, and orientation correlation approaches `exp(-2*Dr*t)` as
the substep shrinks. Movement uses a midpoint heading approximation. Initial
headings are isotropic on the sphere. At each tumble exit, one reorientation is
drawn from a gamma distribution, conditioned on [0,180] degrees, with uniform
azimuth and conversion by pi/180. There is no rotational swim during a tumble.

[Sneddon et al. (2012), Methods](https://emonet.biology.yale.edu/assets/papers/past/sneddon_emonet_2012_paper.pdf)
gives the adopted 20 um/s speed, 0.062 rad^2/s rotational diffusion, stationary
tumbles and gamma shape 4/scale 18.32 degrees. The legacy -4.60-degree shift is
removed; truncation is an explicit bounded-angle assumption. Their multi-motor,
flagellar-conformation and noisy adaptation model is not implemented here.
The deterministic reduced receptor/adaptation model remains a separate choice.

## Mechanics, layers and parallelism

Swimming reflects at radius-aware walls, including multiple crossings. Contact
relaxation projects growth/division corrections inside the walls. A sweep
searches all voxels within the cell radius plus the largest particle radius,
including diagonal and distant bins. Pair corrections use stable ID order;
coincident centers use a deterministic pair-derived direction. Sweeps repeat
until overlap is below tolerance or fail with the remaining overlap. This is a
geometric exclusion algorithm, not an inertial or calibrated elastic model.

With `layered=true`, biomass height is the greatest `z + total_radius` over all
particles, including EPS and dead cells. Voxel centers at/below that height use
biofilm material; a boundary layer extends above it; higher centers are fixed
bulk reservoirs. Any target reset when layer membership changes is recorded as
reservoir transfer. This is a horizontally uniform, voxel-center layer model,
not a resolved local biofilm surface. Explicit reservoir overrides persist.

Per-agent `mt19937_64` streams are seeded by a SplitMix64-style mix of run seed
and stable ID. Fixed partitioning, disjoint worker results, serial commits and
fixed reduction orders make the complete state independent of worker count
and scheduling on one build/platform. Workers are limited to available work;
creation failures join already-created threads and worker exceptions propagate.
Standard-library probability distributions and floating-point libraries can
differ across platforms, so bitwise cross-platform equality is not promised.

Transactional step copies and explicit solvers prioritize correctness over
large-population throughput. Dense or physically impossible packings may fail
mechanical convergence. Resource limits, tolerances and resolution are recorded;
successful execution alone does not establish numerical convergence or empirical
validity for a new scientific scenario.
