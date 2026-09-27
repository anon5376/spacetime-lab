# Spacetime Lab — Physics reference

Per-module equations, approximations, and limits for the N-body sandbox core
and the analytic relativity modules. Companion to `CONTRACT.md` (API/build
contract), `README.md` (run/controls), and `doc/PLAN-SANDBOX.md` (sandbox
build plan). Verification runs behind §13 were executed by the
integrator on 2026-09-27; the table records their measured output.

## 1. Unit system

SI doubles on every API boundary: meters (m), kilograms (kg), seconds (s),
radians (rad). Angles are radians unless the name says otherwise (`*_deg`,
`*_arcsec`). All physics uses C `double`; `float` and `long double` are
forbidden in sim code (platform variance).

Constants (`include/stc/units.h`):

| Symbol | Value | Unit |
|---|---|---|
| `STC_C` | 299792458 (exact) | m/s |
| `STC_G` | 6.67430e-11 | m^3 kg^-1 s^-2 |
| `STC_M_SUN` | 1.98847e30 (IAU 2015 nominal) | kg |
| `STC_M_EARTH` | 5.97237e24 | kg |
| `STC_R_SUN` | 6.957e8 (IAU 2015 nominal) | m |
| `STC_AU` | 149597870700 (exact, IAU 2012) | m |
| `STC_YEAR` | 31557600 (Julian year, 365.25 d) | s |
| `STC_DAY` | 86400 (exact) | s |
| `STC_PI` | 3.14159265358979323846 | dimensionless |

Helpers: Schwarzschild radius `rs = 2GM/c^2` (`stc_schwarzschild_radius`,
meters); geometric mass `M_geom = GM/c^2` (`stc_geometric_M`, meters);
arcsecond/microarcsecond in radians (`stc_arcsec`, `stc_microarcsec`).
Kerr/embedding/geodesic lengths are meters (geometric meters where noted);
frame-dragging rates are SI rad/s (`geometric 1/m × c`).

## 2. Determinism statement

Same argv + seed = byte-identical CSV/PPM. Rules (all builders):

- Seeded `stc/rng.h` only: mulberry32 (`StcRng32`) and splitmix64
  (`StcRng64`), uniform `[0,1)` plus Box-Muller normals. No `rand`,
  `random`, `time`, `getpid`, pointer hashing, wall-clock, ASLR-derived, or
  uninitialized values in sim code. `NULL` RNG maps to a safe default.
- Fixed iteration counts only. No `while (err > eps)` without a documented
  cap and no error-driven adaptivity anywhere (Verlet, RK4, geodesic, and
  Barnes-Hut all use fixed steps/counts).
- Fixed evaluation/insertion/traversal order; no unstable `qsort` on ties.
- Libs never `exit`/`abort`; they return NaN/identity/zero or an error code
  documented per function. App/verify own the exit code.

## 3. Integrated core: N-body sandbox (`stc/body.h`, `src/body/body.c`)

Direct Newtonian gravity integrated live per step. State: `StcSim`
(`b[STC_BODY_MAX]`, `n`, `t`, `G`, `soft2`, `restitution`, `merge`, `vmax`);
`StcBody` (`name`, `m > 0` kg, `r >= 0` m, `p` m, `v` m/s, `fixed`, `alive`).

### 3.1 Equations

- Plummer-softened acceleration (`stc_sim_accel`): for body `i`,
  `a_i = Σ_{j≠i} G m_j (p_j − p_i) / (|p_j − p_i|^2 + soft2)^(3/2)`,
  with `soft2` in m^2 from `stc_sim_init`. Dead (`alive == 0`) slots are
  skipped; `fixed` bodies exert gravity normally.
- Velocity Verlet, kick-drift-kick (`stc_sim_step`): half-kick with `a(p0)`,
  drift `p`, recompute `a(p1)`, half-kick with `a(p1)`; `t += dt`; then
  `stc_sim_poststep` — one collision pass, then the `vmax` cap (Gate 0:
  the post-step is exported so RK4 callers run the identical pass).
  Fixed bodies never move; dead slots skipped. `n` is kept compacted.
- Collisions, pairwise `dist < r1 + r2`:
  - Merge mode (`merge != 0`): momentum-conserving merge. Mass adds,
    radius is the volume sum `cbrt(r1^3 + r2^3)`, lower index survives
    (other `alive = 0`, then `n` compacted). If either body is fixed, the
    survivor is fixed at the fixed position with `v = 0`.
  - Bounce mode (`merge == 0`): impulse along the line of centers with
    coefficient `s->restitution` (`e ∈ [0,1]`); fixed bodies act as
    infinite mass.
- Speed cap: `vmax > 0` rescales every free-body velocity to `|v| ≤ vmax`
  post-step; `0` disables.
- Total energy (`stc_sim_energy`): kinetic of free bodies only
  `Σ ½ m v^2` plus pairwise Plummer potential
  `−Σ_{i<j} G m_i m_j / sqrt(r^2 + soft2)` over all live pairs.
- Kepler elements (`stc_sim_elements`, about `center`): vis-viva
  semi-major axis `a` from `v^2 = μ(2/r − 1/a)` (`μ = G m_center`),
  eccentricity-vector magnitude `e`, period `T = 2π√(a^3/μ)`; `T = 0`
  when unbound/invalid. Any out-pointer may be `NULL`.
- Helpers: gravitational redshift proper rate
  `stc_body_proper_rate = √(1 − rs/r)`, clamped to `[0,1]`, `0` for
  `r ≤ rs`; Lorentz factor `stc_gamma(v) = 1/√(1−(v/c)^2)`, `-1` for
  `v ≥ c`.

### 3.2 Presets (`stc_preset_load`, seed-deterministic)

SI, XZ plane with y-up, prograde toward +Z. `stc_preset_count()` is 4;
bad index fails without touching `s`.

| # | Name | Content |
|---|---|---|
| 0 | solar | Free Sun `STC_M_SUN` (1.98847e30 kg, Gate-0: was 1.989e30) plus 8 planets with real m/r/a, circular `v`, Sun velocity momentum-compensated |
| 1 | binary | 2 × 1e30 kg at ±1e10 m with `v = √(GM/4d)` plus a 6e24 kg planet at 1e11 m, circular about total mass |
| 2 | hole_disk | FIXED 10 Msun hole + 24 seeded ring particles (seeded radii 1e6–3e6 m, circular) + one eccentric probe |
| 3 | slingshot | FIXED Jupiter-mass planet at origin + probe at (−6e9, 0, 3e9) m with (1e4, 0, 0) m/s |

### 3.3 Approximations and limits (normative)

- Verlet order: velocity Verlet is 2nd-order global (`O(dt^2)`), symplectic
  (bounded energy error on orbits, no secular drift at fixed `dt`), with
  fixed steps and no adaptivity. Accuracy is resolution-limited: the
  `nbody_orbit` check (dt = T/2000, one period) requires radius and energy
  drift `< 1e-3`; eccentric/close encounters need smaller `dt`.
- Softening bias: Plummer `soft2` removes the `1/r^2` singularity but
  systematically weakens close-range gravity (`a → G m d / soft^3` as
  `d → 0`). Orbits with periapsis ≲ `soft` are biased (slower, less bound);
  energy includes the softened potential, so it is self-consistent but not
  the Newtonian value. Use `soft = 0` for exact Kepler tests.
- Fixed-body momentum caveat: `fixed` bodies exert gravity but never move
  (infinite inertia). Total momentum and the free-body kinetic energy are
  therefore not conserved whenever fixed bodies participate: the fixed
  mass acts as an external field/source. Merge with a fixed parent pins
  the survivor (`v = 0`) and drops the incoming momentum by construction.
- Restricted slingshot: preset 3 is a restricted model — the planet is a
  FIXED point mass and the probe is a test particle. There is no back
  reaction, no third body (no Sun), no relativity, and no finite-radius /
  atmosphere effects. It demonstrates bending/energy exchange with a moving
  frame only in the sense of a static bend; quantitative flyby numbers are
  illustrative.
- Capacity/failure: `STC_BODY_MAX` is 64; `stc_sim_add` rejects `m ≤ 0`,
  `r < 0`, non-finite state, and add 65 (returns −1). `stc_sim_remove`
  swap-with-last compacts `n`.

## 4. Integrator option: RK4 (`stc/rk4.h`, `src/body/rk4.c`)

Classic 4th-order Runge-Kutta (Kutta 1901) on the full state `y = (p, v)`,
`dy/dt = (v, a(p))`: `k1 = f(y0)`, `k2 = f(y0 + dt/2·k1)`,
`k3 = f(y0 + dt/2·k2)`, `k4 = f(y0 + dt·k3)`,
`y1 = y0 + dt/6·(k1 + 2k2 + 2k3 + k4)`. Each stage acceleration is the same
Plummer-softened N-body gravity as the Verlet core (via `stc_sim_accel` on
stage-shifted positions), so the force model matches exactly. Respects
fixed/alive; `t += dt`. No collisions and no `vmax` cap inside — the caller
pairs it with the shared post-integration pass: since Gate 0 the app runs
`stc_sim_poststep` (the same merge/bounce collision resolution and `vmax`
cap `stc_sim_step` runs internally) after every `stc_sim_step_rk4`, so
`--integrator rk4` resolves identical collision/cap policy to Verlet and
differs only in stepping scheme. 4th-order global (`O(dt^4)`),
non-symplectic: more accurate per step than Verlet at small `dt` but with
secular energy drift on long integrations. No RNG, no allocation.

## 5. Approximate gravity: Barnes-Hut (`stc/bhut.h`, `src/bhut/bhut.c`)

Barnes-Hut hierarchical octree (Barnes & Hut 1986, Nature 324): far cells
are replaced by total mass at center of mass, near cells opened recursively
down to direct Plummer-softened summation matching `stc_sim_accel`.
`stc_bhut_accel(s, theta, out)` fills `out[n]` (m/s^2); fixed bodies exert
and receive like direct; dead slots get zero. Guidance `theta ≈ 0.7`
(`STC_BHUT_THETA_DEFAULT`); `theta ≤ 0` opens every cell (exact,
direct-equivalent). Target: `< 1%` relative vs direct on a 24-body ring
config. Deterministic fixed insertion/traversal order; per-call malloc/free
(node pool + links), freed before return, with exact-direct fallback on
allocation failure. Approximation error grows with `theta` and with
clustered/close configurations; tree rebuild is per call (no incremental
update).

## 6. Fragmentation (`stc/frag.h`, `src/frag/frag.c`)

`stc_frag_impact(s, a, b, seed)`: destroys parents `a, b` (larger index
removed first) and spawns `N = 2 + ⌊3·log10(1 + E_imp/E_REF)⌋` debris,
`E_REF = 1e20 J`, `E_imp = ½μ|vrel|^2` (center-of-mass impact energy),
capped at free slots after removing the parents (min 2). Mass `M` is split
by weights in `[0.5, 1.5)` normalized to sum exactly to `M`; velocities are
`v_cm` plus isotropic Gaussian kicks with RMS
`σ = 0.1·|vrel| + 0.1·v_esc` (`v_esc = √(2GM/(ra+rb))`), mean-corrected so
the mass-weighted mean kick is zero; positions are `p_cm` plus
uniform-in-sphere offsets of radius `(ra+rb)/2`, likewise mean-corrected;
radii conserve volume `ri = cbrt((ra^3+rb^3)·mi/M)`. All debris free,
named `frag0…`. Deterministic: one splitmix64 stream, fixed draw order.
Momentum is conserved to floating-point roundoff ("approximately": fixed
parents contribute stored velocity though fixed means infinite inertia,
plus roundoff). Returns debris count or −1 (invalid input, `s` untouched).
No I/O.

## 7. Tides (`stc/tides.h`, `src/tides/tides.c`)

Pure functions, no sim mutation, no I/O.

- Tidal acceleration (Newtonian differential gravity, first order in
  `R/r`): `a_tidal = 2 G M R / r^3` along `+r̂` (center-to-center),
  near-side-minus-center to leading order. Valid for `R ≪ r` (drops
  `O((R/r)^2)`); test-mass limit `M_central ≫ m_body`.
- Roche limit (E. Roche 1848): fluid
  `d = 2.44·R_primary·(ρ_primary/ρ_satellite)^(1/3)` (strengthless
  synchronously rotating satellite; the conventional quoted value,
  `stc_roche_limit` == fluid) and rigid
  `d = 1.26·R_primary·(ρ_primary/ρ_satellite)^(1/3)`. Self-gravity only;
  material strength, rotation state, and eccentricity ignored.
- Eccentricity tidal heating (Peale-Cassen-Segatz mechanism; Io/Jupiter
  type example): `dE/dt = (21/2)·(k2/Q)·G·M^2·R^5·n·e^2/a^6`,
  `n = √(GM/a^3)`. ORDER OF MAGNITUDE only: assumes `M ≫ m`, synchronous
  rotation, zero obliquity, small `e` (leading `e^2` term), and
  frequency-independent `k2/Q`. Rocky `k2/Q ~ 1e-3–1e-2`,
  gas giants `~ 1e-5–1e-4`. Circular orbits (`e = 0`) give zero.
- Guards: invalid/non-finite inputs yield NaN (scalars) or NaN-filled /
  zero vectors as documented per function; `NULL out` is a no-op.

## 8. Drag and radiation pressure (`stc/drag.h`, `src/drag/drag.c`)

Pure functions, no sim mutation, no I/O.

- Atmosphere density (barometric formula, isothermal exponential):
  `ρ(h) = ρ0·exp(−h/H)`. Isothermal, constant-g, single-species
  approximation; negative `h` extrapolates above `ρ0`.
- Aerodynamic drag (quadratic drag equation):
  `a = −½·ρ·|v|·Cd·A/m·v` (opposes relative velocity). Zero velocity or
  zero `ρ/Cd/A` gives zero; subsonic/supersonic and rarefied regimes
  share one `Cd`.
- Radiation pressure (radial inverse-square):
  `a = Q·A·L/(4π·r^2·c·m)·r̂` away from the source. `Q = 1` perfect
  absorber, `Q = 2` perfect normal reflector; eclipse/shadowing and
  non-radial reflection ignored.
- Guards: `NULL out` no-op; `NULL` vector input, zero `r`, or
  invalid/non-finite scalars fill `out` with NaN.

## 9. Orbital elements helpers (`stc/elements.h`, `src/elements/elements.c`)

Pure closed forms + truncated series (no iteration, no RNG), SI.

- Hill radius (Hill 1878): `r_H = a·(m/3M)^(1/3)`. Circular orbit,
  `m ≪ M`; true stability boundary is fuzzy (orbits beyond
  `~r_H/2–r_H/3` typically unstable).
- Laplace sphere of influence: `r_SOI = a·(m/M)^(2/5)`.
  Patched-conic heuristic (perturbation-ratio boundary), not a sharp
  dynamical edge.
- Escape velocity: `v_esc = √(2GM/r)`. Newtonian point mass from rest at
  infinity; ignores atmosphere, rotation, other bodies.
- Circular velocity: `v_circ = √(μ/r)`. Newtonian two-body circular.
- Lagrange points (CR3BP, Lagrange 1772): frame `m1` at origin, `m2` at
  `(d,0,0)`, plane `z = 0`; L4 (`+y`) leads, L5 (`−y`) trails. Collinear
  L1/L2/L3 use a truncated series in `α = (μ/3)^(1/3)` through
  `O(α^3) = O(μ)` (L1/L2) and `O(μ)` (L3); triangular L4/L5 are the exact
  equilateral points. Error `< 0.1%` for `μ ≲ 1e-3`, `~0.3%` at
  Earth-Moon `μ = 0.012`, percent-level at `μ = 0.1`; exact only as
  `μ → 0` (L4/L5 exact for all `μ`).
- Guards: `m = 0` yields `0` radius (Hill/SOI); invalid/non-finite inputs
  yield NaN; Lagrange `NULL` outs skipped individually.

## 10. Analytic relativity modules (closed-form / fixed-step)

These modules evaluate textbook GR/SR formulas directly (analytic or
fixed-step integration). They do not step the N-body `StcSim` and are not
coupled back into it; the sandbox core (§3–§6) is the only live N-body
integrator, while §7–§9 helpers are pure functions the caller may apply.

- SR (`stc/sr.h`): Lorentz factor `γ = 1/√(1−β^2)`; x-boost
  `t′ = γ(t − vx/c^2)`, `x′ = γ(x − vt)`; interval
  `s^2 = (c·dt)^2 − dx^2`; Einstein collinear addition
  `w = (u+v)/(1+uv/c^2)`; longitudinal Doppler
  `f_obs/f_emit = √((1−β)/(1+β))`, transverse `f_obs = f0/γ`, redshift
  `z = √((1+β)/(1−β)) − 1`; time dilation, length contraction,
  `E0 = mc^2`, `E = γmc^2`, `K = (γ−1)mc^2` (cancellation-free),
  `p = γmv`, rapidity `φ = artanh β`; twin and barn-pole/ladder paradox
  solvers; caller-bounded worldline `τ` sum. Flat spacetime only;
  velocities guarded subluminal (never NaN; saturates documented limits).
- Geodesic (`stc/geodesic.h`): Schwarzschild equatorial orbit equation in
  `u = 1/r`: `u″ + u = c0 + 1.5·rs·u^2` (photon `c0 = 0`, timelike
  `c0 = GM/L^2`), fixed-step classical RK4 in `φ` (documented caps, no
  adaptivity); photon first integral `(u′)^2 = 1/b^2 − u^2 + rs·u^3`;
  conserved energy/angular momentum/impact parameter; finite-start
  corrected photon deflection; light-cone tilt `(π/4)(rs/r)` rendering
  aid. Exact Schwarzschild (all orders) within step resolution; equatorial
  slice only; capture/escape by radius checks.
- Embedding (`stc/embedding.h`): Flamm paraboloid
  `w = 2√(rs(r−rs))`, `r > rs` (clamped `w = 0` at/inside throat), slope
  `dw/dr = √(rs/(r−rs))`; polar/Cartesian mesh samplers into caller
  buffers (no allocation). Multi-well `H = −Σw_i` is an ILLUSTRATIVE
  superposition for visualization only — not a GR solution.
- Lensing (`stc/lensing.h`): weak-field deflection `α = 4GM/(c^2·b)`;
  Einstein angle `θE = √(4GM/c^2·Dls/(Dl·Ds))`, radius `RE = θE·Dl`;
  thin-lens mapping `β = θ − θE^2/θ`, image positions, signed/total
  magnifications, aligned-ring points. Weak-field thin-lens only; point
  source on caustic gives formally infinite magnification.
- Waves (`stc/waves.h`): Newtonian quadrupole, leading-order (0PN),
  circular adiabatic inspiral: chirp mass, strain
  `h0 = (4/d)(GMc/c^2)^(5/3)(πf/c)^(2/3)` and separation form
  `h = 4G^2·m1·m2/(c^4·a·R)` (identical via Kepler), chirp `f(τ)`,
  `τ(f)`, phase `Φ(τ)`, TT-gauge ring response, inclination-projected
  waveform; seeded `Φ_c` from splitmix64. Valid for `v ≪ c`, `a ≫ rs`;
  `f/Φ` diverge as `τ → 0` (cut off at ISCO in the scenario layer); no
  PN/spin/eccentricity/merger-ringdown.
- Kerr (`stc/kerr.h`): horizons `r± = M ± √(M^2−a^2)`, ergosurface,
  Bardeen-Press-Teukolsky ISCO (6M/1M/9M at `a = 0/+M/−M`), ZAMO
  frame-dragging `w = −g_tφ/g_φφ` in SI rad/s, mean shadow-radius
  estimator `3√3·M(1 − 0.06(a/M)^2)` (parabolic fit, few-percent,
  inclination-averaged; not ray-traced). Geometric meters; `|a| ≤ M`
  by clamp.
- Clocks (`stc/clocks.h`): proper rate `√(1−rs/r)`; redshift; circular
  `√(GM/r)`; SR transverse-Doppler and Schwarzschild GR fractions; GPS
  SR/GR/net μs/day COMPUTED from orbit mechanics (never hardcoded,
  netting ≈ +38); Shapiro PPN one-way and grazing-log delays; Einstein
  perihelion `6πGM/(c^2·a(1−e^2))` with Mercury ≈ 43″/cy via Kepler-III;
  crude Hafele-Keating equatorial model (signs/order right, not the
  observed values).

## 11. Analytic-scenario vs integrated-core distinction

Two independent physics paths share constants/RNG but no state:

- Integrated core (§3–§6): time-steps a live `StcSim` (Verlet or RK4,
  optional Barnes-Hut forces, fragmentation events). State evolves by
  numerical integration; accuracy depends on `dt`, `soft`, `theta`, and
  step count. Owns presets 0–3 above.
- Analytic scenarios (`stc/scenario.h`, 10 presets: twin, pole, doppler,
  mercury, gps, bend, ring, flamm, chirp, shadow): closed-form steppers
  with fixed seeds and frame counts evaluating the §10 formulas per frame
  into `obs[8]` slots. They illustrate textbook results (e.g. Mercury
  43″/cy, GPS +38 μs/day, 1.75″ bending) without integrating an N-body
  state; their numbers are formula evaluations, not emergent dynamics.
- Pure helpers (§7–§9: tides/drag/elements) sit between: exact
  closed-form physics the app may query alongside either path, but they
  never mutate a sim on their own.

Do not compare scenario observables against sandbox trajectories as if
they shared a model: e.g. the `mercury` scenario's Einstein shift is
imposed analytically, while a Mercury-like sandbox orbit is Newtonian +
softening + integrator error. Cross-checks must name both models.

## 12. Non-physics support (no equations of motion)

- Trails (`stc/trails.h`): per-body fixed-capacity ring buffers
  (default depth 128), oldest-to-newest push/iterate/get/count/clear;
  pure data, no render calls, no sim mutation; `init`/`free` paired,
  no leaks on failure.
- Formatting (`stc/fmt.h`): `stc_fmt_time` auto-units (`<120s`,
  `<7200s`, `<172800s`, `<63072000s`, else Julian-year) and `stc_fmt_sci`
  (`%.3g` + NULL-safe unit); `snprintf` into caller buffer, no I/O.
- Render/scenario/app/CSV layers carry or plot physics but add none.

## 13. Verification results

Filled by the integrator from their own run on 2026-09-27
(`make clean && make all`, then the commands below; exit 0 everywhere).

| Check / binary | Command | Computed | Expected | Tol | Result |
|---|---|---|---|---|---|
| sr_invariance (stverify) | `./build/stverify` | γ(0.6c) 1.250000000000, interval relerr 1.604e-16, twin 8.000000/10.000000 s | γ(0.6c)=1.25, interval relerr 0, twin 8/10 s | 1e-12 / 1e-9 / 1e-12 | PASS |
| photon_bending (stverify) | `./build/stverify` | lens 1.751243″, geo 1.751255″, sphere 2/1/0 | 1.75″ (lens + geodesic) | 1% | PASS |
| mercury_precession (stverify) | `./build/stverify` | clocks 42.9821″/cy, geo 42.9821″/cy | 43″/cy | 5% | PASS |
| gps_offset (stverify) | `./build/stverify` | SR −7.1101, GR +45.7198, net +38.6097 μs/day | SR −7.2, GR +45.7, net +38 μs/day | 5% (+ net 1.0 μs) | PASS |
| shapiro_delay (stverify) | `./build/stverify` | EV 1-way 116.284 μs, 2-way 232.569 μs, Mars 2-way 247.248 μs | EV 1-way ~120 μs, 2-way ~233 μs | 10% | PASS |
| gw_chirp (stverify) | `./build/stverify` | 1/d scaling 0.500000000000, ident relerr 4.913e-15, roundtrip 100.000000000 | 1/d scaling 0.5, τ roundtrip 100 s | 1e-12 / 1e-6 / 1e-9 | PASS |
| kerr_isco (stverify) | `./build/stverify` | 6.000000000M / 1.000000000M / 9.000000000M | 6M / 1M / 9M | 1e-6 | PASS |
| nbody_orbit (stverify) | `./build/stverify` | r 1.000000e+07 m, rerr 7.451e-15, eerr 1.968e-15 | rerr 0, eerr 0 (r0 = 1e7 m) | 1e-3 / 1e-3 | PASS |
| momentum_merge (stverify) | `./build/stverify` | speed 1.920e-15, merr 0.000e+00, n 1 | \|v\| 0, merr 0, n 1 | 1e-9 / 1e-12 | PASS |
| rk4_poststep (stverify, Gate 0) | `./build/stverify` | speed 1.870e-15, merr 0.000e+00, n 1 | \|v\| 0, merr 0, n 1 | 1e-9 / 1e-12 | PASS |
| stverify summary | `./build/stverify` | stverify: 10/10 PASS, exit 0 | 10/10 PASS, exit 0 | exact | PASS |
| nbody API/determinism (sttest) | `./build/sttest` | 6/6 PASS: n=64 idx63=63 add65=-1; invalid rejected; remove n=2 moved 3e+10; presets count=4 n=[9 3 26 2]; 200 steps maxdiff 0.000e+00 memcmp identical; energies finite | add-64 ok / 65th fails; invalid rejected; remove compacts; presets 0–3 ok; 200-step redeterminism identical; energies finite; exit 0 | exact/memcmp | PASS |
| solar year (stsolar) | `./build/stsolar` | 8/8 drift < 2% (worst Mars 0.2407%); Earth swept 6.268824 rad err 0.2286%; dE/E 0.0000%; exit 0 | per-planet drift < 2%, Earth angle within 2% of 2π, \|dE/E\| < 0.5%; exit 0 | 2% / 2% / 0.5% | PASS |
| fuzz (stfuzz) | `./build/stfuzz` | ALL_CLEAN cases=200 steps_each=50; exit 0 | ALL_CLEAN (200 seeded sims, finite state); exit 0 | exact | PASS |
| bench (stbench) | `./build/stbench` | N=8: 3508771.9 steps/s; N=16: 884955.8; N=32: 231213.9; N=64: 52015.6 steps/s, 104863459.0 pairs/s; exit 0 | info-only steps/s + pairs/s table; exit 0 | n/a | PASS (info) |
| Barnes-Hut accuracy (integrator /tmp probe vs direct) | /tmp throwaway (hub + 23 ring sats, theta 0.7) | max relerr 4.110e-16 over 23 sats; hub excluded (exact accel 1.17e-23 by symmetry, abs diff 9.87e-24) | max relative error vs direct < 1% on 24-body ring | 1% | PASS |

How to reproduce (integrator/verifier, after `make clean && make all`):

```sh
./build/stverify          # expect: stverify: 10/10 PASS, exit 0
./build/sttest            # expect: exit 0, all API/determinism lines pass
./build/stsolar           # expect: exit 0, per-planet table within tolerances
./build/stfuzz            # expect: ALL_CLEAN, exit 0
./build/stbench           # expect: exit 0, throughput table (info only)
```

Build flags: `cc -std=c11 -Wall -Wextra -Werror -O2 -I include`, link
`-lm`. The §13 table above was filled from the integrator's own runs of
these exact commands (all exit 0); re-run them to reproduce.
