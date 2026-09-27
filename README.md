# Spacetime Lab

Pure-C11 relativity simulator + N-body gravity sandbox. Zero dependencies
(C stdlib + libm only), SI doubles, seeded determinism: same argv =
byte-identical CSV/PPM.

One binary (`build/stlab`), four modes: the 10 analytic relativity
scenarios; `--sandbox`, a full 3D N-body universe sandbox with 4D
spacetime readouts (character-grid viewport or a macOS native window
via `--visual`; spawn any body, 4 seeded presets, time warp, CSV + PPM
export); `--space`, a real Solar System simulator built on JPL
Horizons ephemerides; and `--merger`, a cinematic binary black hole
coalescence view. Double-click `stlab.command` to
build (if needed) and launch the sandbox.

## Run

```sh
make all            # builds all eight binaries (see below)
make verify         # stverify 11/11 + sttest + stsolar + stfuzz +
                    # stspace + stephem + stlab smokes
make test           # sttest only (N-body API/determinism)
make solar          # stsolar only (1-year solar accuracy)
make fuzz           # stfuzz only (seeded fuzz sweep)
make bench          # stbench only (scaling benchmark, info only)
make ephem-accuracy # stephem --accuracy: J2000->today vs Horizons (slow,
                    # not part of make verify)
make smoke          # headless 60 frames x 10 scenarios, CSV+PPM to /tmp
make sandbox-smoke  # headless 60 frames x 4 sandbox presets to /tmp/sbx
make run            # headless 60-frame twin demo, no output files
make sandbox        # interactive sandbox, preset 0 (solar)
make merger         # BBH merger view, windowed (macOS)
make merger-shot    # headless merger still -> build/merger/shot.ppm
make space          # space simulator, windowed (macOS)
make space-shot     # headless space still -> build/space/shot.ppm
```

Binaries: `build/stlab` (scenarios + sandbox + space + merger),
`build/stverify` (11 acceptance checks), `build/sttest` (N-body
API/determinism), `build/stsolar` (1-year solar accuracy),
`build/stfuzz` (seeded fuzz), `build/stbench` (scaling benchmark, info
only), `build/stspace` (space-core dynamics checks), `build/stephem`
(ephemeris catalog/calendar checks; `--accuracy` runs the slow JPL
comparison).

Direct runs:

```sh
./build/stlab --list
./build/stlab --scenario all --frames 60 --ppm /tmp --csv /tmp --smoke --no-adapt
./build/stlab --scenario twin --frames 600   # interactive when stdin is a TTY
./build/stlab --sandbox --preset 0           # interactive sandbox (solar)
./build/stlab --sandbox --preset 0 --visual  # native window (macOS only)
./build/stlab --sandbox --preset 0 --frames 60 --visual-shot /tmp/vis.ppm  # headless 1280x800 PPM
./build/stlab --sandbox --preset 2 --frames 60 --ppm /tmp/h.ppm --csv /tmp/h.csv --smoke --no-adapt
./build/stverify && ./build/sttest && ./build/stsolar && ./build/stfuzz && \
    ./build/stspace && ./build/stephem && ./build/stbench
```

Build flags: `cc -std=c11 -Wall -Wextra -Werror -O2 -I include`, link
`-lm` (`-lpthread` off macOS).

### Space simulator (--space)

```sh
./build/stlab --space                    # windowed (macOS), preset 0
./build/stlab --space --preset 1         # today's Solar System
./build/stlab --space --preset 0 --frames 300 --visual-shot /tmp/s.ppm
./build/stlab --space --preset 1 --frames 1200 --movie /tmp/frames --movie-every 4
```

A real Solar System simulator: all catalogued bodies start from NASA
JPL Horizons state vectors and evolve under Newtonian N-body gravity
with a 1PN correction, J2/J4 zonal harmonics, collisions, tidal
disruption, and black-hole capture -- plus a sandbox console (replace
the Sun with a black hole, inject a rogue, spawn Earth swarms) and
checkpoint-exact time travel in both directions. Six presets (0-5);
windowed on macOS by default, `--visual-shot`/`--movie` run headless.
Full documentation: `doc/SPACE.md`.

### Black hole merger view (--merger)

```sh
./build/stlab --merger                   # windowed (macOS), preset 0
./build/stlab --merger --preset 0 --visual-shot /tmp/m.ppm
./build/stlab --merger --preset 0 --frames 2400 --movie /tmp/f --movie-every 2
tools/merger_movie.sh 0 2400 2           # frames + MP4 + GIF via ffmpeg
```

A cinematic edge-on view of a binary black-hole coalescence: a Peters
0PN inspiral timeline, shaped plunge, damped-sinusoid ringdown, and
NR-calibrated remnant/QNM fits drive a rendered spacetime sheet whose
ripples are the retarded-time quadrupole strain. Four presets (0-3);
windowed on macOS by default, `--visual-shot`/`--movie` run headless.
Full documentation: `doc/MERGER.md`.

## Controls

CLI flags (`--help` documents all): `--scenario <name|all>`,
`--frames N` (0..1000000, default 600, smoke default 1),
`--ppm <path>` (file for one scenario, directory when scenario=all;
sandbox always takes a file), `--csv <path>` (same rule),
`--sandbox` (N-body sandbox mode), `--merger` (BBH merger view),
`--space` (Solar System simulator), `--preset N` (sandbox preset 0..3:
solar, binary, hole_disk, slingshot; merger preset 0..3; space
scenario 0..5; default 0), `--seed N` (sandbox
RNG seed, default 20260927), `--integrator verlet|rk4` (sandbox
stepper, default verlet; both run the same collision/speed-cap
post-step, only the stepping scheme differs), `--visual` (sandbox in a
native window, macOS only; accepted but not required by --merger and
--space, which window by default), `--visual-shot <ppm>` (headless visual
render: one 1280x800 P6 PPM, no window), `--movie <dir>`
(--merger/--space headless PPM frame sequence), `--movie-every N`
(movie frame stride, default 2), `--smoke` (headless),
`--list` (the 10 preset names), `--help`, `--no-adapt` (disable 60fps
pacing).

Scenario keys (TTY, raw mode): WASD/arrows orbit camera, +/-/R/F zoom,
1-9,0 switch scenario, P pause, E export current frame, Q quit
(Ctrl-C/Esc also quit).

Sandbox keys: N spawn form (name, mass_kg, radius_m, px/py/pz_m,
vx/vy/vz_m_s, fixed; numbers accept `1.989e30` syntax; form prefills a
circular orbit about body 0; Up/Down move, Enter on SPAWN> commits, Esc
cancels), K remove followed body, Tab cycle follow (body list shows
mass/speed/proper-time rate `taurate`; followed body shows a/e/T),
G cycle preset,
Z/X halve/double dt (1e-6..1e9 s), C warp x1/x10/x100, T trails
on/off, P pause, E export CSV+PPM, ? key guide, Q quit. Camera:
WASD/arrows orbit, +/- zoom. Header shows the stepper as
integ=verlet|rk4.

- Visual modes add D: an ILLUSTRATIVE rubber-sheet grid in the orbital
  plane extending toward the horizon (potential-shaped view-scale
  dent, geometric line spacing far out, not the spacetime metric)
  plus an x/y/z axis triad and a 100 px scale bar (m/km/AU
  auto-units). In visual mode D takes the WASD yaw-right binding; the
  Right arrow still yaws.

User guide with spawn recipes: `doc/GUIDE.md`. Physics equations:
`doc/PHYSICS.md`. Man page: `doc/stlab.1`.

## Physics reference (equation names)

- body: velocity-Verlet N-body core (StcSim/StcBody, up to 64 bodies),
  Plummer softening, merge-or-bounce collisions, optional speed cap,
  fixed anchors, total energy, vis-viva Kepler elements (a/e/T from
  state vectors), per-body proper rate + Lorentz gamma. 4 seeded
  presets (solar, binary + circumbinary planet, black-hole disk,
  slingshot flyby).
- rk4: classic 4th-order Runge-Kutta stepper alternative for the core
  (`--integrator rk4`; the stepper itself holds no collision logic --
  the app runs the shared post-step merge/cap pass after each rk4
  step, identical policy to Verlet).
- bhut: Barnes-Hut octree approximate gravity (opening angle theta,
  default 0.7; theta <= 0 is direct-equivalent).
- frag: energy-scaled fragmentation debris on destructive collisions.
- tides: tidal acceleration, Roche limits (fluid + rigid), eccentricity
  heating order-of-magnitude dE/dt.
- drag: exponential-atmosphere drag + radiation pressure.
- elements: Hill radius, Laplace sphere of influence, escape/circular
  velocity, Lagrange points (collinear series + exact triangular).
- trails/fmt/help: per-body trajectory ring buffers (depth 128);
  auto-unit time + "%.3g unit" formatting; static key-guide/splash
  strings for the sandbox screen.
- core-math: column-major Mat4 lookAt/perspective, Schwarzschild radius
  rs=2GM/c^2, geometric mass GM/c^2, mulberry32/splitmix64 RNG.
- sr: Lorentz factor gamma=1/sqrt(1-v^2/c^2); x-boost Lorentz transform;
  interval s^2=(c dt)^2-dx^2; Einstein velocity addition; longitudinal and
  transverse Doppler; redshift z; time dilation; length contraction;
  E0=mc^2, E=gamma mc^2, K=(gamma-1)mc^2, p=gamma m v; rapidity artanh(v/c);
  twin paradox; barn-pole/ladder paradox; worldline proper time.
- geodesic: relativistic Binet orbit equation u''+u=c0+1.5 rs u^2
  (photon c0=0, timelike c0=GM/L^2); photon first integral;
  classical RK4 in phi; Schwarzschild energy equation; light-cone tilt
  (pi/4)(rs/r).
- embedding: Flamm paraboloid w=2 sqrt(rs (r-rs)); Flamm slope;
  illustrative multi-well superposition H=-sum w_i (visualization only,
  not a GR solution); polar and Cartesian mesh sampling.
- lensing: Einstein weak-field deflection alpha=4GM/(c^2 b); Einstein
  angle thetaE; thin-lens mapping beta=theta-thetaE^2/theta; image
  positions; signed and total magnifications; Einstein ring.
- waves: chirp mass Mc; circular-binary quadrupole strain amplitude;
  separation-form strain; chirp frequency f(tau); time to coalescence
  tau(f); inspiral phase Phi(tau); TT-gauge ring response; chirp waveform
  with inclination projection. Newtonian 0PN circular adiabatic inspiral.
- kerr: outer/inner horizon r+-=M+-sqrt(M^2-a^2); ergosurface;
  Bardeen-Press-Teukolsky ISCO (6M/1M/9M limits); ZAMO frame-dragging
  rate w=-g_tphi/g_phiphi; parabolic mean shadow-radius estimator.
- clocks: Schwarzschild proper rate sqrt(1-rs/r); gravitational redshift;
  circular velocity sqrt(GM/r); SR transverse-Doppler fraction; GR rate
  ratio; GPS SR/GR/net us/day (computed); Shapiro PPN one-way and
  grazing-log delay; Einstein perihelion shift; Hafele-Keating model.
- scenario: 10 presets (twin, pole, doppler, mercury, gps, bend, ring,
  flamm, chirp, shadow) with fixed seeds and frame counts; closed-form
  steppers per preset.
- render: pinhole camera VP=perspective*lookAt; clip/NDC/viewport project;
  Bresenham lines; edge-function flat triangles; flat Lambert shading;
  PPM P6 raster; ANSI SGR escapes; phosphor/amber ramps.
- par: persistent POSIX worker pool behind `stc_par_for` (dynamic
  chunks, caller participates, per-index slots only so results are
  bitwise thread-count independent); backs the space core's force and
  substep loops and the space view's row passes.
- space: large-N dynamics core (--space): SoA massive + tracer
  particles; fixed-substep 6th-order Yoshida composition of KDK
  leapfrog; optional 1PN Schwarzschild test-particle term from
  STAR/BLACKHOLE sources via implicit-midpoint kicks; J2/J4 zonal
  harmonics with momentum-conserving reaction; merge/absorb/capture/
  tidal-disruption event rules; Barnes-Hut octree above a body-count
  threshold; checkpoint ring + deterministic seek for exact rewind.
- ephem: JPL Horizons Solar System catalog (34 entries x 2 epochs,
  barycentric J2000 ecliptic, EXTERNAL data generated into
  src/space/ephem_data.c), planet-vs-system-barycentre loading rule,
  seeded main-belt tracers, TDB <-> proleptic-Gregorian calendar.
- merger: closed-form BBH inspiral-merger-ringdown timeline in
  geometric units (Peters 0PN inspiral, shaped plunge, damped-sinusoid
  ringdown, NR-fit remnant mass/spin + QNM), retarded-time strain
  accessors; checked by stverify's bbh_timeline.
- vis_space / vis_merger: front ends for --space and --merger (macOS
  window via vis_mac, headless via --visual-shot/--movie); the merger
  view renders the retarded GW strain onto a sheet + thin-lens
  background (ILLUSTRATIVE visuals over the merger timeline).

## Module map

| Module | Header | Source |
|---|---|---|
| core-math | include/stc/vec.h, units.h, rng.h, fmt.h | src/core/vec.c, units.c, rng.c, fmt.c |
| sr | include/stc/sr.h | src/sr/sr.c |
| geodesic | include/stc/geodesic.h | src/gr/geodesic.c (see note) |
| embedding | include/stc/embedding.h | src/embedding/embedding.c |
| lensing | include/stc/lensing.h | src/lensing/lensing.c |
| waves | include/stc/waves.h | src/waves/waves.c |
| kerr | include/stc/kerr.h | src/kerr/kerr.c |
| clocks | include/stc/clocks.h | src/clocks/clocks.c |
| scenario | include/stc/scenario.h | src/scenario/scenario.c |
| render | include/stc/render.h | src/render/render.c |
| body | include/stc/body.h | src/body/body.c |
| rk4 | include/stc/rk4.h | src/body/rk4.c |
| bhut | include/stc/bhut.h | src/bhut/bhut.c |
| tides | include/stc/tides.h | src/tides/tides.c |
| drag | include/stc/drag.h | src/drag/drag.c |
| frag | include/stc/frag.h | src/frag/frag.c |
| elements | include/stc/elements.h | src/elements/elements.c |
| par | include/stc/par.h | src/core/par.c |
| space | include/stc/space.h | src/space/space.c |
| ephem | include/stc/ephem.h | src/space/ephem.c, src/space/ephem_data.c (generated by tools/gen_ephem.py from data/horizons/) |
| merger | include/stc/merger.h | src/merger/merger.c |
| spawn | include/stc/spawn.h | src/app/spawn.c |
| vis | include/stc/vis.h | src/vis/vis_render.c, vis_mac.c |
| vis_space | include/stc/vis_space.h | src/vis/vis_space.c |
| vis_merger | include/stc/vis_merger.h | src/vis/vis_merger.c |
| app | include/stc/app.h, ui.h, csv.h, sandbox.h, trails.h, help.h | src/app/main.c, ui.c, csv.c, sandbox.c, trails.c, help.c -> build/stlab (links all modules incl. vis_*) |
| verify | (all stc headers) | src/verify/verify.c -> build/stverify; nbody_test.c -> build/sttest; solar_check.c -> build/stsolar; fuzz.c -> build/stfuzz; bench.c -> build/stbench; space_check.c -> build/stspace; ephem_check.c -> build/stephem |

Integration notes: canonical sources only (src/gr/geodesic.c per
CONTRACT Sec 8, src/kerr/kerr.c, src/scenario/scenario.c); byte-identical
duplicates were removed on 2026-09-27. The app drives the scenario layer
directly (setup/step per frame with deterministic preset looping) and
renders through the renderer (ASCII live view, PPM export); CSV carries
frame, sim time, head position, and all 8 observable slots per preset.
Sandbox smoke CSV carries frame, sim time, total energy, and body count
per frame; sandbox PPM is a P6 320x200 projection of the live viewport.

## Verification (from `make verify`, this run)

| Check | Computed | Expected | Tol | Result |
|---|---|---|---|---|
| sr_invariance | gamma(0.6c) 1.25, interval relerr 1.6e-16, twin 8/10 s | 1.25, 0.0, 8/10 s | 1e-12/1e-9/1e-12 | PASS |
| photon_bending | lens 1.751243", geo 1.751255" | 1.75" | 1% | PASS |
| mercury_precession | 42.9821"/cy | 43"/cy | 5% | PASS |
| gps_offset | SR -7.1101, GR +45.7198, net +38.6097 us/day | -7.2/+45.7/+38 | 5% | PASS |
| shapiro_delay | EV 1-way 116.284 us, 2-way 232.569 us | ~120/~233 us | 10% | PASS |
| gw_chirp | 1/d scaling 0.5, roundtrip 100.0 s | 0.5, 100 s | 1e-12/1e-9 | PASS |
| kerr_isco | 6M / 1M / 9M | 6M/1M/9M | 1e-6 | PASS |
| nbody_orbit | r 1e7 m, rerr 7.5e-15, eerr 2.0e-15 | r0, 0.0, 0.0 | 1e-3/1e-3 | PASS |
| momentum_merge | speed 1.9e-15, merr 0.0, n 1 | 0.0, 0.0, 1 | 1e-9/1e-12 | PASS |
| rk4_poststep | rk4 step + shared poststep merge: speed 1.9e-15, merr 0.0, n 1 | 0.0, 0.0, 1 | 1e-9/1e-12 | PASS |
| bbh_timeline | Peters a(t) relerr 0, omega relerr 1.2e-10, mf 0.9546, af 0.6865, phase-seam jumps 3.2e-10/5.5e-10 | 0, 0, mf in (0.94,0.96), af in (0.66,0.71), 0, 0 | 1e-12/1e-6/range/range/1e-6/1e-6 | PASS |
| sttest (API/determinism) | 6/6: 64-cap, invalid rejected, remove compacts, 4 presets, 200-step identical, energies finite | exact/memcmp | exact | PASS |
| stsolar (1-year) | 8/8 drift < 2%, Earth angle err 0.2286%, dE/E 0.0% | 2%/2%/0.5% | 2%/2%/0.5% | PASS |
| stfuzz (seeded fuzz) | ALL_CLEAN, 200 cases x 50 steps | finite state | exact | PASS |
| stspace (space core) | 15 checks: kepler 1.2e-11 AU, order 5.98, reversible 2.6e-12, Mercury 42.98"/cy (Newtonian -1.7e-4), merge, BH disrupt+capture, seek bitwise both ways, J2 node -2.608e-8 vs -2.605e-8 rad/s, J2 momentum 3e-15, 1PN substep-indep 0.27 km/10 yr, serial-vs-parallel identical state hash, tree force err median 6.5e-4 | per-check | per-check | PASS |
| stephem (ephemeris) | 11/11: calendar x5 (2000-trial roundtrip worst 0.0 s), catalog 34 entries, loads 9/27/34, barycentre residual 1.5e-16, Earth J2000 sample, two-epoch distance, 500 belt tracers in-band | per-check | per-check | PASS |
| stbench (benchmark) | N=64: 52016 steps/s, 1.05e8 pairs/s | info only | n/a | PASS (info) |

`stverify: 11/11 PASS` (`bbh_timeline` added with the merger model),
all test binaries exit 0. Full measured table:
`doc/PHYSICS.md` §13.

`make ephem-accuracy` (not in `make verify` — too slow) ran
`./build/stephem --accuracy` this run: 12/13 PASS — `accuracy_core`
9/9 and `accuracy_all` 33/34, with Titania at 20,783 km just over its
20,000 km gate (residual Uranus-system gravity error; tolerance not
loosened, see `doc/SPACE.md` "Measured accuracy").
