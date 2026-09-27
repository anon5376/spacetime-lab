# SPACE.md — the space simulator (`stlab --space`)

`--space` is a real Solar System simulator: every catalogued body starts
from an actual NASA JPL Horizons state vector, and the engine then
integrates Newtonian N-body gravity with a 1PN correction, J2/J4 zonal
harmonics, collisions, tidal disruption, and black-hole capture. It is
also a sandbox: a command console can replace the Sun with a black hole,
inject a rogue interloper, or spawn a million test Earths, and the clock
runs both directions with checkpoint-exact rewind.

Run it:

```sh
make space                                   # windowed (macOS)
./build/stlab --space --preset 1             # today's Solar System
./build/stlab --space --preset 0 --frames 300 --visual-shot /tmp/s.ppm
./build/stlab --space --preset 1 --frames 1200 --movie /tmp/frames \
    --movie-every 4                          # PPM frame sequence
```

`--preset N` selects the scenario (0–5). Without `--visual-shot` or
`--movie` the simulator opens a window via `stc_svis_run`, which needs
macOS (AppKit); on other platforms use the headless flags. On-screen
dates are TDB, not UTC (see "What this is not").

## Scenarios

| key | `--preset` | name (as shown in the HUD) |
|-----|-----------|-----------------------------|
| 1 | 0 | `SOLAR SYSTEM  J2000.0` |
| 2 | 1 | `SOLAR SYSTEM  2026-09-27` |
| 3 | 2 | `SUN REPLACED BY 1 MSUN BLACK HOLE` |
| 4 | 3 | `ROGUE 10 MSUN BLACK HOLE FLYBY` |
| 5 | 4 | `EARTH SWARM  2000 EARTHS AT 1 AU` |
| 6 | 5 | `BLACK HOLE + 1,000,000 TEST EARTHS` |

- **0 — Solar System at J2000.0.** Full catalog (34 massive bodies: Sun,
  8 planets, Pluto + Ceres, 4 asteroids, 19 moons) plus 5,000 sampled
  main-belt tracers, at catalog epoch 0 (JD 2451545.0 TDB). Substep
  h = 1800 s, warp 10 d/s.
- **1 — Solar System today.** Same loadout at catalog epoch 1
  (2026-09-27 00:00 TDB, JD 2461310.5).
- **2 — Sun replaced by a black hole.** Scenario 1 with the Sun swapped
  for a 1 M_sun black hole (Schwarzschild radius ≈ 2.95 km). Orbits are
  unchanged to leading order — same exterior monopole (Birkhoff); only
  the 1PN source strength and the capture/disruption event rules differ.
- **3 — Rogue flyby.** CORE planets + dwarf planets (no moons), 3,000
  belt tracers, and a 10 M_sun black hole inbound from
  (−60, +15, +22) AU at (45, 0, −14) km/s. h = 6 h, warp 0.5 yr/s,
  camera pulled back to 90 AU.
- **4 — Earth swarm.** CORE planets plus 2,000 real Earth-mass bodies in
  a near-circular ring at 1 AU (r ±10%, i < 2°, seeded). h = 4 h,
  warp 5 d/s.
- **5 — Million test Earths.** Empty space: a 10 M_sun black hole at the
  origin, 1,000,000 massless tracers on a disk around 1 AU, and 12 real
  Earths in a ring at 0.6 AU (which interact, collide, and get captured
  or disrupted). h = 3600 s, warp 2 d/s. The tracers are ILLUSTRATIVE —
  a sampled population, not named bodies.

## Keys

Windowed mode (`stc_svis_key`, plus the mouse handlers in
`src/vis/vis_mac.c`):

| input | effect |
|-------|--------|
| `1`–`6` | load scenario 0–5 |
| `Space` | pause / resume |
| `,` `.` | step down/up the warp ladder: 1, 60, 3600, 21600, 1 d, 7 d, 30 d, 1 yr, 10 yr, 100 yr per displayed second |
| `B` | reverse time (negative warp) |
| `[` `]` | jump one warp interval back/forward in time (checkpoint restore when possible, otherwise integration) |
| `Tab` | cycle camera focus through stars, planets, dwarf planets, black holes |
| `F` | focus the camera on the selected body |
| `0` | free the camera (back to the barycentre) |
| `O` | toggle osculating-orbit traces |
| `L` | toggle name labels |
| `G` | toggle the ecliptic reference grid + scale bar |
| `T` | toggle tracers |
| `H` | toggle the HUD |
| `W` `S` `A` `D` or arrow keys | orbit camera pitch / yaw |
| `+` `=` / `-` `_` | zoom in / out |
| `/` or `:` | open the command console (`Esc` closes, `Enter` executes, `Backspace` edits; sim pauses while open) |
| left-drag | orbit camera |
| scroll / pinch | zoom |
| click (< 3 px travel) | select the nearest body (readout panel + reticle) |
| `Q`, `Esc`, `Ctrl-C` | quit |

The HUD shows the TDB calendar date and JD, warp rate, integrator
substep, live massive/tracer counts, 1PN on/off, and the relative energy
drift dE/E since scenario load. When the requested warp needs more
substeps than the per-frame budget, the status shows `LIMITED` in amber
— sim time then advances slower than the warp implies, never by
skipping physics.

## Console commands

`sv_exec` parses these (case-insensitive verbs; body names match
case-insensitively):

| command | effect |
|---------|--------|
| `HELP` or `?` | print the command summary in the status line |
| `PRESET n` / `SCENARIO n` | load scenario 0–5 |
| `DATE YYYY-MM-DD [HH:MM]` / `GOTO` | travel to that TDB calendar date: restores the nearest checkpoint if the time is inside the checkpointed span, else integrates forward or backward |
| `WARP <dur>` | set the warp directly; `<dur>` is a number with suffix `y` (Julian years), `d`, `h`, `m` (minutes), or bare seconds; negative runs backward (`WARP -1y`) |
| `BH <name> [msun]` | turn a body into a black hole — radius becomes 2GM/c²; optional argument sets the mass in solar masses (default: keep the body's GM) |
| `MASS <name> <x>` | multiply the body's GM by `x` |
| `SPAWN EARTH <n> [RING\|DISK\|CLOUD\|SHELL] [au]` | spawn `n` Earths about the focused body (or the most massive body): `RING` near-circular at `au`±10%, i<2°; `DISK` flat 0.3–1.5 au; `CLOUD` isotropic inclinations, 0.5–1.5 au; `SHELL` isotropic at `au`. Up to 3,000 spawn as fully interacting Earth-mass bodies; above that they are massless tracers. Limit 5,000,000 |
| `ROGUE [msun] [km/s] [au]` | inject an interloper at (−150, b, 5) AU moving +x at `km/s`, aiming `b` AU off the barycentre; mass < 3 M_sun spawns a star (`ROGUE STAR`, radius ∝ m^0.8), ≥ 3 M_sun a black hole. Defaults 1 M_sun, 30 km/s, 10 AU |
| `FOCUS <name>` | camera follows that body and selects it |
| `DEL <name>` / `DELETE` | remove a body from the simulation |
| `PN ON` / `PN OFF` | toggle the 1PN correction |
| `STEP <dur>` | set the integrator substep `h` (e.g. `STEP 30m`) |

## Data provenance

Initial conditions are `EXTERNAL` data; nothing is typed in by hand.

- **State vectors, GM, radii.** NASA JPL Horizons
  (`ssd.jpl.nasa.gov/horizons`), fetched 2026-09-27 by
  `tools/fetch_horizons.sh` — every query URL and UTC fetch time is in
  `data/horizons/MANIFEST.md`, and the raw responses are kept under
  `data/horizons/`. Horizons sources seen in the fetched set: DE441
  (planets/Sun), JPL#129/179/36/48/74 (asteroids/dwarfs), and the
  satellite-ephemeris kernels jup365, mar099, nep098, plu060, sat441l,
  ura184. `tools/gen_ephem.py` compiles them into
  `src/space/ephem_data.c` (generated; do not hand-edit). Two epochs per
  body: **e0 = J2000.0** (JD 2451545.0 TDB) and **e1 = 2026-09-27
  00:00 TDB** (JD 2461310.5), barycentric J2000 ecliptic, km→m.
- **System barycentres.** When moons are not loaded, each planet is
  replaced by its Horizons *system* barycentre vector (ids 1–9), so the
  outer dynamics still see planet+moons mass. Horizons reports no GM for
  a barycentre; `gm_system` is the sum of the fetched component GMs
  (major moons only — a documented approximation).
- **Zonal harmonics (J2, J4) and spin poles.** From the NAIF SPICE
  satellite-kernel comment files matching the kernels that produced each
  planet's moon vectors — jup365 (Jupiter), sat441 (Saturn), ura184
  (Uranus), nep098_part-1 (Neptune) — plus Earth's J2 (IERS 2010) from
  the fetched Horizons response. Pole RA/Dec (constant term only) from
  the generic PCK `pck00011.tpc`, rotated ICRF-equatorial → J2000
  ecliptic via the IAU76 obliquity 84381.448″. Exact source lines are
  recorded in `data/naif/MANIFEST.md`. Zonal terms are applied only when
  the MOONS group is loaded; a system-barycentre planet stays a point
  mass.
- **The Uranus choice — URA182.21, not URA184.** `ura184.cmt` bundles
  several historical solutions. Its headline block gives J2 = 4.036e-3,
  but integrating this model (point-mass + J2/J4 core, J2000 →
  2026-09-27) against the fetched epoch-1 vectors misses Miranda by
  247,000 km. The older `URA182.21` block in the same file (J2 =
  3.509e-3) lands the same moons at 2.6–4.6 thousand km except Titania
  (see the accuracy table). Horizons' Uranian-moon vectors are evidently
  consistent with ≈3.51e-3, so `gen_ephem.py` parses the URA182.21
  block. The full comparison table and reasoning are in
  `data/naif/MANIFEST.md` §EPHEM-3 and the `ephem_data.c` header.
- **Dropped bodies.** Eris, Makemake, Haumea, and Davida were requested
  for the catalog but dropped: the Horizons small-body lookup returned
  `GM= n.a.` for all four (Eris/Makemake/Haumea also `RAD= n.a.`), and
  rule R8 forbids substituting literature values that were not actually
  fetched. The DWARFS group therefore contains Pluto (with Charon as a
  moon) and Ceres; ASTEROIDS contains Vesta, Pallas, Hygiea, Interamnia.

## Dynamics

All positions/velocities are barycentric, J2000 ecliptic axes; time is
TDB seconds since J2000.0. SI doubles throughout; same inputs and thread
count produce bitwise-identical state.

- **Integrator.** Fixed-substep 6th-order Yoshida composition (solution
  A) of kick-drift-kick leapfrog — 8 kicks + 7 drifts per substep,
  force-evaluated 7 times per substep in steady Newtonian operation via
  kick caching. `h` is the max substep (default 3600 s; scenarios set
  1800 s–6 h). Negative `dt` integrates backward.
- **1PN correction** (`pn`, default on): for every STAR or BLACKHOLE,
  the harmonic-gauge Schwarzschild test-particle term
  `a += GM/(c²r³)[(4GM/r − v²) r + 4(r·v) v]` applies to every other
  particle, in the source's frame. The velocity dependence is handled by
  an implicit-midpoint fixed-point solve inside each kick (so the
  composition keeps its order); reaction on the source is neglected.
  Measured: Mercury's perihelion advance comes out at 42.98″/century.
- **Zonal J2/J4** per body with oblateness set (Earth, Jupiter, Saturn,
  Uranus, Neptune): axisymmetric multipole acceleration about the body's
  spin pole, felt by massive bodies and tracers alike; the reaction on
  the oblate body is applied for massive targets, so total momentum is
  conserved (checked to 3e-15).
- **Gravity backend.** Direct pairwise summation while the live massive
  count is ≤ `tree_threshold` (default 1024); above it, a monopole
  Barnes–Hut octree rebuilt each force evaluation, opening angle
  `theta` = 0.5 against the nearest-point test (more conservative than
  classic size/d). Tracers are direct-over-massive below the threshold
  and tree-sampled above it. Measured on a 5,000-body Plummer sphere:
  median relative force error 6.5e-4, p99 3.7e-3.
- **Events** (`collisions`, default on; applied after every substep):
  massive-massive contact → perfectly inelastic merge conserving mass,
  momentum, centre of mass (radius = summed-volume sphere, or
  Schwarzschild radius if either is a black hole); a body within
  2 r_s of a black hole → CAPTURE; a massive body inside the tidal
  radius R(M_bh/m)^{1/3} but outside capture → DISRUPT, replaced by
  `disrupt_tracers` (default 400) tracers spread through its volume with
  its velocity; tracers are absorbed on contact or captured inside
  2 r_s. All events are logged (`MERGE`/`IMPACT`/`CAPTURE`/`TIDAL
  DISRUPTION`, shown in the HUD event log).
- **Checkpoints and seek.** Full-state snapshots on an absolute time
  grid (`ckpt_interval`, default 30 d; ring of `ckpt_max`, default 64),
  taken at exact substep boundaries. `stc_space_seek(t)` restores the
  newest checkpoint ≤ t and re-integrates deterministically — so
  scrubbing backward across a merge or capture reproduces the pre-event
  state bitwise. This is the "exact time travel" path; it exists because
  the velocity-dependent 1PN term breaks exact backward integration at
  ~1e-8 relative per orbit.
- **Tracers** are massless test particles: they feel massive bodies,
  never exert force, and never interact with each other.

## Measured accuracy

`./build/stephem --accuracy` (`make ephem-accuracy`) loads the catalog
at J2000, integrates with the production `stc_space_advance` (pn=1) to
epoch 1 = 2026-09-27, and compares every body to its fetched Horizons
vector. Real output from this build (M5 Pro, 2026-09-27; dt = 26.736 yr):

```text
-- accuracy_core: groups=0x1 h=21600s pn=1 bodies=9 dt=26.736yr --
body             err_km     tol_km   arcsec
Sun               129.0     1000.0      n/a
Mercury           126.8     1000.0    0.374
Venus             129.2     1000.0    0.245
Earth            3133.6     5000.0    4.310
Mars             121.5      1000.0    0.108
Jupiter           133.6     1924.0    0.035
Saturn            127.7     3421.8    0.019
Uranus            112.8     7050.4    0.008
Neptune           105.2    10834.9    0.005
PASS accuracy_core: 9/9 bodies within tolerance

-- accuracy_all: groups=0xf h=1800s pn=1 bodies=34 dt=26.736yr --
body             err_km     tol_km   arcsec
Sun               115.8     1000.0      n/a
Mercury           114.2     1000.0    0.337
Venus             116.4     1000.0    0.221
Earth             115.8     1000.0    0.159
Mars              120.7     1000.0    0.107
Jupiter            49.6     1924.0    0.013
Saturn            248.2     3421.8    0.036
Uranus           2609.5     7050.4    0.185
Neptune          1998.9    10834.9    0.092
Pluto             268.5     2000.0    0.010
Ceres             164.0     2000.0    0.084
Vesta             122.0     2000.0    0.069
Pallas            115.6     2000.0    0.056
Hygiea             83.1     2000.0    0.036
Interamnia        108.2     2000.0    0.045
Moon              166.4     2000.0   90.748
Io               3746.5    20000.0 1838.495
Europa           1830.3    20000.0  567.573
Ganymede          686.6    20000.0  131.988
Callisto           35.5    20000.0    3.873
Mimas           11577.7    20000.0 13082.930
Enceladus        1804.6    20000.0 1567.230
Tethys           3644.8    20000.0 2551.206
Dione             508.7    20000.0  277.810
Rhea              387.1    20000.0  151.482
Titan             519.9    20000.0   89.192
Iapetus           268.1    20000.0   15.306
Miranda          3053.2    20000.0 4843.570
Ariel            2577.2    20000.0 2782.716
Umbriel          3258.8    20000.0 2536.015
Titania         20782.6    20000.0 9821.256 <-- FAIL
Oberon           4557.5    20000.0 1612.712
Triton           9288.7    20000.0 5400.101
Charon            711.7     2000.0 7490.501
FAIL accuracy_all: 33/34 bodies within tolerance
```

Planets land at ≈0.1–0.3″ over 26.7 years. **Titania (20,783 km)
exceeds its 20,000 km gate and the check reports FAIL.** The
20,000 km moon bar is a first gate on the oblateness upgrade; Titania
sits ≈4% over it — consistent with residual Uranus-gravity error (the
URA182.21 compromise, plus unmodelled higher zonal terms and pole
motion), not an integrator defect. The tolerance is not loosened. The
CORE-mode Earth row (3,134 km) reflects the Earth–Moon-barycentre
point-mass limitation — a point EMB cannot feel the real lunar tidal
coupling; with a real Moon loaded (ALL) Earth is at 116 km. Mimas's
11,578 km is within its gate and reflects how little of Saturn's inner
system structure survives 26.7 years in this model.

Tolerances are fixed in `src/verify/ephem_check.c` before the run:
Sun/planets `max(1000 km, 0.5″ × heliocentric distance)`; CORE-mode
Earth 5,000 km; Moon/Charon/dwarfs/asteroids 2,000 km; other moons
20,000 km.

## Performance

Integration budget is deterministic, not wall-clock: each displayed
frame gets a substep budget computed from live particle counts
(`sv_budget`, calibrated at ~1e7 cost units ≈ 11 ms per frame on the
M5 Pro). If the warp needs more, the HUD shows `LIMITED` and sim time
runs behind — never at reduced accuracy.

Measured substep costs this run (`stspace` perf lines, M5 Pro):
45 massive + 20,000 tracers (tree): 2.81 ms; 2,000 massive direct:
8.13 ms; 20,000 massive tree: 418 ms; 20 massive + 1,000,000 tracers:
33.8 ms (pn off), 61.9 ms (pn on). In the windowed build the full Solar
System scenario runs at ≈60 fps; the million-tracer scenario 5 runs at
≈10 fps — dominated by the tracer force loop and the per-particle
screen projection.

## Verification

`stspace` (`src/verify/space_check.c`, in `make verify`) covers the
core: Kepler closure (1.2e-11 AU over one orbit), observed convergence
order 5.98, bitwise reversibility with pn=0, Mercury 42.98″/cy with pn=1
vs −0.0002″/cy without, merge conservation, black-hole
disrupt/capture events, checkpoint seek bitwise both directions, J2
nodal regression vs the analytic rate (1.3e-3 rel), zonal momentum
conservation, 1PN substep independence (0.27 km over 10 yr at h = 6 h
vs 1 h), thread-count-independence (identical state hash), and the tree
force error above. `stephem` (in `make verify`) covers the catalog:
calendar round-trips, catalog integrity (34 entries), load counts
(9 CORE / 27 CORE+MOONS / 34 ALL), barycentre residual 1.5e-16, the
Earth J2000 sample tie-in, and the belt-tracer envelope. The
`--accuracy` suite above is too slow for `make verify`; it is
`make ephem-accuracy`.

## Claim classes

| what | class | basis |
|------|-------|-------|
| catalog states, GM, radii at both epochs | `EXTERNAL` | JPL Horizons fetch, `data/horizons/` |
| J2/J4, rref, spin poles | `EXTERNAL` | NAIF kernel comments + PCK, `data/naif/` |
| evolved positions, velocities, events, energy/momentum, osculating readouts | `MEASURED` | `stc_space_advance` integration |
| belt tracers, `SPAWN` swarms, rogue bodies | `ILLUSTRATIVE` | seeded sampled populations, not real catalogued objects |
| orbit-ellipse overlay | the drawn curve is the *instantaneous* osculating two-body orbit computed from measured state, not a prediction |
| starfield, grid, labels, lit spheres, BH shadow glyph | `ILLUSTRATIVE` | presentation only (R6: the renderer never feeds back into state) |

## What this is not

- **Not a DE440/441-grade ephemeris.** No J6 or higher zonals, no
  tesseral harmonics, no figure (torque/spin) dynamics — poles are fixed
  vectors. No tidal forces (tidal *disruption* is an event rule, not a
  stress model). No non-gravitational forces: no radiation pressure, no
  Poynting–Robertson, no atmospheric drag, no outgassing; `drag` exists
  in the library but is not wired in.
- **1PN is a test-particle correction**, not a full post-Newtonian
  N-body formulation: reaction on the sources is neglected, there is no
  1PN source-to-source coupling, and the core shows residual phase
  drift consistent with that.
- **CORE mode's "Earth" is the Earth–Moon barycentre as a point mass**;
  it cannot reproduce the lunar tidal coupling, which is why that row
  carries its own tolerance. Load MOONS (scenarios 0–2) for the real
  Earth+Moon pair.
- **Times are TDB.** The TDB−UTC offset (~69 s) is not applied and no
  leap seconds exist in the calendar math (proleptic Gregorian). Dates
  you read are TDB calendar dates.
- **Tracers are one-way coupled** — belts and swarms show structure,
  not self-gravity.
- **Events are instantaneous bookkeeping** (merge at contact, capture
  inside 2 r_s, disruption at the Hill-style tidal radius), not
  hydrodynamics or a Roche-limit simulation.
