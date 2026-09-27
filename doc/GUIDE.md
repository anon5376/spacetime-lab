# Sandbox user guide

Launch the sandbox with `./build/stlab --sandbox --preset 0`, or
double-click `stlab.command`. You get a 3D gravity universe on a
character grid: a live viewport with orbit trails, a body list, and full
time controls. Press N to spawn any body you can describe in SI units.

## Presets

Four seeded starting points, selected with `--preset N` or cycled live
with G (cycling resets the simulation):

| N | Name | Bodies | Contents |
|---|---|---|---|
| 0 | solar | 9 | Sun plus 8 planets, real masses and distances |
| 1 | binary | 3 | Twin stars with one circumbinary planet |
| 2 | hole_disk | 26 | Black hole with a 25-particle accretion disk |
| 3 | slingshot | 2 | Planet plus an inbound probe set up for a flyby |

`--seed N` reseeds the randomized parts (default 20260927). Same
arguments always produce byte-identical output.

## Keys

| Key | Action |
|---|---|
| N | Open the spawn form |
| K | Remove the followed body |
| Tab | Cycle follow target (body list marker `>`) |
| G | Next preset (resets the sim) |
| Z / X | Halve / double the timestep, clamped to 1e-6..1e9 s |
| C | Cycle warp x1 / x10 / x100 (steps per frame) |
| T | Toggle orbit trails |
| D | Toggle spacetime-grid sheet + orientation aids (visual mode only) |
| P | Pause |
| E | Export CSV + PPM (to `--csv`/`--ppm` paths, else `sandbox.csv`/`sandbox.ppm`) |
| WASD / arrows, + / - | Orbit camera, zoom |
| ? | Toggle the key guide (any key returns) |
| Q or Esc | Quit |

The header shows sim time, timestep, warp, body count, and total
energy. Each body row shows mass, radius, speed, and `taurate` -- the
instantaneous proper-time rate d(tau)/dt relative to body 0 (1.0 means
no measurable dilation); it is a rate factor, not accumulated proper
time. FIX marks anchored bodies. The follow line shows the followed body's Kepler
elements about body 0: semi-major axis a, eccentricity e, and period T
(or "unbound").

## Spawning bodies

The form has 11 fields: name, mass_kg, radius_m, three positions in
meters, three velocities in m/s, fixed, and a SPAWN> commit line. All
numbers accept scientific notation such as `1.989e30`. Mass must be
positive, radius zero or more. The fixed field takes `y` or `1` for an
immovable anchor, `n`, `0`, or empty for a free body. An empty name
becomes "body". Up and Down move between fields, Enter on SPAWN>
commits, Esc cancels.

The form opens prefilled with a probe on a circular orbit about body 0
through the middle of the view, so pressing Enter immediately gives you
a working satellite. At most 64 bodies fit; the 65th spawn is refused.

Bodies that touch merge into one, conserving momentum; the lower slot
number survives and absorbs the other's mass. (Elastic bounce exists in
the core API; the sandbox always merges.) Spawned bodies join the live
simulation at once, and the follow target jumps to the newcomer.

## Time controls

The default step is 1000 s. Earth's orbit takes ~31558 such steps, so use
warp x100 (C) to watch a year pass in about 316 frames. Shorten dt with
Z when bodies pass close together: a step should stay well under a
tenth of the fastest orbit, or orbits smear and energy drifts. Lengthen
it with X for deep-time coasting. Pause with P before delicate spawns;
the sim holds still while the form stays live.

## Spawn recipes

All recipes use SI units and assume preset 0 (Sun, 1.98847e30 kg, at the
origin) unless stated. Circular-orbit speed at distance r from mass M
is v = sqrt(GM/r) with G = 6.6743e-11.

Earth twin. A second Earth sharing our orbit, placed opposite the Sun:
position (1.496e11, 0, 0), velocity (0, 29780, 0), mass 5.972e24,
radius 6.371e6. It laps the Sun in a year. Offset the start angle to
watch conjunctions.

Retrograde comet. Position (5e11, 0, 0), velocity (0, -8000, 0), mass
1e14, radius 5e3. Negative vy sends it against planetary traffic; the
Sun bends it into a steep ellipse. Raise dt and warp to follow the
full loop.

Sundiver. Position (1.496e11, 0, 0), velocity (0, 0, 0), mass 1e10.
With no sideways motion it falls straight in. Lower dt to 100 s (Z a
few times) before impact so the final plunge resolves instead of
tunneling through. `--integrator rk4` swaps the Verlet stepper for
classic 4th-order Runge-Kutta; collisions merge and the speed cap
applies exactly as under Verlet -- the same post-step pass runs after
every step under both integrators, so only the stepping scheme
differs. The header names the active stepper.

Escape probe. Circular speed at 1 AU is 29.8 km/s; anything past
sqrt(2) times that, about 42.1 km/s, leaves forever. Spawn at
(1.496e11, 0, 0) with (0, 45000, 0) and watch the follow line flip to
"unbound" as it exits.

Binary suns. Start from an empty-feeling scene (preset 3 has only 2
bodies, or remove down to none and build up): two 1.989e30 kg stars at
(-1e11, 0, 0) and (1e11, 0, 0) with velocities (0, -1.82e4, 0) and (0,
1.82e4, 0). Each circles the midpoint with a period near 68 days. Then
spawn a circumbinary planet far out at (1e12, 0, 0) with (0, 13000, 0)
and check it survives.

Fixed anchor plus moon. Spawn a fixed Earth: mass 5.972e24 at the
origin, fixed `y`. A low orbiter at (6.771e6, 0, 0) with (0, 7660, 0)
circles in 92 minutes; dt 10 s with warp x10 shows it smoothly. A
geostationary satellite sits at 4.2164e7 m with 3075 m/s and hangs
over one spot (there is no rotating frame, so "stationary" means
matching a 24 h period).

Slingshot tuning (preset 3). The preset probe is already inbound. Pause
at closest approach, note its speed gain on the way out, then K to
remove it, tweak, and respawn with slightly different aim to see how
flyby distance trades against deflection.

## Headless runs

`--smoke` runs without interaction and exits 0 on success, writing a
time-series CSV (frame, sim time, total energy, body count) and one PPM
frame of the final viewport:

```sh
./build/stlab --sandbox --preset 0 --frames 60 \
  --csv /tmp/run.csv --ppm /tmp/run.ppm --smoke --no-adapt
```

`make sandbox-smoke` runs all four presets into `/tmp/sbx`. Repeating
any command reproduces its files byte for byte.

## Visual mode

`--sandbox --visual` runs the same sandbox in a native window instead
of the character grid -- macOS only (AppKit driven through the
Objective-C runtime, system frameworks, no third-party dependencies).
On other platforms it prints `stlab: --visual needs macOS (AppKit)`
and exits 1. The window shows a software-rasterized 1280x800
framebuffer: shaded body discs, fading trails, a starfield, and the
same HUD body list (the `taurate` column) plus the same keys and spawn
form. `--frames N` caps the session; `--smoke` is ignored.

`--sandbox --visual-shot out.ppm` renders the visual view headlessly:
it steps the sim `--frames N` times (default 60), then writes the
final framebuffer as a 1280x800 P6 PPM. No window, no pacing, works on
any platform, and identical arguments produce a byte-identical file.
If both visual flags are given, `--visual-shot` wins. Without
`--sandbox`, both visual flags are ignored (they are sandbox modes).

### Spacetime grid (visual mode)

D toggles a wireframe lattice lying in the orbital plane -- the
classic rubber-sheet picture. Each node is displaced by a
potential-shaped dent under the live bodies,
`h(x,z) = -K * sum_i w_i * s^2 / (d_i^2 + s^2)`, with the kernel width
`s` fixed in view units (a couple of lattice spacings, so every funnel
is a few nodes wide at any zoom) and compressed mass weights
`w_i = (m_i / M_max)^0.25` (the heaviest body digs the deep funnel,
planets still leave dimples). K is set so the deepest lattice node
dips by 18% of the scene half-span -- an arbitrary view-scale depth.
The lattice extends toward the horizon: a uniform 33x33 core around
the fitted scene continues with geometrically spaced lines -- the
node coordinate doubles each step out to ~16384x the view half-span,
far enough that the finite edge compresses into the horizon line and
the far lattice packs to sub-pixel density -- so the sheet reads as
unbounded and flattens to y=0 as the dent tails to zero. This is an ILLUSTRATIVE
exaggeration for orientation only: it draws
the Newtonian potential as a sheet, it is not the spacetime metric, GR
plays no part in it, and nothing about it feeds back into the
physics. The same key also shows a small x/y/z axis triad (top right)
and a 100 px scale bar (bottom right, auto-units m/km/AU). Note that
in visual mode D takes over the WASD yaw-right binding; the Right
arrow still orbits the camera.

## Limits

Newtonian gravity only; relativity does not bend the trajectories
(the `taurate` column is an instantaneous d(tau)/dt readout, not
accumulated proper time; gamma is a compiled-in helper). Touching
bodies merge; the bounce, fragmentation, tide, drag, and
radiation-pressure modules are API helpers, not live forces in the
sandbox stepper. Details and equations: `doc/PHYSICS.md`.
