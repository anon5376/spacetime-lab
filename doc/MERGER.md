# MERGER.md — the binary black-hole merger view (`stlab --merger`)

`--merger` is a cinematic edge-on view of a binary black-hole coalescence:
two black holes inspiral on a gravitational-wave-driven decay, plunge,
merge into a common horizon, and ring down, while the surrounding
spacetime "sheet" ripples and the background lensing bends around the
pair. Time, radii, chirp frequency, remnant mass and spin, and the QNM
parameters on the HUD are real numbers evaluated from the shipped model;
the drawn sheet depth, ripple, lensing, and shadows are illustrative
visual choices (see "Claim classes").

Run it:

```sh
make merger                          # windowed (macOS)
./build/stlab --merger --preset 0
make merger-shot                     # headless still -> build/merger.png
./build/stlab --merger --preset 0 --frames 2400 --movie /tmp/f \
    --movie-every 2                  # PPM frames for tools/merger_movie.sh
tools/merger_movie.sh 0 2400 2       # frames + MP4 + GIF via ffmpeg
```

`--preset N` selects the binary (0–3). The windowed path
(`stc_mvis_run`) needs macOS; elsewhere use `--visual-shot` or
`--movie`. The renderer is deterministic and thread-count independent.

| preset | name (`kPre[]`, shown uppercased in the HUD) | q | a0 | inspiral ends (t_m) | M_tot |
|--------|--------------------------------------------|---|-----|---------|-------|
| 0 | `GW150914-like  36 + 29 Msun` | 1.24 | 11.0 M | 862 M | 65 M_sun |
| 1 | `equal mass  30 + 30 Msun` | 1.00 | 11.0 M | 854 M | 60 M_sun |
| 2 | `unequal  q = 3  30 + 10 Msun` | 3.00 | 10.0 M | 689 M | 40 M_sun |
| 3 | `extreme  q = 8  24 + 3 Msun` | 8.00 | 8.5 M | 513 M | 27 M_sun |

The presets still pass `a_merge = 3.4 M` to `stc_bbh_init`; the argument
is validated but superseded — the inspiral now always ends at the ISCO,
a = 6 M (see below).

## Keys

Windowed mode (`stc_mvis_key`, plus the `stc_mvis_run` mouse handlers):

| input | effect |
|-------|--------|
| `1`–`4` | load preset 0–3 |
| `Space` / `P` | pause / resume |
| `R` | restart the timeline at t = 0 |
| `+` `=` / `-` `_` | playback speed ×2 / ÷2 (0.0625–16) |
| `A` `D` / left/right arrows | orbit yaw |
| `S` / down arrow, up arrow | pitch down / up (`W` is taken by wave gain) |
| `G` | toggle the spacetime sheet |
| `W` | cycle gravitational-wave ripple gain (×0.5, ×1, ×2, ×4) |
| `L` | toggle gravitational lensing |
| `T` | toggle orbital trails |
| `C` | toggle the guided camera (auto-dolly toward merger) |
| `H` | toggle the HUD |
| left-drag | orbit camera |
| scroll / pinch | zoom |
| `Q`, `Esc`, `Ctrl-C` | quit |

The HUD shows the stage badge (`INSPIRAL`, `PLUNGE / COMMON HORIZON`,
`RINGDOWN`, then `REMNANT KERR BLACK HOLE` once
t − t_peak > 6 τ_qnm), component masses, separation in M and km, GW
frequency in Hz, time-to-peak, and (after merger) remnant mass, spin,
radiated energy, and the QNM220 frequency and damping time — plus a
live chirp trace of the face-on strain r h_+ / peak vs T − T_PEAK.

## The model (as implemented in `src/merger/merger.c`)

Geometric units G = c = 1, total initial mass M = 1; all lengths and
times are in M. `stc_bbh_init(q, a0, a_merge)` integrates the inspiral
once into a heap table (`stc_bbh_free` releases it); `stc_bbh_eval(t)`
evaluates the timeline deterministically.

- **Inspiral — `ANALYTIC` (3.5PN).** TaylorT4 for non-spinning circular
  binaries: `dv/dt` from Buonanno, Iyer, Ochsner, Pan & Sathyaprakash
  2009 (arXiv:0907.0700, Eq. 6), with `v = (Mω)^{1/3}` and `dφ/dt = v³`,
  integrated at init by fixed-step RK4 from `v0 = a0^{−1/2}`
  (leading-order Kepler) forward to the ISCO `v = 6^{−1/2}` and backward
  to t = −400 M for retarded-time lookups; evaluated by cubic Hermite
  interpolation (node derivatives are exact). Separation for display
  `a = v^{−2}` and strain scale `r h = 4η v²` are leading-order.
- **Plunge — `ILLUSTRATIVE`.** From the ISCO (a = 6 M) to coalescence over
  `STC_BBH_PLUNGE_T = 14 M` (u = 0..1): separation falls parabolically,
  the orbital frequency ramps linearly ω_ISCO → ω_qnm/2 (GW frequency
  2ω_ISCO → ω_qnm), the strain amplitude smoothsteps up 40% to the peak,
  and a smoothstep `fuse` blend sweeps the mass budget 1 → m_f. Phase
  and frequency are continuous at both seams. A shaped transition, not a
  solve: post-ISCO merger physics is outside post-Newtonian theory.
- **Ringdown — `ILLUSTRATIVE`.** A single l = m = 2 damped sinusoid:
  `amp ∝ exp(−(t − t_peak)/τ_qnm)`, GW frequency ω_qnm.
- **Remnant and QNM — `ILLUSTRATIVE`.** Phenomenological polynomial fits
  in η, as coded: `e_rad = η(1 − √(8/9)) + 0.498 η²`, `m_f = 1 − e_rad`;
  `a_f = 2√3 η − 3.5171 η² + 2.5763 η³`; `ω_qnm = (1.5251 −
  1.1568(1−a_f)^0.1292)/m_f`, `τ_qnm = 2(0.7 + 1.4187(1−a_f)^−0.4990)/ω_qnm`.
  At q = 1: m_f ≈ 0.9546, a_f ≈ 0.6865, inside the numerical-relativity
  ranges the check gates on (0.94–0.96, 0.66–0.71). No source is cited
  for the coefficient sets; treat them as illustrative.
- **Strain accessors.** `stc_bbh_rh_plus(b, t_ret, az)` returns the
  edge-on plus-polarization `r h_+ = amp cos(2φ − 2 az)/2`;
  `stc_bbh_rh_faceon(b, t_ret)` returns `amp cos 2φ`. Both take the
  already-retarded time. `rh_faceon` backs the HUD strain plot; the
  sheet's ripple table is filled from `stc_bbh_eval`.

### Checks (`stverify`, 12/12)

- `bbh_timeline`: dφ/dt vs ω at five interior times (relative error
  5.3e-11), a(t_m) = 6 M, phase and frequency continuity at both seams,
  remnant ranges, bitwise determinism.
- `bbh_t2_crosscheck`: for q = 1.24 (65 M_sun), time from f_GW = 35 Hz to
  the ISCO from the production TaylorT4 table vs the independent TaylorT2
  closed form (arXiv:0907.0700, Eq. 8b) evaluated in the harness:
  0.1193 s vs 0.1134 s, 5.2 % apart (tolerance 10 %; the two differ
  beyond 3.5PN order).
- `bbh_gw150914_context` (information only): 35 Hz → ISCO 0.119 s,
  35 Hz → strain peak 0.124 s, 6.1 GW cycles. LIGO's GW150914 paper
  (arXiv:1602.03837) describes ~0.2 s and ~8 cycles from 35 to 150 Hz;
  that interval includes the post-ISCO merger, which PN cannot model,
  so it is context rather than a pass/fail target. The previous 0PN
  model gave 0.191 s / ~10.5 cycles — closer to that number by
  coincidence, not by better physics.

## Rendering and retarded time

The renderer (`src/vis/vis_merger.c`) keeps a ring of GW strain samples
`r h_+` vs retarded radius (0.25 M spacing, gain adjustable with `W`).
Each sheet vertex and background point is displaced by the strain that
was emitted `t − r` — a real retarded-time bookkeeping, applied to an
illustrative displacement. Rendering is CPU-only and scanline-parallel
on a private pthread row pool inside `vis_merger.c` (MV_MAXT = 16
workers max; every pixel is computed independently, so the image is
identical for any thread count — it does not use `stc_par_for`). Each
`stc_mvis_step` advances the timeline by `mv_rate · speed` in M per
frame, where `mv_rate = 1.5(0.24 + 0.76·smooth(12, 170, |t − t_peak|))`
is the slow-mo pacing — deterministic per frame index, never driven by
wall-clock time. The 60 fps pacing lives only in the windowed event
loop, so frame 2400 of a `--movie` sequence is the same state as frame
2400 in the window.

## Claim classes

| what | class | basis |
|------|-------|-------|
| orbital phase/frequency during inspiral | `ANALYTIC` | TaylorT4 3.5PN (arXiv:0907.0700 Eq. 6), checked by `bbh_timeline` and `bbh_t2_crosscheck` |
| separation a = v⁻², strain scale r h = 4ηv² and HUD readouts derived from them | `ANALYTIC` | leading-order expressions |
| remnant mass, spin, QNM ω/τ | `ILLUSTRATIVE` | phenomenological fits, uncited, gated to NR ranges |
| plunge trajectory/amplitude, ringdown waveform | `ILLUSTRATIVE` | shaped blend + single damped sinusoid |
| sheet depth (softened Newtonian dip, −30 units, ε = 2.6, half-extent 170 M) | `ILLUSTRATIVE` | visual proxy for curvature, not a metric |
| ripple amplitude/falloff (75× base gain, near-zone taper, (r+12)^{−0.45} envelope) | `ILLUSTRATIVE` | exaggerated so the GW is visible |
| black-hole shadows (screen-space metaballs ~ photon-sphere sized) and the peanut common-horizon at merge | `ILLUSTRATIVE` | sized to read on screen, not a raytraced Kerr shadow |
| thin-lens background deflection (Einstein radius² = 2.1 × shadow²) | `ILLUSTRATIVE` | plausible qualitative bending, not integrated geodesics |
| photon-ring glow, bloom, ACES tone map, vignette, generated starfield sky, playback pacing (slow-mo near merger), camera dolly | `ILLUSTRATIVE` | cinematic treatment |
