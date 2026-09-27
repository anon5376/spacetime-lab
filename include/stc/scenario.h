#ifndef STC_SCENARIO_H
#define STC_SCENARIO_H

/* SPACETIME LAB — scenario presets (builder 9).
 * Deterministic demo table: 10 presets, fixed seeds, fixed frame counts.
 * Plain-data context (no rendering, no dynamic allocation). SI doubles.
 * Depends only on stc/core (vec/units/rng) + C stdlib. No sibling includes.
 *
 * Canonical names (CONTRACT Sec 5, CLI/CSV stable):
 *   twin pole doppler mercury gps bend ring flamm chirp shadow
 * Accepted aliases (same physics, alternate lookup names):
 *   minkowski->doppler ladder->pole photon->bend lensing->ring
 *   waves->chirp kerr->shadow
 *
 * Observable slots obs[8] (SI unless noted), set every step:
 *   twin:    [0]=earth time s [1]=ship proper time s [2]=gamma [3]=speed m/s
 *   pole:    [0]=contracted length m [1]=rest length m [2]=gamma [3]=barn m
 *   doppler: [0]=redshift z [1]=f_obs Hz [2]=f_emit Hz [3]=beta
 *            [4]=pulses emitted [5]=pulses received
 *   mercury: [0]=accumulated precession arcsec [1]=orbits [2]=a m [3]=e
 *   gps:     [0]=SR offset us/day [1]=GR offset us/day [2]=net us/day
 *            [3]=accumulated offset ns
 *   bend:    [0]=deflection arcsec [1]=impact parameter m
 *            [2]=photon x m [3]=photon y m
 *   ring:    [0]=theta_E arcsec [1]=theta_plus arcsec [2]=theta_minus arcsec
 *            [3]=beta arcsec [4]=seeded beta offset arcsec
 *   flamm:   [0]=rs m [1]=w(r) m [2]=r m [3]=w_max m [4..7]=azimuths rad
 *   chirp:   [0]=f_gw Hz [1]=strain h [2]=h_plus [3]=h_cross [4]=Mc kg
 *   shadow:  [0]=spin chi a/M [1]=prograde ISCO radius m [2]=r_plus m
 *            [3]=shadow-radius estimate m
 */

#include <stddef.h>
#include <stdint.h>

#include "stc/rng.h"
#include "stc/vec.h"

#define STC_SCENARIO_COUNT 10
#define STC_SCENARIO_NPOS 4
#define STC_SCENARIO_NOBS 8
#define STC_SCENARIO_NAME_LEN 32

typedef struct {
    const char *name;     /* canonical preset name (stable) */
    const char *alias;    /* alternate lookup name, NULL if none */
    const char *desc;     /* one-line description */
    const char *equation; /* defining equation name */
    uint64_t seed;        /* fixed seed for this preset */
    int frames;           /* fixed frame count */
    double dt;            /* seconds per frame */
} StcScenarioInfo;

typedef struct {
    char name[STC_SCENARIO_NAME_LEN]; /* canonical preset name */
    uint64_t seed;                   /* rng seed in force */
    int index;                       /* table index, -1 when unset */
    int frame;                       /* frames completed, 0..frames */
    int frames;                      /* fixed frame count */
    double dt;                       /* seconds per frame */
    double t;                        /* elapsed sim time, s */
    double phase;                    /* oscillator/orbital phase, rad */
    StcVec3 pos[STC_SCENARIO_NPOS];  /* body/feature positions, m */
    StcVec3 vel[STC_SCENARIO_NPOS];  /* body/feature velocities, m/s */
    double obs[STC_SCENARIO_NOBS];   /* observables, SI (see table above) */
    StcRng64 rng;                    /* deterministic stream, seeded at setup */
    int done;                        /* 1 when frame == frames */
} StcScenarioCtx;

/* Preset-table size identity: count = 10 always. Returns STC_SCENARIO_COUNT. */
size_t stc_scenario_count(void);

/* Preset-table accessor: list[i], i in [0, count). Returns static table. */
const StcScenarioInfo *stc_scenario_list(void);

/* Name lookup (exact match on name or alias): info for `name`, NULL if
 * unknown or if `name` is NULL. */
const StcScenarioInfo *stc_scenario_find(const char *name);

/* Setup with the preset fixed seed: ctx = initial state(scenario, seed0).
 * Returns 0 on success, -1 on NULL ctx or unknown name. */
int stc_scenario_setup(StcScenarioCtx *ctx, const char *name);

/* Setup with explicit seed: ctx = initial state(scenario, seed).
 * Returns 0 on success, -1 on NULL ctx or unknown name. */
int stc_scenario_setup_with_seed(StcScenarioCtx *ctx, const char *name,
                                 uint64_t seed);

/* Single-frame advance: state(f+1) = step(state(f)). Returns 0 when more
 * frames remain, 1 when this step completed the run (or run was already
 * complete; then a no-op), -1 on NULL ctx or unset ctx. */
int stc_scenario_step(StcScenarioCtx *ctx);

/* Multi-frame advance: run up to max_frames steps (negative = to
 * completion). Returns frames advanced, -1 on NULL ctx or unset ctx. */
int stc_scenario_run(StcScenarioCtx *ctx, int max_frames);

#endif /* STC_SCENARIO_H */
