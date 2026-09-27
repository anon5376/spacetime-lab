#ifndef STC_UI_H
#define STC_UI_H

/* SPACETIME LAB — app: termios raw-mode UI + orbit camera.
 * This header is pure C11 (no termios here; POSIX is confined to
 * src/app/ui.c). SI doubles: meters for distance, radians for angles.
 * Camera model: eye orbits `target` at `dist_m` with yaw/pitch.
 */

#include "stc/vec.h"

/* Key codes from stc_ui_poll_key. Plain ASCII passes through unchanged
 * (e.g. 'W', 'q', '1', 27 for lone ESC); special keys use codes >= 256. */
#define STC_UI_KEY_NONE 0
#define STC_UI_KEY_ERROR (-1)
#define STC_UI_KEY_CTRL_C 3
#define STC_UI_KEY_ESC 27
#define STC_UI_KEY_UP 256
#define STC_UI_KEY_DOWN 257
#define STC_UI_KEY_LEFT 258
#define STC_UI_KEY_RIGHT 259

/* Orbit camera: `eye`/`up` are derived from yaw/pitch/dist by update. */
typedef struct {
    StcVec3 target;   /* look-at center, meters */
    double yaw_rad;   /* azimuth around target, radians */
    double pitch_rad; /* elevation above horizon, radians */
    double dist_m;    /* eye distance from target, meters */
    StcVec3 eye;      /* derived eye position, meters */
    StcVec3 up;       /* derived up vector (unit) */
} StcUiCamera;

/* Equation: raw-mode enter — tcgetattr snapshot + cfmakeraw-equivalent flags.
 * Installs SIGINT/SIGTERM handlers that restore the snapshot, then re-raise.
 * Returns 0 on success (or if already raw), -1 on tcgetattr/tcsetattr failure,
 * -2 if stdin is not a TTY. Safe to call twice (second call is a no-op). */
int stc_ui_enter_raw(void);

/* Equation: raw-mode exit — tcsetattr restore + SIG_DFL handlers.
 * Idempotent and safe to call without a matching enter (also via atexit). */
void stc_ui_exit_raw(void);

/* Equation: raw-state query — nonzero iff raw mode is currently active. */
int stc_ui_is_raw(void);

/* Equation: key poll — nonblocking read + ANSI CSI/SS3 arrow decode.
 * Returns an ASCII byte, STC_UI_KEY_UP/DOWN/LEFT/RIGHT, STC_UI_KEY_NONE
 * when no key is waiting (or raw mode is off), STC_UI_KEY_ERROR on read
 * failure. Lone ESC returns STC_UI_KEY_ESC; Ctrl-C arrives as STC_UI_KEY_CTRL_C. */
int stc_ui_poll_key(void);

/* Equation: camera init — default target/yaw/pitch/dist, then update.
 * NULL-safe no-op. */
void stc_ui_camera_init(StcUiCamera *cam);

/* Equation: orbit camera — eye = target + dist*(cos p*cos y, sin p, cos p*sin y).
 * Clamps pitch to [-1.45, 1.45] rad and dist to [0.5, 1e6] m; up = +Y.
 * NULL-safe no-op. */
void stc_ui_camera_update(StcUiCamera *cam);

/* Equation: camera key step — WASD/arrows adjust yaw/pitch, +/- zoom dist.
 * Accepts upper- and lowercase WASD; arrow codes match poll_key.
 * Returns 1 if the key moved the camera, 0 otherwise (incl. NULL cam). */
int stc_ui_camera_key(StcUiCamera *cam, int key);

/* Equation: key name — static display string for a poll_key code.
 * Never returns NULL ("printable" for ASCII 32..126, "unknown" otherwise). */
const char *stc_ui_key_name(int key);

#endif /* STC_UI_H */
