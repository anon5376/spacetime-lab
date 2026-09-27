/* SPACETIME LAB — app: termios raw-mode UI + orbit camera implementation.
 * POSIX (termios/unistd/signal) is confined to this file per contract Sec 4.
 */

#define _POSIX_C_SOURCE 200809L

#include "stc/ui.h"
#include "stc/vec.h"

#include <errno.h>
#include <math.h>
#include <signal.h>
#include <stddef.h>
#include <termios.h>
#include <unistd.h>

/* Camera step sizes: 0.08 rad per keypress, x1.1 zoom per +/-. */
#define STC_UI_YAW_STEP 0.08
#define STC_UI_PITCH_STEP 0.08
#define STC_UI_ZOOM_RATIO 1.1
#define STC_UI_PITCH_MAX 1.45
#define STC_UI_DIST_MIN 0.5
#define STC_UI_DIST_MAX 1000000.0

static struct termios g_saved;
static int g_raw = 0;
static int g_have_saved = 0;

/* Equation: signal-safe restore — best-effort tcsetattr, then re-raise. */
static void stc_ui_sig_restore(int sig)
{
    if (g_raw && g_have_saved) {
        tcsetattr(STDIN_FILENO, TCSANOW, &g_saved);
        g_raw = 0;
    }
    signal(sig, SIG_DFL);
    raise(sig);
}

int stc_ui_enter_raw(void)
{
    struct termios raw;

    if (g_raw) {
        return 0;
    }
    if (!isatty(STDIN_FILENO)) {
        return -2;
    }
    if (tcgetattr(STDIN_FILENO, &g_saved) != 0) {
        return -1;
    }
    g_have_saved = 1;
    raw = g_saved;
    raw.c_iflag &= (tcflag_t) ~(BRKINT | ICRNL | INPCK | ISTRIP | IXON);
    raw.c_oflag &= (tcflag_t) ~(OPOST);
    raw.c_cflag |= (tcflag_t)CS8;
    raw.c_lflag &= (tcflag_t) ~(ECHO | ICANON | IEXTEN | ISIG);
    raw.c_cc[VMIN] = 0;
    raw.c_cc[VTIME] = 0;
    if (tcsetattr(STDIN_FILENO, TCSANOW, &raw) != 0) {
        return -1;
    }
    g_raw = 1;
    signal(SIGINT, stc_ui_sig_restore);
    signal(SIGTERM, stc_ui_sig_restore);
    return 0;
}

void stc_ui_exit_raw(void)
{
    if (g_raw && g_have_saved) {
        tcsetattr(STDIN_FILENO, TCSANOW, &g_saved);
    }
    g_raw = 0;
    signal(SIGINT, SIG_DFL);
    signal(SIGTERM, SIG_DFL);
}

int stc_ui_is_raw(void)
{
    return g_raw;
}

int stc_ui_poll_key(void)
{
    unsigned char c;
    unsigned char seq0;
    unsigned char seq1;
    ssize_t n;

    if (!g_raw) {
        return STC_UI_KEY_NONE;
    }
    n = read(STDIN_FILENO, &c, 1);
    if (n == 0) {
        return STC_UI_KEY_NONE;
    }
    if (n < 0) {
        if (errno == EAGAIN || errno == EINTR) {
            return STC_UI_KEY_NONE;
        }
        return STC_UI_KEY_ERROR;
    }
    if (c != 0x1b) {
        return (int)c;
    }
    /* ESC: decode CSI "ESC [ X" and SS3 "ESC O X" arrows; lone ESC -> 27. */
    n = read(STDIN_FILENO, &seq0, 1);
    if (n <= 0) {
        return STC_UI_KEY_ESC;
    }
    if (seq0 != '[' && seq0 != 'O') {
        return STC_UI_KEY_ESC;
    }
    n = read(STDIN_FILENO, &seq1, 1);
    if (n <= 0) {
        return STC_UI_KEY_ESC;
    }
    switch (seq1) {
    case 'A':
        return STC_UI_KEY_UP;
    case 'B':
        return STC_UI_KEY_DOWN;
    case 'C':
        return STC_UI_KEY_RIGHT;
    case 'D':
        return STC_UI_KEY_LEFT;
    default:
        return STC_UI_KEY_ESC;
    }
}

void stc_ui_camera_init(StcUiCamera *cam)
{
    if (cam == NULL) {
        return;
    }
    cam->target.x = 0.0;
    cam->target.y = 0.0;
    cam->target.z = 0.0;
    cam->yaw_rad = 0.7;
    cam->pitch_rad = 0.35;
    cam->dist_m = 10.0;
    stc_ui_camera_update(cam);
}

void stc_ui_camera_update(StcUiCamera *cam)
{
    double cp;

    if (cam == NULL) {
        return;
    }
    cam->pitch_rad = stc_clamp(cam->pitch_rad, -STC_UI_PITCH_MAX, STC_UI_PITCH_MAX);
    cam->dist_m = stc_clamp(cam->dist_m, STC_UI_DIST_MIN, STC_UI_DIST_MAX);
    cp = cos(cam->pitch_rad);
    cam->eye.x = cam->target.x + cam->dist_m * cp * cos(cam->yaw_rad);
    cam->eye.y = cam->target.y + cam->dist_m * sin(cam->pitch_rad);
    cam->eye.z = cam->target.z + cam->dist_m * cp * sin(cam->yaw_rad);
    cam->up.x = 0.0;
    cam->up.y = 1.0;
    cam->up.z = 0.0;
}

int stc_ui_camera_key(StcUiCamera *cam, int key)
{
    if (cam == NULL) {
        return 0;
    }
    switch (key) {
    case 'W':
    case 'w':
    case STC_UI_KEY_UP:
        cam->pitch_rad += STC_UI_PITCH_STEP;
        break;
    case 'S':
    case 's':
    case STC_UI_KEY_DOWN:
        cam->pitch_rad -= STC_UI_PITCH_STEP;
        break;
    case 'A':
    case 'a':
    case STC_UI_KEY_LEFT:
        cam->yaw_rad -= STC_UI_YAW_STEP;
        break;
    case 'D':
    case 'd':
    case STC_UI_KEY_RIGHT:
        cam->yaw_rad += STC_UI_YAW_STEP;
        break;
    case '+':
    case '=':
    case 'R':
    case 'r':
        cam->dist_m /= STC_UI_ZOOM_RATIO;
        break;
    case '-':
    case '_':
    case 'F':
    case 'f':
        cam->dist_m *= STC_UI_ZOOM_RATIO;
        break;
    default:
        return 0;
    }
    stc_ui_camera_update(cam);
    return 1;
}

const char *stc_ui_key_name(int key)
{
    switch (key) {
    case STC_UI_KEY_NONE:
        return "none";
    case STC_UI_KEY_ERROR:
        return "error";
    case STC_UI_KEY_CTRL_C:
        return "ctrl-c";
    case STC_UI_KEY_ESC:
        return "esc";
    case STC_UI_KEY_UP:
        return "up";
    case STC_UI_KEY_DOWN:
        return "down";
    case STC_UI_KEY_LEFT:
        return "left";
    case STC_UI_KEY_RIGHT:
        return "right";
    default:
        break;
    }
    if (key >= 32 && key < 127) {
        return "printable";
    }
    return "unknown";
}
