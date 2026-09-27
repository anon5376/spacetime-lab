/* SPACETIME LAB — vis_mac: native macOS window for the visual sandbox.
 * Pure C over the Objective-C runtime plus AppKit/Foundation (system
 * frameworks only, no third-party dependencies, no ObjC source).
 * Own event loop at 60fps: poll keys, step the sim, repaint the
 * software framebuffer, blit through an NSBitmapImageRep.
 */

/* Must precede every header: exposes clock_gettime/CLOCK_MONOTONIC under
 * -std=c11 on Linux (glibc gates them on _POSIX_C_SOURCE >= 199309L when
 * __STRICT_ANSI__ is set); harmless on macOS, which does not gate them. */
#define _POSIX_C_SOURCE 199309L

#include "stc/vis.h"
#include "stc/vis_merger.h"
#include "stc/vis_space.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifdef __APPLE__
#include <objc/message.h>
#include <objc/runtime.h>

/* Minimal Foundation/AppKit surface (no framework headers needed).
 * CGFloat is double on 64-bit; NSInteger long, NSUInteger unsigned long.
 * BOOL arrives as int on every shim: narrow stack-passed integer args
 * leave the upper bits undefined on ARM64, which once flipped isPlanar
 * on and banded the display into R/G/B thirds. int is exact-safe for a
 * callee reading BOOL or int. */
typedef double CGFloat_;
typedef struct {
    CGFloat_ x;
    CGFloat_ y;
} NSPoint_;
typedef struct {
    CGFloat_ w;
    CGFloat_ h;
} NSSize_;
typedef struct {
    NSPoint_ o;
    NSSize_ s;
} NSRect_;

extern void *NSDefaultRunLoopMode;
extern void *NSDeviceRGBColorSpace;
extern void *objc_autoreleasePoolPush(void);
extern void objc_autoreleasePoolPop(void *);

static NSRect_ NSMakeRect_(double x, double y, double w, double h)
{
    NSRect_ r;

    r.o.x = x;
    r.o.y = y;
    r.s.w = w;
    r.s.h = h;
    return r;
}

/* Typed objc_msgSend shims. Every call goes through an exact prototype,
 * which keeps argument passing ABI-correct on arm64. No method that
 * returns a struct is ever called (drawRect receives its rect; nothing
 * reads back frames, bounds, or points). */
static id vm_id(id o, SEL s)
{
    id (*f)(id, SEL) = (id(*)(id, SEL))objc_msgSend;

    return f(o, s);
}

static void vm_v_id(id o, SEL s, id a)
{
    void (*f)(id, SEL, id) = (void (*)(id, SEL, id))objc_msgSend;

    f(o, s, a);
}

static void vm_v_bool(id o, SEL s, int b)
{
    void (*f)(id, SEL, int) = (void (*)(id, SEL, int))objc_msgSend;

    f(o, s, b);
}

static void vm_v_long(id o, SEL s, long v)
{
    void (*f)(id, SEL, long) = (void (*)(id, SEL, long))objc_msgSend;

    f(o, s, v);
}

static long vm_long(id o, SEL s)
{
    long (*f)(id, SEL) = (long (*)(id, SEL))objc_msgSend;

    return f(o, s);
}

static signed char vm_bool(id o, SEL s)
{
    signed char (*f)(id, SEL) =
        (signed char (*)(id, SEL))objc_msgSend;

    return f(o, s);
}

static double vm_dbl(id o, SEL s)
{
    double (*f)(id, SEL) = (double (*)(id, SEL))objc_msgSend;

    return f(o, s);
}

static const char *vm_cstr(id o, SEL s)
{
    const char *(*f)(id, SEL) =
        (const char *(*)(id, SEL))objc_msgSend;

    return f(o, s);
}

/* locationInWindow returns NSPoint (two CGFloat/double): on both arm64
 * and x86_64 that's a 16-byte all-floating-point aggregate returned in
 * registers, so a plain typed objc_msgSend cast (no objc_msgSend_stret)
 * is correct here. */
static NSPoint_ vm_point(id o, SEL s)
{
    NSPoint_ (*f)(id, SEL) = (NSPoint_ (*)(id, SEL))objc_msgSend;

    return f(o, s);
}

static id vm_id_rect(id o, SEL s, NSRect_ r)
{
    id (*f)(id, SEL, NSRect_) = (id(*)(id, SEL, NSRect_))objc_msgSend;

    return f(o, s, r);
}

static void vm_v_rect(id o, SEL s, NSRect_ r)
{
    void (*f)(id, SEL, NSRect_) =
        (void (*)(id, SEL, NSRect_))objc_msgSend;

    f(o, s, r);
}

static id vm_window(id o, SEL s, NSRect_ r, unsigned long style,
                    unsigned long backing, int defer)
{
    id (*f)(id, SEL, NSRect_, unsigned long, unsigned long, int) =
        (id(*)(id, SEL, NSRect_, unsigned long, unsigned long,
               int))objc_msgSend;

    return f(o, s, r, style, backing, defer);
}

static id vm_next_event(id o, SEL s, unsigned long mask, id date, id mode,
                        int dequeue)
{
    id (*f)(id, SEL, unsigned long, id, id, int) =
        (id(*)(id, SEL, unsigned long, id, id, int))objc_msgSend;

    return f(o, s, mask, date, mode, dequeue);
}

static id vm_bitmap(id o, SEL s, unsigned char **planes, long w, long h,
                    long bps, long spp, int alpha, int opaque,
                    id cs, long bpr, long bpp)
{
    id (*f)(id, SEL, unsigned char **, long, long, long, long, int,
            int, id, long, long) =
        (id(*)(id, SEL, unsigned char **, long, long, long, long, int,
               int, id, long, long))objc_msgSend;

    return f(o, s, planes, w, h, bps, spp, alpha, opaque, cs, bpr, bpp);
}

static unsigned char *vm_bytes(id o, SEL s)
{
    unsigned char *(*f)(id, SEL) =
        (unsigned char *(*)(id, SEL))objc_msgSend;

    return f(o, s);
}

/* Display state (one window). */
static id g_bitmap = 0;
static NSRect_ g_view_rect;
static long g_pix_w = 0;
static long g_pix_h = 0;
static Class g_view_super;
static int g_drag_on = 0;

static void VisDrawRect(id self, SEL cmd, NSRect_ r)
{
    (void)self;
    (void)cmd;
    (void)r;
    if (g_bitmap != 0) {
        vm_v_rect(g_bitmap, sel_registerName("drawInRect:"), g_view_rect);
    }
}

/* AppKit resizes the content view with its window. Intercepting the
 * size change (points) into g_view_rect lets the frame loop re-derive
 * backing pixels without ever calling a struct-returning method —
 * NSView must still run to commit the frame, hence msgSendSuper. */
static void VisSetFrameSize(id self, SEL cmd, NSSize_ s)
{
    struct objc_super sup;
    void (*f)(struct objc_super *, SEL, NSSize_) =
        (void (*)(struct objc_super *, SEL, NSSize_))objc_msgSendSuper;

    g_view_rect.s = s;
    sup.receiver = self;
    sup.super_class = g_view_super;
    f(&sup, cmd, s);
}

/* Fresh meshed RGB bitmap rep at w*h pixels (isPlanar:0 — the BOOL/int
 * args are the banding fix; keep them exact). Caller owns the rep. */
static id vm_make_bitmap(long w, long h)
{
    id rep;

    rep = vm_bitmap(vm_id((id)objc_getClass("NSBitmapImageRep"),
                          sel_registerName("alloc")),
                    sel_registerName("initWithBitmapDataPlanes:"
                                     "pixelsWide:pixelsHigh:"
                                     "bitsPerSample:samplesPerPixel:"
                                     "hasAlpha:isPlanar:colorSpaceName:"
                                     "bytesPerRow:bitsPerPixel:"),
                    NULL, w, h, 8L, 3L, 0, 0,
                    (id)NSDeviceRGBColorSpace, w * 3L, 24L);
    if (rep == 0 || vm_bytes(rep, sel_registerName("bitmapData")) == NULL) {
        if (rep != 0) {
            vm_id(rep, sel_registerName("release"));
        }
        return 0;
    }
    return rep;
}

/* If the view's backing size moved, rebuild the bitmap rep and the sim
 * framebuffer at the new pixel size. Keeps the old rep on allocation
 * failure (retried on the next frame). g_fb_* track the framebuffer's
 * real dims separately so a failed stc_vis_resize is retried each frame
 * instead of leaving the window permanently blank. */
static long g_fb_w = 0;
static long g_fb_h = 0;

static int vm_check_resize(StcVis *v, id win)
{
    double scale;
    long nw;
    long nh;
    id nb;

    scale = vm_dbl(win, sel_registerName("backingScaleFactor"));
    if (!(scale >= 1.0)) {
        scale = 1.0;
    }
    nw = (long)(g_view_rect.s.w * scale + 0.5);
    nh = (long)(g_view_rect.s.h * scale + 0.5);
    if (nw <= 0 || nh <= 0) {
        return 0;
    }
    if (nw == g_pix_w && nh == g_pix_h && nw == g_fb_w && nh == g_fb_h) {
        return 0;
    }
    if (nw != g_pix_w || nh != g_pix_h) {
        nb = vm_make_bitmap(nw, nh);
        if (nb == 0) {
            return 0;
        }
        vm_id(g_bitmap, sel_registerName("release"));
        g_bitmap = nb;
        g_pix_w = nw;
        g_pix_h = nh;
    }
    if (stc_vis_resize(v, (int)nw, (int)nh) == 0) {
        g_fb_w = nw;
        g_fb_h = nh;
    }
    return 1;
}

/* The merger renderer is CPU-heavy, so its framebuffer (and therefore the
 * bitmap rep, which drawInRect: scales to fill the view regardless of its
 * pixel size) is capped: the LONGER edge is held to this many pixels (a
 * tall narrow window would otherwise cap only the width and leave the
 * height, e.g. 200 x 8000, unbounded), the other edge scaled down by the
 * same factor to keep the backing store's aspect ratio. Reuses the same
 * bitmap and pixel-size globals as the visual sandbox (only one of the
 * two runs per process). */
#define STC_MVIS_MAX_W 1920

static int vm_check_resize_mvis(StcMVis *v, id win)
{
    double scale;
    long nw;
    long nh;
    long fw;
    long fh;
    long longer;
    id nb;

    scale = vm_dbl(win, sel_registerName("backingScaleFactor"));
    if (!(scale >= 1.0)) {
        scale = 1.0;
    }
    nw = (long)(g_view_rect.s.w * scale + 0.5);
    nh = (long)(g_view_rect.s.h * scale + 0.5);
    if (nw <= 0 || nh <= 0) {
        return 0;
    }
    fw = nw;
    fh = nh;
    longer = fw > fh ? fw : fh;
    if (longer > STC_MVIS_MAX_W) {
        double k = (double)STC_MVIS_MAX_W / (double)longer;

        fw = (long)((double)fw * k + 0.5);
        fh = (long)((double)fh * k + 0.5);
        if (fw < 1) {
            fw = 1;
        }
        if (fh < 1) {
            fh = 1;
        }
    }
    if (fw == g_fb_w && fh == g_fb_h) {
        return 0;
    }
    /* Resize the framebuffer before touching the bitmap: while the
     * window sits at a size stc_mvis_resize rejects (e.g. below the
     * renderer's minimum, which live-resize dragging passes through
     * continuously), this returns 0 every frame without ever
     * allocating/releasing an NSBitmapImageRep. g_pix_w/h (the bitmap's
     * real size, which the blit trusts) and g_fb_w/h only ever advance
     * together, once both the framebuffer resize and the matching
     * bitmap allocation have actually succeeded — never leaving the fb
     * and bitmap sizes mismatched (see the blit's exact bytesPerRow
     * check in stc_mvis_run). */
    if (stc_mvis_resize(v, (int)fw, (int)fh) != 0) {
        return 0;
    }
    nb = vm_make_bitmap(fw, fh);
    if (nb == 0) {
        return 0;
    }
    if (g_bitmap != 0) {
        vm_id(g_bitmap, sel_registerName("release"));
    }
    g_bitmap = nb;
    g_pix_w = fw;
    g_pix_h = fh;
    g_fb_w = fw;
    g_fb_h = fh;
    return 1;
}

/* Same framebuffer-capping rule as the merger view (the space simulator
 * is at least as CPU-heavy: N-body + Barnes-Hut + 1PN), and the same
 * longer-edge cap and resize-before-bitmap-swap ordering as
 * vm_check_resize_mvis (see its comments). Reuses the same bitmap/
 * pixel-size globals (only one of stc_vis_run/stc_mvis_run/stc_svis_run
 * runs per process). */
#define STC_SVIS_MAX_W 1920

static int vm_check_resize_svis(StcSVis *v, id win)
{
    double scale;
    long nw;
    long nh;
    long fw;
    long fh;
    long longer;
    id nb;

    scale = vm_dbl(win, sel_registerName("backingScaleFactor"));
    if (!(scale >= 1.0)) {
        scale = 1.0;
    }
    nw = (long)(g_view_rect.s.w * scale + 0.5);
    nh = (long)(g_view_rect.s.h * scale + 0.5);
    if (nw <= 0 || nh <= 0) {
        return 0;
    }
    fw = nw;
    fh = nh;
    longer = fw > fh ? fw : fh;
    if (longer > STC_SVIS_MAX_W) {
        double k = (double)STC_SVIS_MAX_W / (double)longer;

        fw = (long)((double)fw * k + 0.5);
        fh = (long)((double)fh * k + 0.5);
        if (fw < 1) {
            fw = 1;
        }
        if (fh < 1) {
            fh = 1;
        }
    }
    if (fw == g_fb_w && fh == g_fb_h) {
        return 0;
    }
    if (stc_svis_resize(v, (int)fw, (int)fh) != 0) {
        return 0;
    }
    nb = vm_make_bitmap(fw, fh);
    if (nb == 0) {
        return 0;
    }
    if (g_bitmap != 0) {
        vm_id(g_bitmap, sel_registerName("release"));
    }
    g_bitmap = nb;
    g_pix_w = fw;
    g_pix_h = fh;
    g_fb_w = fw;
    g_fb_h = fh;
    return 1;
}

/* NSEvent types consumed by the frame loop. */
#define VM_KEY_DOWN 10
#define VM_LMOUSE_DOWN 1
#define VM_LMOUSE_UP 2
#define VM_LMOUSE_DRAGGED 6
#define VM_SCROLL_WHEEL 22
#define VM_MAGNIFY 30

/* Orbit drag sensitivity, radians per point of mouse travel. */
#define VM_ORBIT_SENS 0.005

/* Map NSEvent key characters (UTF-8) to stc/ui.h key codes. */
static int vm_map_key(const char *utf8)
{
    unsigned char b0;
    unsigned char b1;
    unsigned char b2;

    if (utf8 == NULL || utf8[0] == '\0') {
        return STC_UI_KEY_NONE;
    }
    b0 = (unsigned char)utf8[0];
    if (b0 == 0x7Fu) {
        return 127;
    }
    if (b0 < 0x80u) {
        if (b0 == 0x09u) {
            return '\t';
        }
        if (b0 == 0x0Du) {
            return '\r';
        }
        if (b0 == 0x1Bu) {
            return STC_UI_KEY_ESC;
        }
        return (int)b0;
    }
    /* Arrow keys arrive as U+F700..U+F703 (EF 9C 80..83). */
    if (utf8[1] == '\0' || utf8[2] == '\0') {
        return STC_UI_KEY_NONE;
    }
    b1 = (unsigned char)utf8[1];
    b2 = (unsigned char)utf8[2];
    if (b0 == 0xEFu && b1 == 0x9Cu) {
        if (b2 == 0x80u) {
            return STC_UI_KEY_UP;
        }
        if (b2 == 0x81u) {
            return STC_UI_KEY_DOWN;
        }
        if (b2 == 0x82u) {
            return STC_UI_KEY_LEFT;
        }
        if (b2 == 0x83u) {
            return STC_UI_KEY_RIGHT;
        }
        /* Forward-delete (fn+Delete) arrives as U+F728 (EF 9C A8). No
         * cursor-position API exists here (single-line console editing
         * only truncates from the end), so it is treated the same as
         * Backspace/DEL. */
        if (b2 == 0xA8u) {
            return 127;
        }
    }
    return STC_UI_KEY_NONE;
}

/* Wall-clock seconds (monotonic), or -1.0 on failure. clock() measures
 * process CPU time summed across threads, so with the worker-pool
 * renderers (stc/par.h) a frame that keeps several cores busy for
 * 1/60 s of wall time reports several times that much "elapsed" CPU
 * time: pacing never engaged (a windowed --frames run finished in far
 * less than N/60 s of wall time). CLOCK_MONOTONIC is wall time and
 * immune to that. */
static double vm_now_s(void)
{
    struct timespec ts;

    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
        return -1.0;
    }
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1.0e-9;
}

static void vm_pace(double start, int no_adapt)
{
    double elapsed;

    if (no_adapt || start < 0.0) {
        return;
    }
    do {
        double now = vm_now_s();

        if (now < 0.0) {
            return;
        }
        elapsed = now - start;
    } while (elapsed < STC_APP_FRAME_DT);
}

/* Shared window/app/view setup for both the visual sandbox and the
 * merger view: NSApplication + the (once-registered) StcVisView class +
 * an NSWindow/NSView pair sized w_pt x h_pt points, titled `title`.
 * g_view_rect is set here. Returns 0 on success (app/win/view filled),
 * 1 on failure (stderr already printed; caller pops its autorelease pool
 * and returns). Does not touch the bitmap or pixel-size globals, and
 * does not activate the app: callers size and allocate the bitmap
 * themselves since the visual sandbox and the merger view use different
 * framebuffer-sizing rules. */
static int vm_open_window(const char *title, double w_pt, double h_pt,
                          id *app_out, id *win_out, id *view_out)
{
    id app;
    id win;
    id view;
    id strcls;
    id titleid;
    Class cls;

    app = vm_id((id)objc_getClass("NSApplication"),
                sel_registerName("sharedApplication"));
    if (app == 0) {
        fprintf(stderr, "stlab: visual unavailable (no NSApplication)\n");
        return 1;
    }
    vm_v_long(app, sel_registerName("setActivationPolicy:"), 0L);
    vm_id(app, sel_registerName("finishLaunching"));

    /* View subclass with a drawRect: that blits our bitmap. Registered
     * once per process: objc_allocateClassPair fails (returns 0) if the
     * name already exists, so fall back to looking it up. */
    cls = (Class)objc_getClass("StcVisView");
    if (cls == 0) {
        Class super = objc_getClass("NSView");
        Class newcls = objc_allocateClassPair(super, "StcVisView", 0);

        if (newcls == 0) {
            fprintf(stderr, "stlab: visual unavailable (no view class)\n");
            return 1;
        }
        class_addMethod(newcls, sel_registerName("drawRect:"),
                        (IMP)VisDrawRect,
                        "v@:{CGRect={CGPoint=dd}{CGSize=dd}}");
        class_addMethod(newcls, sel_registerName("setFrameSize:"),
                        (IMP)VisSetFrameSize, "v@:{CGSize=dd}");
        g_view_super = super;
        objc_registerClassPair(newcls);
        cls = newcls;
    } else {
        g_view_super = objc_getClass("NSView");
    }

    g_view_rect = NSMakeRect_(0.0, 0.0, w_pt, h_pt);
    win = vm_window(vm_id((id)objc_getClass("NSWindow"),
                          sel_registerName("alloc")),
                    sel_registerName("initWithContentRect:styleMask:"
                                     "backing:defer:"),
                    NSMakeRect_(200.0, 200.0, w_pt, h_pt),
                    15UL, 2UL, 0);
    if (win == 0) {
        fprintf(stderr, "stlab: visual unavailable (no window)\n");
        return 1;
    }
    strcls = (id)objc_getClass("NSString");
    /* stringWithUTF8String: takes (const char *): send fully typed. */
    {
        id (*f)(id, SEL, const char *) =
            (id(*)(id, SEL, const char *))objc_msgSend;

        titleid = f(strcls, sel_registerName("stringWithUTF8String:"),
                   title);
    }
    if (titleid != 0) {
        vm_v_id(win, sel_registerName("setTitle:"), titleid);
    }
    /* Flat chrome: no window drop shadow (BOOL arg passed as int). */
    vm_v_bool(win, sel_registerName("setHasShadow:"), 0);
    view = vm_id_rect(vm_id((id)cls, sel_registerName("alloc")),
                      sel_registerName("initWithFrame:"), g_view_rect);
    if (view == 0) {
        fprintf(stderr, "stlab: visual unavailable (no view)\n");
        return 1;
    }
    vm_v_id(win, sel_registerName("setContentView:"), view);
    *app_out = app;
    *win_out = win;
    *view_out = view;
    return 0;
}

int stc_vis_run(const StcAppOpts *o, int preset, unsigned long long seed)
{
    StcAppOpts def;
    StcVis v;
    id app;
    id win;
    id view;
    id datecls;
    id past;
    double scale = 1.0;
    long cap;
    int quit = 0;
    void *setup_pool = NULL;

    if (o == NULL) {
        stc_app_opts_init(&def);
        o = &def;
    }
    if (preset < 0 || preset >= stc_preset_count()) {
        fprintf(stderr, "stlab: bad --preset %d (0..%d)\n", preset,
                stc_preset_count() - 1);
        return 2;
    }
    cap = o->frames_given ? o->frames : -1L;
    setup_pool = objc_autoreleasePoolPush();

    if (vm_open_window("Spacetime Lab", (double)STC_VIS_DEF_W,
                       (double)STC_VIS_DEF_H, &app, &win, &view) != 0) {
        objc_autoreleasePoolPop(setup_pool);
        return 1;
    }
    scale = vm_dbl(win, sel_registerName("backingScaleFactor"));
    if (!(scale >= 1.0)) {
        scale = 1.0;
    }
    g_pix_w = (long)((double)STC_VIS_DEF_W * scale + 0.5);
    g_pix_h = (long)((double)STC_VIS_DEF_H * scale + 0.5);

    g_bitmap = vm_make_bitmap(g_pix_w, g_pix_h);
    if (g_bitmap == 0) {
        fprintf(stderr, "stlab: visual unavailable (no bitmap)\n");
        objc_autoreleasePoolPop(setup_pool);
        return 1;
    }

    if (stc_vis_init(&v, preset, seed, o->integrator, (int)g_pix_w,
                     (int)g_pix_h) != 0) {
        fprintf(stderr, "stlab: visual unavailable (init failed)\n");
        objc_autoreleasePoolPop(setup_pool);
        return 1;
    }
    v.csv_path = o->csv;
    v.ppm_path = o->ppm;

    vm_v_id(win, sel_registerName("makeKeyAndOrderFront:"), 0);
    vm_v_bool(app, sel_registerName("activateIgnoringOtherApps:"), 1);

    datecls = (id)objc_getClass("NSDate");
    for (;;) {
        void *pool = objc_autoreleasePoolPush();
        id ev;

        past = vm_id(datecls, sel_registerName("distantPast"));
        while ((ev = vm_next_event(app,
                                   sel_registerName("nextEventMatchingMask:"
                                                    "untilDate:inMode:"
                                                    "dequeue:"),
                                   ~0UL, past,
                                   (id)NSDefaultRunLoopMode, 1)) != 0) {
            long etype = vm_long(ev, sel_registerName("type"));

            if (etype == VM_KEY_DOWN) {
                id chars = vm_id(ev, sel_registerName("characters"));
                const char *utf8 = chars != 0
                    ? vm_cstr(chars, sel_registerName("UTF8String"))
                    : NULL;
                int key = vm_map_key(utf8);

                if (key != STC_UI_KEY_NONE && stc_vis_key(&v, key)) {
                    quit = 1;
                }
            } else if (etype == VM_LMOUSE_DOWN) {
                g_drag_on = 1;
            } else if (etype == VM_LMOUSE_UP) {
                g_drag_on = 0;
            } else if (etype == VM_LMOUSE_DRAGGED && g_drag_on) {
                double dx = vm_dbl(ev, sel_registerName("deltaX"));
                double dy = vm_dbl(ev, sel_registerName("deltaY"));

                if (isfinite(dx) && isfinite(dy)) {
                    v.cam.yaw_rad += dx * VM_ORBIT_SENS;
                    v.cam.pitch_rad += dy * VM_ORBIT_SENS;
                    stc_ui_camera_update(&v.cam);
                }
            } else if (etype == VM_SCROLL_WHEEL) {
                double dy = vm_dbl(ev,
                                   sel_registerName("scrollingDeltaY"));

                if (isfinite(dy) && dy != 0.0) {
                    v.cam.dist_m *= pow(0.95, dy);
                    stc_ui_camera_update(&v.cam);
                }
            } else if (etype == VM_MAGNIFY) {
                double mag = vm_dbl(ev, sel_registerName("magnification"));

                if (isfinite(mag) && mag > -0.9) {
                    v.cam.dist_m /= (1.0 + mag);
                    stc_ui_camera_update(&v.cam);
                }
            }
            /* keyDown events are consumed by stc_vis_key; forwarding them
             * to sendEvent would walk the responder chain and NSBeep on
             * every press — we have no menu key equivalents to preserve. */
            if (etype != VM_KEY_DOWN) {
                vm_v_id(app, sel_registerName("sendEvent:"), ev);
            }
            if (quit) {
                break;
            }
        }
        objc_autoreleasePoolPop(pool);
        if (quit || !vm_bool(win, sel_registerName("isVisible"))) {
            break;
        }
        /* Resize is honored before the first draw too: setFrameSize:
         * may land during makeKeyAndOrderFront or while events drain. */
        (void)vm_check_resize(&v, win);
        {
            double t0 = vm_now_s();
            unsigned char *dst;

            stc_vis_step(&v);
            stc_vis_render(&v);
            dst = vm_bytes(g_bitmap, sel_registerName("bitmapData"));
            if (dst != NULL && v.fb.px != NULL) {
                long bpr = vm_long(g_bitmap,
                                   sel_registerName("bytesPerRow"));
                if (bpr == g_pix_w * 3L) {
                    memcpy(dst, v.fb.px,
                           (size_t)g_pix_w * (size_t)g_pix_h * 3u);
                } else if (bpr > g_pix_w * 3L) {
                    int row;
                    for (row = 0; row < g_pix_h; row++) {
                        memcpy(dst + (size_t)row * (size_t)bpr,
                               v.fb.px + (size_t)row * (size_t)g_pix_w * 3u,
                               (size_t)g_pix_w * 3u);
                    }
                }
            }
            vm_v_bool(view, sel_registerName("setNeedsDisplay:"), 1);
            vm_id(app, sel_registerName("updateWindows"));
            if (cap >= 0 && v.frames >= cap) {
                break;
            }
            vm_pace(t0, o->no_adapt);
        }
    }
    printf("\nbye: t=%g s steps=%ld bodies=%d\n", v.sim.t, v.frames,
           v.sim.n);
    stc_vis_free(&v);
    objc_autoreleasePoolPop(setup_pool);
    return 0;
}

/* Windowed cinematic binary black hole merger session: same event loop
 * shape as stc_vis_run (drain events, orbit/zoom on drag/scroll/magnify,
 * resize, step, render, blit), against the stc/vis_merger.h API instead
 * of stc/vis.h. Framebuffer capped at STC_MVIS_MAX_W (see
 * vm_check_resize_mvis). */
int stc_mvis_run(const StcAppOpts *o, int preset)
{
    StcAppOpts def;
    StcMVis v;
    id app;
    id win;
    id view;
    id datecls;
    id past;
    double scale = 1.0;
    long cap;
    long nw;
    long nh;
    long fw;
    long fh;
    int quit = 0;
    void *setup_pool = NULL;

    if (o == NULL) {
        stc_app_opts_init(&def);
        o = &def;
    }
    if (preset < 0 || preset >= STC_MVIS_NPRESET) {
        fprintf(stderr, "stlab: bad --merger preset %d (0..%d)\n", preset,
                STC_MVIS_NPRESET - 1);
        return 2;
    }
    cap = o->frames_given ? o->frames : -1L;
    setup_pool = objc_autoreleasePoolPush();

    if (vm_open_window("Spacetime Lab \xE2\x80\x94 Black Hole Merger",
                       (double)STC_VIS_DEF_W, (double)STC_VIS_DEF_H,
                       &app, &win, &view) != 0) {
        objc_autoreleasePoolPop(setup_pool);
        return 1;
    }
    scale = vm_dbl(win, sel_registerName("backingScaleFactor"));
    if (!(scale >= 1.0)) {
        scale = 1.0;
    }
    nw = (long)((double)STC_VIS_DEF_W * scale + 0.5);
    nh = (long)((double)STC_VIS_DEF_H * scale + 0.5);
    fw = nw;
    fh = nh;
    if (fw > STC_MVIS_MAX_W) {
        fh = (long)((double)fh * (double)STC_MVIS_MAX_W / (double)fw + 0.5);
        fw = STC_MVIS_MAX_W;
        if (fh < 1) {
            fh = 1;
        }
    }

    g_bitmap = vm_make_bitmap(fw, fh);
    if (g_bitmap == 0) {
        fprintf(stderr, "stlab: visual unavailable (no bitmap)\n");
        objc_autoreleasePoolPop(setup_pool);
        return 1;
    }
    g_pix_w = fw;
    g_pix_h = fh;
    g_fb_w = fw;
    g_fb_h = fh;

    if (stc_mvis_init(&v, preset, (int)fw, (int)fh) != 0) {
        fprintf(stderr, "stlab: visual unavailable (merger init failed)\n");
        objc_autoreleasePoolPop(setup_pool);
        return 1;
    }

    vm_v_id(win, sel_registerName("makeKeyAndOrderFront:"), 0);
    vm_v_bool(app, sel_registerName("activateIgnoringOtherApps:"), 1);

    datecls = (id)objc_getClass("NSDate");
    for (;;) {
        void *pool = objc_autoreleasePoolPush();
        id ev;

        past = vm_id(datecls, sel_registerName("distantPast"));
        while ((ev = vm_next_event(app,
                                   sel_registerName("nextEventMatchingMask:"
                                                    "untilDate:inMode:"
                                                    "dequeue:"),
                                   ~0UL, past,
                                   (id)NSDefaultRunLoopMode, 1)) != 0) {
            long etype = vm_long(ev, sel_registerName("type"));

            if (etype == VM_KEY_DOWN) {
                id chars = vm_id(ev, sel_registerName("characters"));
                const char *utf8 = chars != 0
                    ? vm_cstr(chars, sel_registerName("UTF8String"))
                    : NULL;
                int key = vm_map_key(utf8);

                if (key != STC_UI_KEY_NONE && stc_mvis_key(&v, key)) {
                    quit = 1;
                }
            } else if (etype == VM_LMOUSE_DOWN) {
                g_drag_on = 1;
            } else if (etype == VM_LMOUSE_UP) {
                g_drag_on = 0;
            } else if (etype == VM_LMOUSE_DRAGGED && g_drag_on) {
                double dx = vm_dbl(ev, sel_registerName("deltaX"));
                double dy = vm_dbl(ev, sel_registerName("deltaY"));

                if (isfinite(dx) && isfinite(dy)) {
                    stc_mvis_orbit(&v, dx * 0.005, dy * 0.005);
                }
            } else if (etype == VM_SCROLL_WHEEL) {
                double dy = vm_dbl(ev,
                                   sel_registerName("scrollingDeltaY"));

                if (isfinite(dy) && dy != 0.0) {
                    stc_mvis_zoom(&v, pow(0.95, dy));
                }
            } else if (etype == VM_MAGNIFY) {
                double mag = vm_dbl(ev, sel_registerName("magnification"));

                if (isfinite(mag) && mag > -0.9) {
                    stc_mvis_zoom(&v, 1.0 / (1.0 + mag));
                }
            }
            /* keyDown events are consumed by stc_mvis_key; forwarding them
             * to sendEvent would walk the responder chain and NSBeep on
             * every press — we have no menu key equivalents to preserve. */
            if (etype != VM_KEY_DOWN) {
                vm_v_id(app, sel_registerName("sendEvent:"), ev);
            }
            if (quit) {
                break;
            }
        }
        objc_autoreleasePoolPop(pool);
        if (quit || !vm_bool(win, sel_registerName("isVisible"))) {
            break;
        }
        /* Resize is honored before the first draw too: setFrameSize:
         * may land during makeKeyAndOrderFront or while events drain. */
        (void)vm_check_resize_mvis(&v, win);
        {
            double t0 = vm_now_s();
            unsigned char *dst;

            stc_mvis_step(&v);
            stc_mvis_render(&v);
            dst = vm_bytes(g_bitmap, sel_registerName("bitmapData"));
            if (dst != NULL && v.fb.px != NULL) {
                long bpr = vm_long(g_bitmap,
                                   sel_registerName("bytesPerRow"));
                if (bpr == g_pix_w * 3L) {
                    memcpy(dst, v.fb.px,
                           (size_t)g_pix_w * (size_t)g_pix_h * 3u);
                } else if (bpr > g_pix_w * 3L) {
                    int row;
                    for (row = 0; row < g_pix_h; row++) {
                        memcpy(dst + (size_t)row * (size_t)bpr,
                               v.fb.px + (size_t)row * (size_t)g_pix_w * 3u,
                               (size_t)g_pix_w * 3u);
                    }
                }
            }
            vm_v_bool(view, sel_registerName("setNeedsDisplay:"), 1);
            vm_id(app, sel_registerName("updateWindows"));
            if (cap >= 0 && v.frames >= cap) {
                break;
            }
            vm_pace(t0, o->no_adapt);
        }
    }
    printf("\nbye: frames=%ld t=%g M\n", v.frames, v.t);
    stc_mvis_free(&v);
    objc_autoreleasePoolPop(setup_pool);
    return 0;
}

/* Windowed space-simulator session: same event-loop shape as
 * stc_mvis_run (drain events, keys, orbit/zoom on drag/scroll/magnify,
 * resize, step, render, blit) against the stc/vis_space.h API, plus
 * click-to-pick: a left mouse-up whose drag traveled under 3 px (view
 * points) since mouse-down is treated as a click and resolves to a
 * framebuffer pixel via the event's locationInWindow (the content view
 * fills the window with origin (0,0), so window-space points equal
 * view-space points here). */
int stc_svis_run(const StcAppOpts *o, int scenario)
{
    StcAppOpts def;
    StcSVis v;
    id app;
    id win;
    id view;
    id datecls;
    id past;
    double scale = 1.0;
    long cap;
    long nw;
    long nh;
    long fw;
    long fh;
    double drag_dist = 0.0;
    int quit = 0;
    void *setup_pool = NULL;

    if (o == NULL) {
        stc_app_opts_init(&def);
        o = &def;
    }
    if (scenario < 0 || scenario >= STC_SVIS_NSCEN) {
        fprintf(stderr, "stlab: bad --space preset %d (0..%d)\n", scenario,
                STC_SVIS_NSCEN - 1);
        return 2;
    }
    cap = o->frames_given ? o->frames : -1L;
    setup_pool = objc_autoreleasePoolPush();

    if (vm_open_window("Spacetime Lab \xE2\x80\x94 Space Simulator",
                       (double)STC_VIS_DEF_W, (double)STC_VIS_DEF_H,
                       &app, &win, &view) != 0) {
        objc_autoreleasePoolPop(setup_pool);
        return 1;
    }
    scale = vm_dbl(win, sel_registerName("backingScaleFactor"));
    if (!(scale >= 1.0)) {
        scale = 1.0;
    }
    nw = (long)((double)STC_VIS_DEF_W * scale + 0.5);
    nh = (long)((double)STC_VIS_DEF_H * scale + 0.5);
    fw = nw;
    fh = nh;
    if (fw > STC_SVIS_MAX_W) {
        fh = (long)((double)fh * (double)STC_SVIS_MAX_W / (double)fw + 0.5);
        fw = STC_SVIS_MAX_W;
        if (fh < 1) {
            fh = 1;
        }
    }

    g_bitmap = vm_make_bitmap(fw, fh);
    if (g_bitmap == 0) {
        fprintf(stderr, "stlab: visual unavailable (no bitmap)\n");
        objc_autoreleasePoolPop(setup_pool);
        return 1;
    }
    g_pix_w = fw;
    g_pix_h = fh;
    g_fb_w = fw;
    g_fb_h = fh;

    if (stc_svis_init(&v, scenario, (int)fw, (int)fh) != 0) {
        fprintf(stderr, "stlab: visual unavailable (space init failed)\n");
        objc_autoreleasePoolPop(setup_pool);
        return 1;
    }

    vm_v_id(win, sel_registerName("makeKeyAndOrderFront:"), 0);
    vm_v_bool(app, sel_registerName("activateIgnoringOtherApps:"), 1);

    datecls = (id)objc_getClass("NSDate");
    for (;;) {
        void *pool = objc_autoreleasePoolPush();
        id ev;

        past = vm_id(datecls, sel_registerName("distantPast"));
        while ((ev = vm_next_event(app,
                                   sel_registerName("nextEventMatchingMask:"
                                                    "untilDate:inMode:"
                                                    "dequeue:"),
                                   ~0UL, past,
                                   (id)NSDefaultRunLoopMode, 1)) != 0) {
            long etype = vm_long(ev, sel_registerName("type"));

            if (etype == VM_KEY_DOWN) {
                id chars = vm_id(ev, sel_registerName("characters"));
                const char *utf8 = chars != 0
                    ? vm_cstr(chars, sel_registerName("UTF8String"))
                    : NULL;
                int key = vm_map_key(utf8);

                if (key != STC_UI_KEY_NONE && stc_svis_key(&v, key)) {
                    quit = 1;
                }
            } else if (etype == VM_LMOUSE_DOWN) {
                g_drag_on = 1;
                drag_dist = 0.0;
            } else if (etype == VM_LMOUSE_UP) {
                if (g_drag_on && drag_dist < 3.0) {
                    NSPoint_ pt = vm_point(ev,
                                          sel_registerName("locationInWindow"));
                    double vw = g_view_rect.s.w;
                    double vh = g_view_rect.s.h;

                    if (vw > 0.0 && vh > 0.0 && isfinite(pt.x) &&
                        isfinite(pt.y)) {
                        int px = (int)(pt.x / vw * (double)v.w);
                        int py = (int)((vh - pt.y) / vh * (double)v.h);

                        stc_svis_pick(&v, px, py);
                    }
                }
                g_drag_on = 0;
            } else if (etype == VM_LMOUSE_DRAGGED && g_drag_on) {
                double dx = vm_dbl(ev, sel_registerName("deltaX"));
                double dy = vm_dbl(ev, sel_registerName("deltaY"));

                if (isfinite(dx) && isfinite(dy)) {
                    drag_dist += fabs(dx) + fabs(dy);
                    stc_svis_orbit(&v, dx * 0.005, dy * 0.005);
                }
            } else if (etype == VM_SCROLL_WHEEL) {
                double dy = vm_dbl(ev,
                                   sel_registerName("scrollingDeltaY"));

                if (isfinite(dy) && dy != 0.0) {
                    stc_svis_zoom(&v, pow(0.93, dy));
                }
            } else if (etype == VM_MAGNIFY) {
                double mag = vm_dbl(ev, sel_registerName("magnification"));

                if (isfinite(mag) && mag > -0.9) {
                    stc_svis_zoom(&v, 1.0 / (1.0 + mag));
                }
            }
            /* keyDown events are consumed by stc_svis_key (including the
             * command console's printable/Enter/Backspace/Esc/'/' keys
             * via vm_map_key); forwarding them to sendEvent would walk
             * the responder chain and NSBeep on every press. */
            if (etype != VM_KEY_DOWN) {
                vm_v_id(app, sel_registerName("sendEvent:"), ev);
            }
            if (quit) {
                break;
            }
        }
        objc_autoreleasePoolPop(pool);
        if (quit || !vm_bool(win, sel_registerName("isVisible"))) {
            break;
        }
        (void)vm_check_resize_svis(&v, win);
        {
            double t0 = vm_now_s();
            unsigned char *dst;

            stc_svis_step(&v);
            stc_svis_render(&v);
            dst = vm_bytes(g_bitmap, sel_registerName("bitmapData"));
            if (dst != NULL && v.fb.px != NULL) {
                long bpr = vm_long(g_bitmap,
                                   sel_registerName("bytesPerRow"));
                if (bpr == g_pix_w * 3L) {
                    memcpy(dst, v.fb.px,
                           (size_t)g_pix_w * (size_t)g_pix_h * 3u);
                } else if (bpr > g_pix_w * 3L) {
                    int row;
                    for (row = 0; row < g_pix_h; row++) {
                        memcpy(dst + (size_t)row * (size_t)bpr,
                               v.fb.px + (size_t)row * (size_t)g_pix_w * 3u,
                               (size_t)g_pix_w * 3u);
                    }
                }
            }
            vm_v_bool(view, sel_registerName("setNeedsDisplay:"), 1);
            vm_id(app, sel_registerName("updateWindows"));
            if (cap >= 0 && v.frames >= cap) {
                break;
            }
            vm_pace(t0, o->no_adapt);
        }
    }
    printf("\nbye: frames=%ld t=%g s\n", v.frames, v.sim.t);
    stc_svis_free(&v);
    objc_autoreleasePoolPop(setup_pool);
    return 0;
}

#else /* !__APPLE__ */

int stc_vis_run(const StcAppOpts *o, int preset, unsigned long long seed)
{
    (void)o;
    (void)preset;
    (void)seed;
    fprintf(stderr, "stlab: --visual needs macOS (AppKit)\n");
    return 1;
}

int stc_mvis_run(const StcAppOpts *o, int preset)
{
    (void)o;
    (void)preset;
    fprintf(stderr, "stlab: --merger --visual needs macOS (AppKit)\n");
    return 1;
}

int stc_svis_run(const StcAppOpts *o, int scenario)
{
    (void)o;
    (void)scenario;
    fprintf(stderr, "stlab: --space --visual needs macOS (AppKit)\n");
    return 1;
}

#endif /* __APPLE__ */
