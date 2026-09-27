CC = cc
CFLAGS = -std=c11 -Wall -Wextra -Werror -O2 -I include
LDLIBS = -lm
# Visual sandbox links macOS system frameworks only (no third-party
# dependencies); elsewhere the visual TU compiles to an unsupported stub.
ifeq ($(shell uname),Darwin)
FRAMEWORKS = -framework AppKit -framework Foundation
else
FRAMEWORKS =
# src/core/par.c uses POSIX threads directly. macOS resolves pthread_*
# from libSystem with no extra flag; Linux's libc needs -lpthread linked
# explicitly. Harmless everywhere else (every binary pulls in par.c via
# CORE_SRCS, so this applies uniformly, not just to stlab).
LDLIBS += -lpthread
endif

CORE_SRCS = src/core/vec.c src/core/units.c src/core/rng.c src/core/fmt.c \
	src/core/par.c
# Canonical sources only: src/gr/geodesic.c (CONTRACT Sec 8), src/kerr/kerr.c,
# src/scenario/scenario.c. Byte-identical duplicates removed 2026-09-27.
LIB_SRCS = src/sr/sr.c src/gr/geodesic.c src/embedding/embedding.c \
	src/lensing/lensing.c src/waves/waves.c src/kerr/kerr.c \
	src/clocks/clocks.c src/scenario/scenario.c src/render/render.c \
	src/body/body.c src/body/rk4.c src/bhut/bhut.c \
	src/tides/tides.c src/drag/drag.c src/frag/frag.c \
	src/elements/elements.c src/merger/merger.c \
	src/space/space.c src/space/ephem.c src/space/ephem_data.c
APP_SRCS = src/app/main.c src/app/ui.c src/app/csv.c src/app/sandbox.c \
	src/app/trails.c src/app/help.c src/app/spawn.c \
	src/vis/vis_render.c src/vis/vis_mac.c src/vis/vis_merger.c \
	src/vis/vis_space.c
VERIFY_SRCS = src/verify/verify.c
TEST_SRCS = src/verify/nbody_test.c
SOLAR_SRCS = src/verify/solar_check.c
FUZZ_SRCS = src/verify/fuzz.c
BENCH_SRCS = src/verify/bench.c
SPACE_SRCS = src/verify/space_check.c
EPHEM_SRCS = src/verify/ephem_check.c

ALL_BINS = build/stlab build/stverify build/sttest build/stsolar \
	build/stfuzz build/stbench build/stspace build/stephem
DEPFILES = $(ALL_BINS:%=%.d)

all: $(ALL_BINS)

# Header dependency tracking: binaries are built by one cc invocation over
# all sources, so -MMD in the compile line is unusable (each TU's depfile
# overwrites the previous; only the last TU's deps survive). Instead each
# recipe runs a separate `cc -MM -MP` pass and rewrites every `*.o:` stanza
# target to the binary path, emitting $@.d; -MP's phony header stubs keep a
# deleted header from breaking the build. -include $(DEPFILES) below turns
# those into real prerequisites on the next make run.
build/stlab: $(CORE_SRCS) $(LIB_SRCS) $(APP_SRCS) | build
	$(CC) $(CFLAGS) -MM -MP $(CORE_SRCS) $(LIB_SRCS) $(APP_SRCS) | sed 's|^[^ ]*\.o:|$@:|' > $@.d
	$(CC) $(CFLAGS) $(CORE_SRCS) $(LIB_SRCS) $(APP_SRCS) -o $@ $(LDLIBS) $(FRAMEWORKS)

build/stverify: $(CORE_SRCS) $(LIB_SRCS) $(VERIFY_SRCS) | build
	$(CC) $(CFLAGS) -MM -MP $(CORE_SRCS) $(LIB_SRCS) $(VERIFY_SRCS) | sed 's|^[^ ]*\.o:|$@:|' > $@.d
	$(CC) $(CFLAGS) $(CORE_SRCS) $(LIB_SRCS) $(VERIFY_SRCS) -o $@ $(LDLIBS)

build/sttest: $(CORE_SRCS) $(LIB_SRCS) $(TEST_SRCS) | build
	$(CC) $(CFLAGS) -MM -MP $(CORE_SRCS) $(LIB_SRCS) $(TEST_SRCS) | sed 's|^[^ ]*\.o:|$@:|' > $@.d
	$(CC) $(CFLAGS) $(CORE_SRCS) $(LIB_SRCS) $(TEST_SRCS) -o $@ $(LDLIBS)

build/stsolar: $(CORE_SRCS) $(LIB_SRCS) $(SOLAR_SRCS) | build
	$(CC) $(CFLAGS) -MM -MP $(CORE_SRCS) $(LIB_SRCS) $(SOLAR_SRCS) | sed 's|^[^ ]*\.o:|$@:|' > $@.d
	$(CC) $(CFLAGS) $(CORE_SRCS) $(LIB_SRCS) $(SOLAR_SRCS) -o $@ $(LDLIBS)

build/stfuzz: $(CORE_SRCS) $(LIB_SRCS) $(FUZZ_SRCS) | build
	$(CC) $(CFLAGS) -MM -MP $(CORE_SRCS) $(LIB_SRCS) $(FUZZ_SRCS) | sed 's|^[^ ]*\.o:|$@:|' > $@.d
	$(CC) $(CFLAGS) $(CORE_SRCS) $(LIB_SRCS) $(FUZZ_SRCS) -o $@ $(LDLIBS)

build/stbench: $(CORE_SRCS) $(LIB_SRCS) $(BENCH_SRCS) | build
	$(CC) $(CFLAGS) -MM -MP $(CORE_SRCS) $(LIB_SRCS) $(BENCH_SRCS) | sed 's|^[^ ]*\.o:|$@:|' > $@.d
	$(CC) $(CFLAGS) $(CORE_SRCS) $(LIB_SRCS) $(BENCH_SRCS) -o $@ $(LDLIBS)

build/stspace: $(CORE_SRCS) $(LIB_SRCS) $(SPACE_SRCS) | build
	$(CC) $(CFLAGS) -MM -MP $(CORE_SRCS) $(LIB_SRCS) $(SPACE_SRCS) | sed 's|^[^ ]*\.o:|$@:|' > $@.d
	$(CC) $(CFLAGS) $(CORE_SRCS) $(LIB_SRCS) $(SPACE_SRCS) -o $@ $(LDLIBS)

build/stephem: $(CORE_SRCS) $(LIB_SRCS) $(EPHEM_SRCS) | build
	$(CC) $(CFLAGS) -MM -MP $(CORE_SRCS) $(LIB_SRCS) $(EPHEM_SRCS) | sed 's|^[^ ]*\.o:|$@:|' > $@.d
	$(CC) $(CFLAGS) $(CORE_SRCS) $(LIB_SRCS) $(EPHEM_SRCS) -o $@ $(LDLIBS)

build:
	mkdir -p build

run: build/stlab
	./build/stlab --scenario twin --frames 60 --smoke --no-adapt

sandbox: build/stlab
	./build/stlab --sandbox --preset 0

merger: build/stlab
	./build/stlab --merger --visual

merger-shot: build/stlab
	mkdir -p build/merger && ./build/stlab --merger --frames 700 --visual-shot build/merger/shot.ppm

space: build/stlab
	./build/stlab --space

space-shot: build/stlab
	mkdir -p build/space && ./build/stlab --space --frames 300 --visual-shot build/space/shot.ppm

# Slow accuracy suite: integrate J2000 -> 2026-09-27 and compare every
# body with JPL Horizons (not part of `make verify`, which stays fast).
ephem-accuracy: build/stephem
	./build/stephem --accuracy

smoke: build/stlab
	./build/stlab --scenario all --frames 60 --ppm /tmp --csv /tmp --smoke --no-adapt

sandbox-smoke: build/stlab
	mkdir -p /tmp/sbx && for p in 0 1 2 3; do \
		./build/stlab --sandbox --preset $$p --frames 60 --ppm /tmp/sbx/p$$p.ppm --csv /tmp/sbx/p$$p.csv --smoke --no-adapt || exit 1; \
	done

# verify exercises the shipped stlab too: bare --smoke runs one
# deterministic frame per scenario with no output files; the sandbox smoke
# steps preset 0 headless for 2 frames, also writing nothing. stbench stays
# a manual target (`make bench`), not a pass/fail check (CONTRACT Sec 7).
verify: build/stverify build/sttest build/stsolar build/stfuzz build/stspace \
		build/stephem build/stlab
	./build/stverify && ./build/sttest && ./build/stsolar && ./build/stfuzz && \
	./build/stspace && ./build/stephem && \
	./build/stlab --smoke --no-adapt && \
	./build/stlab --sandbox --preset 0 --frames 2 --smoke --no-adapt

test: build/sttest
	./build/sttest

solar: build/stsolar
	./build/stsolar

fuzz: build/stfuzz
	./build/stfuzz

bench: build/stbench
	./build/stbench

clean:
	rm -f $(ALL_BINS) $(DEPFILES)

-include $(DEPFILES)

.PHONY: ephem-accuracy all run sandbox merger merger-shot space space-shot smoke \
	sandbox-smoke verify test solar fuzz bench clean
