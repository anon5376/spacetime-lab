#!/bin/sh
# Double-click to build (if needed) and run Spacetime Lab's space simulator.
cd "$(dirname "$0")" || exit 1
[ -x build/stlab ] || make -s all
exec ./build/stlab --space
