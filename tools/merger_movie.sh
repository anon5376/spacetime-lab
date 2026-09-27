#!/bin/sh
# merger_movie.sh [preset] [frames] [every] [outdir]
#
# One-command MP4/GIF preview pipeline for `stlab --merger`:
# renders headless 1280x800 P6 PPM frames via
#   ./build/stlab --merger --preset P --frames F --movie <outdir>/frames \
#                 --movie-every K
# then encodes <outdir>/merger_p<P>.mp4 (libx264, yuv420p, crf 18,
# +faststart, framerate 60/K to match the live 60fps window) and a small
# 640px 20fps palette GIF.
set -eu

PRESET=${1:-0}
NFRAMES=${2:-2400}
EVERY=${3:-2}
OUTDIR=${4:-build/merger/movie}

die() {
    printf 'merger_movie: error: %s\n' "$*" >&2
    exit 1
}

case $PRESET in
    ''|*[!0-9]*) die "preset must be an integer 0..3, got '$PRESET'" ;;
esac
[ "$PRESET" -le 3 ] || die "preset must be 0..3, got '$PRESET'"
case $NFRAMES in
    ''|*[!0-9]*) die "frames must be a positive integer, got '$NFRAMES'" ;;
esac
[ "$NFRAMES" -ge 1 ] || die "frames must be >= 1, got '$NFRAMES'"
case $EVERY in
    ''|*[!0-9]*) die "every must be a positive integer, got '$EVERY'" ;;
esac
[ "$EVERY" -ge 1 ] || die "every must be >= 1, got '$EVERY'"

# Run from the repo root so `make` and ./build/stlab work regardless of
# the caller's cwd; a relative OUTDIR is then anchored at the root.
SCRIPT_DIR=$(CDPATH= cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR/.." || die "cannot cd to repo root"

FFMPEG=${FFMPEG:-/opt/homebrew/bin/ffmpeg}
if [ ! -x "$FFMPEG" ]; then
    FFMPEG=$(command -v ffmpeg 2>/dev/null || true)
fi
[ -n "$FFMPEG" ] && [ -x "$FFMPEG" ] \
    || die "ffmpeg not found (checked /opt/homebrew/bin/ffmpeg and PATH)"

STLAB=./build/stlab
if [ ! -x "$STLAB" ]; then
    printf 'merger_movie: %s missing, building\n' "$STLAB"
    make build/stlab || die "make build/stlab failed"
fi
if ! "$STLAB" --help 2>&1 | grep -q -- '--merger'; then
    # Binary may predate the --merger sources; make rebuilds if they are
    # newer, otherwise report the missing flag.
    make build/stlab || die "make build/stlab failed"
    "$STLAB" --help 2>&1 | grep -q -- '--merger' \
        || die "$STLAB has no --merger support (is src/vis/vis_merger.c wired into the Makefile?)"
fi

FRAMES_DIR=$OUTDIR/frames
mkdir -p "$FRAMES_DIR" || die "cannot create $FRAMES_DIR"
rm -f "$FRAMES_DIR"/merger_*.ppm

printf 'merger_movie: rendering preset %s, %s frames (every %s) -> %s\n' \
    "$PRESET" "$NFRAMES" "$EVERY" "$FRAMES_DIR"
"$STLAB" --merger --preset "$PRESET" --frames "$NFRAMES" \
    --movie "$FRAMES_DIR" --movie-every "$EVERY" \
    || die "stlab --merger run failed"

NWRITTEN=$(ls -1 "$FRAMES_DIR"/merger_*.ppm 2>/dev/null | wc -l | tr -d ' ')
[ "$NWRITTEN" -gt 0 ] \
    || die "stlab wrote no frames to $FRAMES_DIR"
printf 'merger_movie: %s frames written\n' "$NWRITTEN"

# 60fps live window sampled every K-th frame plays back at 60/K fps.
FPS=$(awk -v k="$EVERY" 'BEGIN { printf "%.10g", 60 / k }')

MP4=$OUTDIR/merger_p$PRESET.mp4
GIF=$OUTDIR/merger_p$PRESET.gif

printf 'merger_movie: encoding %s (%sfps)\n' "$MP4" "$FPS"
"$FFMPEG" -hide_banner -loglevel error -y \
    -framerate "$FPS" -pattern_type glob -i "$FRAMES_DIR/merger_*.ppm" \
    -c:v libx264 -pix_fmt yuv420p -crf 18 -movflags +faststart "$MP4" \
    || die "ffmpeg mp4 encode failed"

printf 'merger_movie: encoding %s (640px, 20fps)\n' "$GIF"
"$FFMPEG" -hide_banner -loglevel error -y \
    -framerate "$FPS" -pattern_type glob -i "$FRAMES_DIR/merger_*.ppm" \
    -filter_complex \
    '[0:v]fps=20,scale=640:-2:flags=lanczos,split[a][b];[a]palettegen[p];[b][p]paletteuse' \
    -loop 0 "$GIF" \
    || die "ffmpeg gif encode failed"

printf 'mp4: %s\n' "$MP4"
printf 'gif: %s\n' "$GIF"
printf 'frames: %s ppm in %s\n' "$NWRITTEN" "$FRAMES_DIR"
