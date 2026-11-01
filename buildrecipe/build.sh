#!/bin/bash
# buildrecipe/build.sh -- ngspice_augmentation private fork build
# (card #455 plan T3.1: the unmodified-base recipe; card #456 plan T3.2/T3.3:
# the --patches extension below, to build OUR changes -- patches/*.cir under
# this repo -- on top of the pin)
#
# Builds ngspice from the pinned pre-master-48 commit (see ../ngspice-base.txt)
# into a private, versioned prefix under ~/Work/apps/ngspice_augmentation/.
# Never touches /usr/local. See buildrecipe/README.md for the full story
# (why these flags, how this was matched to the system 45.2 build, the
# equivalence results, and how to roll back).
#
# Usage: buildrecipe/build.sh [--patches[=DIR]] [--prefix-name=NAME] [source_dir] [jobs]
#   source_dir     default: ~/Work/dev/ngspice_augmentation/ngspice-src
#   jobs           default: nproc
#   --patches[=DIR]  apply this repo's patch series on top of the pin before
#                    building (default DIR: <repo>/patches, which must hold
#                    numbered `git format-patch` files plus a `series` file
#                    listing them in order, one per line). Without this flag,
#                    the build is the plain, unmodified-pin recipe (card
#                    #455's original path) -- unchanged.
#   --prefix-name=NAME  install into ~/Work/apps/ngspice_augmentation/NAME
#                    instead of the default (the short commit hash of
#                    whatever HEAD is *after* --patches, if given, is
#                    applied). Required in practice with --patches, since
#                    the post-patch commit hash is otherwise the only name
#                    and is not the pin's hash.
#
# --patches checks out the pin FRESH (git checkout -f, detached HEAD -- any
# uncommitted local state in $SRC is discarded, same as this script already
# implied by building "the pin") and applies the series with `git am`, so
# the build is reproducible from the pin plus the patches alone, with no
# dependency on a local `aug` branch existing.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

PATCH_DIR=""
PREFIX_NAME=""
POSITIONAL=()
for arg in "$@"; do
    case "$arg" in
        --patches)
            PATCH_DIR="$REPO_ROOT/patches"
            ;;
        --patches=*)
            PATCH_DIR="${arg#--patches=}"
            ;;
        --prefix-name=*)
            PREFIX_NAME="${arg#--prefix-name=}"
            ;;
        *)
            POSITIONAL+=("$arg")
            ;;
    esac
done
set -- "${POSITIONAL[@]+"${POSITIONAL[@]}"}"

SRC="${1:-$HOME/Work/dev/ngspice_augmentation/ngspice-src}"
JOBS="${2:-$(nproc)}"

if [ ! -d "$SRC" ]; then
    echo "error: source dir '$SRC' not found." >&2
    echo "Clone it first (see README.md step 1) -- it is gitignored and not" >&2
    echo "part of this repo." >&2
    exit 1
fi

cd "$SRC"

PIN_FILE="$REPO_ROOT/ngspice-base.txt"
PINNED_COMMIT=$(awk -F'= ?' '/^commit *=/{print $2; exit}' "$PIN_FILE")

if [ -n "$PATCH_DIR" ]; then
    if [ ! -f "$PATCH_DIR/series" ]; then
        echo "error: $PATCH_DIR/series not found (need a patch series to apply)." >&2
        exit 1
    fi
    if [ -z "$PINNED_COMMIT" ]; then
        echo "error: no pinned commit in $PIN_FILE; --patches needs a pin to check out." >&2
        exit 1
    fi
    echo "=== checking out the pin $PINNED_COMMIT (fresh, detached) ==="
    git am --abort >/dev/null 2>&1 || true
    git checkout -f "$PINNED_COMMIT"
    echo "=== applying patch series from $PATCH_DIR ==="
    while IFS= read -r p; do
        [ -z "$p" ] && continue
        case "$p" in \#*) continue ;; esac
        echo "  git am $PATCH_DIR/$p"
        # committer identity is set here so the recipe works on a machine with no git
        # identity (the development workstation, AMD Ryzen AI Max+ 395 / Linux, has none); the patch authorship inside each file is kept
        git -c user.name="ngspice_augmentation build" -c user.email="build@localhost" am "$PATCH_DIR/$p"
    done < "$PATCH_DIR/series"
    echo
fi

HEAD_COMMIT=$(git rev-parse HEAD)

if [ -z "$PATCH_DIR" ] && [ -n "$PINNED_COMMIT" ] && [ "$HEAD_COMMIT" != "$PINNED_COMMIT" ]; then
    echo "warning: $SRC is at $HEAD_COMMIT, not the pinned $PINNED_COMMIT" >&2
    echo "         (ngspice-base.txt). Continuing anyway -- pass a clean" >&2
    echo "         checkout of the pin if you want a reproducible build." >&2
fi

if [ -n "$PREFIX_NAME" ]; then
    PREFIX="$HOME/Work/apps/ngspice_augmentation/$PREFIX_NAME"
else
    SHORT=$(git rev-parse --short=12 HEAD)
    PREFIX="$HOME/Work/apps/ngspice_augmentation/$SHORT"
fi

echo "=== ngspice_augmentation build ==="
echo "source:  $SRC"
echo "commit:  $HEAD_COMMIT$([ -n "$PATCH_DIR" ] && echo "  (pin $PINNED_COMMIT + $(wc -l < "$PATCH_DIR/series" | tr -d ' ') patch(es) from $PATCH_DIR)")"
echo "prefix:  $PREFIX"
echo "jobs:    $JOBS"
echo

# --- autogen (regenerate configure/Makefile.in from configure.ac/Makefile.am) ---
./autogen.sh

# --- out-of-tree build dir ---
mkdir -p release
cd release

# --- configure: matches the system /usr/local/bin/ngspice (45.2) build ---
# Feature flags found present on the installed 45.2 (via `strings` on the
# ELF: KLU/XSPICE/OSDI banners, "OpenMP multithreading ... enabled", and
# ldd showing libgomp/libreadline/libX* but no libngspice.so, i.e. the
# default non-shared build, not --enable-cider (no cider strings found),
# not --adms). PSS is default-on in pre-master-48 (Enhancement-117; see
# configure.ac's `--disable-pss` help text) so no flag is needed for it.
../configure \
    --with-x \
    --enable-klu \
    --enable-xspice \
    --enable-openmp \
    --prefix="$PREFIX" \
    CFLAGS="-m64 -O2" \
    LDFLAGS="-m64 -s"

echo
echo "=== configure feature summary ==="
grep -iE "KLU solver enabled|PSS .* enabled|XSPICE features included|OSDI features included|OpenMP feature enabled|checking for readline" config.log || true
echo

make -j"$JOBS"
make install

PREFIX_BIN="$PREFIX/bin"

# --- wrapper mirroring the system wrapper's OMP_WAIT_POLICY=passive ---
cat > "$PREFIX_BIN/ngspice-aug" <<WRAP
#!/bin/sh
# ngspice-aug wrapper -- private pre-master-48 fork build (card #455 / T3.1).
# Mirrors the OMP_WAIT_POLICY=passive behaviour of the system wrapper at
# /usr/local/bin/ngspice (kanban #89): idle OpenMP threads busy-spin by
# default, which burns CPU on multi-hour parameter sweeps for no wall-clock
# benefit.
# REVERT: delete this file, or the whole prefix ($PREFIX) to remove this
# build entirely. Nothing here touches /usr/local.
HERE=\$(CDPATH= cd -- "\$(dirname -- "\$0")" && pwd)
[ -n "\$OMP_WAIT_POLICY" ] || OMP_WAIT_POLICY=passive
export OMP_WAIT_POLICY
exec "\$HERE/ngspice" "\$@"
WRAP
chmod 755 "$PREFIX_BIN/ngspice-aug"

echo
echo "=== done ==="
echo "prefix: $PREFIX"
"$PREFIX_BIN/ngspice-aug" -v
