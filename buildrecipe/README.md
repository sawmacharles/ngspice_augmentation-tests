# buildrecipe — private pre-master-48 ngspice build (card #455, plan T3.1)

How to reproduce this project's private ngspice fork build from scratch, on
the development workstation (AMD Ryzen AI Max+ 395, Linux), without touching `/usr/local`. This is the recipe; `ngspice-base.txt`
(repo root) is the pin.

## What this is / is not

This is an **unmodified upstream checkout**, built to match the feature set
of the system `/usr/local/bin/ngspice` (45.2), at a newer upstream commit
(`pre-master-48`). It is not yet a fork with our own changes — project ruling
#446 B (2026-09-27) authorizes the private fork; this card (T3.1) is the
first build of the unmodified base, so equivalence with 45.2 can be
characterized before anything of ours is added on top.

**Never installs into or modifies `/usr/local`.** Never uses `sudo`. If a
build dependency is missing, this recipe reports it — it does not install
system packages.

## 1. Clone the source (gitignored, not part of this repo)

```sh
mkdir -p ~/Work/dev/ngspice_augmentation
cd ~/Work/dev/ngspice_augmentation
git clone https://git.code.sf.net/p/ngspice/ngspice ngspice-src
cd ngspice-src
git checkout pre-master-48
```

`ngspice-src/` is listed in this repo's `.gitignore`. The pinned commit is
recorded in `../ngspice-base.txt` at the repo root — check that file's
`commit=` line against `git rev-parse HEAD` before building if you want the
exact commit this was validated against; check out a clean tree at that
commit for a reproducible build (`build.sh` warns, but does not refuse, if
they differ).

## 2. Dependencies

Found present on the development workstation (Ubuntu 24.04, gcc 13.3.0): gcc, g++, make, bison,
flex, autoconf, automake, libtool, pkg-config, libx11-dev, libxaw7-dev,
libreadline-dev, libsuitesparseconfig7 (KLU/AMD/COLAMD runtime libs).

**Found missing** (not installed, per the no-`sudo`/no-package-install
rule): `libsndfile-dev` (`sndfile.h`). This only matters if you bisect or
otherwise check out a source commit older than the point where
`src/frontend/sndprint.c`'s `#include <sndfile.h>` is guarded by
`#if defined(HAVE_LIBSNDFILE) && defined(HAVE_LIBSAMPLERATE)` (guarded as
of the pinned commit aaba50ba — a normal build from the pin is unaffected).
If you ever need libsndfile-backed audio features, that is a build-package
separate decision, not something this recipe installs.

## 3. Build

```sh
buildrecipe/build.sh                    # source: ../ngspice-src, jobs: nproc
buildrecipe/build.sh /path/to/src 16    # override source dir / job count
```

This runs, in an out-of-tree `ngspice-src/release` directory:

```sh
./autogen.sh
mkdir release && cd release
../configure --with-x --enable-klu --enable-xspice --enable-openmp \
    --prefix="$HOME/Work/apps/ngspice_augmentation/<short-commit>" \
    CFLAGS="-m64 -O2" LDFLAGS="-m64 -s"
make -jN
make install
```

then writes a small wrapper (`ngspice-aug`) into the prefix's `bin/`.

### Why these flags (matched to the installed 45.2)

Determined by `strings`/`ldd`/`file` on `/usr/local/ngspice-45.2/bin/ngspice`
(read-only inspection, no install dir has a `config.log` to read from):

| feature | evidence in the installed 45.2 binary | flag used here |
|---|---|---|
| KLU direct solver | `** Compiled with KLU Direct Linear Solver`, `Using KLU as Direct Linear Solver` | `--enable-klu` (also configure.ac default-on) |
| XSPICE | `** XSPICE extensions included`, `xspice_enabled` | `--enable-xspice` (also default-on) |
| OpenMP | `** OpenMP multithreading for BSIM3, BSIM4 enabled`, `ldd` shows `libgomp.so.1` | `--enable-openmp` (also default-on in this configure.ac: `: ${enable_openmp:=yes}`) |
| OSDI | `osdi_enabled`, `OSDI_VERSION_MAJOR` etc. strings | not passed explicitly — default-on, confirmed in `configure.log`'s "OSDI features included" |
| X11 (nutmeg plotting) | `ldd` shows `libXaw`/`libXt`/`libX11`/... | `--with-x` |
| readline | `ldd` shows `libreadline.so.8` | default (readline found via `pkg-config`/headers) |
| CIDER | **absent** — no `cider` strings anywhere in the 45.2 binary | not passed (CIDER is opt-in via `--enable-cider`; the upstream `compile_linux.sh` example script does pass `--enable-cider`, but the installed binary shows it was **not** used to build 45.2, so this recipe does not use it either) |
| ADMS | absent (no `adms` strings) | not passed |
| shared library (`libngspice.so`) | absent from `/usr/local/ngspice-45.2/lib` | default (non-shared); use `compile_linux_shared.sh`'s flags instead if a shared build is ever wanted |
| compiler | `GCC: (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0` embedded in the ELF | same gcc 13.3.0 on the development workstation (`gcc --version`) |
| PSS (periodic steady state) | not present in 45.2 at all (45.2 predates it) | default-on in `pre-master-48` (Enhancement-117: "PSS is now a shipped, default-on analysis"); no flag needed, confirmed via `--disable-pss` help text in `configure.ac` and "PSS (periodic steady state) analysis enabled" in `configure.log` |

`CFLAGS="-m64 -O2" LDFLAGS="-m64 -s"` matches the upstream `compile_linux.sh`
release-build flags (`-s` strips, matching that the installed 45.2 binary is
reported `stripped` by `file`).

## 4. Prefix layout

```
~/Work/apps/ngspice_augmentation/<short-commit>/
  bin/
    ngspice        # the real binary (not a libtool wrapper -- non-shared build)
    ngspice-aug    # our wrapper: sets OMP_WAIT_POLICY=passive unless already set
  lib/ngspice/      # XSPICE code models (analog.cm, digital.cm, ...)
  share/ngspice/    # scripts, help
  share/man/
```

Same layout as `/usr/local/ngspice-45.2` (no `libngspice.so`, same set of
`.cm` code models).

## 5. Use it

```sh
NGSPICE=~/Work/apps/ngspice_augmentation/<short-commit>/bin/ngspice-aug \
    python3 tests/test_stb.py
```

(`tools/stb.py` honours the `NGSPICE` env var; default is `ngspice` on
`$PATH`, which resolves to the system `/usr/local/bin/ngspice` wrapper.)

## 6. Roll back

Delete the prefix:

```sh
rm -rf ~/Work/apps/ngspice_augmentation/<short-commit>
```

and, if you also want the source gone:

```sh
rm -rf ~/Work/dev/ngspice_augmentation/ngspice-src
```

Nothing under `/usr/local` is ever touched by this recipe — there is
nothing there to roll back.

## 7. Building OUR changes (card #456, plan T3.2/T3.3)

Everything above builds the **unmodified** pin. This project's own changes
(so far: the native `.stb` loop-stability analysis) live as a numbered
`git format-patch` series in `patches/` (repo root), with `patches/series`
listing the files in apply order — not as commits in `ngspice-src` itself
(that directory is gitignored, local-only, and gets rebuilt from the pin +
patches on demand).

```sh
buildrecipe/build.sh --patches --prefix-name=<name> [source_dir] [jobs]
```

- `--patches` (optionally `--patches=<dir>` for a different patch set)
  makes the script check out the pinned commit **fresh** (`git checkout -f`,
  detached HEAD — any local edits or branches in `source_dir` are
  discarded) and `git am` the series before configuring/building. This
  means the augmented build is reproducible from the pin plus the patches
  alone, with no dependency on any local branch (e.g. `aug`) existing.
- `--prefix-name=<name>` installs into
  `~/Work/apps/ngspice_augmentation/<name>` instead of the default
  (a short hash of whatever HEAD ends up being after the patches apply,
  which is not the pin's hash and not a meaningful name on its own).
- Everything else (configure flags, prefix layout, the `ngspice-aug`
  wrapper, roll-back) is identical to the unmodified-pin path above.
- The **original, no-`--patches` invocation is unchanged** — it still
  builds the unmodified pin at whatever commit `source_dir` is currently
  checked out to, exactly as card #455 left it. `--patches` and
  `--prefix-name` are both opt-in.

Example (what produced the `aaba50ba3f49-stb1` prefix used in the card #456
report):

```sh
buildrecipe/build.sh --patches --prefix-name=aaba50ba3f49-stb1
```

Validated (card #456): a full `--patches` run against a **fresh clone** of
the pin (not the working `ngspice-src` checkout) applied both patches
cleanly with `git am` and reproduced the same `.stb` results as the
in-place build.

## Equivalence results (card #455, 2026-09-28)

Full numbers, the two upstream commits identified as the cause of the
small (non-build-config) differences against system 45.2, the PSS/PAC
smoke test, and everything skipped, are in the card #455 report (private). Summary:

- `tests/test_stb.py`: 22/22 with both the system `ngspice` and this build
  (`NGSPICE=.../ngspice-aug`), identical numbers.
- `refs/r2_switched_rc.cir` (native `SW` switch, no B-sources, open-loop):
  matches the system build to 7 significant figures on every printed
  measurement.
- `refs/r3_vm_buck.cir`, `r4_pcm_buck.cir`, `r5_hysteretic_buck.cir`
  (closed-loop switching converters built from `tanh()` B-source
  comparators/switches): small (~1e-3 to ~3e-2 relative) differences in
  switching-edge-adjacent measurements (ripple, duty), bisected to two
  genuine upstream changes between 45.2 and `pre-master-48`:
    - `f3871c3ea2c03f55765802588f835baf4177727a` — `dctran.c`, equalises
      the last two transient time steps before a breakpoint (affects any
      circuit with switching-edge breakpoints).
    - `0638aaa1604c3f94221f0a3e489f956ed870043a` — `asrcconv.c`/`niconv.c`,
      a genuine convergence bug fix for B-source (`ASRC`) devices
      (`ASRCconvTest` was comparing old-vs-old current instead of
      old-vs-new, so B-source non-convergence was silently never enforced
      in 45.2).
- `refs/r6_vm_buck_ideal.cir` (native `sw` + XSPICE `d_pwm`, no B-sources):
  steady-state `vout`/ripple match exactly; `duty_meas`/`duty_meas_final`
  (window-averaged over the switching edge) drift slightly, consistent
  with the same breakpoint-timestep-equalising commit above (not
  separately bisected for R6 — see "what was not done").
- `.pss`/`.pac` smoke test on the R2 circuit body: both ran to completion
  (exit 0, well under a second). Not validated numerically — that is a
  later phase's job.

### What was not done

- R6's divergent `duty_meas` was not bisected separately (inferred, not
  confirmed, to share `f3871c3ea2c03f55765802588f835baf4177727a`'s cause).
- The 522-commit range was not exhaustively audited; only the two commits
  that bisection isolated for R3/R5 were examined in detail.
- PSS/PAC numeric correctness was not validated (explicitly out of scope
  for this card).
- No dependency install was attempted for the missing `libsndfile-dev`
  (not needed for the pinned commit; reported, not installed).

## AddressSanitizer build (2026-09-29)

`tests/test_pac.py` checks multi-point `.pac`/`.pstb` runs under ASan. Build it from the same series into its own prefix:

    git clone <private repo>/ngspice-src <dir>   # branch aug = pin + patches/
    cd <dir> && ./autogen.sh && mkdir asan && cd asan
    ../configure --prefix=$HOME/Work/apps/ngspice_augmentation/<pin>-s<N>-asan --with-x --enable-klu --enable-xspice --enable-openmp \
      CFLAGS="-m64 -O1 -g -fsanitize=address -fno-omit-frame-pointer" LDFLAGS="-m64 -fsanitize=address"
    ASAN_OPTIONS=detect_leaks=0 make -j8 && ASAN_OPTIONS=detect_leaks=0 make install

`ASAN_OPTIONS=detect_leaks=0` is required during `make`: a build-time helper (cmpp) leaks and LeakSanitizer would fail the build. Point `ASAN_NGSPICE` at `<prefix>/bin/ngspice` (default: `aaba50ba3f49-s37-asan`).

## Source of the patch series: the public history (card #724, 2026-10-23)

Since 2026-10-23 `patches/` is generated from the public repository
(ngspice_augmentation-src on GitHub), not from a local tree.
Before that, each release was made in a local tree (`ngspice-src-743` and
earlier) whose history was rewritten for publication on 2026-10-17 (#701), and
every release was cherry-picked into the public repository by hand. The two
histories had different hashes, and their comments differed (internal names
scrubbed in the public one).

- **Tree:** `ngspice-src-gh`, a clone of the public repository with the remote
  named `github` (there is no `origin`). The local identity is
  `ngspice-augmentation <noreply@ngspice-augmentation.invalid>`.
- **Regenerate the series** (code commits only; the README / CHANGES-FORK docs
  commits are left out):

  ```sh
  cd ngspice-src-gh; n=0
  for c in $(git rev-list --reverse upstream-base..sNNN -- . ':!README.md' ':!CHANGES-FORK.md'); do
    n=$((n+1)); git format-patch -q -1 $c --start-number $n -o <dir>; done
  (cd <dir> && ls *.patch > series)
  ```

- **New work:** branch from public `main` in a worktree of `ngspice-src-gh`
  (`git -C ngspice-src-gh fetch github && git -C ngspice-src-gh worktree add ../ngspice-src-NNN -b wNNN github/main`).
  Branch from `main`, not from the release tag: `main` carries the README /
  CHANGES-FORK docs commit on top of the tag. The commits are then published as
  they are, with no cherry-pick and the same hashes (#714, `s150`, was branched
  from the tag and still needed one cherry-pick). Before
  pushing, check the added lines for internal names:
  `git diff sPREV wNNN | grep '^+' | grep -n -i -E '<internal-name list>'`.
- **Proof of the switch:** the series regenerated from public `s149` was built
  as `aaba50ba3f49-s149g`. Its `src` tree equals the public tag. `bin/ngspice`
  disassembles identically to the `s149` release build made from the old series
  (sections 8110306 / 390248 / 254608 in both), and every `.cm` code model
  disassembles identically. `ivlng.so` differs only in data offsets, because the
  install path is one character longer. The handbook also ran on the new build
  (results on card #724).
- **Old trees:** `ngspice-src`, `-559`, `-642`, `-644`, `-649` (with the `aug-*`
  branches, never pushed, #640 B), `-682` and `-7xx` are kept untouched. A
  release older than `s149` is rebuilt from the series in this repo's history
  (`git show <commit>:patches/...`).
