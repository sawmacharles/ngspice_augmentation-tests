# ngspice_augmentation-tests

Test suites, tools, reference decks and documentation for ngspice_augmentation,
a fork of the ngspice circuit simulator. The simulator source and the Linux
release binaries are at https://github.com/sawmacharles/ngspice_augmentation-src.
This repository is a snapshot taken at release s152.

## Layout

- `tests/`: test suites (`test_*.py`), their fixture decks (`tests/fixtures/`) and bug reproductions (`tests/bugs/`)
- `tools/`: Python drivers the tests and docs use (FRA, stb, PAC references, envelope, campaign, loss, Monte Carlo, ...)
- `refs/`: reference converter and amplifier decks (R1..R9) and their README
- `bench/decks/`: generated benchmark decks (synthetic)
- `stb/`: the loop-gain probe subcircuits (`lgprobe.lib`, `diffprobe.lib`)
- `docs/`: one page per feature, worked examples (`docs/examples/`) and the handbook examples (`docs/examples/handbook/`, run by `run_all.sh`)
- `buildrecipe/`: how the fork is built (`build.sh`, README)

## Running the tests

1. Download the s152 release archive from the src repository's releases and
   extract it. The binary is `<extracted>/ngspice-augmentation-s152/bin/ngspice-aug`.
2. Point the tests at it and run a suite from the repository root:

       export NGSPICE=<extracted>/ngspice-augmentation-s152/bin/ngspice-aug
       python3 tests/test_stb.py

Each suite prints one line per check and a `N/M pass` summary.

Exit status (most suites): 0 = every check passed, 1 = a check failed,
2 = could not run (binary, fixture or reference missing). Exceptions:
`tests/test_accel.py` exits 1 for could-not-run as well, and
`tests/test_tranbench.py` exits 0 or 1 only. Each file's docstring states its own
convention, its references and its run time.

### Environment variables

Common:

- `NGSPICE`: the binary under test. When unset, most suites use `ngspice` from
  `PATH`; some default to `~/Work/apps/ngspice_augmentation/<build>/bin/ngspice`
  (the development layout), so set `NGSPICE` explicitly.
- `NGSPICE_BASE`: an earlier build, for the checks that compare against it
  (`test_digital_pac.py`, `test_envelope.py`, `test_pss_digclock.py`,
  `test_pstb_probe.py`, `test_pac_dcm.py`, `test_pss_krylov.py`); without it
  those checks are skipped or report could-not-run.
- `NGSPICE_ASAN`: an AddressSanitizer build (see `buildrecipe/README.md`) for
  the ASan checks (`test_pac.py`, `test_pac_charge.py`, `test_pxf.py`,
  `test_pss_klu.py`, `test_pss_krylov.py`, `test_xspice.py`, `test_accel17.py`,
  `test_handbook.py`).
- `NGSPICE_NUMPARAM_NOPATCH`: a build without the numparam patch, for the red
  proof of `test_numparam.py` (could-not-run without it).
- `NGSPICE_AUG`, `NGSPICE_PREV`, `NGSPICE_STOCK`: other builds read by a few
  suites (see their docstrings).
- `PAC_EXACT_PYTHON`: a Python interpreter with numpy and scipy, for
  `tests/test_pac_exact.py` (else `python3` if it imports scipy, else that check
  is could-not-run).
- `OMP_NUM_THREADS`: the suites were run with 2.

Output folders: outputs never go into the repository. By default they go under
`~/Work/data/ngspice_augmentation/<suite folder>/` (several reference decks in
`refs/` and `tests/fixtures/` write there with `write`/`wrdata`, ngspice expands
the `~`). Some suites take an override:
`ENVELOPE_OUT`, `DIG714_OUT`, `TLINE_PSS_OUT`, `PAC_CHARGE_OUT`, `PAC_EXACT_OUT`,
`PSS_KRYLOV_TEST_OUT`, `PNOISE_NQS_DATA`, `PNOISE_BSIM4_DATA`, `MC_BSIM4_DATA`.

### Handbook examples

`tests/test_handbook.py` runs `docs/examples/handbook/run_all.sh` against an
install prefix (`--prefix <dir>`, or the prefix two levels above `$NGSPICE`);
the extracted release folder `ngspice-augmentation-s152` is such a prefix.
`run_all.sh` calls `<prefix>/bin/ngspice` directly in places; in the release
folder that is the bare binary, which finds `spinit` and the XSPICE code models
only through the variables the `ngspice-aug` wrapper sets, so export them first:

    P=<extracted>/ngspice-augmentation-s152
    export SPICE_LIB_DIR=$P/share/ngspice NGSPICE_AUG_PKGLIB=$P/lib/ngspice OMP_WAIT_POLICY=passive
    python3 tests/test_handbook.py --prefix $P

Example 43 checks that the build label contains `sNNN`. Some checks also need
`NGSPICE_ASAN`, earlier release builds next to the prefix, a stock ngspice 45.2
(`STOCK_NGSPICE`, else `/usr/local/bin/ngspice`), or converged FRA curves under
`~/Work/data/ngspice_augmentation/fra/`; without them those checks report
could-not-run.

## Requirements

- Linux, Python 3 (run with 3.12). The test suites use the standard library only.
- numpy and scipy: only `tools/pac_exact.py` (and through it one check of
  `tests/test_pac_exact.py`).
- `systemd-run` and `/usr/bin/time` are used by a few suites when present.

## Not public

Some files refer to material that is not published: `research/*.md`,
`patches/` and `patches-wip/` of the development repository, `aug-*` branch
names, `bench/*.py` generators (only their generated decks are in
`bench/decks/`), the handbook HTML page, and tracker card numbers (`#482`,
`card #719`, ...). The same card numbers appear in the src repository's history.
`tools/s4view/build_r4d.py` needs a schematic editor core that is not
published (its `r4d_devices.lib` is used by `tests/test_floquet_subckt.py`).

## License

BSD 3-Clause, see `LICENSE`.
