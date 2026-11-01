#!/usr/bin/env bash
# run_all.sh -- runs every handbook example deck/command (docs/analog/
# ngspice-augmentation-handbook.html) against the private-fork build and
# prints the key output lines each one actually produced. Card #507.
#
# Usage:
#   docs/examples/handbook/run_all.sh [NGSPICE_PREFIX]
#
# NGSPICE_PREFIX defaults to the pin+116-patches build (s115) this handbook was
# written against (card #560; s108 -> s115; examples 22-31 from card #534, 32-37 added in #543, 38-44 in #552, 45-50 in #560). Output (raw logs, wrdata traces, tool manifests) is
# written under ~/Work/data/ngspice_augmentation/handbook/run_all/ --
# never into this repository. Runs everything SERIALLY (never more than
# one ngspice process at a time), matching this project's own
# at-most-2-concurrent-ngspice rule with margin.
#
# Exit status: 0 if every step ran to completion (ngspice exit 0 / tool
# exit 0), 1 if any step failed to run. This script does not itself check
# the NUMBERS each step prints match the handbook page -- that is
# tests/test_handbook.py's job.
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
PREFIX="${1:-$HOME/Work/apps/ngspice_augmentation/aaba50ba3f49-s152}"
BIN="$PREFIX/bin/ngspice"
BIN_AUG="$PREFIX/bin/ngspice-aug"
# stock ngspice 45.2 (examples 32tran): $STOCK_NGSPICE, else /usr/local/bin/ngspice, else a copy of
# /usr/local/ngspice-45.2 synced next to the builds as <apps>/stock-45.2 (card #550)
STOCK="${STOCK_NGSPICE:-/usr/local/bin/ngspice}"
[ -x "$STOCK" ] || STOCK="$(dirname "$PREFIX")/stock-45.2/bin/ngspice"
# every tool (loss.py, fra.py, convmeas.py, mc.py, compopt.py, pac_perturbed.py)
# reads $NGSPICE; without this they fell back to the stock ngspice on PATH.
export NGSPICE="$BIN" NGSPICE_AUG="$BIN_AUG"
OUT="$HOME/Work/data/ngspice_augmentation/handbook/run_all"
mkdir -p "$OUT"

FAIL=0
step() {
  local name="$1"; shift
  echo "=== $name ==="
  if "$@"; then
    echo "--- $name: ok ---"
  else
    echo "--- $name: FAILED (exit $?) ---"
    FAIL=1
  fi
}

run_deck() {
  # run_deck NAME DECK [extra args to ngspice]
  local name="$1" deck="$2"; shift 2
  local d="$OUT/$name"; mkdir -p "$d"
  (cd "$d" && "$BIN" -b "$HERE/$deck" "$@") > "$OUT/$name.log" 2>&1
  local rc=$?
  tail -n 25 "$OUT/$name.log"
  return $rc
}

run_deck_aug() {
  local name="$1" deck="$2"; shift 2
  local d="$OUT/$name"; mkdir -p "$d"
  (cd "$d" && "$BIN_AUG" -b "$HERE/$deck" "$@") > "$OUT/$name.log" 2>&1
  local rc=$?
  tail -n 25 "$OUT/$name.log"
  return $rc
}

step "01 stb LDO loop"        run_deck  01 01_stb_ldo_loop.cir
step "02 stb diffamp DM/CM"   run_deck  02 02_stb_diffamp_dmcm.cir
step "03 stb device mode"     run_deck  03 03_stb_device_mode.cir
step "04 stb sens="           run_deck  04 04_stb_sens.cir
step "05 stb instant"         run_deck  05 05_stb_instant.cir
step "06 pss method=shoot"    run_deck  06 06_pss_shoot.cir
step "07 pac solver=td"       run_deck  07 07_pac_td.cir
step "08 pstb"                run_deck  08 08_pstb.cir
step "09 pxf"                 run_deck  09 09_pxf.cir
step "10 saltation ideal sw"  run_deck  10 10_saltation_ideal_sw.cir

echo "=== 11 pac_perturbed ==="
mkdir -p "$OUT/11run"
python3 "$ROOT/tools/pac_perturbed.py" sweep \
    --deck "$HERE/11_pac_perturbed_deck.cir" \
    --fs 100e3 --inject source:Vin --outputs out \
    --freqs 10000,30000 --amplitude 0.01 \
    --backend pss --stabtime 0 \
    --ngspice "$BIN_AUG" --max-procs 2 \
    --out "$OUT/11run" > "$OUT/11.log" 2>&1
rc11=$?
tail -n 15 "$OUT/11.log"
[ $rc11 -eq 0 ] || { echo "--- 11: FAILED ---"; FAIL=1; }

echo "=== 12 fra --adaptive-settle ==="
mkdir -p "$OUT/12run"
python3 "$ROOT/tools/fra.py" "$HERE/12_fra_adaptive_settle_deck.cir" \
    --mode transfer --source vin --output "v(vx)" \
    --freqs 10000 --amplitude 1e-3 --adaptive-settle \
    --settle-period 10e-6 --residual-bound 0.01 --settle-cap 3e-4 \
    --settle-periods 6 --measure-periods 4 --points-per-period 1000 \
    --out "$OUT/12run" > "$OUT/12.log" 2>&1
rc12=$?
tail -n 15 "$OUT/12.log"
[ $rc12 -eq 0 ] || { echo "--- 12: FAILED ---"; FAIL=1; }

echo "=== 13 convmeas ==="
mkdir -p "$OUT/13run"
python3 "$ROOT/tools/convmeas.py" "$HERE/13_convmeas_deck.cir" \
    --vin Vin --vout-node vout --rload-ohms 1e9 \
    --period 198.69e-6 --periods 1 --tstop 1.01e-3 \
    --load-step-time 10e-6 --out "$OUT/13run" > "$OUT/13.log" 2>&1
rc13=$?
tail -n 15 "$OUT/13.log"
[ $rc13 -eq 0 ] || { echo "--- 13: FAILED ---"; FAIL=1; }

echo "=== 14 loss.py ==="
mkdir -p "$OUT/14run"
python3 "$ROOT/tools/loss.py" "$HERE/14_loss_deck.cir" \
    --vin Vin --load Rload --period 2e-6 --periods 5 --tstop 1.5e-3 \
    --dcr Rdcr --esr Resr --diode Ddiode --switch Bhs \
    --out "$OUT/14run" > "$OUT/14.log" 2>&1
rc14=$?
tail -n 15 "$OUT/14.log"
[ $rc14 -eq 0 ] || { echo "--- 14: FAILED ---"; FAIL=1; }

echo "=== 15 mc.py ==="
mkdir -p "$OUT/15run"
(cd "$HERE" && python3 "$ROOT/tools/mc.py" run 15_mc_spec.json --out "$OUT/15run") \
    > "$OUT/15.log" 2>&1
rc15=$?
tail -n 20 "$OUT/15.log"
[ $rc15 -eq 0 ] || { echo "--- 15: FAILED ---"; FAIL=1; }

echo "=== 16 yieldest.py ==="
(cd "$HERE" && NGSPICE="$BIN" python3 16_yield_run.py) > "$OUT/16.log" 2>&1
rc16=$?
tail -n 15 "$OUT/16.log"
[ $rc16 -eq 0 ] || { echo "--- 16: FAILED ---"; FAIL=1; }

echo "=== 17 compopt.py ==="
mkdir -p "$OUT/17run"
python3 "$ROOT/tools/compopt.py" "$HERE/17_compopt_deck.cir" \
    --lgprobe-inst XP \
    --param cc=1e-13:3e-11:log --param c1=1e-14:1e-11:log \
    --start cc=1e-12 --start c1=0.2e-12 \
    --target fc=7.6e6 --target pm=65 \
    --sweep "dec 30 10 1e9" --restarts 2 --max-iter 80 --verify stb \
    --ngspice "$BIN_AUG" --out "$OUT/17run" > "$OUT/17.log" 2>&1
rc17=$?
tail -n 15 "$OUT/17.log"
[ $rc17 -eq 0 ] || { echo "--- 17: FAILED ---"; FAIL=1; }

echo "=== 18 checkpoint save/restore ==="
mkdir -p "$OUT/18run"
(
  cd "$OUT/18run"
  "$BIN_AUG" -b "$HERE/18_checkpoint_uninterrupted.cir"
  echo "--- save ---"
  "$BIN_AUG" -b "$HERE/18_checkpoint_save.cir"
  echo "--- restore ---"
  "$BIN_AUG" -b "$HERE/18_checkpoint_restore.cir"
) > "$OUT/18.log" 2>&1
rc18=$?
grep -E "vx_at_|Checkpoint written|Restored checkpoint" "$OUT/18.log" || true
[ $rc18 -eq 0 ] || { echo "--- 18: FAILED ---"; FAIL=1; }

echo "=== 19 once() vs without ==="
mkdir -p "$OUT/19run"
(
  cd "$OUT/19run"
  echo "--- without once() ---"
  "$BIN_AUG" -b "$HERE/19_mc_once_without.cir"
  echo "--- with once() ---"
  "$BIN_AUG" -b "$HERE/19_mc_once_with.cir"
) > "$OUT/19.log" 2>&1
rc19=$?
grep -E "resistance" "$OUT/19.log" || true
[ $rc19 -eq 0 ] || { echo "--- 19: FAILED ---"; FAIL=1; }

echo "=== 20 pss driven-detection red-proof (expected REFUSAL) ==="
mkdir -p "$OUT/20run"
(cd "$OUT/20run" && "$BIN" -b "$HERE/20_pss_driven_detection_badclock.cir") \
    > "$OUT/20.log" 2>&1
# This step is EXPECTED to exit nonzero (a red-proof) -- pass iff the
# specific documented refusal message appears, not on exit code.
if grep -q "no V/I source with a time function was found" "$OUT/20.log"; then
  echo "--- 20: ok (refused as documented) ---"
else
  echo "--- 20: FAILED (expected refusal message not found) ---"
  cat "$OUT/20.log"
  FAIL=1
fi

echo "=== 21 saltation: switch node without a capacitor (card #510, fixed) ==="
mkdir -p "$OUT/21run"
(cd "$OUT/21run" && "$BIN" -b "$HERE/21_saltation_direct_control_zero.cir") \
    > "$OUT/21.log" 2>&1
rc21=$?
grep -A3 "Index" "$OUT/21.log" || true
[ $rc21 -eq 0 ] || { echo "--- 21: FAILED ---"; FAIL=1; }

echo "=== 22 pss closed-loop ideal-switch buck (card #498) ==="
run_deck 22 22_pss_closed_loop_sw.cir > "$OUT/22.tail" 2>&1; rc22=$?
grep -E "CONVERGED|FLOQUET" "$OUT/22.log" || true
[ $rc22 -eq 0 ] || { echo "--- 22: FAILED ---"; FAIL=1; }

echo "=== 22off red-proof: same deck with PSS_NO_EVENT_LOCATE=1 (expected NON-convergence) ==="
mkdir -p "$OUT/22off"
(cd "$OUT/22off" && PSS_NO_EVENT_LOCATE=1 "$BIN" -b "$HERE/22_pss_closed_loop_sw.cir") > "$OUT/22off.log" 2>&1
if grep -q "FAILED to converge within" "$OUT/22off.log"; then
  grep "FAILED to converge within" "$OUT/22off.log"
  echo "--- 22off: ok (refused as documented) ---"
else
  echo "--- 22off: FAILED (expected non-convergence not seen) ---"
  FAIL=1
fi

echo "=== 23 pstb closed-loop ideal-switch buck + independent pac_perturbed --loop-gain ==="
run_deck 23 23_pstb_closed_loop_sw.cir > "$OUT/23.tail" 2>&1; rc23=$?
grep -E "CONVERGED" "$OUT/23.log" || true
[ $rc23 -eq 0 ] || { echo "--- 23: FAILED ---"; FAIL=1; }
mkdir -p "$OUT/23ref"
python3 "$ROOT/tools/pac_perturbed.py" sweep \
    --deck "$HERE/23_pstb_closed_loop_sw_ref_deck.cir" --fs 500e3 --loop-gain --probe-nodes fbtop:vout \
    --freqs 30000,100000 --amplitude 0.01 --backend tran --tran-stabtime-mult 30 \
    --ngspice "$BIN_AUG" --max-procs 2 --out "$OUT/23ref" > "$OUT/23ref.log" 2>&1
rc23r=$?
tail -n 4 "$OUT/23ref.log"
[ $rc23r -eq 0 ] || { echo "--- 23ref: FAILED ---"; FAIL=1; }

step "24 pnoise switched RC"   run_deck  24 24_pnoise_switched_rc.cir
step "25 pnoise track-and-hold" run_deck  25 25_pnoise_track_hold.cir
step "26 pss R3d dead-time"    run_deck  26 26_pss_r3d_deadtime.cir
step "27 pss R4d dead-time"    run_deck  27 27_pss_r4d_deadtime.cir
step "28 pstb R3d dead-time"   run_deck  28 28_pstb_r3d_deadtime.cir

step "30 pac free-running R5"  run_deck  30 30_pac_freerun_r5.cir
echo "=== 30t1 same deck at num_threads=1 (card #517: result must not depend on the solver's thread/step path) ==="
mkdir -p "$OUT/30t1"
sed 's/^set filetype=ascii/set filetype=ascii\nset num_threads=1/' "$HERE/30_pac_freerun_r5.cir" > "$OUT/30t1/deck.cir"
(cd "$OUT/30t1" && "$BIN" -b deck.cir) > "$OUT/30t1.log" 2>&1
rc30t=$?
python3 - "$OUT/30/pac_r5.txt" "$OUT/30t1/pac_r5.txt" <<'PYEOF' || FAIL=1
import cmath, math, sys
def rd(p):
    rows = [ln.split() for ln in open(p) if ln.strip() and not ln.strip().startswith("freq")]
    return [(float(r[0]), complex(float(r[1]), float(r[2]))) for r in rows]
a, b = rd(sys.argv[1]), rd(sys.argv[2])
for (f, x), (_, y) in zip(a, b):
    print(f"30t1 f={f:.0f} Hz: 2 threads {20*math.log10(abs(x)):.4f} dB/{math.degrees(cmath.phase(x)):.3f} deg, "
          f"1 thread {20*math.log10(abs(y)):.4f} dB/{math.degrees(cmath.phase(y)):.3f} deg, "
          f"delta {20*math.log10(abs(y)/abs(x)):+.5f} dB/{math.degrees(cmath.phase(y/x)):+.5f} deg")
PYEOF
[ $rc30t -eq 0 ] || { echo "--- 30t1: FAILED ---"; FAIL=1; }

echo "=== 29 speed switches (rounds 13-15; interleaved, same binary) ==="
python3 "$HERE/29_speed_switches.py" "$BIN" "$OUT/29run" 7 > "$OUT/29.log" 2>&1
rc29=$?
cat "$OUT/29.log"
[ $rc29 -eq 0 ] || { echo "--- 29: FAILED ---"; FAIL=1; }

echo "=== 31 campaign/mc thread rule (tools/ompthreads.py, card #527) ==="
python3 - "$ROOT" > "$OUT/31.log" 2>&1 <<'PYEOF'
import sys
sys.path.insert(0, sys.argv[1] + "/tools")
import ompthreads as o
for cores in (32, 8):
    print(f"THREADRULE cores={cores}: " + "  ".join(f"N={n}:{o.auto_threads(n, cores)}" for n in (1, 2, 4, 6, 8, 12, 16, 32)))
print("THREADRULE settings: " + "  ".join(f"{x!r}->{o.resolve_setting(x, 12, 32)}" for x in ("auto", 3, "default", 0)))
PYEOF
rc31=$?
cat "$OUT/31.log"
[ $rc31 -eq 0 ] || { echo "--- 31: FAILED ---"; FAIL=1; }

# ------------------------------------------------------------------------------------------- card #543 (s96)
echo "=== 32 R7 transistor-level (BSIM4) buck: .pss / .pstb / .pac / .pxf (decks from 32_r7_decks.py) ==="
mkdir -p "$OUT/32run"
python3 "$HERE/32_r7_decks.py" "$OUT/32run" || { echo "--- 32: deck generation FAILED ---"; FAIL=1; }
for n in r7_pss r7_pss_reltol3e-4 r7_pss_reltol1e-3 r7_pstb r7_pstb_reltol1e-3 r7_pac r7_pxf_fwd r7_pxf_adj; do
  (cd "$OUT/32run" && "$BIN" -b "$n.cir") > "$OUT/32_$n.log" 2>&1
  rc=$?
  tr '\r' '\n' < "$OUT/32_$n.log" | grep -a -E "PSS-SHOOT: CONVERGED|FAILED to converge|^pstb: |PAC-TD: generic|PXF-TD: generic|PSTB: generic|PXF-TD: adjoint" | head -4
  [ $rc -eq 0 ] || { echo "--- 32 $n: FAILED (exit $rc) ---"; FAIL=1; }
done
echo "=== 32old R7 .pstb on the build BEFORE the charge-history fix (s91, 0092 not applied): silently wrong ==="
OLDBIN="$(dirname "$PREFIX")/aaba50ba3f49-s91/bin/ngspice"
mkdir -p "$OUT/32old"
if [ -x "$OLDBIN" ]; then
  (cd "$OUT/32old" && "$OLDBIN" -b "$OUT/32run/r7_pstb.cir") > "$OUT/32_r7_pstb_s91.log" 2>&1
  tr '\r' '\n' < "$OUT/32_r7_pstb_s91.log" | grep -a -E "^pstb: |PAC-TD: generic" | head -3
else
  echo "s91 build not found at $OLDBIN: red-proof skipped" | tee "$OUT/32_r7_pstb_s91.log"
fi
echo "=== 32tran R7 ramp-PULSE overlap trap: stock ngspice and the fork, old ramp line vs the fixed one ==="
mkdir -p "$OUT/32tran"
python3 "$HERE/32_r7_decks.py" --tran-trap "$OUT/32tran" || { echo "--- 32tran: deck generation FAILED ---"; FAIL=1; }
for v in old fixed; do
  for who in stock fork; do
    if [ $who = stock ]; then B="$STOCK"; else B="$BIN"; fi
    (cd "$OUT/32tran" && "$B" -b "r7_tran_$v.cir") > "$OUT/32tran_${v}_$who.log" 2>&1
    echo "32tran $v $who: $(tr '\r' '\n' < "$OUT/32tran_${v}_$who.log" | grep -a -E 'Timestep too small|^vo ' | head -1)"
  done
done

echo "=== 33 pac_perturbed on the ideal-switch deck of example 10: exact frequency, defaults; and the pre-#536 settings (expected REFUSED) ==="
mkdir -p "$OUT/33run" "$OUT/33old"
python3 "$ROOT/tools/pac_perturbed.py" sweep \
    --deck "$HERE/10_saltation_ideal_sw.cir" --fs 500e3 --inject source:Vctrl --outputs vout \
    --freqs 1000 --exact-freq --backend tran --stabtime 100e-6 --tran-stabtime-mult 30 --residual-bound 0.1 \
    --ngspice "$BIN_AUG" --max-procs 2 --out "$OUT/33run" > "$OUT/33.log" 2>&1
rc33=$?
tail -n 6 "$OUT/33.log"
[ $rc33 -eq 0 ] || { echo "--- 33: FAILED (exit $rc33) ---"; FAIL=1; }
python3 "$ROOT/tools/pac_perturbed.py" sweep \
    --deck "$HERE/10_saltation_ideal_sw.cir" --fs 500e3 --inject source:Vctrl --outputs vout \
    --freqs 1000 --backend tran --stabtime 100e-6 --tran-stabtime-mult 30 --residual-bound 0.1 \
    --qmax 500 --amplitude 0.01 --points-per-cycle 200 --reltol 1e-3 --no-halve \
    --ngspice "$BIN_AUG" --max-procs 2 --out "$OUT/33old" > "$OUT/33old.log" 2>&1
tail -n 4 "$OUT/33old.log"     # exit status is not checked: a REFUSED point is the documented outcome

echo "=== 34 pac_perturbed: requested vs measured frequency, and the PULSE-overlap note (example 11's deck, ramp clock) ==="
mkdir -p "$OUT/34run" "$OUT/34runfix"
python3 "$ROOT/tools/pac_perturbed.py" sweep \
    --deck "$HERE/34_pac_perturbed_overlap_deck.cir" --fs 100e3 --inject source:Vin --outputs out \
    --freqs 1000 --amplitude 0.01 --backend pss --stabtime 0 \
    --ngspice "$BIN_AUG" --max-procs 2 --out "$OUT/34run" > "$OUT/34.log" 2>&1
tail -n 6 "$OUT/34.log"
PACP_FIX_PULSE=1 python3 "$ROOT/tools/pac_perturbed.py" sweep \
    --deck "$HERE/34_pac_perturbed_overlap_deck.cir" --fs 100e3 --inject source:Vin --outputs out \
    --freqs 1000 --exact-freq --amplitude 0.01 --backend pss --stabtime 0 \
    --ngspice "$BIN_AUG" --max-procs 2 --out "$OUT/34runfix" > "$OUT/34fix.log" 2>&1
tail -n 6 "$OUT/34fix.log"

echo "=== 35 free-running stability verdict (card #540): R5 + 2 extracted RC trees (10 states), matrix-free vs dense, normal vs forced-noise pass ==="
mkdir -p "$OUT/35run"
python3 "$HERE/35_freerun_verdict.py" "$BIN" "$OUT/35run" > "$OUT/35.log" 2>&1
rc35=$?
cat "$OUT/35.log"
[ $rc35 -eq 0 ] || { echo "--- 35: FAILED ---"; FAIL=1; }
echo "=== 35b example 30's own deck, dense engine (NGSPICE_PSS_KRYLOV=0) ==="
mkdir -p "$OUT/35d"
(cd "$OUT/35d" && NGSPICE_PSS_KRYLOV=0 "$BIN" -b "$HERE/30_pac_freerun_r5.cir") > "$OUT/35d.log" 2>&1
rc35d=$?
tr '\r' '\n' < "$OUT/35d.log" | grep -a -E "PSS-FREERUN: CONVERGED"
[ $rc35d -eq 0 ] || { echo "--- 35d: FAILED ---"; FAIL=1; }

echo "=== 36 PEX-size .pss (490 states): SPARSE vs KLU, dense monodromy vs matrix-free (n=2, rotated) ==="
mkdir -p "$OUT/36run"
python3 "$HERE/36_pex_krylov.py" "$BIN" "$OUT/36run" 2 > "$OUT/36.log" 2>&1
rc36=$?
cat "$OUT/36.log"
[ $rc36 -eq 0 ] || { echo "--- 36: FAILED ---"; FAIL=1; }

echo "=== 37 .pss on a PWL clock written 'PWL(...) r=0' (patch 0096, card #541); the same deck on s91 is refused ==="
mkdir -p "$OUT/37run" "$OUT/37old"
(cd "$OUT/37run" && "$BIN" -b "$HERE/37_pss_pwl_clock.cir") > "$OUT/37.log" 2>&1
rc37=$?
tr '\r' '\n' < "$OUT/37.log" | grep -a -E "PSS-SHOOT: CONVERGED|^vmax|^vmin"
[ $rc37 -eq 0 ] || { echo "--- 37: FAILED ---"; FAIL=1; }
if [ -x "$OLDBIN" ]; then
  (cd "$OUT/37old" && "$OLDBIN" -b "$HERE/37_pss_pwl_clock.cir") > "$OUT/37_s91.log" 2>&1
  grep -a "no V/I source with a time function" "$OUT/37_s91.log" | head -1
fi

# ------------------------------------------------------------------------------------------- card #552 (s108)
echo "=== 38 .pnoise: noise at the operating point of each step (diode / BJT, LTI bias and 100 kHz LO bias; s96 and PAC_NOISE_STALE_STATE=1 as the before) ==="
mkdir -p "$OUT/38run"
python3 "$HERE/38_pnoise_opoint.py" "$BIN" "$OUT/38run" > "$OUT/38.log" 2>&1
rc38=$?
cat "$OUT/38.log" | cut -c1-400
[ $rc38 -eq 0 ] || { echo "--- 38: FAILED ---"; FAIL=1; }

echo "=== 39 .pnoise of a MOSFET track-and-hold (MOS1 and BSIM4 switch) against the closed-form kT/C; red proof PAC_NATIVE_HISTORY=1 ==="
for k in mos1 bsim4; do
  for v in normal native; do
    mkdir -p "$OUT/39_${k}_${v}"
    if [ $v = native ]; then E="PAC_NATIVE_HISTORY=1"; else E="X39=1"; fi
    (cd "$OUT/39_${k}_${v}" && env $E "$BIN" -b "$HERE/39_pnoise_th_${k}.cir") > "$OUT/39_${k}_${v}.log" 2>&1
    rc=$?
    echo "39 $k $v:"; cat "$OUT/39_${k}_${v}/pnoise_th_${k}.txt"
    [ $rc -eq 0 ] || { echo "--- 39 $k $v: FAILED (exit $rc) ---"; FAIL=1; }
  done
done

echo "=== 40 R7 after cards #546/#549: hidden-state continuation (orbit vs settled tran), pacreltol=, method=gear ==="
mkdir -p "$OUT/40run"
python3 "$HERE/40_r7_hidden_state.py" "$BIN" "$OUT/40run" > "$OUT/40.log" 2>&1
rc40=$?
cut -c1-420 "$OUT/40.log"
[ $rc40 -eq 0 ] || { echo "--- 40: FAILED ---"; FAIL=1; }

echo "=== 41 R8 free-running BSIM4 hysteretic buck: .pss trigger=cx (hidden state, and PSS_NO_HIDDEN=1), .pac ==="
mkdir -p "$OUT/41run"
python3 "$HERE/41_r8_freerun.py" "$BIN" "$OUT/41run" > "$OUT/41.log" 2>&1
rc41=$?
cut -c1-420 "$OUT/41.log"
[ $rc41 -eq 0 ] || { echo "--- 41: FAILED ---"; FAIL=1; }

echo "=== 42 ideal-switch events with device charges (R6o with a BSIM4 follower in the control path); red proofs ==="
mkdir -p "$OUT/42run"
python3 "$HERE/42_salt_charge.py" "$BIN" "$OUT/42run" > "$OUT/42.log" 2>&1
rc42=$?
cut -c1-330 "$OUT/42.log"
[ $rc42 -eq 0 ] || { echo "--- 42: FAILED ---"; FAIL=1; }

echo "=== 43 XSPICE code models under .pss and td .pac (s108; s105 as the before; red proof PSS_XS_NO_STATE=1; refusals) ==="
mkdir -p "$OUT/43run"
python3 "$HERE/43_xspice.py" "$BIN" "$OUT/43run" > "$OUT/43.log" 2>&1
rc43=$?
cut -c1-330 "$OUT/43.log"
[ $rc43 -eq 0 ] || { echo "--- 43: FAILED ---"; FAIL=1; }

echo "=== 44 R9 digitally controlled buck: .pss (T = 1.2 ms) vs settled tran, .pstb refused, 7-bit limit-cycling variant refused (about 10 minutes) ==="
mkdir -p "$OUT/44run"
python3 "$HERE/44_r9.py" "$BIN" "$OUT/44run" > "$OUT/44.log" 2>&1
rc44=$?
cut -c1-420 "$OUT/44.log"
[ $rc44 -eq 0 ] || { echo "--- 44: FAILED ---"; FAIL=1; }

# ------------------------------------------------------------------------------------------- card #560 (s115)
echo "=== 45 .pnoise folded flicker: chopper vs the closed form; red proof PNOISE_FLICKER_OLD=1 ==="
mkdir -p "$OUT/45run"
python3 "$HERE/45_pnoise_fold.py" "$BIN" "$OUT/45run" > "$OUT/45.log" 2>&1
rc45=$?
cut -c1-330 "$OUT/45.log"
[ $rc45 -eq 0 ] || { echo "--- 45: FAILED ---"; FAIL=1; }

echo "=== 46 PSS memory: peak RSS on the PEX decks, s108 vs s115 vs PSS_KEEP_RUN_DATA=1; plots listed after .pss ==="
mkdir -p "$OUT/46run"
python3 "$HERE/46_pss_memory.py" "$BIN" "$OUT/46run" > "$OUT/46.log" 2>&1
rc46=$?
cut -c1-330 "$OUT/46.log"
[ $rc46 -eq 0 ] || { echo "--- 46: FAILED ---"; FAIL=1; }

echo "=== 47 speed rounds 16-17: automatic KLU, PSS_FORK_COLS, PAC_MULTIRHS, PSS_SEED_STAGE1, Arnoldi stop (interleaved, n=3) ==="
mkdir -p "$OUT/47run"
python3 "$HERE/47_speed_r16_r17.py" "$BIN" "$OUT/47run" 3 > "$OUT/47.log" 2>&1
rc47=$?
cut -c1-420 "$OUT/47.log"
[ $rc47 -eq 0 ] || { echo "--- 47: FAILED ---"; FAIL=1; }

echo "=== 48 free-running R5: .pxf adjoint=1 = forward, .pnoise in band, band guard ==="
mkdir -p "$OUT/48run"
python3 "$HERE/48_freerun_noise.py" "$BIN" "$OUT/48run" > "$OUT/48.log" 2>&1
rc48=$?
cut -c1-420 "$OUT/48.log"
[ $rc48 -eq 0 ] || { echo "--- 48: FAILED ---"; FAIL=1; }

echo "=== 49 OMP_NUM_THREADS does not limit the OpenMP threads; set num_threads does ==="
mkdir -p "$OUT/49run"
python3 "$HERE/49_omp_threads.py" "$BIN" "$OUT/49run" > "$OUT/49.log" 2>&1
rc49=$?
cut -c1-330 "$OUT/49.log"
[ $rc49 -eq 0 ] || { echo "--- 49: FAILED ---"; FAIL=1; }

echo "=== 50 td .pac through an s_xfer compensator (R3 buck) vs the cited FRA; linear s_xfer vs .ac; refusals ==="
mkdir -p "$OUT/50run"
python3 "$HERE/50_xspice_pac_through.py" "$BIN" "$OUT/50run" > "$OUT/50.log" 2>&1
rc50=$?
cut -c1-420 "$OUT/50.log"
[ $rc50 -eq 0 ] || { echo "--- 50: FAILED ---"; FAIL=1; }

echo
if [ "$FAIL" -eq 0 ]; then
  echo "run_all.sh: all steps completed"
else
  echo "run_all.sh: at least one step FAILED (see above)"
fi
exit $FAIL
