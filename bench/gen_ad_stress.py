#!/usr/bin/env python3
"""gen_ad_stress.py -- generate bench/decks/g_ad_stress.cir (ngspice_augmentation
card #486, plan T9.4, round 4): a B-source deck exercising every PT_*/PTF_*
node type and operator src/spicelib/parser/inpptree.c's PTdifferentiate()
handles, argument-swept through PWL ramp voltage sources (not literal
constants -- see the round-4 report for why: mkb()/mkf() constant-fold a
literal-only expression at PARSE time, before it ever reaches PTeval()/
PTevalDual() at all, so a constant-only deck like
tests/regression/parser/bxpressn-1.cir exercises the PARSER's folding, not
the runtime evaluator this round's candidate (0106, forward-mode AD)
changes).

This is a VALUE-level regression/crash/coverage check (every B-source's
v(node) output, sampled at several times along its input ramp, compared
between the baseline and 0106-forward-ad candidate via tests/test_accel.py)
-- not by itself proof that the AD candidate's *derivatives* agree with the
symbolic path (value computation is decoupled from derivative computation in
both code paths, so a wrong derivative would not necessarily show up in
v(node) here, especially without a feedback loop -- see research/
transient-acceleration.md section 9 for the direct derivative-vs-derivative
comparison, done separately with a throwaway self-check instrumented build,
never committed). This deck's job is coverage (does every function/operator
parse, run, and produce the same VALUE) and to be a crash/regression guard
going forward, not the primary derivative-correctness evidence.

Do not hand-edit bench/decks/g_ad_stress.cir -- regenerate it from this
script instead (same convention as every other bench/gen_*.py).
"""
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "decks" / "g_ad_stress.cir"

# Sample times (ms) along every ramp's 0..4ms sweep, used for every .meas.
SAMPLE_TIMES_MS = [0.4, 1.2, 2.0, 2.8, 3.6]

# (name, PT_*/PTF_* coverage note, B-source RHS expression)
# vw   : wide ramp, -3.7 .. 3.7           (abs/sgn/cos/cosh/sinh/sin/tanh/
#                                           exp/ceil/floor/nint/uramp/u/u2/
#                                           uminus/asinh/atan/ddt)
# vtan : narrow ramp, -1.2 .. 1.2         (tan -- avoids the +-pi/2 asymptote)
# vu   : unit ramp, -0.95 .. 0.95         (acos/asin/atanh)
# vh   : >=1 ramp, 1.05 .. 4.0            (acosh)
# vp   : positive ramp, 0.05 .. 5.0       (sqrt/ln/log/log10)
# va   : positive ramp, 0.5 .. 3.0        (pow/pwr/^ base, min/max operand)
# vb   : signed ramp, -2.0 .. 3.0         (pow/pwr/^ exponent, min/max operand)
# vc   : signed ramp, -3.0 .. 3.0         (eq0/ne0/gt0/lt0/ge0/le0, comparisons,
#                                           booleans, ternary condition)
# vpwlx: -12 .. 22                        (pwl/pwl_derivative, matches the
#                                           breakpoint table below)
FUNCS = [
    ("abs", "PTF_ABS", "abs(V(vw))"),
    ("sgn", "PTF_SGN", "sgn(V(vw))"),
    ("acos", "PTF_ACOS", "acos(V(vu))"),
    ("acosh", "PTF_ACOSH", "acosh(V(vh))"),
    ("asin", "PTF_ASIN", "asin(V(vu))"),
    ("asinh", "PTF_ASINH", "asinh(V(vw))"),
    ("atan", "PTF_ATAN", "atan(V(vw))"),
    ("atanh", "PTF_ATANH", "atanh(V(vu))"),
    ("cos", "PTF_COS", "cos(V(vw))"),
    ("cosh", "PTF_COSH", "cosh(V(vw))"),
    ("exp", "PTF_EXP", "exp(V(vw))"),
    ("ln", "PTF_LOG (via 'ln')", "ln(V(vp))"),
    ("log", "PTF_LOG (via 'log')", "log(V(vp))"),
    ("log10", "PTF_LOG10", "log10(V(vp))"),
    ("sin", "PTF_SIN", "sin(V(vw))"),
    ("sinh", "PTF_SINH", "sinh(V(vw))"),
    ("sqrt", "PTF_SQRT", "sqrt(V(vp))"),
    ("tan", "PTF_TAN", "tan(V(vtan))"),
    ("tanh", "PTF_TANH", "tanh(V(vw))"),
    ("ustep", "PTF_USTEP (via 'u')", "u(V(vw))"),
    ("uramp", "PTF_URAMP", "uramp(V(vw))"),
    ("ceil", "PTF_CEIL", "ceil(V(vw))"),
    ("floor", "PTF_FLOOR", "floor(V(vw))"),
    ("nint", "PTF_NINT", "nint(V(vw))"),
    ("uminus", "PTF_UMINUS (unary '-')", "-V(vw)"),
    ("ustep2", "PTF_USTEP2 (via 'u2')", "u2(V(vw))"),
    ("pwl", "PTF_PWL (its own derivative -- PTF_PWL_DERIVATIVE -- is "
             "exercised internally computing this entry's gradient; a "
             "DIRECT user call to pwl_derivative(...) is not a supported "
             "form in the baseline either -- PT_mkfnode() only special-"
             "cases PTF_PWL/PTF_DDT's own 'prepare_PTF_*' comma-list-to-"
             "data step, never PTF_PWL_DERIVATIVE's -- confirmed pre-"
             "existing on aaba50ba3f49-s37 before this round's patch, not "
             "attempted here, see the round-4 report)",
     "pwl(V(vpwlx), -10,-8, 2,4, 12,9, 22,19)"),
    ("eq0", "PTF_EQ0", "eq0(V(vc))"),
    ("ne0", "PTF_NE0", "ne0(V(vc))"),
    ("gt0", "PTF_GT0", "gt0(V(vc))"),
    ("lt0", "PTF_LT0", "lt0(V(vc))"),
    ("ge0", "PTF_GE0", "ge0(V(vc))"),
    ("le0", "PTF_LE0", "le0(V(vc))"),
    ("pow_bconst", "PTF_POW, b-constant branch", "pow(V(va),2)"),
    ("pow_aconst", "PTF_POW, a-constant branch", "pow(2,V(vb))"),
    ("pow_general", "PTF_POW, general branch", "pow(V(va),V(vb))"),
    ("pwr_bconst", "PTF_PWR, b-constant branch", "pwr(V(va),3)"),
    ("pwr_general", "PTF_PWR, general branch", "pwr(V(va),V(vb))"),
    ("min", "PTF_MIN", "min(V(va),V(vb))"),
    ("max", "PTF_MAX", "max(V(va),V(vb))"),
    ("ddt", "PTF_DDT", "ddt(V(vw))"),
    ("caret_bconst", "PT_POWER '^', b-constant branch", "V(va)^2"),
    ("caret_aconst", "PT_POWER '^', a-constant branch", "2^V(vb)"),
    ("caret_general", "PT_POWER '^', general branch", "V(va)^V(vb)"),
    ("tern", "PT_TERN ('?:')", "V(vc) > 0 ? sin(V(vw)) : cos(V(vw))"),
    ("ternfcn", "PT_TERN (explicit ternary_fcn)",
     "ternary_fcn(V(vc), V(va), V(vb))"),
    ("cmp_eq", "eq0(a-b) via '=='", "V(va) == V(vb)"),
    ("cmp_ne", "ne0(a-b) via '!='", "V(va) != V(vb)"),
    ("cmp_gt", "gt0(a-b) via '>'", "V(va) > V(vb)"),
    ("cmp_lt", "lt0(a-b) via '<'", "V(va) < V(vb)"),
    ("cmp_ge", "ge0(a-b) via '>='", "V(va) >= V(vb)"),
    ("cmp_le", "le0(a-b) via '<='", "V(va) <= V(vb)"),
    ("booland", "eq0(eq0+eq0) via '&&'", "V(va) && V(vb)"),
    ("boolor", "ne0(ne0+ne0) via '||'", "V(va) || V(vb)"),
    ("boolnot", "eq0(x) via '!'", "!V(vc)"),
    ("ptime", "PT_TIME", "time*1000"),
    ("ptemper", "PT_TEMPERATURE", "temper"),
    ("compound", "mixed PT_PLUS/MINUS/TIMES/DIVIDE/POWER + nested funcs, "
                  "4 controlling variables",
     "tanh(2*(V(va)-V(vb)))*sin(V(vw)) + sqrt(abs(V(vc))) "
     "- V(va)^2/exp(0.1*V(vb))"),
]

HEADER = """\
* g_ad_stress -- forward-mode AD coverage stress deck (ngspice_augmentation
* card #486, plan T9.4, round 4; generated by bench/gen_ad_stress.py -- do
* not hand-edit, regenerate instead). One B-source per PT_*/PTF_* node type
* and operator src/spicelib/parser/inpptree.c's PTdifferentiate() handles
* (see that script's own header for why literal-constant expressions, e.g.
* tests/regression/parser/bxpressn-1.cir, do NOT exercise this: they fold
* away at parse time before ever reaching the runtime evaluator). Every
* argument is a PWL-ramped voltage (not a literal), so every B-source is a
* real PT_VAR-rooted tree with a live gradient. Synthetic, not derived
* from any real design -- synthetic coverage only.
* run: ngspice -b <this file> (writes nothing outside
*  ~/Work/data/ngspice_augmentation/accel/bench/g_ad_stress/, never into the
*  repo)

Vw   vw   0 PWL(0 -3.7    4m 3.7)
Vtan vtan 0 PWL(0 -1.2    4m 1.2)
Vu   vu   0 PWL(0 -0.95   4m 0.95)
Vh   vh   0 PWL(0 1.05    4m 4.0)
Vp   vp   0 PWL(0 0.05    4m 5.0)
Va   va   0 PWL(0 0.5     4m 3.0)
Vb   vb   0 PWL(0 -2.0    4m 3.0)
Vc   vc   0 PWL(0 -3.0    4m 3.0)
Vpwlx vpwlx 0 PWL(0 -12   4m 22)

* nodeset: ngspice's default Newton initial guess is 0V for every node,
* which is OUT OF DOMAIN for acosh (x>=1) -- without a nodeset, the very
* first Newton trial at t=0 evaluates acosh's derivative formula
* (1/sqrt(u^2-1)) at u=0, hitting sqrt(-1), which forces gmin/source-
* stepping fallback and (found by actually running this deck, not assumed)
* still fails to find an operating point at all on THIS fork's unmodified
* baseline (aaba50ba3f49-s37) even with that fallback -- a pre-existing
* ngspice limitation for any acosh(V(node))-shaped B-source with a poor
* initial guess, not something this round's AD patch introduced (confirmed
* by reproducing it with the same expression on the baseline binary
* directly). nodeset every ramp's t=0 value so Newton starts in-domain and
* this stays a clean, converging correctness deck rather than an
* accidental op-point-finder stress test (that is a different, real
* research question -- see the round-4 report -- not this deck's job).
.nodeset v(vw)=-3.7 v(vtan)=-1.2 v(vu)=-0.95 v(vh)=1.05 v(vp)=0.05
+ v(va)=0.5 v(vb)=-2.0 v(vc)=-3.0 v(vpwlx)=-12
"""

FOOTER = """
.tran 10u 4m 0
.options reltol=1e-6

.control
set filetype=binary
rusage cputime
run
rusage cputime
"""


def main():
    lines = [HEADER]
    for name, note, expr in FUNCS:
        lines.append(f"* {name}: {note}")
        lines.append(f"B{name} out_{name} 0 V = {expr}")
        lines.append("")
    lines.append(FOOTER)
    for name, _note, _expr in FUNCS:
        for t in SAMPLE_TIMES_MS:
            tname = str(t).replace(".", "p")
            lines.append(
                f"meas tran {name}_t{tname} find v(out_{name}) at={t}m"
            )
    lines.append("rusage cputime")
    lines.append("quit")
    lines.append(".endc")
    lines.append(".end")

    OUT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT} ({len(FUNCS)} functions/operators, "
          f"{len(FUNCS) * len(SAMPLE_TIMES_MS)} .meas points)")


if __name__ == "__main__":
    main()
