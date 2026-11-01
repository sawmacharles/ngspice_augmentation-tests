"""Loop-gain maths for the stb testbench (card #451, 2026-09-27).

Pure standard library. Formulas: Tian et al. 2001 eq. 30 as mapped to the lgprobe
subcircuit in research/cadence-stb-pstb.md §2.4; Middlebrook 1975 from the same runs.
T is Bode's return ratio (textbook sign: positive at low frequency for negative
feedback, -1 at oscillation). Spectre's loopGain is -T.
"""
import cmath
import math


def tian(v1, i1, v2, i2):
    """T from run 1 (voltage injection) and run 2 (current injection) probe readings."""
    s = 2 * (i1 * v2 - v1 * i2) + v1 + i2
    return s / (1 - s)


def middlebrook(v1, va1, i2, ib2):
    """Middlebrook 1975 double injection. va1 = v(a) in run 1, ib2 = i(v.<inst>.vb) in run 2.

    Orientation dependent: the probe must have a = downstream, b = upstream (see stb/lgprobe.lib).
    Comparison only.
    """
    tv = -v1 / va1
    ti = -ib2 / i2
    return (tv * ti - 1) / (tv + ti + 2)


def tian_from_y(y11, y12, y21, y22):
    """T from the two-port short-circuit Y parameters of the broken loop (independent reference)."""
    return (y12 + y21) / (y11 + y22)


def db(z):
    return 20 * math.log10(abs(z)) if z != 0 else -math.inf


def unwrapped_phase_deg(values):
    """Phase in degrees, unwrapped along the list, first point in (-180, 180]."""
    out = []
    prev = None
    for z in values:
        p = math.degrees(cmath.phase(z))
        if prev is not None:
            while p - prev > 180:
                p -= 360
            while p - prev < -180:
                p += 360
        out.append(p)
        prev = p
    return out


def _interp_logf(f0, f1, y0, y1, target):
    """Frequency where a quantity linear in log f crosses target between two points."""
    if y1 == y0:
        return f0, 0.0
    t = (target - y0) / (y1 - y0)
    return f0 * (f1 / f0) ** t, t


POSITIVE_FEEDBACK_NOTE = (
    "positive feedback at low frequency: loop is unstable if |T|>1 there; "
    "phase/gain margins are not meaningful"
)

# card #476: card #464's original low_f_positive_feedback trigger (Re(T) < 0 at
# the lowest swept frequency) false-triggered on healthy integrator-type loops
# (e.g. a type-III compensated buck, tests/fixtures/compopt_r3_averaged_detuned.cir
# at fstart=10 Hz): T is genuinely close to -90 deg there, and a few degrees of
# extra lag from a real pole makes Re(T) slightly negative even though the loop
# is fine. DECISION (option A, card #476): only the two conditions below
# fire the WARNING; a loop that is merely near an integrator gets a neutral
# label instead, and a loop with |arg T|>150 deg but |T|<=1 gets a note, not a
# warning (it cannot be read off Bode margins either way, but it is not
# necessarily unstable the way |T|>1 with that phase is).
POSITIVE_FEEDBACK_LOW_GAIN_NOTE = "positive feedback at low frequency, loop gain < 1"
INTEGRATOR_LABEL = "integrator-like low-frequency behaviour"
LOW_F_POSFB_ARG_DEG = 150.0    # |arg T(fmin)| must exceed this for either positive-feedback case
LOW_F_INTEGRATOR_LO_DEG = -120.0  # arg T(fmin) in [LO, HI] -> integrator-like label
LOW_F_INTEGRATOR_HI_DEG = -60.0   # (-90 +/- 30 deg)


def margins(freqs, t_values):
    """Phase and gain margins of T (Bode sign).

    Returns a dict with every crossing found, and the worst margin of each kind.
    Phase margin at |T| = 1: PM = 180 + phase(T) (unwrapped from low frequency).
    Gain margin where phase(T) = -180 - 360k: GM = -20 log10|T|.
    No crossing is reported as None, never as a number.

    Each 0 dB crossing (in pm_crossings, and the selected "pm") carries a
    "falling" flag: True when |T| is decreasing through that crossing (the
    usual stability-margin crossing), False when |T| is rising through it (a
    PM computed there, even one under 180 deg, is not the textbook margin --
    see card #464).

    Low-frequency behaviour of T at the lowest swept frequency (card #464,
    narrowed by card #476 -- see POSITIVE_FEEDBACK_LOW_GAIN_NOTE/INTEGRATOR_LABEL
    above for why): let arg0 = arg(T(fmin)) in degrees, mag0 = |T(fmin)|.

    - "low_f_positive_feedback" is True only when |arg0| > 150 deg AND
      mag0 > 1 -- positive feedback WITH loop gain, so the loop cannot be
      stable in the usual sense; "note" then carries POSITIVE_FEEDBACK_NOTE.
      This is the WARNING case (found 2026-09-29 on refs/r1c_fully_diff.cir
      before its card #464 fix: T(1 Hz) = -12.9/-3333 read back as 274/186
      deg "phase margins", which are crossings of a loop that is actually
      unstable at DC for any |T|>1 -- not real margins).
    - "low_f_positive_feedback_low_gain" is True when |arg0| > 150 deg but
      mag0 <= 1 -- positive feedback, but loop gain under 1 there, so this is
      reported as a note (low_f_positive_feedback_low_gain_note), not the
      warning above; "low_f_positive_feedback" is False in this case.
    - "low_f_integrator" is True when arg0 is within -90 +/- 30 deg (i.e. in
      [-120, -60]) -- ordinary integrator-like behaviour at low frequency (a
      type-III compensator's near-integrator branch, for instance); reported
      as the neutral informational low_f_integrator_note, never a warning.

    These three are mutually exclusive (disjoint arg0 ranges) and each may be
    False/None on a loop whose T(fmin) doesn't fall in any of them (the
    ordinary case: Re(T(fmin)) > 0, |arg0| well under 150 deg).
    """
    mags = [db(z) for z in t_values]
    phs = unwrapped_phase_deg(t_values)
    pm_list, gm_list = [], []
    for k in range(len(freqs) - 1):
        f0, f1 = freqs[k], freqs[k + 1]
        m0, m1 = mags[k], mags[k + 1]
        if (m0 - 0.0) * (m1 - 0.0) < 0 or (m1 == 0.0 and m0 != 0.0):
            fc, t = _interp_logf(f0, f1, m0, m1, 0.0)
            ph = phs[k] + t * (phs[k + 1] - phs[k])
            pm_list.append({"freq": fc, "pm_deg": 180.0 + ph, "falling": m1 < m0})
        # phase crossings of -180 + 360n; a point sitting exactly on a target counts once,
        # in the interval it ends
        p0, p1 = phs[k], phs[k + 1]
        lo, hi = min(p0, p1), max(p0, p1)
        n = math.ceil((lo + 180) / 360)
        while -180 + 360 * n <= hi:
            target = -180 + 360 * n
            if target != p0:
                fc, t = _interp_logf(f0, f1, p0, p1, target)
                mag = m0 + t * (m1 - m0)
                gm_list.append({"freq": fc, "gm_db": -mag, "phase_deg": target})
            n += 1
    pm = min(pm_list, key=lambda d: d["pm_deg"]) if pm_list else None
    gm = min(gm_list, key=lambda d: d["gm_db"]) if gm_list else None

    # card #476: arg0/mag0 of T at the lowest swept frequency, and the three
    # (mutually exclusive) low-frequency classifications derived from them.
    if t_values:
        arg0 = phs[0]  # cmath.phase(t_values[0]) in degrees -- phs[0] is unwrapped
                        # from nothing yet, so it IS the raw phase of the first point
        mag0 = abs(t_values[0])
    else:
        arg0 = mag0 = None

    low_f_positive_feedback = bool(t_values) and abs(arg0) > LOW_F_POSFB_ARG_DEG and mag0 > 1.0
    low_f_positive_feedback_low_gain = (bool(t_values) and abs(arg0) > LOW_F_POSFB_ARG_DEG
                                         and mag0 <= 1.0)
    low_f_integrator = (bool(t_values)
                         and LOW_F_INTEGRATOR_LO_DEG <= arg0 <= LOW_F_INTEGRATOR_HI_DEG)

    note = POSITIVE_FEEDBACK_NOTE if low_f_positive_feedback else None
    low_f_positive_feedback_low_gain_note = (
        POSITIVE_FEEDBACK_LOW_GAIN_NOTE if low_f_positive_feedback_low_gain else None)
    low_f_integrator_note = INTEGRATOR_LABEL if low_f_integrator else None

    return {"pm": pm, "gm": gm, "pm_crossings": pm_list, "gm_crossings": gm_list,
            "low_f_positive_feedback": low_f_positive_feedback, "note": note,
            "low_f_positive_feedback_low_gain": low_f_positive_feedback_low_gain,
            "low_f_positive_feedback_low_gain_note": low_f_positive_feedback_low_gain_note,
            "low_f_integrator": low_f_integrator, "low_f_integrator_note": low_f_integrator_note}


def max_rel_err(a, b):
    """max |a-b| / |b| over paired complex lists."""
    worst = 0.0
    for x, y in zip(a, b, strict=True):
        worst = max(worst, abs(x - y) / abs(y))
    return worst
