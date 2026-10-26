#!/usr/bin/env python3
"""
power_calc.py: closed-form power / sample-size / implied-effect calculations
for the trial-powering skill.

Standard library only (uses math.erf for the normal CDF). Prints JSON to
stdout so the calling agent can parse results deterministically.

Subcommands
-----------
forward   Compute power given effect + N (and the rest).
n         Compute N per arm given effect + target power.
implied   Solve for the effect size given N + target power.
grid      Sweep one or more parameters and emit a power table.
cp        Conditional power at an interim: P(final significant | interim).
          Flags: --z-interim Z (observed interim z, sign = direction of
          benefit), --info-frac T (information fraction 0-1, events-based for
          TTE), plus --alpha/--sided. Reports CP under three drifts: current
          trend, original hypothesis (pass expected final z via --drift-z,
          e.g. from `forward`: z_drift = z_alpha + z_power), and null.
assurance Probability of success: power averaged over a prior on the effect
          (O'Hagan; "Bayesian expected power"), counting only a significant
          result in the direction of benefit. Flags: --prior-mean and
          --prior-sd (HR for tte, rate ratio for count, delta for continuous,
          p2 for binary), plus the design flags (--events for tte, otherwise
          --n1; --alpha,
          --sided, --ratio, --sd for continuous, --p1 for binary, etc.). For
          tte and count the prior is lognormal, and --prior-mean/--prior-sd
          are its mean and sd on the ratio scale (the mean, not the median).
          For continuous and binary the prior is normal. Build the prior from
          the Phase 3 assumptions audit: prior mean = precedent-informed
          effect, prior sd = spread of the reference class.

Non-inferiority: pass --ni-margin M to forward/n, for binary, continuous or
tte endpoints only; anything else stops with an error. Binary: M is the absolute
risk-difference margin (treatment worse by M is the null). Continuous: M on
the delta scale. TTE: M is the HR margin (e.g. 1.3). NI is always analyzed
one-sided; pass --sided 1 --alpha 0.025 unless the SAP says otherwise. The
true effect defaults to zero difference (HR 1.0) unless you supply
--p2/--delta/--hr.

Conventions (documented because they matter at the margins):
- Binary two-sample uses POOLED variance under H0 (matches the standard SAP
  z-test). Tools using unpooled H0 variance differ by ~1pp at unequal
  allocation.
- Continuous uses the normal approximation, not t. At n<30/arm this is
  optimistic by roughly one patient per arm; add ~z^2/4 ≈ 1 per arm or use
  exact t software when N is small.
- grid --vary-dropout models ATTRITION (loss of information), not ITT
  dilution. Effect dilution from noncompliance must be put into the effect
  parameters per references/02-stats-methodology.md.

Endpoint types: binary | continuous | tte | count

Common flags
------------
--alpha A           significance level (default 0.05)
--sided {1,2}       sidedness (default 2)
--power P           target power for `n` and `implied` (default 0.80)
--ratio R           allocation ratio treat:control, e.g. 1 or 2 (default 1)
--dropout D         fraction lost; inflates required N by 1/(1-D) (default 0)
--alpha-spent X     alpha already spent at interims; uses (alpha - X) at final.
                    This is conservative: if the SAP gives the final nominal
                    alpha, pass that as --alpha instead.
--direction D       which way benefit runs: higher or lower. Pass it for binary
                    and continuous endpoints, so an effect in the harmful
                    direction doesn't count as power. Without it, forward and
                    grid count either direction, `assurance` assumes the way
                    the prior mean points and says so, and `implied --type
                    binary` reports both solutions. For tte and count, benefit
                    is always a ratio below 1.

Output is always JSON. --dropout applies to forward, implied, assurance and grid as
attrition (N or events multiplied by 1-D), and inflates required N in `n`.
So with --dropout, --n1 means patients enrolled, not patients evaluable. An
event-driven trial runs until it reaches its target events whatever the
dropout, so for tte pass --dropout only for a calendar-driven analysis; the
output carries a warning when events were reduced. `n` also reports each size
rounded up to a whole number (n1_ceil, n2_ceil, n_total_ceil,
events_required_ceil); quote those.

binary
------
--p1 P              control proportion
--p2 P              treatment proportion (forward/n)  OR
--delta D           absolute risk difference (alt to p2)
--n1 N              control arm N (forward/implied)
--single-arm --p_null P   one-sample test vs historical p_null: forward (--p2,
                    --n1), n (--p2) and implied (--n1). Not supported by grid.

continuous
----------
--delta D           mean difference
--sd S              common SD
--ancova-r2 R2      baseline-covariate r^2 (reduces effective variance)
--n1 N              per-arm N (forward/implied; ratio applied for arm 2)

tte (time-to-event, log-rank / Schoenfeld)
-----------------------------------------
--hr H              hazard ratio (treat vs control)
--events E          total events (forward/implied operate on events)
--n-from-events     also emit N given --med-ctrl, --accrual, --followup
--med-ctrl M        control median (same units as accrual/followup)
--accrual A         uniform accrual duration
--followup F        minimum follow-up after last patient in

count (negative binomial rate ratio)
------------------------------------
--rate-ctrl R       control event rate per unit time
--rr X              rate ratio (treat/control)
--exposure T        mean exposure time per patient
--dispersion K      NB overdispersion (Var = mu + k*mu^2); k=0 -> Poisson
--n1 N              per-arm N (forward/implied)

grid
----
--vary VALUES        comma list of effect values to sweep, as effect=a,b,c or a,b,c
                    (meaning depends on type: binary->p2, continuous->delta,
                    tte->hr, count->rr). control=, sd= and dropout= are also
                    accepted and behave like the flags below, and events=a,b,c
                    sweeps the event count for tte. Pass --vary once.
--vary-control      comma list of control-side values (binary->p1, count->rate-ctrl;
                    an error for continuous and tte, where power depends only on
                    delta or events: sweep the effect instead)
--vary-sd           comma list of SDs (continuous) or dispersion (count)
--vary-dropout      comma list of dropout fractions
All other params fixed via the flags above. Output: list of {inputs, power}.

Examples
--------
# Implied HR for a 330-event TTE design at 90% power, 1:1, two-sided 0.05
python power_calc.py implied --type tte --events 330 --power 0.90

# Power for two proportions 0.30 vs 0.45, n=120/arm
python power_calc.py forward --type binary --p1 0.30 --p2 0.45 --n1 120

# Sensitivity grid on HR
python power_calc.py grid --type tte --events 330 --vary effect=0.65,0.70,0.75,0.80,0.85
"""
import argparse, json, math, sys

SQRT2 = math.sqrt(2.0)

def norm_cdf(x):
    return 0.5 * (1.0 + math.erf(x / SQRT2))

def norm_ppf(p):
    # Acklam's rational approximation for the inverse normal CDF.
    if p <= 0.0 or p >= 1.0:
        raise ValueError("norm_ppf domain (0,1)")
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
          1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
          6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2*math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > phigh:
        q = math.sqrt(-2*math.log(1-p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5; r = q*q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)

def z_alpha(alpha, sided):
    return norm_ppf(1 - alpha / sided)

# ---------- binary (two independent proportions) ----------
def binary_power(p1, p2, n1, ratio, alpha, sided, sign=0):
    # sign: 0 = benefit is whichever way p2 differs; +1/-1 = only higher/lower counts
    n2 = n1 * ratio
    pbar = (n1*p1 + n2*p2) / (n1 + n2)
    se0 = math.sqrt(pbar*(1-pbar)*(1/n1 + 1/n2))          # under H0 (pooled)
    se1 = math.sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2)          # under H1
    za = z_alpha(alpha, sided)
    z = ((sign*(p2 - p1) if sign else abs(p2 - p1)) - za*se0) / se1
    return norm_cdf(z)

def binary_n(p1, p2, power, ratio, alpha, sided):
    za = z_alpha(alpha, sided); zb = norm_ppf(power)
    pbar = (p1 + ratio*p2) / (1 + ratio)
    q1, q2, qbar = 1-p1, 1-p2, 1-pbar
    num = (za*math.sqrt((1+1/ratio)*pbar*qbar) + zb*math.sqrt(p1*q1 + p2*q2/ratio))**2
    return num / (p2 - p1)**2  # n1 (control arm)

def binary_single_arm_power(p_null, p_alt, n, alpha, sided):
    se0 = math.sqrt(p_null*(1-p_null)/n)
    se1 = math.sqrt(p_alt*(1-p_alt)/n)
    za = z_alpha(alpha, sided)
    return norm_cdf((abs(p_alt - p_null) - za*se0) / se1)

def binary_single_arm_n(p_null, p_alt, power, alpha, sided):
    za = z_alpha(alpha, sided); zb = norm_ppf(power)
    num = za*math.sqrt(p_null*(1-p_null)) + zb*math.sqrt(p_alt*(1-p_alt))
    return (num / (p_alt - p_null))**2

# ---------- continuous (two means) ----------
def continuous_power(delta, sd, n1, ratio, alpha, sided, r2=0.0, sign=0):
    sd_eff = sd * math.sqrt(max(1 - r2, 1e-12))
    n2 = n1 * ratio
    se = sd_eff * math.sqrt(1/n1 + 1/n2)
    za = z_alpha(alpha, sided)
    return norm_cdf((sign*delta if sign else abs(delta))/se - za)

def continuous_n(delta, sd, power, ratio, alpha, sided, r2=0.0):
    sd_eff = sd * math.sqrt(max(1 - r2, 1e-12))
    za = z_alpha(alpha, sided); zb = norm_ppf(power)
    return ((za + zb)**2) * (sd_eff**2) * (1 + 1/ratio) / (delta**2)  # n1

# ---------- time-to-event (Schoenfeld) ----------
def tte_power(hr, events, ratio, alpha, sided, sign=0):
    za = z_alpha(alpha, sided)
    pi = ratio / (1 + ratio)  # fraction in treatment arm
    se = 1.0 / math.sqrt(events * pi * (1 - pi))
    return norm_cdf((sign*math.log(hr) if sign else abs(math.log(hr))) / se - za)

def tte_events(hr, power, ratio, alpha, sided):
    za = z_alpha(alpha, sided); zb = norm_ppf(power)
    pi = ratio / (1 + ratio)
    return ((za + zb)**2) / (math.log(hr)**2 * pi * (1 - pi))

def tte_prob_event(med, accrual, followup):
    """P(event observed) for one arm: exponential survival, uniform accrual.
    Simpson's rule on S(followup), S(followup+accrual/2), S(followup+accrual)."""
    if med <= 0:
        return 1.0
    lam = math.log(2) / med
    def S(t): return math.exp(-lam * t)
    if accrual <= 0:
        return 1 - S(followup)
    return 1 - (S(followup) + 4*S(followup + accrual/2) + S(followup + accrual)) / 6.0

def tte_n_from_events(events, med_ctrl, hr, accrual, followup, ratio):
    med_trt = med_ctrl / hr if hr > 0 else float('inf')
    d_c = tte_prob_event(med_ctrl, accrual, followup)
    d_t = tte_prob_event(med_trt, accrual, followup)
    pi = ratio / (1 + ratio)
    d_bar = (1 - pi) * d_c + pi * d_t
    return events / d_bar  # total N

# ---------- count (negative binomial rate ratio) ----------
def count_power(rate_ctrl, rr, exposure, k, n1, ratio, alpha, sided, sign=0):
    r1 = rate_ctrl; r2 = rate_ctrl * rr
    v1 = (1.0 / (r1 * exposure)) + k
    v2 = (1.0 / (r2 * exposure)) + k
    se = math.sqrt(v1 / n1 + v2 / (n1 * ratio))
    za = z_alpha(alpha, sided)
    return norm_cdf((sign*math.log(rr) if sign else abs(math.log(rr))) / se - za)

def count_n(rate_ctrl, rr, exposure, k, power, ratio, alpha, sided):
    r1 = rate_ctrl; r2 = rate_ctrl * rr
    v1 = (1.0 / (r1 * exposure)) + k
    v2 = (1.0 / (r2 * exposure)) + k
    za = z_alpha(alpha, sided); zb = norm_ppf(power)
    return ((za + zb)**2) * (v1 + v2 / ratio) / (math.log(rr)**2)  # n1

# ---------- implied-effect solvers (bisection on power - target) ----------
def bisect(f, lo, hi, tol=1e-6, maxit=200):
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None
    for _ in range(maxit):
        mid = (lo + hi) / 2; fm = f(mid)
        if abs(fm) < tol or (hi - lo) < tol:
            return mid
        if flo * fm < 0: hi, fhi = mid, fm
        else: lo, flo = mid, fm
    return (lo + hi) / 2

def implied_binary(p1, n1, ratio, alpha, sided, power, direction="higher"):
    # solve for p2 on the chosen side of p1
    f = lambda p2: binary_power(p1, p2, n1, ratio, alpha, sided) - power
    if direction == "lower":
        return bisect(f, max(1e-4, p1*0.01), p1 - 1e-6)
    hi = min(0.9999, p1 + 0.99*(1-p1))
    return bisect(f, p1 + 1e-6, hi)

def implied_continuous(sd, n1, ratio, alpha, sided, power, r2=0.0):
    sd_eff = sd * math.sqrt(max(1 - r2, 1e-12))
    za = z_alpha(alpha, sided); zb = norm_ppf(power)
    n2 = n1 * ratio
    return (za + zb) * sd_eff * math.sqrt(1/n1 + 1/n2)  # delta

def implied_tte(events, ratio, alpha, sided, power):
    za = z_alpha(alpha, sided); zb = norm_ppf(power)
    pi = ratio / (1 + ratio)
    lhr = (za + zb) / math.sqrt(events * pi * (1 - pi))
    return math.exp(-lhr)  # HR < 1

def implied_count(rate_ctrl, exposure, k, n1, ratio, alpha, sided, power):
    # power isn't monotone near rr=0, so step down from 1 to the first crossing
    f = lambda rr: count_power(rate_ctrl, rr, exposure, k, n1, ratio, alpha, sided) - power
    hi = 0.999
    while hi > 0.01:
        lo = hi - 0.01
        if f(hi) * f(lo) <= 0:
            return bisect(f, lo, hi)
        hi = lo
    return None

# ---------- CLI ----------
def parse_kv_list(s):
    # "effect=0.7,0.75,0.8"  or  "0.7,0.75,0.8"
    if s is None: return None, []
    if "=" in s: key, vals = s.split("=", 1)
    else: key, vals = None, s
    try:
        return key, [float(x) for x in vals.split(",") if x.strip()]
    except ValueError:
        raise SystemExit(f"--vary expects numbers, as in effect=0.7,0.8 or dropout=0.1,0.2; got {s!r}")

# ---------- assurance (power averaged over a normal prior) ----------
def assurance(power_fn, prior_mean, prior_sd, steps=4001):
    """E[power(effect)] for effect ~ N(prior_mean, prior_sd^2). Power is close
    to a step function for large trials, so integrate on a dense grid over
    +/- 6 sd rather than with a few quadrature nodes."""
    if prior_sd <= 0:
        return power_fn(prior_mean)
    total = 0.0; h = 12.0 / steps
    for i in range(steps):
        z = -6.0 + (i + 0.5) * h
        total += math.exp(-0.5*z*z) * power_fn(prior_mean + prior_sd*z) * h
    return total / math.sqrt(2.0*math.pi)

def lognormal_params(mean, sd):
    """Log-scale mean and sd of a lognormal with this mean and sd."""
    if mean <= 0:
        raise SystemExit("assurance: --prior-mean must be positive for tte and count (it's a ratio)")
    s2 = math.log(1 + (sd/mean)**2)
    return math.log(mean) - s2/2, math.sqrt(s2)

def conditional_power(z_interim, t, z_alpha_final, drift):
    """B-value formulation. drift = E[Z(1)] under the assumed truth.
    CP = Phi( (B(t) + drift*(1-t) - z_alpha_final) / sqrt(1-t) )"""
    b = z_interim * math.sqrt(t)
    return norm_cdf((b + drift * (1 - t) - z_alpha_final) / math.sqrt(1 - t))

TTE_DROPOUT_NOTE = ("--dropout reduced the event count; an event-driven trial runs until it reaches its "
                    "target events, so use this only for a calendar-driven analysis, otherwise vary --events")

DIRECTION_NOTE = ("--direction not given, so benefit was assumed to run the way the prior mean points; "
                  "a sceptical prior can reverse this, so pass --direction")

def finish(out, a):
    """Echo the flags used if the branch did not, add whole-number sizes, print."""
    if not out["inputs"]:
        skip = {"cmd", "type"}
        if a.type != "count": skip |= {"exposure", "dispersion"}
        if a.type != "continuous": skip |= {"ancova_r2"}
        out["inputs"] = {k: v for k, v in vars(a).items()
                         if k not in skip and v is not None and v is not False}
    for k in ("n1", "n2", "n", "events_required"):
        if k in out: out[k + "_ceil"] = math.ceil(out[k] - 1e-9)
    if "n1_ceil" in out and "n2_ceil" in out:
        out["n_total_ceil"] = out["n1_ceil"] + out["n2_ceil"]
    elif "n_total" in out:
        out["n_total_ceil"] = math.ceil(out["n_total"] - 1e-9)
    json.dump(out, sys.stdout, indent=2); print()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["forward", "n", "implied", "grid", "cp", "assurance"])
    ap.add_argument("--type", required=False, choices=["binary", "continuous", "tte", "count"])
    ap.add_argument("--z-interim", type=float)
    ap.add_argument("--info-frac", type=float)
    ap.add_argument("--drift-z", type=float, help="expected final z under the design hypothesis (z_alpha + z_power)")
    ap.add_argument("--ni-margin", type=float, help="non-inferiority margin: abs risk diff (binary), delta (continuous), or HR bound like 1.3 (tte)")
    ap.add_argument("--prior-mean", type=float, help="assurance: prior mean of the effect (HR for tte, rate ratio for count)")
    ap.add_argument("--prior-sd", type=float, help="assurance: prior sd of the effect, on the same scale as --prior-mean")
    ap.add_argument("--direction", choices=["higher", "lower"], help="which way benefit runs, for binary and continuous endpoints")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--sided", type=int, default=2, choices=[1, 2])
    ap.add_argument("--power", type=float, default=0.80)
    ap.add_argument("--ratio", type=float, default=1.0)
    ap.add_argument("--dropout", type=float, default=0.0)
    ap.add_argument("--alpha-spent", type=float, default=0.0)
    # binary
    ap.add_argument("--p1", type=float); ap.add_argument("--p2", type=float)
    ap.add_argument("--delta", type=float)
    ap.add_argument("--single-arm", action="store_true"); ap.add_argument("--p_null", type=float)
    # continuous
    ap.add_argument("--sd", type=float); ap.add_argument("--ancova-r2", type=float, default=0.0)
    # shared n
    ap.add_argument("--n1", type=float)
    # tte
    ap.add_argument("--hr", type=float); ap.add_argument("--events", type=float)
    ap.add_argument("--med-ctrl", type=float); ap.add_argument("--accrual", type=float)
    ap.add_argument("--followup", type=float); ap.add_argument("--n-from-events", action="store_true")
    # count
    ap.add_argument("--rate-ctrl", type=float); ap.add_argument("--rr", type=float)
    ap.add_argument("--exposure", type=float, default=1.0); ap.add_argument("--dispersion", type=float, default=0.0)
    # grid
    ap.add_argument("--vary", action="append"); ap.add_argument("--vary-control"); ap.add_argument("--vary-sd"); ap.add_argument("--vary-dropout")
    a = ap.parse_args()
    if a.cmd != "cp" and a.type is None:
        raise SystemExit(f"{a.cmd} requires --type")
    if a.ni_margin is not None and (a.cmd not in ("forward", "n") or a.type == "count"):
        raise SystemExit("--ni-margin is only supported by forward and n, for binary, continuous and tte endpoints")
    if a.single_arm and a.cmd in ("grid", "assurance"):
        raise SystemExit(f"--single-arm is not supported by {a.cmd}; call forward once per value")

    alpha = max(a.alpha - a.alpha_spent, 1e-6)
    out = {"cmd": a.cmd, "type": a.type, "inputs": {}}

    if a.cmd == "cp":
        if a.z_interim is None or a.info_frac is None:
            raise SystemExit("cp requires --z-interim and --info-frac")
        t = a.info_frac
        za = z_alpha(alpha, a.sided)
        current_trend_drift = a.z_interim / math.sqrt(t)
        out["inputs"] = {"z_interim": a.z_interim, "info_frac": t, "alpha_final": alpha, "sided": a.sided}
        out["cp_current_trend"] = conditional_power(a.z_interim, t, za, current_trend_drift)
        out["cp_null"] = conditional_power(a.z_interim, t, za, 0.0)
        if a.drift_z is not None:
            out["cp_design_hypothesis"] = conditional_power(a.z_interim, t, za, a.drift_z)
        json.dump(out, sys.stdout, indent=2); print(); return

    def inflate(n):  # dropout inflation for required-N outputs
        return n / max(1 - a.dropout, 1e-9)

    if a.cmd in ("forward", "implied", "assurance") and a.dropout:  # attrition, as in grid
        if a.n1: a.n1 *= 1 - a.dropout
        if a.events:
            a.events *= 1 - a.dropout
            out["warnings"] = [TTE_DROPOUT_NOTE]
        out["dropout_applied"] = a.dropout

    t = a.type
    sgn = {"higher": 1, "lower": -1}.get(a.direction, 0)  # 0: benefit is whichever way the effect points

    if t == "tte" and a.cmd != "n" and not a.events and not (a.vary and a.vary[0].startswith("events=")):
        raise SystemExit(f"{a.cmd} --type tte needs --events; estimate them from N first (see references/02-stats-methodology.md)")
    if a.n_from_events and not (a.med_ctrl and a.accrual is not None and a.followup is not None):
        raise SystemExit("--n-from-events needs --med-ctrl, --accrual and --followup")

    if a.cmd == "assurance":
        if a.prior_mean is None or a.prior_sd is None:
            raise SystemExit("assurance requires --prior-mean and --prior-sd")
        if t == "tte":
            # prior given on HR scale; convert to log scale (lognormal approx)
            m_log, s_log = lognormal_params(a.prior_mean, a.prior_sd)
            fn = lambda lhr: tte_power(math.exp(lhr), a.events, a.ratio, alpha, a.sided, sign=-1)
            val = assurance(fn, m_log, s_log)
            out["inputs"] = {"prior_hr_mean": a.prior_mean, "prior_hr_sd": a.prior_sd, "events": a.events}
        elif t == "continuous":
            sign = {"higher": 1, "lower": -1}.get(a.direction, 1 if a.prior_mean >= 0 else -1)
            fn = lambda d: continuous_power(d, a.sd, a.n1, a.ratio, alpha, a.sided, a.ancova_r2, sign=sign)
            val = assurance(fn, a.prior_mean, a.prior_sd)
            if not a.direction:
                out.setdefault("warnings", []).append(DIRECTION_NOTE)
            out["inputs"] = {"prior_delta_mean": a.prior_mean, "prior_delta_sd": a.prior_sd, "n1": a.n1, "sd": a.sd,
                             "benefit_direction": "higher" if sign > 0 else "lower"}
        elif t == "binary":
            sign = {"higher": 1, "lower": -1}.get(a.direction, 1 if a.prior_mean >= a.p1 else -1)
            fn = lambda p2: binary_power(a.p1, min(max(p2, 1e-4), 0.9999), a.n1, a.ratio, alpha, a.sided, sign=sign)
            val = assurance(fn, a.prior_mean, a.prior_sd)
            if not a.direction:
                out.setdefault("warnings", []).append(DIRECTION_NOTE)
            out["inputs"] = {"p1": a.p1, "prior_p2_mean": a.prior_mean, "prior_p2_sd": a.prior_sd, "n1": a.n1,
                             "benefit_direction": "higher" if sign > 0 else "lower"}
        elif t == "count":
            m_log, s_log = lognormal_params(a.prior_mean, a.prior_sd)
            fn = lambda lrr: count_power(a.rate_ctrl, math.exp(lrr), a.exposure, a.dispersion, a.n1, a.ratio, alpha, a.sided, sign=-1)
            val = assurance(fn, m_log, s_log)
            out["inputs"] = {"prior_rr_mean": a.prior_mean, "prior_rr_sd": a.prior_sd, "n1": a.n1}
        out["assurance"] = val
        out["note"] = "probability of success averaged over the prior; compare against conditional power at the point assumption"
        json.dump(out, sys.stdout, indent=2); print(); return

    # Non-inferiority: test against the margin instead of zero/one. dist is the
    # signed distance from the NI null on the test's scale; a true effect worse
    # than the margin gives dist <= 0, so power falls below alpha.
    if a.ni_margin is not None:
        sign = -1 if a.direction == "lower" else 1
        za = z_alpha(alpha, a.sided)
        if a.sided == 2:
            out.setdefault("warnings", []).append("NI is conventionally one-sided (alpha 0.025); you passed --sided 2")
        elif alpha > 0.025:
            out.setdefault("warnings", []).append("NI is conventionally one-sided at alpha 0.025; you passed a larger alpha")
        if t == "binary":
            # Direct NI formulas with variances at the TRUE rates (Blackwelder).
            p2 = a.p2 if a.p2 is not None else a.p1 + (a.delta or 0.0)   # default: truly equal
            dist = sign*(p2 - a.p1) + a.ni_margin
            var1 = a.p1*(1-a.p1); var2 = p2*(1-p2)
            out["ni"] = {"margin_abs_diff": a.ni_margin, "true_p2_assumed": p2, "benefit_direction": a.direction or "higher"}
        elif t == "continuous":
            d = a.delta if a.delta is not None else 0.0
            dist = sign*d + a.ni_margin
            var1 = var2 = (a.sd**2) * max(1 - a.ancova_r2, 1e-12)
            out["ni"] = {"margin_delta": a.ni_margin, "true_delta_assumed": d, "benefit_direction": a.direction or "higher"}
        else:  # tte
            hr = a.hr if a.hr is not None else 1.0     # default: truly equal hazards
            dist = math.log(a.ni_margin) - math.log(hr)
            pi = a.ratio / (1 + a.ratio)
            out["ni"] = {"margin_hr": a.ni_margin, "true_hr_assumed": hr}
        if a.cmd == "n":
            if dist <= 0:
                raise SystemExit("the assumed true effect is at or beyond the NI margin, so no sample size can demonstrate non-inferiority")
            zb = norm_ppf(a.power)
            if t == "tte":
                ev = (za + zb)**2 / (dist**2 * pi * (1 - pi))
                out["events_required"] = ev
                if a.n_from_events and a.med_ctrl and a.accrual is not None and a.followup is not None:
                    out["n_total"] = inflate(tte_n_from_events(ev, a.med_ctrl, hr, a.accrual, a.followup, a.ratio))
            else:
                n1 = (za + zb)**2 * (var1 + var2/a.ratio) / dist**2
                out.update({"n1": inflate(n1), "n2": inflate(n1*a.ratio), "n_total": inflate(n1*(1+a.ratio))})
        else:
            if t == "tte":
                se = 1.0 / math.sqrt(a.events * pi * (1 - pi))
            else:
                se = math.sqrt(var1/a.n1 + var2/(a.n1*a.ratio))
            out["power"] = norm_cdf(dist/se - za)
        finish(out, a); return

    if a.cmd in ("forward", "n", "implied"):
        if t == "binary":
            if a.single_arm:
                if a.cmd == "forward":
                    out["inputs"] = {"p_null": a.p_null, "p_alt": a.p2, "n": a.n1}
                    out["power"] = binary_single_arm_power(a.p_null, a.p2, a.n1, alpha, a.sided)
                elif a.cmd == "n":
                    out["inputs"] = {"p_null": a.p_null, "p_alt": a.p2}
                    out["n"] = inflate(binary_single_arm_n(a.p_null, a.p2, a.power, alpha, a.sided))
                else:
                    f = lambda p: binary_single_arm_power(a.p_null, p, a.n1, alpha, a.sided) - a.power
                    if a.direction == "lower":
                        out["implied_p_alt"] = bisect(f, 1e-4, a.p_null - 1e-6)
                    else:
                        out["implied_p_alt"] = bisect(f, a.p_null + 1e-6, 0.9999)
            elif a.cmd == "forward":
                p2 = a.p2 if a.p2 is not None else a.p1 + a.delta
                out["inputs"] = {"p1": a.p1, "p2": p2, "n1": a.n1, "ratio": a.ratio}
                out["power"] = binary_power(a.p1, p2, a.n1, a.ratio, alpha, a.sided, sign=sgn)
            elif a.cmd == "n":
                p2 = a.p2 if a.p2 is not None else a.p1 + a.delta
                n1 = binary_n(a.p1, p2, a.power, a.ratio, alpha, a.sided)
                out.update({"n1": inflate(n1), "n2": inflate(n1*a.ratio), "n_total": inflate(n1*(1+a.ratio))})
            else:  # implied
                if a.direction:
                    p2 = implied_binary(a.p1, a.n1, a.ratio, alpha, a.sided, a.power, a.direction)
                    out.update({"implied_p2": p2, "implied_abs_diff": (p2 - a.p1) if p2 else None})
                else:
                    out["implied_p2_if_higher"] = implied_binary(a.p1, a.n1, a.ratio, alpha, a.sided, a.power, "higher")
                    out["implied_p2_if_lower"] = implied_binary(a.p1, a.n1, a.ratio, alpha, a.sided, a.power, "lower")
                    out["note"] = "pass --direction, or use the root that matches this endpoint's direction of benefit"
        elif t == "continuous":
            if a.cmd == "forward":
                out["inputs"] = {"delta": a.delta, "sd": a.sd, "n1": a.n1, "ratio": a.ratio, "r2": a.ancova_r2}
                out["power"] = continuous_power(a.delta, a.sd, a.n1, a.ratio, alpha, a.sided, a.ancova_r2, sign=sgn)
            elif a.cmd == "n":
                n1 = continuous_n(a.delta, a.sd, a.power, a.ratio, alpha, a.sided, a.ancova_r2)
                out.update({"n1": inflate(n1), "n2": inflate(n1*a.ratio), "n_total": inflate(n1*(1+a.ratio))})
            else:
                if a.sd is not None:
                    d = implied_continuous(a.sd, a.n1, a.ratio, alpha, a.sided, a.power, a.ancova_r2)
                    out.update({"implied_delta": d, "implied_cohens_d": d / a.sd})
                else:
                    d = implied_continuous(1.0, a.n1, a.ratio, alpha, a.sided, a.power, a.ancova_r2)
                    out.update({"implied_cohens_d": d, "note": "sd not supplied; reporting standardized effect"})
        elif t == "tte":
            if a.cmd == "forward":
                out["inputs"] = {"hr": a.hr, "events": a.events, "ratio": a.ratio}
                out["power"] = tte_power(a.hr, a.events, a.ratio, alpha, a.sided, sign=-1)
            elif a.cmd == "n":
                ev = tte_events(a.hr, a.power, a.ratio, alpha, a.sided)
                out["events_required"] = ev
                if a.n_from_events and a.med_ctrl and a.accrual is not None and a.followup is not None:
                    out["n_total"] = inflate(tte_n_from_events(ev, a.med_ctrl, a.hr, a.accrual, a.followup, a.ratio))
            else:
                hr = implied_tte(a.events, a.ratio, alpha, a.sided, a.power)
                out.update({"implied_hr": hr, "implied_reduction_pct": (1-hr)*100})
        elif t == "count":
            if a.cmd == "forward":
                out["inputs"] = {"rate_ctrl": a.rate_ctrl, "rr": a.rr, "exposure": a.exposure, "k": a.dispersion, "n1": a.n1}
                out["power"] = count_power(a.rate_ctrl, a.rr, a.exposure, a.dispersion, a.n1, a.ratio, alpha, a.sided, sign=-1)
            elif a.cmd == "n":
                n1 = count_n(a.rate_ctrl, a.rr, a.exposure, a.dispersion, a.power, a.ratio, alpha, a.sided)
                out.update({"n1": inflate(n1), "n2": inflate(n1*a.ratio), "n_total": inflate(n1*(1+a.ratio))})
            else:
                rr = implied_count(a.rate_ctrl, a.exposure, a.dispersion, a.n1, a.ratio, alpha, a.sided, a.power)
                out.update({"implied_rr": rr, "implied_reduction_pct": (1-rr)*100 if rr else None})

    elif a.cmd == "grid":
        if a.vary and len(a.vary) > 1:
            raise SystemExit("pass --vary once; use --vary-control, --vary-sd and --vary-dropout for other axes")
        key, vary_vals = parse_kv_list(a.vary[0] if a.vary else None)
        if key not in (None, "effect", "control", "sd", "dropout", "events"):
            raise SystemExit(f"--vary: unknown parameter {key!r}; use effect, control, sd, dropout or events")
        if key == "events" and t != "tte":
            raise SystemExit("--vary events applies to tte only")
        eff_vals = vary_vals if key in (None, "effect") else []
        event_vals = vary_vals if key == "events" else [a.events]
        _, ctrl_vals = parse_kv_list(a.vary_control)
        _, sd_vals = parse_kv_list(a.vary_sd)
        _, drop_vals = parse_kv_list(a.vary_dropout)
        if key == "control": ctrl_vals = vary_vals
        if key == "sd": sd_vals = vary_vals
        if key == "dropout": drop_vals = vary_vals
        if sd_vals and t in ("binary", "tte"):
            raise SystemExit(f"--vary-sd has no effect for {t}")
        if ctrl_vals and t in ("continuous", "tte"):
            raise SystemExit(f"--vary-control has no effect for {t}: power depends on the effect, so sweep that instead")
        eff_vals = eff_vals or [None]; ctrl_vals = ctrl_vals or [None]
        sd_vals = sd_vals or [None]; drop_vals = drop_vals or [a.dropout]
        rows = []
        for e in eff_vals:
            for c in ctrl_vals:
                for s in sd_vals:
                    for dr in drop_vals:
                        n1_eff = (a.n1 or 0) * (1 - dr) if a.n1 else None
                        ev_eff = None
                        if t == "binary":
                            p1 = c if c is not None else a.p1
                            p2 = e if e is not None else (a.p2 if a.p2 is not None else p1 + a.delta)
                            pw = binary_power(p1, p2, n1_eff, a.ratio, alpha, a.sided, sign=sgn)
                            rows.append({"p1": p1, "p2": p2, "dropout": dr, "n1_eff": n1_eff, "power": pw})
                        elif t == "continuous":
                            sd = s if s is not None else a.sd
                            d = e if e is not None else a.delta
                            pw = continuous_power(d, sd, n1_eff, a.ratio, alpha, a.sided, a.ancova_r2, sign=sgn)
                            rows.append({"delta": d, "sd": sd, "dropout": dr, "n1_eff": n1_eff, "power": pw})
                        elif t == "tte":
                            hr = e if e is not None else a.hr
                            for ev in event_vals:
                                ev_eff = ev * (1 - dr)
                                pw = tte_power(hr, ev_eff, a.ratio, alpha, a.sided, sign=-1)
                                rows.append({"hr": hr, "events_eff": ev_eff, "dropout": dr, "power": pw})
                        elif t == "count":
                            rc = c if c is not None else a.rate_ctrl
                            rr = e if e is not None else a.rr
                            k = s if s is not None else a.dispersion
                            pw = count_power(rc, rr, a.exposure, k, n1_eff, a.ratio, alpha, a.sided, sign=-1)
                            rows.append({"rate_ctrl": rc, "rr": rr, "k": k, "dropout": dr, "n1_eff": n1_eff, "power": pw})
        out["grid"] = rows
        if t == "tte" and a.events and any(drop_vals):
            out["warnings"] = [TTE_DROPOUT_NOTE]

    finish(out, a)

if __name__ == "__main__":
    main()
