"""Deterministic checks for trial-powering's power_calc.py.

Pinned values come from textbook results (Schoenfeld events, Blackwelder NI
sizing, O'Hagan assurance) and from direct numerical checks.

Stdlib only: run with `python3 test_power_calc.py` from this directory.
"""

import json
import pathlib
import subprocess
import sys
import unittest

SCRIPT = pathlib.Path(__file__).resolve().parent / "power_calc.py"


def run(*args):
    out = subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True, text=True, check=True,
    )
    return json.loads(out.stdout)


class TestBinary(unittest.TestCase):
    def test_forward_two_proportions(self):
        # Two-proportion pooled-variance z-test: 0.30 vs 0.45 at 200/arm.
        r = run("forward", "--type", "binary",
                "--p1", "0.30", "--p2", "0.45", "--n1", "200")
        self.assertAlmostEqual(r["power"], 0.875, delta=0.005)

    def test_ni_sizing_at_expected_rates(self):
        # Eval case 10a: cure rate 0.70 both arms, 10pp NI margin,
        # one-sided 0.025, 80% power -> ~330/arm (band 320-340).
        r = run("n", "--type", "binary", "--p1", "0.70", "--p2", "0.70",
                "--ni-margin", "0.10", "--sided", "1", "--alpha", "0.025",
                "--power", "0.80")
        n = r.get("n1") or r.get("n_per_arm") or r["inputs"].get("n1")
        self.assertTrue(320 <= n <= 340, f"n/arm {n} outside 320-340")

    def test_ni_sizing_rate_drift_increases_n(self):
        # Eval case 10c: at 0.60 the variance grows -> ~377/arm (370-385).
        r = run("n", "--type", "binary", "--p1", "0.60", "--p2", "0.60",
                "--ni-margin", "0.10", "--sided", "1", "--alpha", "0.025",
                "--power", "0.80")
        n = r.get("n1") or r.get("n_per_arm") or r["inputs"].get("n1")
        self.assertTrue(370 <= n <= 385, f"n/arm {n} outside 370-385")


class TestContinuous(unittest.TestCase):
    def test_forward_full_trial(self):
        # Eval case 8: delta 4, SD 12, 200/arm -> power ~0.92 (0.90-0.93).
        r = run("forward", "--type", "continuous",
                "--delta", "4", "--sd", "12", "--n1", "200")
        self.assertTrue(0.90 <= r["power"] <= 0.93, r["power"])

    def test_forward_half_subgroup(self):
        # Eval case 8: biomarker subgroup at 100/arm -> ~0.65 (0.60-0.70).
        r = run("forward", "--type", "continuous",
                "--delta", "4", "--sd", "12", "--n1", "100")
        self.assertTrue(0.60 <= r["power"] <= 0.70, r["power"])

    def test_mde_at_80(self):
        # Eval case 8: MDE at 80% power, 200/arm, SD 12 -> 3.3-3.4 points.
        r = run("implied", "--type", "continuous",
                "--sd", "12", "--n1", "200", "--power", "0.80")
        delta = r["implied_delta"]
        self.assertTrue(3.3 <= delta <= 3.4, f"MDE {delta} outside 3.3-3.4")

    def test_sizing_with_dropout(self):
        # Eval case 2: delta 4.2, SD 10, 15% dropout -> ~210 total at 80%.
        r = run("n", "--type", "continuous", "--delta", "4.2", "--sd", "10",
                "--dropout", "0.15", "--power", "0.80")
        n_total = r.get("n_total") or 2 * (r.get("n1") or 0)
        self.assertTrue(205 <= n_total <= 214, f"n_total {n_total}")


class TestTTE(unittest.TestCase):
    def test_schoenfeld_events(self):
        # Eval case 11 forward check: HR 0.72 at 90% power -> ~389 events
        # (band 385-395).
        r = run("n", "--type", "tte", "--hr", "0.72", "--power", "0.90")
        e = r["events_required"]
        self.assertTrue(385 <= e <= 395, f"events {e} outside 385-395")

    def test_power_at_reference_class_center(self):
        # Eval case 11: with 390 events, power at the class-median HR 0.80
        # is ~0.60 (band 0.57-0.62).
        r = run("forward", "--type", "tte", "--hr", "0.80", "--events", "390")
        self.assertTrue(0.57 <= r["power"] <= 0.62, r["power"])

    def test_assurance_under_precedent_prior(self):
        # Eval case 11: assurance with prior N(HR 0.80, sd 0.10) at 390
        # events -> ~0.56 (band 0.52-0.60).
        r = run("assurance", "--type", "tte", "--events", "390",
                "--prior-mean", "0.80", "--prior-sd", "0.10")
        a = r.get("assurance") or r.get("probability_of_success")
        self.assertTrue(0.52 <= a <= 0.60, f"assurance {a} outside 0.52-0.60")


def fails(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True).returncode != 0


class Regressions(unittest.TestCase):
    def test_ni_unsupported_paths_stop(self):
        self.assertTrue(fails("grid", "--type", "tte", "--events", "300", "--ni-margin", "1.3", "--vary", "effect=1.0"))
        self.assertTrue(fails("implied", "--type", "tte", "--events", "300", "--ni-margin", "1.3"))
        self.assertTrue(fails("forward", "--type", "count", "--rate-ctrl", ".5", "--rr", "1", "--n1", "100", "--ni-margin", "1.2"))

    def test_ni_true_effect_worse_than_margin(self):
        p = run("forward", "--type", "tte", "--events", "300", "--hr", "1.6", "--ni-margin", "1.3", "--sided", "1", "--alpha", ".025")
        self.assertLess(p["power"], 0.01)
        p = run("forward", "--type", "continuous", "--delta", "-0.8", "--sd", "1", "--n1", "100", "--ni-margin", "0.3", "--sided", "1", "--alpha", ".025")
        self.assertLess(p["power"], 0.01)
        self.assertTrue(fails("n", "--type", "tte", "--hr", "1.6", "--ni-margin", "1.3"))

    def test_ni_uses_true_hr_for_n_from_events(self):
        p = run("n", "--type", "tte", "--hr", "1.0", "--ni-margin", "1.3", "--sided", "1", "--alpha", ".025",
                "--n-from-events", "--med-ctrl", "12", "--accrual", "12", "--followup", "12")
        self.assertAlmostEqual(p["n_total"], 713.4, delta=1)  # 456 events / 0.639 event probability

    def test_ni_binary_delta_and_direction(self):
        a = run("forward", "--type", "binary", "--p1", ".8", "--delta", ".05", "--n1", "200", "--ni-margin", ".1", "--sided", "1", "--alpha", ".025")
        b = run("forward", "--type", "binary", "--p1", ".8", "--p2", ".85", "--n1", "200", "--ni-margin", ".1", "--sided", "1", "--alpha", ".025")
        self.assertAlmostEqual(a["power"], b["power"], places=6)
        m = run("n", "--type", "binary", "--p1", ".20", "--p2", ".22", "--ni-margin", ".10", "--sided", "1", "--alpha", ".025", "--direction", "lower")
        self.assertAlmostEqual(m["n1"], 407, delta=3)

    def test_assurance_lower_is_better(self):
        c = run("assurance", "--type", "continuous", "--prior-mean", "-3", "--prior-sd", "1", "--sd", "8", "--n1", "200")
        self.assertAlmostEqual(c["assurance"], 0.868, delta=0.01)
        b = run("assurance", "--type", "binary", "--p1", ".30", "--prior-mean", ".18", "--prior-sd", ".02", "--n1", "250")
        self.assertGreater(b["assurance"], 0.7)

    def test_assurance_accuracy_for_wide_priors(self):
        c = run("assurance", "--type", "continuous", "--prior-mean", "2", "--prior-sd", "3", "--sd", "10", "--n1", "1000")
        self.assertAlmostEqual(c["assurance"], 0.650, delta=0.005)

    def test_implied_binary_reports_both_roots(self):
        p = run("implied", "--type", "binary", "--p1", ".30", "--n1", "250", "--power", ".9")
        self.assertAlmostEqual(p["implied_p2_if_lower"], 0.177, delta=0.003)
        self.assertAlmostEqual(p["implied_p2_if_higher"], 0.439, delta=0.003)

    def test_grid_vary_keys(self):
        g = run("grid", "--type", "continuous", "--delta", "3", "--sd", "10", "--n1", "200", "--vary", "sd=8,10")
        self.assertEqual([r["delta"] for r in g["grid"]], [3, 3])
        self.assertEqual([r["sd"] for r in g["grid"]], [8, 10])
        d = run("grid", "--type", "tte", "--hr", ".7", "--events", "300", "--vary", "dropout=0.1,0.2")
        self.assertEqual([r["hr"] for r in d["grid"]], [.7, .7])
        self.assertTrue(fails("grid", "--type", "tte", "--events", "330", "--vary", "effect,control,dropout"))
        self.assertTrue(fails("grid", "--type", "tte", "--events", "330", "--vary-control", "5,6"))

    def test_single_arm(self):
        n = run("n", "--type", "binary", "--single-arm", "--p_null", ".2", "--p2", ".35", "--power", ".9")
        p = run("forward", "--type", "binary", "--single-arm", "--p_null", ".2", "--p2", ".35", "--n1", str(n["n"]))
        self.assertAlmostEqual(p["power"], 0.9, delta=0.001)
        i = run("implied", "--type", "binary", "--single-arm", "--p_null", ".2", "--n1", "50", "--power", ".9")
        self.assertGreater(i["implied_p_alt"], 0.2)
        self.assertTrue(fails("grid", "--type", "binary", "--single-arm", "--p_null", ".2", "--n1", "50", "--vary", "effect=.3"))

    def test_dropout_consistent_across_commands(self):
        f = run("forward", "--type", "continuous", "--delta", "3", "--sd", "10", "--n1", "200", "--dropout", ".3")
        g = run("grid", "--type", "continuous", "--delta", "3", "--sd", "10", "--n1", "200", "--dropout", ".3")
        self.assertAlmostEqual(f["power"], g["grid"][0]["power"], places=6)

    def test_implied_count_finds_root(self):
        p = run("implied", "--type", "count", "--rate-ctrl", ".3", "--dispersion", "1", "--n1", "100")
        self.assertIsNotNone(p["implied_rr"])

    def test_missing_type_stops(self):
        self.assertTrue(fails("forward"))


if __name__ == "__main__":
    unittest.main()
