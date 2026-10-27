"""Unit tests for sample_size_calculator.py. They need scipy, like the calculator."""

import json
import math
import pathlib
import subprocess
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from sample_size_calculator import (  # noqa: E402
    calculate_binary_sample_size,
    calculate_continuous_sample_size,
)

SCRIPT = pathlib.Path(__file__).resolve().parent / "sample_size_calculator.py"


def continuous_arms(**kwargs):
    sizes = calculate_continuous_sample_size(dropout_rate=0.0, **kwargs)["sample_size"]
    return sizes["treatment_arm"], sizes["control_arm"]


def binary_arms(**kwargs):
    sizes = calculate_binary_sample_size(dropout_rate=0.0, **kwargs)["sample_size"]
    return sizes["treatment_arm"], sizes["control_arm"]


class TestContinuous(unittest.TestCase):
    def test_equal_allocation(self):
        # Chow, Shao & Wang, Sample Size Calculations in Clinical Research, sec. 3.2, worked example:
        # SD 0.10, difference 0.05, two-sided alpha 0.05, power 0.80 gives 63 per arm.
        self.assertEqual(continuous_arms(effect_size=0.05, std_dev=0.10), (63, 63))

    def test_five_to_one(self):
        # M. Bland, "Allocation ratios in clinical trials" (University of York): a trial that
        # needs 1000 people at 1:1 (0.205 SD, 90% power) needs 1800 at 5:1.
        treatment, control = continuous_arms(
            effect_size=0.205, std_dev=1.0, power=0.90, allocation_ratio=5.0
        )
        self.assertAlmostEqual(treatment + control, 1800, delta=18)
        self.assertEqual(treatment, 5 * control)

    def test_other_ratios(self):
        # Julious (2004), Statistics in Medicine 23:1921: with r = treatment / control, the control
        # arm needs (r + 1) (z_alpha + z_beta)^2 SD^2 / (r difference^2), here 31.4 (r + 1) / r.
        # 1:2 is not quite 2:1 reversed, because the control arm is the one rounded up first.
        for ratio, expected in ((0.5, (48, 95)), (1.5, (80, 53)), (2.0, (96, 48))):
            with self.subTest(ratio=ratio):
                arms = continuous_arms(effect_size=0.5, std_dev=1.0, allocation_ratio=ratio)
                self.assertEqual(arms, expected)

    def test_sensitivity_rows(self):
        rows = calculate_continuous_sample_size(5.0, 15.0, allocation_ratio=2.0)["sensitivity_analysis"]
        self.assertEqual(rows["power_90_percent"]["total"], 426)
        self.assertEqual(rows["effect_size_reduced_10_percent"]["total"], 393)


class TestBinary(unittest.TestCase):
    def test_sizes_for_60_vs_75_percent(self):
        # Fleiss, Levin & Paik, Statistical Methods for Rates and Proportions, 3rd ed., ch. 4
        # (unequal sample sizes, no continuity correction), two-sided alpha 0.05, power 0.80.
        for ratio, expected in (
            (0.5, (116, 231)), (1.0, (152, 152)), (1.5, (189, 126)), (2.0, (224, 112)), (3.0, (297, 99)),
        ):
            with self.subTest(ratio=ratio):
                self.assertEqual(binary_arms(p1=0.60, p2=0.75, allocation_ratio=ratio), expected)

    def test_sensitivity_rows(self):
        rows = calculate_binary_sample_size(0.60, 0.75)["sensitivity_analysis"]
        self.assertEqual(rows["power_90_percent"]["total"], 406)
        self.assertEqual(rows["effect_size_reduced_10_percent"]["total"], 380)


class TestInputs(unittest.TestCase):
    def test_allocation_must_be_a_positive_number(self):
        for ratio in (0.0, -1.0, math.nan, math.inf):
            with self.subTest(ratio=ratio):
                for result in (
                    calculate_continuous_sample_size(1.0, 1.0, allocation_ratio=ratio),
                    calculate_binary_sample_size(0.60, 0.75, allocation_ratio=ratio),
                ):
                    self.assertIn("must be a positive number", result["error"])

    def test_no_difference_and_no_spread_are_refused(self):
        self.assertIn("must be a nonzero number", calculate_continuous_sample_size(0.0, 1.0)["error"])
        self.assertIn("must be a positive number", calculate_continuous_sample_size(1.0, 0.0)["error"])
        self.assertIn("must differ", calculate_binary_sample_size(0.60, 0.60)["error"])

    def test_both_spellings_of_non_inferiority_give_the_same_one_sided_result(self):
        outputs = [
            subprocess.run(
                [sys.executable, str(SCRIPT), "--type", "continuous", "--effect-size", "5",
                 "--std-dev", "15", "--design", design],
                capture_output=True, text=True, check=True,
            ).stdout
            for design in ("non-inferiority", "noninferiority")
        ]
        self.assertEqual(outputs[0], outputs[1])
        result = json.loads(outputs[0])
        self.assertEqual(result["study_design"], "non-inferiority")
        # 2 (1.645 + 0.842)^2 / (5 / 15)^2 = 111.3
        self.assertEqual(result["sample_size"]["treatment_arm"], 112)


if __name__ == "__main__":
    unittest.main()
