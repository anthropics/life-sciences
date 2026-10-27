#!/usr/bin/env python3
"""
FDA Device Clinical Trial Sample Size Calculator

Simple, accurate statistical power calculations for medical device trials.
Supports continuous and binary primary endpoints.
"""

import argparse
import json
import sys
import math
from typing import Dict, Any

ALLOCATION_ERROR = (
    "Allocation ratio (--allocation) must be a positive number: treatment "
    "participants per control participant, e.g. 2 for 2:1"
)


def _continuous_arms(z_alpha: float, z_beta: float, cohens_d: float, allocation_ratio: float):
    """Return (n_treatment, n_control) for comparing two means, where r = n_treatment / n_control.

    Both are rounded up, and n_treatment is r times the rounded n_control, so a whole-number r gives
    groups exactly in that ratio. Formula: Julious 2004, Stat Med 23:1921-1986.
    """
    r = allocation_ratio
    n_control = math.ceil((1 + 1 / r) * ((z_alpha + z_beta) ** 2) / (cohens_d ** 2))
    # round() so floating-point error (1.1 * 50 = 55.00000000000001) cannot add a participant
    return math.ceil(round(r * n_control, 9)), n_control


def _settings_error(alpha: float, power: float, dropout_rate: float, design: str):
    """Return an error message for settings the calculator cannot use, or None if they are fine."""
    if design != "superiority":
        return (
            "This calculator cannot size non-inferiority trials because it has no option for a "
            "non-inferiority margin."
        )
    if not 0 < alpha < 1:
        return "Alpha (--alpha) must be between 0 and 1 (exclusive), e.g. 0.05"
    if not 0 < power < 1:
        return "Power (--power) must be between 0 and 1 (exclusive), e.g. 0.80"
    if not 0 <= dropout_rate < 1:
        return "Dropout rate (--dropout) must be from 0 up to but not including 1, e.g. 0.15 for 15%"
    return None


def calculate_continuous_sample_size(
    effect_size: float,
    std_dev: float,
    alpha: float = 0.05,
    power: float = 0.80,
    allocation_ratio: float = 1.0,
    dropout_rate: float = 0.15,
    design: str = "superiority"
) -> Dict[str, Any]:
    """
    Calculate sample size for continuous primary endpoint.

    Uses two-sample t-test power calculation.

    Args:
        effect_size: Expected mean difference between groups
        std_dev: Standard deviation (pooled or common)
        alpha: Type I error rate (default 0.05 for two-sided test)
        power: Statistical power (1 - Type II error rate)
        allocation_ratio: Ratio of treatment to control (1.0 = equal allocation)
        dropout_rate: Expected dropout rate (0.15 = 15%)
        design: "superiority" or "non-inferiority"

    Returns:
        Dictionary with sample size results
    """
    try:
        from scipy import stats
    except ImportError:
        return {
            "error": "scipy not installed. Run: pip install scipy statsmodels numpy"
        }

    settings_error = _settings_error(alpha, power, dropout_rate, design)
    if settings_error:
        return {"error": settings_error}
    if not (math.isfinite(allocation_ratio) and allocation_ratio > 0):
        return {"error": ALLOCATION_ERROR}
    if not (math.isfinite(effect_size) and effect_size != 0):
        return {
            "error": "Effect size (--effect-size) must be a nonzero number: the difference in "
            "means the trial should detect. This calculator cannot size a trial that expects no "
            "difference, such as some non-inferiority trials, because it has no option for a "
            "non-inferiority margin."
        }
    if not (math.isfinite(std_dev) and std_dev > 0):
        return {"error": "Standard deviation (--std-dev) must be a positive number"}

    # Calculate effect size (Cohen's d)
    cohens_d = abs(effect_size) / std_dev

    # Determine sidedness
    if design == "superiority":
        alpha_adj = alpha / 2  # Two-sided test
        sidedness = "two-sided"
    else:  # non-inferiority
        alpha_adj = alpha  # One-sided test
        sidedness = "one-sided"

    # Z-scores for alpha and beta
    z_alpha = stats.norm.ppf(1 - alpha_adj)
    z_beta = stats.norm.ppf(power)

    n_treatment, n_control = _continuous_arms(z_alpha, z_beta, cohens_d, allocation_ratio)

    total_n = n_treatment + n_control

    # Adjust for dropout
    total_with_dropout = math.ceil(total_n / (1 - dropout_rate))
    n_treatment_with_dropout = math.ceil(n_treatment / (1 - dropout_rate))
    n_control_with_dropout = math.ceil(n_control / (1 - dropout_rate))

    # Sensitivity analysis: 90% power
    z_beta_90 = stats.norm.ppf(0.90)
    total_90 = sum(_continuous_arms(z_alpha, z_beta_90, cohens_d, allocation_ratio))
    total_90_dropout = math.ceil(total_90 / (1 - dropout_rate))

    # Sensitivity analysis: Effect size reduced by 10%
    cohens_d_reduced = cohens_d * 0.9
    total_reduced = sum(_continuous_arms(z_alpha, z_beta, cohens_d_reduced, allocation_ratio))
    total_reduced_dropout = math.ceil(total_reduced / (1 - dropout_rate))

    return {
        "endpoint_type": "continuous",
        "study_design": design,
        "statistical_test": f"Two-sample t-test ({sidedness})",
        "effect_size": effect_size,
        "standard_deviation": std_dev,
        "cohens_d": round(cohens_d, 3),
        "alpha": alpha,
        "power": power,
        "allocation_ratio": f"{allocation_ratio}:1" if allocation_ratio != 1.0 else "1:1",
        "dropout_rate": dropout_rate,
        "sample_size": {
            "treatment_arm": n_treatment,
            "control_arm": n_control,
            "total": total_n,
            "total_with_dropout": total_with_dropout,
            "treatment_with_dropout": n_treatment_with_dropout,
            "control_with_dropout": n_control_with_dropout
        },
        "sensitivity_analysis": {
            "power_90_percent": {
                "total": total_90,
                "total_with_dropout": total_90_dropout
            },
            "effect_size_reduced_10_percent": {
                "total": total_reduced,
                "total_with_dropout": total_reduced_dropout
            }
        },
        "assumptions": [
            "Normally distributed continuous outcome",
            "Equal variances between groups",
            "Independent observations",
            f"Effect size of {effect_size} is clinically meaningful and realistic"
        ],
        "disclaimers": [
            "Sample size calculation is preliminary and based on assumptions",
            "Final sample size requires biostatistician review and validation",
            "FDA Pre-Submission meeting may require adjustments",
            "Consider pilot study to validate assumptions"
        ]
    }


def _binary_arms(z_alpha: float, z_beta: float, p1: float, p2: float, allocation_ratio: float):
    """Return (n_treatment, n_control) for a two-proportion z-test; p1 is the control proportion,
    p2 the treatment proportion, and r = n_treatment / n_control.

    The variance is pooled under the null hypothesis and unpooled under the alternative.
    Rounding is as in _continuous_arms.
    """
    r = allocation_ratio
    p_pooled = (p1 + r * p2) / (1 + r)
    n_control = math.ceil((z_alpha * math.sqrt(p_pooled * (1 - p_pooled) * (1 + 1 / r)) +
                           z_beta * math.sqrt(p1 * (1 - p1) + p2 * (1 - p2) / r)) ** 2 / ((p2 - p1) ** 2))
    # round() so floating-point error cannot add a participant, as in _continuous_arms
    return math.ceil(round(r * n_control, 9)), n_control


def calculate_binary_sample_size(
    p1: float,
    p2: float,
    alpha: float = 0.05,
    power: float = 0.80,
    allocation_ratio: float = 1.0,
    dropout_rate: float = 0.15,
    design: str = "superiority"
) -> Dict[str, Any]:
    """
    Calculate sample size for binary primary endpoint (proportions).

    Uses two-proportion z-test power calculation.

    Args:
        p1: Expected proportion in control group
        p2: Expected proportion in treatment group
        alpha: Type I error rate (default 0.05 for two-sided test)
        power: Statistical power (1 - Type II error rate)
        allocation_ratio: Ratio of treatment to control (1.0 = equal allocation)
        dropout_rate: Expected dropout rate (0.15 = 15%)
        design: "superiority" or "non-inferiority"

    Returns:
        Dictionary with sample size results
    """
    try:
        from scipy import stats
    except ImportError:
        return {
            "error": "scipy not installed. Run: pip install scipy statsmodels numpy"
        }

    settings_error = _settings_error(alpha, power, dropout_rate, design)
    if settings_error:
        return {"error": settings_error}

    # Validate proportions
    if not (0 < p1 < 1 and 0 < p2 < 1):
        return {"error": "Proportions must be between 0 and 1 (exclusive)"}
    if p1 == p2:
        return {
            "error": "The control and treatment proportions (--p1 and --p2) must differ: the "
            "trial is sized to detect the difference between them. This calculator cannot size "
            "a trial that expects no difference, such as some non-inferiority trials, because "
            "it has no option for a non-inferiority margin."
        }

    if not (math.isfinite(allocation_ratio) and allocation_ratio > 0):
        return {"error": ALLOCATION_ERROR}

    # Determine sidedness
    if design == "superiority":
        alpha_adj = alpha / 2  # Two-sided test
        sidedness = "two-sided"
    else:  # non-inferiority
        alpha_adj = alpha  # One-sided test
        sidedness = "one-sided"

    # Z-scores
    z_alpha = stats.norm.ppf(1 - alpha_adj)
    z_beta = stats.norm.ppf(power)

    # Effect size
    effect_size = abs(p2 - p1)

    n_treatment, n_control = _binary_arms(z_alpha, z_beta, p1, p2, allocation_ratio)

    total_n = n_treatment + n_control

    # Adjust for dropout
    total_with_dropout = math.ceil(total_n / (1 - dropout_rate))
    n_treatment_with_dropout = math.ceil(n_treatment / (1 - dropout_rate))
    n_control_with_dropout = math.ceil(n_control / (1 - dropout_rate))

    # Sensitivity analysis: 90% power
    z_beta_90 = stats.norm.ppf(0.90)
    total_90 = sum(_binary_arms(z_alpha, z_beta_90, p1, p2, allocation_ratio))
    total_90_dropout = math.ceil(total_90 / (1 - dropout_rate))

    # Sensitivity analysis: Effect size reduced by 10%
    effect_reduced = effect_size * 0.9
    p2_reduced = p1 + effect_reduced if p2 > p1 else p1 - effect_reduced

    total_reduced = sum(_binary_arms(z_alpha, z_beta, p1, p2_reduced, allocation_ratio))
    total_reduced_dropout = math.ceil(total_reduced / (1 - dropout_rate))

    return {
        "endpoint_type": "binary",
        "study_design": design,
        "statistical_test": f"Two-proportion z-test ({sidedness})",
        "control_proportion": p1,
        "treatment_proportion": p2,
        "effect_size": round(effect_size, 4),
        "alpha": alpha,
        "power": power,
        "allocation_ratio": f"{allocation_ratio}:1" if allocation_ratio != 1.0 else "1:1",
        "dropout_rate": dropout_rate,
        "sample_size": {
            "treatment_arm": n_treatment,
            "control_arm": n_control,
            "total": total_n,
            "total_with_dropout": total_with_dropout,
            "treatment_with_dropout": n_treatment_with_dropout,
            "control_with_dropout": n_control_with_dropout
        },
        "sensitivity_analysis": {
            "power_90_percent": {
                "total": total_90,
                "total_with_dropout": total_90_dropout
            },
            "effect_size_reduced_10_percent": {
                "total": total_reduced,
                "total_with_dropout": total_reduced_dropout
            }
        },
        "assumptions": [
            "Binary outcome (success/failure, event/no event)",
            "Independent observations",
            "Large enough sample for normal approximation",
            f"Effect size of {round(effect_size * 100, 1)}% is clinically meaningful and realistic"
        ],
        "disclaimers": [
            "Sample size calculation is preliminary and based on assumptions",
            "Final sample size requires biostatistician review and validation",
            "FDA Pre-Submission meeting may require adjustments",
            "Consider pilot study to validate assumptions"
        ]
    }


def main():
    parser = argparse.ArgumentParser(
        description="FDA Device Clinical Trial Sample Size Calculator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Continuous endpoint (mean difference):
    python sample_size_calculator.py --type continuous --effect-size 5.0 --std-dev 15.0

  Binary endpoint (proportions):
    python sample_size_calculator.py --type binary --p1 0.60 --p2 0.75

  Custom parameters:
    python sample_size_calculator.py --type continuous --effect-size 5.0 --std-dev 15.0 \\
        --alpha 0.05 --power 0.90 --dropout 0.20 --output results.json
        """
    )

    parser.add_argument(
        "--type",
        required=True,
        choices=["continuous", "binary"],
        help="Type of primary endpoint"
    )

    # Continuous endpoint parameters
    parser.add_argument(
        "--effect-size",
        type=float,
        help="Expected mean difference (for continuous endpoints)"
    )
    parser.add_argument(
        "--std-dev",
        type=float,
        help="Standard deviation (for continuous endpoints)"
    )

    # Binary endpoint parameters
    parser.add_argument(
        "--p1",
        type=float,
        help="Expected proportion in control group (for binary endpoints)"
    )
    parser.add_argument(
        "--p2",
        type=float,
        help="Expected proportion in treatment group (for binary endpoints)"
    )

    # Common parameters
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Type I error rate (default: 0.05)"
    )
    parser.add_argument(
        "--power",
        type=float,
        default=0.80,
        help="Statistical power (default: 0.80)"
    )
    parser.add_argument(
        "--dropout",
        type=float,
        default=0.15,
        help="Expected dropout rate (default: 0.15)"
    )
    parser.add_argument(
        "--allocation",
        type=float,
        default=1.0,
        help="Allocation ratio treatment:control (default: 1.0 for equal allocation)"
    )
    parser.add_argument(
        "--design",
        choices=["superiority", "non-inferiority", "noninferiority"],
        default="superiority",
        help="Study design type (default: superiority); either non-inferiority spelling stops with an error"
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output JSON file path (if not specified, prints to stdout)"
    )

    args = parser.parse_args()

    # Validate inputs based on endpoint type
    if args.type == "continuous":
        if args.effect_size is None or args.std_dev is None:
            parser.error("--effect-size and --std-dev required for continuous endpoints")

        result = calculate_continuous_sample_size(
            effect_size=args.effect_size,
            std_dev=args.std_dev,
            alpha=args.alpha,
            power=args.power,
            allocation_ratio=args.allocation,
            dropout_rate=args.dropout,
            design=args.design
        )

    elif args.type == "binary":
        if args.p1 is None or args.p2 is None:
            parser.error("--p1 and --p2 required for binary endpoints")

        result = calculate_binary_sample_size(
            p1=args.p1,
            p2=args.p2,
            alpha=args.alpha,
            power=args.power,
            allocation_ratio=args.allocation,
            dropout_rate=args.dropout,
            design=args.design
        )

    # Check for errors
    if "error" in result:
        print(f"ERROR: {result['error']}", file=sys.stderr)
        sys.exit(1)

    # Output results
    output_json = json.dumps(result, indent=2)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(output_json)
        print(f"Results written to {args.output}")
    else:
        print(output_json)


if __name__ == "__main__":
    main()
