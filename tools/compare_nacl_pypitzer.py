#!/usr/bin/env python3
"""Optional PyPitzer v2 NaCl reference comparison; never used to qualify the C++ solver.

The published PyPI release is PyPitzer 1.x; this uses an explicitly pinned
2.0 source commit. Its _properties entry point is private, so treat it as a
version-locked research dependency, not as the project's supported API.

Run: python tools/compare_nacl_pypitzer.py --aqueous-cli build/aqueous_cli --output nacl.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

REFERENCE_COMMIT = "5a774ae52e41914a7237fee55ef0bbb1e7c76eff"
REFERENCE_LN_GAMMA_NA_AT_ONE_M = -0.41707260461259255
MOLALITIES = (0.001, 0.01, 0.05, 0.09, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0)
T_K = 298.15

def run_cpp_davies(exe: Path, apparent_molarity: float) -> float:
    """Return C++ Davies gamma from a neutral NaCl model at the same numeric c.

    This deliberately treats input molality as an equal numerical MOLARITY;
    it is a sensitivity illustration, NOT a valid molal-to-molar conversion.
    """
    model = (
        "temperature_k 298.15\n"
        "reference_temperature_k 298.15\n"
        "activity_model davies\n"
        f"component Na+ 1 {apparent_molarity:.16g}\n"
        f"component Cl- -1 {apparent_molarity:.16g}\n"
    )
    with tempfile.TemporaryDirectory(prefix="nacl_comparison_") as temp:
        filename = Path(temp) / "nacl.txt"
        filename.write_text(model, encoding="utf-8")
        completed = subprocess.run(
            [str(exe.resolve()), str(filename)],
            check=False, text=True, capture_output=True, timeout=90,
        )
    if completed.returncode:
        raise RuntimeError(
            f"C++ Davies calculation at numerical molarity {apparent_molarity:g} failed: "
            + completed.stderr.strip()
        )
    result = json.loads(completed.stdout)
    if result["activity_model"] != "davies":
        raise RuntimeError("Unexpected C++ activity model")
    gamma_na = result["activities"]["Na+"] / result["concentrations_mol_L"]["Na+"]
    gamma_cl = result["activities"]["Cl-"] / result["concentrations_mol_L"]["Cl-"]
    if not math.isclose(gamma_na, gamma_cl, rel_tol=1e-10):
        raise RuntimeError("Unexpected C++ ion-specific Davies asymmetry for NaCl")
    return math.sqrt(gamma_na * gamma_cl)

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--aqueous-cli", type=Path, required=True)
    ap.add_argument("--output", type=Path, default=Path("nacl-pypitzer-comparison.csv"))
    args = ap.parse_args()

    if not args.aqueous_cli.is_file():
        ap.error(f"C++ solver executable not found: {args.aqueous_cli}")
    try:
        # The property calculator is PRIVATE in PyPitzer 2.0 and intentionally
        # imported only by this optional, source-pinned external benchmark.
        from pypitzer._properties import PropertyCalculator
    except ImportError as exc:
        raise SystemExit(
            "PyPitzer 2.0 is required. Install the pinned GitHub revision from "
            "docs/pypitzer-nacl-benchmark.md, not the PyPI 1.x release."
        ) from exc

    calculator = PropertyCalculator()
    rows = []
    for m in MOLALITIES:
        props = calculator.calculate({"Na+": m, "Cl-": m}, temperature_k=T_K)
        ln_na = props.ln_activity_coefficients["Na+"]
        ln_cl = props.ln_activity_coefficients["Cl-"]
        gamma_mean = math.exp((ln_na + ln_cl) / 2.0)
        if not all(map(math.isfinite, (ln_na, ln_cl, gamma_mean, props.water_activity))):
            raise RuntimeError(f"Nonfinite reference output at m={m}")
        if m == 1.0 and not math.isclose(
            ln_na, REFERENCE_LN_GAMMA_NA_AT_ONE_M, rel_tol=1e-10, abs_tol=1e-11
        ):
            raise RuntimeError(
                "Reference NaCl 1 m result disagrees with PyPitzer v2's own regression test; "
                "check reference commit, database, and dependency versions."
            )
        cpp_gamma = run_cpp_davies(args.aqueous_cli, m) if m <= 0.09 else None
        deviation = (
            100 * (cpp_gamma / gamma_mean - 1.0) if cpp_gamma is not None else None
        )
        rows.append({
            "temperature_K": T_K,
            "nacl_molality_mol_kg_water": m,
            "pypitzer_gamma_na": math.exp(ln_na),
            "pypitzer_gamma_cl": math.exp(ln_cl),
            "pypitzer_mean_gamma": gamma_mean,
            "pypitzer_water_activity": props.water_activity,
            "pypitzer_osmotic_coefficient": props.osmotic_coefficient,
            "cpp_davies_gamma_at_equal_numeric_mol_L": cpp_gamma,
            "davies_relative_difference_percent": deviation,
            "comparison_status": (
                "qualitative_only_molality_vs_molarity"
                if cpp_gamma is not None else "pypitzer_only_beyond_cpp_davies_domain"
            ),
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print("PyPitzer2 reference commit:", REFERENCE_COMMIT)
    print("Source of reference: GitHub 2.0 development code; NOT PyPI 1.x")
    print("25 C NaCl; m is mol/kg WATER. C++ only tested for numerically equal mol/L up to 0.09.")
    print("This is model comparison, NOT externally validated chemistry.")
    print("   m (mol/kg)  Pitzer mean gamma   C++ Davies gamma*  relative difference (%)")
    for row in rows:
        g = row["cpp_davies_gamma_at_equal_numeric_mol_L"]
        pct = row["davies_relative_difference_percent"]
        print(
            f"{row['nacl_molality_mol_kg_water']:12.3g}"
            f"{row['pypitzer_mean_gamma']:19.8f}"
            f"{g:20.8f}{pct:25.4f}" if g is not None else
            f"{row['nacl_molality_mol_kg_water']:12.3g}"
            f"{row['pypitzer_mean_gamma']:19.8f}"
            f"{'outside C++ domain':>20}{'n/a':>25}"
        )
    print("* same numeric amount, not a molality-to-molarity density conversion")
    print("CSV:", args.output)
    return 0

if __name__ == "__main__":
    sys.exit(main())
