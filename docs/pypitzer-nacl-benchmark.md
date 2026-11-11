# 25 °C NaCl activity coefficients: PyPitzer 2 vs C++ Davies

**Classification: reproducible cross-model comparison, not a validation of geochemical predictions.**

## Reference provenance

- PyPitzer [repository](https://github.com/yiping-git/pypitzer), source commit `5a774ae52e41914a7237fee55ef0bbb1e7c76eff`.
- The external benchmark uses PyPitzer **2.0 source code**, **not** the published PyPI 1.1.5 package (they have incompatible APIs).
- PyPitzer 2.0's internal `pypitzer._properties.PropertyCalculator` computes ion activity coefficients, osmotic coefficient, and water activity from a bundled read-only SQLite parameter database (database `2026.1`).
- The dependency is pinned and used only by an optional reference workflow; there is **no change to the production C++ equilibrium or activity models**.
- The reference implementation includes a test asserting `ln(gamma_Na+) = -0.41707260461259255` at **1 mol/kg NaCl and 298.15 K**, which is used as a benchmark self-check.

## Scope and limitations

For 0.001, 0.01, 0.05, 0.09, 0.1, 0.5, 1, 2, 3, and 5 mol/kg NaCl:

1. Compute PyPitzer mean ionic activity coefficient (sqrt(gamma_Na gamma_Cl)) and water activity.
2. For **up to 0.09**, run the *actual compiled C++ aqueous_cli* with a model set to an equal **numerical concentration in mol/L**, and compute its Davies mean ionic coefficient using the returned activities and concentrations.
3. Save and print the cross-model results. Values greater than 0.09 mol/kg are *Pitzer only*; the C++ Davies model restricts ionic strength to approximately 0.1 mol/L and must not be extrapolated beyond its guard.

**Critical:** molality (mol/kg of solvent) is not molarity (mol/L of solution). Using equal numerical values is a *dilute-solution qualitative comparison*, not an exact thermodynamic reference. Replacing it with a defensible molality conversion requires solution density and volume/composition modeling. Moreover, the C++ Davies model uses a molarity-based approximation and has no water-activity treatment. Differences between the curves do **not** prove a C++ bug.

## Run the optional reference benchmark

With Python 3.12, Git, a compiler and CMake:

```bash
python -m pip install "git+https://github.com/yiping-git/pypitzer.git@5a774ae52e41914a7237fee55ef0bbb1e7c76eff"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
python tools/compare_nacl_pypitzer.py --aqueous-cli build/aqueous_cli --output nacl-pypitzer-comparison.csv
```

On Windows the executable may be in `build/Release/aqueous_cli.exe` depending on the chosen generator. The dedicated GitHub Actions reference workflow runs this benchmark and uploads the CSV for inspection.

## Interpretation and next validation gate

At 1 m NaCl, the PyPitzer self-check verifies its published code regression; this is not an independent physical measurement. To establish chemistry fidelity, compare to vetted experimental NaCl mean coefficients and/or a trusted PHREEQC Pitzer database, **with matched Pitzer interaction parameters, standard states, temperature, solute basis and water activity**. After verifying the single-electrolyte case, extend to multicomponent mixtures, then consider an independent C++ Pitzer implementation rather than directly embedding a Python interpreter inside the equilibrium solver.

Do not label the prototype Pitzer-capable or thermodynamically validated solely because this benchmark runs.
