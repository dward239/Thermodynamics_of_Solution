# MVP+1 Verification and Validation Plan

**Current classification:** code builds and synthetic numerical checks pass; independent chemical validation is outstanding.

## Reference conditions

For comparisons, explicitly fix temperature, pressure, reference state, formation constants, enthalpies, chemical-component totals, water model, and activity convention. Unaligned constants or conventions invalidate numerical comparisons.

## Verification cases

| Case | Expected behavior | Status |
| --- | --- | --- |
| 1:1 ideal A+B ⇌ AB | Match closed-form concentration | Automated |
| Carbonate mass balances | Normalized residual under configured tolerance | Automated |
| Water autoionization | Activity product matches input Kw | Automated |
| Charge balance with free H+ | Absolute charge residual small | Automated |
| Davies coefficient | Charged ion activity reduced at dilute I | Automated |
| Temperature | Positive formation ΔH increases K with T | Automated |
| Input validation | Reject charge mismatch, high I, out-of-range T | Automated |
| Nonconvergence | Report error rather than partial result | Intended; add focused case |

## Independent comparison cases required before scientific validation

1. **Ideal carbonate system at 25 °C:** Use a reference solver with identical synthetic constants, basis reactions, totals and ideal activities. Compare pH and species concentrations.
2. **Dilute nonideal electrolyte:** Use a reference solver configured with Davies activity coefficients and the same standard-state/ionic-strength basis. Document any differences due to the current molarity approximation.
3. **Temperature series:** Evaluate 0, 10, 25, 40 and 50 °C with traceable measured or published thermodynamic data, matching van't Hoff treatment.
4. **Adverse numerical conditions:** Very dilute, strongly associated, near singular and high-charge mixtures; characterize the convergence domain rather than claim universal robustness.

## Proposed error thresholds

- Normalized component and charge residuals < 1e-9 for normal synthetic cases.
- Log mass-action residuals < 1e-8 (must be explicitly measured by a follow-up automated report).
- Independent reference pH difference < 0.02 only when equations and constants are identical and both models are valid.
- Independent species concentrations within 1% for non-trace species, with documented exceptions near analytical detection limits.

These are **proposed targets, not achieved independent validation results**.

Do not label the program design-grade or nuclear-grade without an appropriate software quality program and scientific validation.
