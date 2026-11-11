# Restricted Davies activity coupling to general fluid mass action

This is an opt-in **approximate** nonideal activity model for the
generalized, single-charge-balanced aqueous phase in a multiphase model.
It uses `aqueous::daviesLogGamma` from the existing aqueous engine; the
new library links to `aqueous_engine` and `ideal_multiphase_equilibrium`.

## Formulation

At a trial ionic strength I, charged aqueous species have

  mu_s/(RT) = mu_s^0/(RT) + ln(gamma_s(I)) + ln(c_s/c^0_p)

where gamma is the molarity-based Davies estimate for 0–50 C and
I <= 0.1 mol/L. Species standard chemical potentials are shifted by
ln(gamma), and the generalized mass-action solver then enforces
elemental/component conservation and one designated aqueous-phase
electroneutrality equation. The new ion concentrations yield

  I_computed = 0.5 * sum(charged aqueous species z_s^2 * c_s, mol/L).

The outer Picard iteration damps I updates and stops when the computed
I agrees with the trial I within tolerance; otherwise it fails clearly.
Neutral organic/gas species retain gamma=1. This is NOT a gas fugacity,
organic-liquid activity, or true phase-specific multicomponent activity model.

The output includes species gamma, dimensionless activities, charge
balance, component residuals and inferred free counterion total.

## Verification

- Analytic NaCl 1:1 inventory and charge balance at 0.01 mol/L.
- Davies gamma at 25 C around 0.9018, matching the existing qualitative
  PyPitzer reference at equal numeric concentration (not molality conversion).
- Nonideal ion pairing produces the expected activity-product constant.
- Cross-solver regression against the historical coupled Davies
  aqueous equilibrium implementation with identical synthetic data.
- Neutral solute partitioned into aqueous, organic and gas phases.
- Error handling for high ionic strength, unsupported temperature,
  and nonsensical convergence parameters.

## Scientific limitations

This formulation is **not a validated molality-based electrolyte model**.
Davies coefficients use the existing approximate molarity convention, and
solvent chemical potential, ion-specific interactions, temperature-dependent
equilibrium constants, pressure dependence, electrified liquid interfaces,
multiple independently electroneutral liquids, charged solids, gas mixtures
and precipitation are not coupled here.

Temperature is supplied to the Davies coefficient only; the base species
standard chemical potentials are fixed unless the user constructs a new
input model using traceable temperature-dependent thermodynamic data.
Radiolysis/decay and full energy balance remain separate modules.

These are mathematical verification and cross-implementation tests,
not physical or regulatory qualification. Do not extrapolate beyond
the domain guards or substitute for a Pitzer/PHREEQC reference solution.
