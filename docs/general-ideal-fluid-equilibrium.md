# General ideal-dilute multiphase equilibrium prototype

This module solves a thermodynamic-equilibrium **subclass** with arbitrary
numbers of neutral conserved components, neutral fluid species and
fixed-volume liquid or gas phases. Every species has user-supplied
dimensionless standard chemical potential g_s = mu_s^0/(RT).

For each species s with stoichiometric chemical-component counts nu_is:

  n_s = V_phase * c_standard_phase
        * exp( sum_i nu_is * lambda_i/(RT) - g_s )

and all component balances are solved simultaneously:

  b_i = sum_s nu_is * n_s.

A damped Newton method with analytic Jacobian and pivoted linear system
solves for the component chemical potentials. The resulting concentrations
and component residuals are returned, with explicit failure for singular
stoichiometric Jacobians and nonconvergent cases.

The formulation is compatible with ideal dilute species chemical
potentials mu = mu^0 + RT ln(c/c_standard). Liquid reference concentration
could be 1000 mol/m3 (=1 mol/L). For an ideal gas at temperature T, choose
c_standard = p_standard/(R*T) and specify consistent standard mu values.

By choosing standard chemical potentials, the model can handle neutral
association reactions and phase partitioning **simultaneously** with
arbitrary populations of fixed-volume fluid phases.

### Exclusions / scientific qualification

- All modeled species must be neutral. No electrolyte electroneutrality
  constraints or activity-coefficient models have been integrated.
- Solids, liquid creation/disappearance, pressure-dependent liquid
  properties and nonideal mixtures are NOT modeled.
- Chemical potential data and standard states are supplied externally;
  no thermodynamic database is bundled or validated for this module.
- No water/solvent material balance, energy equation or time integration
  occurs in this equilibrium solve.
- The Newton iteration is currently dense and intended for small to
  medium well-posed systems, not a claimed industrial-scale solver.

Tests include two-liquid phase partitioning, analytically solvable
1:1 association with a chosen equilibrium constant, three-phase
partitioning and rank/validation safeguards. These are mathematical
verification tests, not experimental chemical validation.
