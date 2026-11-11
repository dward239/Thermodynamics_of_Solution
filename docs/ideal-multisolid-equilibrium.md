# Multiple pure-solid phase equilibrium (experimental)

The neutral ideal-dilute fluid solver now has an optional **initially absent**
pure-solid candidate extension. For each solid j, user supplies stoichiometric
chemical component counts and dimensionless standard chemical potential
g_j = mu_j^0/(RT). The solver enforces component conservation across
all fluids and precipitated solids, plus thermodynamic complementarity:

  n_j >= 0;  a_j = sum_i nu_ij * lambda_i/(RT) - g_j <= 0;
  n_j * a_j = 0.

Present solids are at saturation (a_j=0); absent solids must be
undersaturated (a_j<=0). The solver enumerates possible active subsets
and performs damped Newton solves of the fluid component balances together
with the active solid saturation conditions. Only candidates satisfying
both inventory positivity and stability are accepted.

The search is exponential; a hard limit of **eight candidate solids**
prevents accidental unbounded runtime. Cases with degenerate solid
chemical potentials, rank deficiency, or poor nonlinear conditioning
may be rejected rather than resolved with a dedicated Gibbs minimizer.

With ideal pure A and B monomers and no neutral complexes, the standard
solid chemical potential g_AB = ln(Ksp) reproduces the analytic
precipitation condition [A][B]/(c_standard)^2 = Ksp.

The new tests cover precipitation, undersaturation, competition between
polymorphs and coexistence of two solids sharing a conserved element.

**Important restrictions**: neutral ideal-dilute fluid species only;
pure solids have no internal mixing, no seed/initial solid inventory or
dissolution kinetics. No nonideal electrolyte equilibrium, real chemical
constants, solid solutions, pressure effects, heat release or phase
equilibrium validation with PHREEQC. This is mathematical verification
and not a qualified radiochemistry predictor.
