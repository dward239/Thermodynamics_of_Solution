# Restricted charge balance in the generalized ideal-fluid kernel

## Supported problem

Optional `electroneutralPhaseId` designates ONE electroneutral fluid
phase. Charged species must belong to this phase. All other phases may
hold neutral species. Exactly one `IdealComponent` must be marked
`fixedTotal=false`, with `totalMol=0`: its total is **inferred** from
charge neutrality, rather than incorrectly imposing an extra independent
equation on a fully specified inventory.

The mathematical unknowns remain the N component potential multipliers.
Of the N equations, N-1 enforce fixed chemical-component conservation;
the free component's equation is replaced by

  sum_s z_s * n_s = 0,

where z_s is the integer charge and n_s the amount in the designated
phase (moles). The logarithmic ideal concentration expression remains

  n_s = V_p * c^0_p * exp(sum_i nu_is * x_i - mu_s^0/(RT)).

The Newton Jacobian uses z_s*nu_js for the charge-balance row and
nu_is*nu_js for component-conservation rows. All residuals and
Jacobian rows are scaled to protect numerical conditioning.

The output reports signed charge imbalance in mol of elementary charge,
the inferred totals for *every* component (including the free one),
and residuals for fixed components only.

## Boundaries

- IDEAL-DILUTE thermodynamics only. No ionic activity corrections,
  Born solvation, pressure dependence, or solvent chemical potentials.
- Only one charge-balanced phase, and no charged-species transfer to
  other phases. Multiple independently electroneutral liquids,
  Donnan potentials, electrochemical potential continuity, and
  interfacial ion transfer are NOT supported.
- Solid-phase enumeration remains restricted to neutral models.
- The free inventory is an inferred counterion reservoir, not a
  conserved specified physical total. The selected model is not a
  closed system with a separately conserved counterion inventory.
- Formal charge assignments and standard chemical potentials require
  externally consistent thermodynamic input. No electron/redox chemistry
  or charge-state-changing reaction bookkeeping is implemented.
- Mass-action and charge-balance tests verify numerical behavior only,
  not experimental accuracy or nuclear material suitability.

Use this extension as a data-model and solver verification checkpoint
before adding rigorous chemical activity models and multiple charged phases.
