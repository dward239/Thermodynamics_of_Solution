# Adaptive time integration — restricted radiolytic vessel

This is a **bounded step-doubling controller** over the existing zero-dimensional
decay, hydrogen production/consumption/transfer, thermal and headspace models.
It is not a general multiphase DAE/ODE solver and does not call the
independent chemical-equilibrium engines.

## Step acceptance and rollback

Given the last accepted state at time `t` and a trial duration `dt`:

1. Evaluate a full trial step from the accepted state.
2. Independently evaluate one half step from the **same** accepted state.
3. Evaluate a second half step from the first half-step state.
4. Compare both candidates across nuclide counts, H2 inventories and
   consumption, temperature and gas pressure.
5. When normalized error `<=1`, audit nucleus/H2 bookkeeping and
   commit **only** the two-half-step state.
6. When error exceeds tolerance, reject without mutating the accepted
   state, reduce dt and retry. Excessive retries and minimum-dt
   failure throw explicitly instead of returning a false success.

The method evaluates a first-order backward-Euler H2 kinetic propagator.
Nuclear decay is independently propagated using uniformization.
Adaptive step changes use an exponent `-1/2` for the first-order H2
local error estimate.

Each accepted history sample includes an absolute timestamp,
last accepted step length, estimated local error, nuclear inventory,
H2 dissolved/gas/consumed inventories, temperature, pressure,
interval emitted/deposited energy, gross radiolytic H2, and conservation
residuals. The test suite covers analytical transfer behavior,
half-life and daughter ingrowth, source-integrated heating, pressure,
rejected-step rollback and failure budgets.

## Model restrictions

- The controller presently operates **only** on the restricted H2
  radiolytic vessel. Other radiolysis products are not coupled.
- User-prescribed deposited-energy fractions and G(H2) are assumed;
  radiation transport and radical-kinetics chemistry are not calculated.
- The first-order H2 consumption sink is an **unmodeled products pool**.
  The full atomic and reaction-enthalpy balances are not closed.
- A single constant-property lumped temperature and rigid ideal-gas
  headspace are assumed; no phase-volume response.
- Other chemical species and phase equilibria are **not recalculated
  at each accepted step**. This is not multiphase Gibbs equilibrium/DAE.
- Local error tolerance measures the internal discretization, not
  physical-model uncertainty or nuclear data accuracy.
- General purpose SUNDIALS/IDA coupling and independently validated
  data are future milestones.
- No end-user adaptive input-format/CLI parser yet; the interface
  is C++17, exercised by CTest.

**Scientific status:** mathematical and conservation verification of
restricted synthetic cases, not qualified for safety, regulatory,
or real radioisotope material prediction.

See the [SysML-style diagram index](../architecture/sysml-v2/views/README.md)
and [requirements mapping](../architecture/sysml-v2/requirements-traceability.md).
