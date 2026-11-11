# PAR-like view — integrated numerical contracts

This view reflects the SysML v1 **Parametric Diagram** concern:
equations and conserved quantities. It is a design-level representation,
not a parsed executable SysML constraint model.

```mermaid
flowchart LR
  Decay["Nuclide decay<br/>dN/dt = A N"] --> Emitted["Emitted energy [J]"]
  Emitted --> Deposit["Deposited energy by phase [J]"]
  Deposit --> Gross["Gross G(H2) source [mol]"]
  Gross --> H2Balance["H2 balance:<br/>L0 + G0 + produced = L1 + G1 + consumed"]
  Deposit --> Thermal["Thermal balance:<br/>C dT/dt = Q_in - UA(T - Tamb)"]
  H2Balance --> P["Ideal gas headspace:<br/>P V = n R T"]
  Thermal --> P
  Decay --> Nucleus["Nucleus count check:<br/>sum N(t+dt) = sum N(t)"]
```

Actual numerical acceptance tests are **scaled**: local trial error
`max_j |x_full,j - x_half,j| / (absTol_j + relTol*max(|x_full,j|, |x_half,j|)) <= 1`.
The controller audits H2 inventory plus the explicit unmodeled consumption
pool, and nuclide **nucleus counts** (not chemical elements, which change
on decay). It does not yet close the entire atom balance across all radiolytic
reaction products or the full thermochemical energy balance.
