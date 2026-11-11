# REQ-like view — time integration and verification

Corresponds to SysML v1 **Requirement Diagram**. Candidate textual
requirement definitions are in
[TimeIntegrationRequirements.sysml](../TimeIntegrationRequirements.sysml);
test-level evidence and finite-precision criteria are in
[requirements-traceability.md](../requirements-traceability.md).

```mermaid
flowchart LR
  R20["«requirement» RT-020 Adaptive trial accuracy"] -->|verify| T1["adaptive_vessel_tests: rejected step and analytic transfer"]
  R21["«requirement» RT-021 Rollback"] -->|verify| T2["adaptive_vessel_tests: input immutability"]
  R22["«requirement» RT-022 Time horizon"] -->|verify| T3["adaptive_vessel_tests: end timestamp"]
  R23["«requirement» RT-023 Nuclear conservation"] -->|verify| T4["adaptive_vessel_tests: atoms closure"]
  R24["«requirement» RT-024 H2 conservation"] -->|verify| T5["adaptive_vessel_tests: H2 closure"]
  R25["«requirement» RT-025 Budget and failure reporting"] -->|verify| T6["adaptive_vessel_tests: step budgets"]
  I["«block» AdaptiveVesselIntegrator"] -->|satisfy intent| R20
  I --> R21
  I --> R22
  I --> R23
  I --> R24
  I --> R25
```

The diagram records software-test evidence, **not** experimental validation
of radioactive chemistry, thermodynamic models, or radiolysis yields.
