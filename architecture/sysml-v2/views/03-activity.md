# ACT-like view — adaptive step execution

This is the **SysML v1 Activity Diagram** analogue. The SysML v2 action
definitions are in [TimeIntegrationBehavior.sysml](../TimeIntegrationBehavior.sysml).

```mermaid
flowchart TD
  Start([Start]) --> Validate["Validate model, state, options and end time"]
  Validate --> Horizon{"Reached end time?"}
  Horizon -->|Yes| Done([Return trajectory])
  Horizon -->|No| Propose["Propose bounded trial dt"]
  Propose --> Full["Calculate full-step trial"]
  Full --> Half1["Calculate first half-step"]
  Half1 --> Half2["Calculate second half-step"]
  Half2 --> Error["Estimate normalized local error"]
  Error --> Valid{"Valid trial and error <= 1?"}
  Valid -->|Yes| Audit["Audit nuclide and H2 closure"]
  Audit --> AuditOk{"Audit passes?"}
  AuditOk -->|No| Fail([Fail explicitly])
  AuditOk -->|Yes| Commit["Commit two-half-step state"]
  Commit --> Record["Record accepted sample and diagnostics"]
  Record --> Grow["Select next dt"]
  Grow --> Horizon
  Valid -->|No| Reject["Log rejection without mutating accepted state"]
  Reject --> Floor{"Retry budget or min step reached?"}
  Floor -->|Yes| Fail
  Floor -->|No| Shrink["Reduce dt"]
  Shrink --> Propose
```

The accepted state is immutable during a trial. Errors from invalid
physical configuration fail immediately. Only local truncation error and
explicitly retryable numerical issues trigger step reduction.
