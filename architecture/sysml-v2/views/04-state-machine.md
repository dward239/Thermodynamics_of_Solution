# STM-like view — adaptive controller lifecycle

The **SysML v1 State Machine Diagram** counterpart maps to
`state def AdaptiveTimeControllerStates` in
[TimeIntegrationBehavior.sysml](../TimeIntegrationBehavior.sysml).

```mermaid
stateDiagram-v2
  [*] --> Ready
  Ready --> Evaluating: start / validated
  Evaluating --> Auditing: trial computed and error <= tolerance
  Evaluating --> Adjusting: trial rejected or numerical retry
  Auditing --> Committed: conservation satisfied
  Auditing --> Failed: conservation failure
  Adjusting --> Evaluating: reduce dt, within retry budget
  Adjusting --> Failed: minimum dt or retry budget
  Committed --> Evaluating: remaining simulation time
  Committed --> Completed: reached requested end time
  Completed --> [*]
  Failed --> [*]
```

| State | Invariant / responsibility |
| --- | --- |
| Ready | Input configuration validated; initial state not modified |
| Evaluating | Full and two half-step candidates calculated from same accepted state |
| Auditing | Nuclide and H2 conservation checked |
| Adjusting | A rejected trial does **not** update accepted state |
| Committed | Two-half-step result becomes the accepted state |
| Completed | Exact requested end-time sample available |
| Failed | Throws; no falsely successful partial trajectory returned |

The state diagram is a behavioral **design view**. The C++ implementation
currently realizes these conditions through its control flow, not a runtime
state-machine engine.
