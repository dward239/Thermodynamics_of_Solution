# Vessel input format

Build with CMake and run:

    ./build/vessel_cli examples/synthetic_vessel.txt synthetic_output.csv

Records (whitespace-delimited, '#' comments):

- nuclide ID ELEMENT HALF_LIFE_SECONDS EMITTED_J_PER_DECAY (zero half-life = stable)
- daughter PARENT_ID DAUGHTER_ID BRANCH_FRACTION (fractions sum to one)
- deposit NUCLIDE_ID PHASE_ID FRACTION (energy absorption; sum no more than one)
- atoms NUCLIDE_ID ATOM_COUNT
- gvalue PHASE_ID PRODUCT_ID MOLECULES_PER_100EV
- thermal MASS_KG CP_J_PER_KG_K UA_W_PER_K AMBIENT_K HEAT_FRACTION
- headspace VOLUME_M3 INITIAL_GAS_MOL INITIAL_TEMPERATURE_K
- simulate TIMESTEP_SECONDS COUNT

Arbitrary declared nuclide counts and branching networks are supported,
subject to the numerical uniformization stiffness limit. Radiation deposition
is user-specified. All source terms are generated from the synthetic or
externally supplied nuclear data, NOT independently looked up in a database.

For fixed-step runs, CSV columns use seconds, kelvin, pascal, moles and joules.
The CSV's 'gross products' column adds gross formation amounts across
configured species and phases for diagnostics only; it is not a net chemical
inventory or headspace gas source. Net gas transfer is fixed to ZERO in this
CLI; pressure changes only due to temperature. No liquid or solid chemical
equilibrium or solution phase redistribution is solved in fixed-step mode.

The example is a deliberately nonphysical mathematical test, not a radionuclide
or radiolysis predictive case. The CLI is unsuitable for safety decisions.

## Adaptive runs and chemical projection

Use the same executable with the bundled adaptive case:

```sh
./build/vessel_cli examples/adaptive_vessel_chemistry.txt adaptive_results.csv
```

Keep the existing nuclide, inventory, deposition, G-value, thermal and
headspace records. Replace `simulate` with `adaptive` and include one
`hydrogen` record. Time starts at zero.

| Record | Values in order |
| --- | --- |
| `adaptive` | End time (s), initial step (s), minimum step (s), maximum step (s), relative tolerance |
| `hydrogen` | Source phase ID, initial dissolved H₂ (mol), initial headspace H₂ (mol), transfer rate (1/s), dissolved consumption rate (1/s), equilibrium dissolved-to-gas mole ratio, kinetic substeps |
| `adaptive_tolerances` | Optional absolute tolerances for temperature (K), pressure (Pa), gas inventory (mol), and nuclide inventory (atoms). Defaults: `1e-6 1e-3 1e-12 1e-3`. |
| `adaptive_limits` | Optional accepted and rejected step limits, each 1–100000. Defaults: `100000 100000`. |

The designated source needs exactly one `gvalue PHASE H2 VALUE`
record. Additional positive G-values are outside the adaptive H₂ model.
The initial headspace gas amount includes any initial headspace H₂.

Add the following records to enable neutral, ideal-dilute chemical projection:

| Record | Values in order |
| --- | --- |
| `chem_component` | Component ID, fixed background inventory (mol). Include every nuclide element ID, even if initially absent. |
| `chem_phase` | Phase ID, kind (`aqueous`, `organic`, or `gas`), fixed volume (m³), standard concentration (mol/m³) |
| `chem_species` | Species ID, phase ID, standard chemical potential divided by RT, then one or more pairs of component ID and positive integer coefficient |
| `chem_solid` | Optional pure solid ID, standard chemical potential divided by RT, then the same component/coefficient pairs. Up to eight candidates. |
| `chem_temperature` | Optional reference, minimum, and maximum temperatures (K) |
| `chem_enthalpy` | Species or solid ID and standard molar enthalpy (J/mol). Required for every species and candidate solid when `chem_temperature` is supplied. |

Adaptive CSV rows include the initial state and each accepted state.
Columns include step diagnostics, gas and thermal histories, nuclide
counts, and optional component, species, and solid amounts. The adaptive
solver completes the trajectory before opening the output file. A failed
integration or chemical projection leaves an existing output untouched.
Chemistry uses vessel temperature as an input. The projected gas and
reaction heat do not feed back into the vessel.

The browser editor currently edits the fixed-step input format. Use the
[adaptive example](../examples/adaptive_vessel_chemistry.txt) as the
starting point for adaptive inputs.
