# User-defined restricted vessel calculation

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

CSV columns use seconds, kelvin, pascal, mole, joules and moles.
The CSV's 'gross products' column adds gross formation amounts across
configured species and phases for diagnostics only; it is not a net chemical
inventory or headspace gas source. Net gas transfer is fixed to ZERO in this
CLI; pressure changes only due to temperature. No liquid or solid chemical
equilibrium or solution phase redistribution is solved by this executable.

The example is a deliberately nonphysical mathematical test, not a radionuclide
or radiolysis predictive case. The CLI is unsuitable for safety decisions.
