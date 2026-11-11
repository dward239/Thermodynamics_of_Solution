#pragma once
#include "hydrogen_kinetics.h"
#include "vessel_simulation.h"
#include <string>
namespace thermo {
// Restriction: only H2 produced in a designated liquid phase is coupled.
// Other gas-forming radiolytic products require separate kinetic models.
struct RadiolyticHydrogenVesselState {
    VesselState vessel;
    HydrogenKineticsState hydrogen;
};
struct RadiolyticHydrogenVesselModel {
    VesselModel vessel;
    HydrogenKineticsParameters kinetics;
    std::string liquidPhaseId="water";
};
struct RadiolyticHydrogenVesselStep {
    RadiolyticHydrogenVesselState state;
    VesselStepResult vesselStep;
    HydrogenKineticsStep hydrogenStep;
};
RadiolyticHydrogenVesselStep advanceRadiolyticHydrogenVessel(
    const RadiolyticHydrogenVesselState& initial,
    double durationSeconds,
    const DecayNetwork& network,
    const RadiolyticHydrogenVesselModel& model);
} // namespace thermo
