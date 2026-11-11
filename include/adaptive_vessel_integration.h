#pragma once
#include "radiolytic_hydrogen_vessel.h"
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace thermo {

// Adaptive, zero-dimensional controller for the RESTRICTED radiolytic-H2
// vessel kernel. It does not couple chemical equilibrium or solve a DAE.
struct AdaptiveVesselOptions {
    double initialStepSeconds = 1.0;
    double minimumStepSeconds = 1e-6;
    double maximumStepSeconds = 1e4;
    double relativeTolerance = 1e-4;
    double absoluteTemperatureToleranceK = 1e-6;
    double absolutePressureTolerancePa = 1e-3;
    double absoluteMolesTolerance = 1e-12;
    double absoluteAtomsTolerance = 1e-3;
    std::size_t maximumAcceptedSteps = 100000;
    std::size_t maximumRejectedSteps = 100000;
};

enum class AdaptiveStepDisposition { Accepted, Rejected };

struct AdaptiveVesselAttempt {
    double startTimeSeconds = 0.0;
    double proposedStepSeconds = 0.0;
    double normalizedError = 0.0; // <=1 satisfies configured local criterion
    AdaptiveStepDisposition disposition = AdaptiveStepDisposition::Rejected;
    std::string reason;
};

struct AdaptiveVesselSample {
    RadiolyticHydrogenVesselState state;
    double acceptedStepSeconds = 0.0;
    double pressurePa = 0.0;
    double estimatedNormalizedError = 0.0;
    double intervalEmittedEnergyJ = 0.0;
    double intervalDepositedEnergyJ = 0.0;
    double intervalGrossH2Mol = 0.0;
    double intervalConsumedH2Mol = 0.0;
    double intervalH2ClosureResidualMol = 0.0;
    double intervalNuclideAtomClosureResidual = 0.0;
};

struct AdaptiveVesselTrajectory {
    RadiolyticHydrogenVesselState finalState;
    std::vector<AdaptiveVesselSample> history;  // includes initial sample
    std::vector<AdaptiveVesselAttempt> attempts;
    std::size_t acceptedSteps = 0;
    std::size_t rejectedSteps = 0;
};

// Integrates to ABSOLUTE endTimeSeconds, not a duration.
// Each trial is computed from the same immutable accepted state: a full
// step and two half steps. Only the two-half-step solution may be accepted.
// A rejected trial never updates the accepted state. On failure, throws
// rather than returning an apparently completed partial trajectory.
// Backward-Euler H2 substeps remain the actual kinetic discretization.
// Optional acceptedObserver is called on the initial and each accepted sample
// BEFORE committing that sample; it is never called for rejected trials.
// Observer exceptions abort integration rather than silently accepting an
// unprojectable chemistry state.
AdaptiveVesselTrajectory integrateRadiolyticHydrogenVessel(
    const RadiolyticHydrogenVesselState& initial,
    double endTimeSeconds,
    const DecayNetwork& network,
    const RadiolyticHydrogenVesselModel& model,
    const AdaptiveVesselOptions& options = {},
    const std::function<void(const AdaptiveVesselSample&)>& acceptedObserver = {});

} // namespace thermo
