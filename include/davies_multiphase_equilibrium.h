#pragma once
#include "ideal_multiphase_equilibrium.h"
#include <map>
#include <string>
namespace thermo {

// Optional approximate Davies aqueous activity correction for the existing
// restricted, single-electroneutral-phase ideal-multiphase formulation.
// Uses aqueous_engine's historical molarity-based Davies approximation.
// No real molal activities, solvent-water activity, Pitzer, charge transfer
// between phases, or charged-solid equilibrium.
struct DaviesMultiphaseOptions {
    double temperatureK=298.15;
    int maxOuterIterations=100;
    int maxInnerIterations=150;
    double ionicStrengthRelativeTolerance=1e-9;
    double innerRelativeTolerance=1e-10;
    double relaxation=0.5; // 0 < relaxation <= 1, damp ionic strength updates
};
struct DaviesMultiphaseResult {
    IdealEquilibriumResult equilibrium;
    double ionicStrengthMolL=0.0;
    double relativeIonicStrengthResidual=0.0;
    int outerIterations=0;
    std::map<std::string,double> activityCoefficients;
    // Dimensionless a_s = gamma_s*(c_s / phase standard concentration).
    std::map<std::string,double> activities;
};
DaviesMultiphaseResult solveDaviesMultiphase(
    const IdealMultiphaseModel& model,
    const DaviesMultiphaseOptions& options={});
} // namespace thermo
