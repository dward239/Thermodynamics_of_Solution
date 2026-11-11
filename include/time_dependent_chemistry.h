#pragma once
#include "adaptive_vessel_integration.h"
#include "ideal_multiphase_equilibrium.h"
#include <map>
#include <string>
#include <vector>

namespace thermo {

// One-way operator-split re-equilibration of nuclear-element inventories
// AFTER each accepted adaptive vessel time step. Chemistry does not feed back
// into decay rates, deposited dose, H2 kinetics, temperature or pressure.
// Neutral, ideal-dilute fluid species and optional <=8 pure-solid candidates
// only. No charged/nonideal phase chemistry, DAE or isotope fractionation.
struct TimeDependentChemistryModel {
    IdealMultiphaseModel equilibriumTemplate;
    // Non-radioactive, time-invariant component inventory (mol).
    // Nuclide atoms add to these totals after conversion using Avogadro.
    std::map<std::string,double> backgroundComponentMol;
    std::vector<IdealPureSolid> candidateSolids;
    int maximumEquilibriumIterations=150;
    double equilibriumRelativeTolerance=1e-10;
    double conservationRelativeTolerance=1e-9;
    double conservationAbsoluteToleranceMol=1e-12;
};

struct TimeDependentChemistrySample {
    AdaptiveVesselSample vessel;
    std::map<std::string,double> inputComponentMol;
    IdealEquilibriumResult fluid;
    std::map<std::string,double> solidMoles;
    // Maximum component-closure error, divided by (atol+rtol*input total).
    double maximumScaledComponentClosure=0.0;
};

struct TimeDependentChemistryTrajectory {
    AdaptiveVesselTrajectory vessel;
    // One chemistry projection for EVERY accepted vessel state, including t0.
    std::vector<TimeDependentChemistrySample> samples;
};

// Run the existing adaptive radiolytic-vessel controller, then project every
// accepted time sample onto chemical equilibrium with changing nuclide-derived
// element totals. An exception prevents returning a partially qualified run.
// The accepted time grid is based on vessel dynamics alone: chemical
// projection errors do not affect step acceptance or adaptivity.
TimeDependentChemistryTrajectory integrateVesselWithChemicalProjection(
    const RadiolyticHydrogenVesselState& initial,
    double endTimeSeconds,
    const DecayNetwork& network,
    const RadiolyticHydrogenVesselModel& vesselModel,
    const TimeDependentChemistryModel& chemistry,
    const AdaptiveVesselOptions& integrationOptions={});

} // namespace thermo
