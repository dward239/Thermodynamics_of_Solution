#pragma once
#include "adaptive_vessel_integration.h"
#include "ideal_multiphase_equilibrium.h"
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace thermo {

struct ChemicalProjectionThermodynamics {
    double referenceTemperatureK=298.15;
    double minimumTemperatureK=273.15;
    double maximumTemperatureK=323.15;
    // One constant standard molar enthalpy per fluid species and pure solid
    // candidate, on the same component basis as the reference potentials.
    std::map<std::string,double> standardEnthalpyJPerMolById;
};

// One-way operator-split re-equilibration of nuclear-element inventories
// at each accepted vessel candidate, before commit. Chemistry does not feed back
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
    // If set, template potentials and gas standard concentrations refer to
    // referenceTemperatureK. Otherwise they are used unchanged at every step.
    std::optional<ChemicalProjectionThermodynamics> thermodynamics;
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

// Project t0 and each numerically accepted candidate onto equilibrium before
// committing the sample. A projection failure aborts the run. The vessel
// dynamics alone determine the time-step error estimate and retry decisions.
TimeDependentChemistryTrajectory integrateVesselWithChemicalProjection(
    const RadiolyticHydrogenVesselState& initial,
    double endTimeSeconds,
    const DecayNetwork& network,
    const RadiolyticHydrogenVesselModel& vesselModel,
    const TimeDependentChemistryModel& chemistry,
    const AdaptiveVesselOptions& integrationOptions={});

} // namespace thermo
