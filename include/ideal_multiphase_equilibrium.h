#pragma once
#include <map>
#include <string>
#include <vector>
namespace thermo {
// Restricted IDEAL-DILUTE multi-fluid equilibrium. Charged species can be
// solved using one specified electroneutral fluid phase and one free component.
// Species may inhabit distinct fixed-volume liquids or a gas headspace.
// Base solve excludes solids; optional bounded pure-solid extension below.
// Nonideal activities, charged interphase transfer, and phase creation unsupported.
struct IdealComponent {
    std::string id;
    double totalMol=0.0;
    bool fixedTotal=true; // false: infer total from electroneutrality, not mass balance
};
struct IdealPhase {
    std::string id;
    double volumeM3=0.0;
    double standardConcentrationMolM3=1000.0;
};
struct IdealSpecies {
    std::string id;
    std::string phaseId;
    std::map<std::string,int> componentStoichiometry;
    double standardMuOverRT=0.0;
    int charge=0; // charge in elementary-charge units; default neutral
};
struct IdealMultiphaseModel {
    std::vector<IdealComponent> components;
    std::vector<IdealPhase> phases;
    std::vector<IdealSpecies> species;
    // Empty: all species must be neutral. Otherwise exactly ONE component
    // must have fixedTotal=false; it is resolved using charge balance.
    // Charged species are restricted to this designated phase.
    std::string electroneutralPhaseId;
};
struct IdealEquilibriumResult {
    std::map<std::string,double> speciesMoles;
    std::map<std::string,double> concentrationsMolM3;
    std::map<std::string,double> componentResidualMol; // fixed totals only
    std::map<std::string,double> calculatedComponentMol; // also includes free total
    double chargeResidualMol=0.0; // signed molar amount of elementary charge
    // Component Lagrange multipliers divided by RT (natural-log standard).
    std::map<std::string,double> componentPotentialOverRT;
    int iterations=0;
    double maxNormalizedResidual=0.0;
};
IdealEquilibriumResult solveIdealMultiphase(
    const IdealMultiphaseModel& model,
    int maxIterations=150,
    double relativeTolerance=1e-10);
// Pure-solid candidates initially absent; no solid solutions or seed inventories.
// At equilibrium: present solid affinity ~0, absent solid affinity <=0.
struct IdealPureSolid {
    std::string id;
    std::map<std::string,int> componentStoichiometry;
    double standardMuOverRT=0.0;
};
struct IdealSolidEquilibriumResult {
    IdealEquilibriumResult fluidSolution;
    std::map<std::string,double> solidMoles;
    // Affinity = sum_i nu_is * lambda_i/(RT) - mu_s^0/(RT).
    // Positive affinity denotes supersaturation; equilibrium requires <=0.
    std::map<std::string,double> solidAffinityOverRT;
};
// Enumerates candidate active-solid subsets. Exponential in solids count,
// so <=8 candidates. Only neutral, ideal-dilute, fixed-volume models;
// charge-balanced input is rejected for the solid extension.
IdealSolidEquilibriumResult solveIdealMultiphaseWithSolids(
    const IdealMultiphaseModel& fluidModel,
    const std::vector<IdealPureSolid>& solidCandidates,
    int maxIterations=150,
    double relativeTolerance=1e-10);
} // namespace thermo
