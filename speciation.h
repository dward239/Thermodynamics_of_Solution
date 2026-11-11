#pragma once
#include <map>
#include <string>
#include <vector>

// Formation equilibrium: complex = product(primary_i ^ stoichiometry_i).
// Under ideal activities, c_complex / c_standard = K * product(c_primary/c_standard)^nu.
struct PrimaryComponent {
    std::string name;
    double total = 0.0;   // mol/L, strictly positive in this solver
    int charge = 0;
};
struct SecondaryComplex {
    std::string name;
    int charge = 0;
    double log10K = 0.0;  // dimensionless formation constant at stated conditions
    std::map<std::string, int> stoichiometry;
};
struct IdealEquilibriumProblem {
    std::vector<PrimaryComponent> components;
    std::vector<SecondaryComplex> complexes;
};
struct EquilibriumOptions {
    int maxIterations = 100;
    double relativeTolerance = 1e-10;
};
struct EquilibriumResult {
    std::map<std::string, double> concentrations; // mol/L
    std::map<std::string, double> balanceResiduals; // mol/L
    double chargeResidual = 0.0; // mol/L of charge; DIAGNOSTIC ONLY
    int iterations = 0;
};

// File grammar is documented in README. Neither charge balance nor pH is independently
// imposed; component totals must specify a physically consistent electrolyte mixture.
IdealEquilibriumProblem loadIdealEquilibriumProblem(const std::string& filename);

// Damped Newton method in log10(primary concentration) coordinates.
// Throws on invalid/inconsistent model or nonconvergence; never returns partial "success".
EquilibriumResult solveIdealEquilibrium(const IdealEquilibriumProblem& problem,
                                        const EquilibriumOptions& options = {});
