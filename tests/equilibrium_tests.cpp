#include "speciation.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

void check(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}
void approx(double actual, double expected, double tolerance, const std::string& label) {
    check(std::abs(actual - expected) <= tolerance, label + " differs from analytical result");
}
int main() {
    // Case 1: A + B = AB, K=100, equal totals 0.01 mol/L.
    auto association = loadIdealEquilibriumProblem("examples/association.txt");
    auto result = solveIdealEquilibrium(association);
    const double T = 0.01, K = 100.0;
    const double free = 2 * T / (1 + std::sqrt(1 + 4 * K * T));
    approx(result.concentrations.at("A+"), free, 1e-11, "A free");
    approx(result.concentrations.at("B-"), free, 1e-11, "B free");
    approx(result.concentrations.at("AB"), K * free * free, 1e-11, "AB");
    check(std::abs(result.chargeResidual) < 1e-12, "Charge residual should vanish for balanced totals");
    check(std::abs(result.balanceResiduals.at("A+")) < 1e-11, "A mass balance");
    check(std::abs(result.balanceResiduals.at("B-")) < 1e-11, "B mass balance");

    // Case 2: 2A = A2, mass balance T = [A] + 2K[A]^2.
    auto dimer = loadIdealEquilibriumProblem("examples/dimer.txt");
    auto dim = solveIdealEquilibrium(dimer);
    const double T2 = 0.02, K2 = 10.0;
    const double freeA = 2 * T2 / (1 + std::sqrt(1 + 8 * K2 * T2));
    approx(dim.concentrations.at("A"), freeA, 1e-11, "A monomer");
    approx(dim.concentrations.at("A2"), K2 * freeA * freeA, 1e-11, "A dimer");

    // Reject duplicate names, unphysical totals, unknown components, charge mismatch.
    bool failed = false;
    try { auto bad = association; bad.components[0].total = 0; solveIdealEquilibrium(bad); }
    catch (const std::invalid_argument&) { failed = true; }
    check(failed, "Zero total must be rejected");
    failed = false;
    try { auto bad = association; bad.complexes[0].stoichiometry["missing"] = 1; solveIdealEquilibrium(bad); }
    catch (const std::invalid_argument&) { failed = true; }
    check(failed, "Unknown basis species must be rejected");
    failed = false;
    try { auto bad = association; bad.complexes[0].charge = 1; solveIdealEquilibrium(bad); }
    catch (const std::invalid_argument&) { failed = true; }
    check(failed, "Charge inconsistency must be rejected");
    failed = false;
    try { EquilibriumOptions o; o.maxIterations = 1; solveIdealEquilibrium(association,o); }
    catch (const std::runtime_error&) { failed = true; }
    check(failed, "Iteration limit must report nonconvergence");

    // Deliberately imbalanced electrolyte totals: solver returns diagnostic rather than hiding it.
    auto unequal = association;
    unequal.components[0].total = 0.02;
    auto nonneutral = solveIdealEquilibrium(unequal);
    approx(nonneutral.chargeResidual, 0.01, 1e-11, "Charge-balance diagnostic");
    std::cout << "Ideal mass-action analytical benchmarks passed\n";
}
