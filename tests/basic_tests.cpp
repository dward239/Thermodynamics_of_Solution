#include "pitzer.h"
#include "solver.h"
#include <cassert>
#include <cmath>
#include <map>
#include <stdexcept>
#include <iostream>

int main() {
    Solution s;
    s.concentrations = {{"Na+", 0.10}, {"Cl-", 0.10}};
    const std::map<std::string, double> charges{{"Na+", 1.0}, {"Cl-", -1.0}};
    assert(std::abs(calculateIonicStrength(s, charges) - 0.10) < 1e-12);

    SolidPhase solid;
    solid.name = "AB";
    solid.logKsp = -4.0;
    solid.stoichiometry = {{"A+", 1}, {"B-", 1}};
    assert(std::abs(calculateSaturationIndex(solid, {{"A+", 0.01}, {"B-", 0.01}})) < 1e-12);
    assert(std::isinf(calculateSaturationIndex(solid, {{"A+", 0.0}, {"B-", 0.01}})));
    bool caught = false;
    try { calculateSaturationIndex(solid, {{"A+", 0.01}}); }
    catch (const std::invalid_argument&) { caught = true; }
    assert(caught);

    caught = false;
    s.concentrations["Na+"] = -1.0;
    try { calculateIonicStrength(s, charges); }
    catch (const std::invalid_argument&) { caught = true; }
    assert(caught);
    std::cout << "Numerical-component tests passed\n";
}
