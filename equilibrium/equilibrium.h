#pragma once
#include "chemistry/model.h"
#include <map>
#include <string>

namespace thermo {
struct SolverOptions {
    int maxIterations = 150;
    double tolerance = 1e-9;
};
struct Result {
    double temperatureC = 25.0;
    double ionicStrength = 0.0;    // mol/kg water
    double pH = 0.0;               // hydrogen activity (dimensionless, molal standard)
    double netCharge = 0.0;        // mol charge/kg water
    double maxScaledResidual = 0.0;
    int iterations = 0;
    std::map<std::string,double> molalities;
    std::map<std::string,double> activities;
    std::map<std::string,double> gammas;
    std::map<std::string,double> componentResiduals; // mol/kg
    std::map<std::string,double> calculatedTotals; // mol/kg
};
Result solve(const Prepared& prepared, const SolverOptions& options = {});
} // namespace thermo
