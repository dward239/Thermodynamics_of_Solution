#include "solver.h"
#include <cmath>
#include <limits>
#include <stdexcept>

double calculateSaturationIndex(const SolidPhase& solid,
    const std::map<std::string, double>& activities) {
    // SI = log10(IAP) - log10(Ksp). Avoid products that overflow/underflow.
    double logIap = 0.0;
    for (const auto& [name, coefficient] : solid.stoichiometry) {
        auto it = activities.find(name);
        if (it == activities.end()) {
            throw std::invalid_argument("Missing activity for " + name);
        }
        if (!std::isfinite(it->second) || it->second < 0.0) {
            throw std::invalid_argument("Invalid activity for " + name);
        }
        if (it->second == 0.0) {
            return -std::numeric_limits<double>::infinity();
        }
        logIap += coefficient * std::log10(it->second);
    }
    return logIap - solid.logKsp;
}

std::map<std::string, double> estimateIdealActivities(const Solution& solution) {
    // Activities approximated as concentrations (gamma = 1), for demonstration only.
    return solution.concentrations;
}
