#include "pitzer.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

std::vector<PitzerPair> loadPitzerParameters(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open Pitzer parameter file: " + filename);
    std::vector<PitzerPair> result;
    std::string line;
    int lineNumber = 0;
    while (std::getline(file, line)) {
        ++lineNumber;
        auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream stream(line);
        PitzerPair p;
        if (!(stream >> p.ion1 >> p.ion2 >> p.beta0 >> p.beta1 >> p.Cphi))
            throw std::runtime_error("Invalid parameter row " + std::to_string(lineNumber));
        if (!std::isfinite(p.beta0) || !std::isfinite(p.beta1) || !std::isfinite(p.Cphi))
            throw std::runtime_error("Nonfinite parameter at row " + std::to_string(lineNumber));
        std::string extra;
        if (stream >> extra && extra[0] != '#')
            throw std::runtime_error("Extra field on row " + std::to_string(lineNumber));
        result.push_back(p);
    }
    return result;
}

double calculateIonicStrength(const Solution& solution,
    const std::map<std::string, double>& charges) {
    double strength = 0.0;
    for (const auto& [name, concentration] : solution.concentrations) {
        if (!std::isfinite(concentration) || concentration < 0.0)
            throw std::invalid_argument("Invalid concentration for " + name);
        auto charge = charges.find(name);
        if (charge == charges.end())
            throw std::invalid_argument("Missing charge for " + name);
        strength += 0.5 * concentration * charge->second * charge->second;
    }
    return strength;
}

void calculateActivityCoefficients(std::vector<AqueousSpecies>& species,
    Solution& solution, const std::vector<PitzerPair>& parameters) {
    std::map<std::string, double> charges;
    for (const auto& s : species) charges[s.name] = s.charge;
    solution.ionicStrength = calculateIonicStrength(solution, charges);
    const double rootI = std::sqrt(solution.ionicStrength);
    for (auto& s : species) {
        // Illustrative approximation, deliberately not labeled a complete Pitzer formulation.
        double lnGamma = -0.392 * s.charge * s.charge * rootI / (1.0 + 1.2 * rootI);
        for (const auto& p : parameters) {
            const std::string* other = nullptr;
            if (p.ion1 == s.name) other = &p.ion2;
            else if (p.ion2 == s.name) other = &p.ion1;
            if (!other) continue;
            auto it = solution.concentrations.find(*other);
            if (it == solution.concentrations.end()) continue;
            lnGamma += it->second * (p.beta0 + p.beta1 * std::exp(-2.0 * rootI)
                + s.charge * charges.at(*other) * p.Cphi * solution.ionicStrength);
        }
        s.activityCoefficient = std::exp(lnGamma);
    }
}
