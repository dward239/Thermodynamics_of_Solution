#include "pitzer.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
// Load Pitzer binary interaction parameters from a file
std::vector<PitzerPair> loadPitzerParameters() {
    std::vector<PitzerPair> params;
    std::ifstream file("pitzer_params.txt");
    std::string line;
    if (!file.is_open()) {
        std::cerr << "Error: Could not open pitzer_params.txt\n";
        return params;
    }
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        PitzerPair p;
        iss >> p.ion1 >> p.ion2 >> p.beta0 >> p.beta1 >> p.Cphi;
        params.push_back(p);
    }
    return params;
}
// Calculate ionic strength I = 0.5 * sum(ci * zi^2)
double calculateIonicStrength(const Solution& solution, const std::map<std::string, double>& charges) {
    double I = 0.0;
    for (const auto& pair : solution.concentrations) {
        const std::string& ion = pair.first;
        double conc = pair.second;
        auto it = charges.find(ion);
        if (it != charges.end()) {
            double z = it->second;
            I += conc * z * z;
        }
    }
    return 0.5 * I;
}
// Calculate the Debye-Hückel term
double debyeHuckelTerm(double z, double I) {
    const double A_phi = 0.392;  // at 25 °C
    const double b = 1.2;
    return -A_phi * z * z * sqrt(I) / (1 + b * sqrt(I));
}
// Apply the Pitzer model to compute activity coefficients
void calculateActivityCoefficients(std::vector<AqueousSpecies>& species, Solution& solution) {
    auto pitzerParams = loadPitzerParameters();
    // Map charges for species
    std::map<std::string, double> charges;
    for (const auto& sp : species) {
        charges[sp.name] = sp.charge;
    }
    // Calculate ionic strength
    solution.ionicStrength = calculateIonicStrength(solution, charges);
    std::cout << "Ionic strength (Pitzer): " << solution.ionicStrength << "\n";
    // Precompute ln(gamma) for each species
    for (auto& sp : species) {
        double z = sp.charge;
        double I = solution.ionicStrength;
        // Debye-Hückel term
        double ln_gamma = debyeHuckelTerm(z, I);
        // Pitzer binary terms
        for (const auto& param : pitzerParams) {
            bool isMatch = false;
            std::string otherIon;
            double beta0 = param.beta0;
            double beta1 = param.beta1;
            double Cphi = param.Cphi;
            if (param.ion1 == sp.name) {
                otherIon = param.ion2;
                isMatch = true;
            } else if (param.ion2 == sp.name) {
                otherIon = param.ion1;
                isMatch = true;
            }
            if (isMatch && solution.concentrations.count(otherIon)) {
                double m_j = solution.concentrations.at(otherIon);
                ln_gamma += m_j * (beta0 + beta1 * exp(-2.0 * sqrt(I)) + z * charges[otherIon] * Cphi * I);
            }
        }
        sp.activityCoefficient = exp(ln_gamma);
    }
}
