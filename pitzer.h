#pragma once
#include "species.h"
#include <map>
#include <string>
#include <vector>

struct PitzerPair {
    std::string ion1;
    std::string ion2;
    double beta0 = 0.0;
    double beta1 = 0.0;
    double Cphi = 0.0;
};

std::vector<PitzerPair> loadPitzerParameters(const std::string& filename);
double calculateIonicStrength(const Solution& solution,
    const std::map<std::string, double>& charges);

// Experimental activity-coefficient approximation; NOT a validated full Pitzer model.
void calculateActivityCoefficients(std::vector<AqueousSpecies>& species,
    Solution& solution, const std::vector<PitzerPair>& parameters);
