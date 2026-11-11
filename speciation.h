#pragma once
#include "species.h"
#include <map>
#include <string>
#include <vector>

struct TotalInput {
    std::map<std::string, double> totals;
};

// Placeholder for future validated mass-action / conservation solver.
// Throws until an explicit set of mass-action reactions and independent
// balance equations is supplied.
void solveSpeciation(std::vector<AqueousSpecies>& species,
    Solution& solution, const TotalInput& input);
