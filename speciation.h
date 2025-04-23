#ifndef SPECIATION_H
#define SPECIATION_H
#include "species.h"
#include <string>
#include <vector>
#include <map>
// Holds total input amounts for conserved elements
struct TotalInput {
   std::map<std::string, double> totals;  // e.g., {"U": 1e-3, "F": 2e-3}
};
// Main speciation solver interface
void solveSpeciation(std::vector<AqueousSpecies>& species,
                    Solution& solution,
                    const TotalInput& input);
#endif
