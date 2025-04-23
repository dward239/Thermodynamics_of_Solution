#include "solver.h"
#include <iostream>
#include <cmath>
#include <limits>
// Calculate the saturation index of a solid given species activities
double calculateSaturationIndex(const SolidPhase& solid, const std::map<std::string, double>& activities) {
   double Q = 1.0;
   for (const auto& pair : solid.stoichiometry) {
       const std::string& species = pair.first;
       int coeff = pair.second;
       auto it = activities.find(species);
       if (it != activities.end() && it->second > 0) {
           Q *= pow(it->second, coeff);
       } else {
           // If species not found or has zero activity, Q becomes 0
           return -std::numeric_limits<double>::infinity();
       }
   }
   return log10(Q) - solid.logKsp;
}
// Fallback: activity = concentration (ideal case)
std::map<std::string, double> estimateActivities(const Solution& solution) {
   std::map<std::string, double> activities;
   for (const auto& pair : solution.concentrations) {
       activities[pair.first] = pair.second; // assume gamma = 1
   }
   return activities;
}
// Main phase stability and equilibrium function
void solveEquilibrium(const std::vector<AqueousSpecies>& species, Solution& solution) {
   std::cout << "\n=== Equilibrium and Phase Stability Check ===\n";
   // Estimate activities (placeholder)
   auto activities = estimateActivities(solution);
   // Display activities
   std::cout << "Estimated activities (assuming gamma = 1):\n";
   for (const auto& act : activities) {
       std::cout << "  " << act.first << ": " << act.second << "\n";
   }
   // Load solids and evaluate their saturation index
   auto solids = loadSolids("solids.txt");
   std::cout << "\nPhase stability (saturation indices):\n";
   for (const auto& solid : solids) {
       double SI = calculateSaturationIndex(solid, activities);
       std::cout << "  " << solid.name << "  SI = " << SI;
       if (SI > 0) std::cout << "  [Supersaturated → may precipitate]";
       else if (SI == 0) std::cout << "  [At equilibrium]";
       else std::cout << "  [Undersaturated]";
       std::cout << "\n";
   }
}
