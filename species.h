#ifndef SPECIES_H
#define SPECIES_H
#include <string>
#include <map>
#include <vector>

// Represents an aqueous complex or ion (e.g., UO2F2, H+, C2O4--, etc.)
struct AqueousSpecies {
   std::string name;                          // Species name
   double charge;                             // Electrical charge
   double logK;                               // log10 equilibrium constant
   std::map<std::string, int> stoichiometry;  // {component -> coefficient}
   double activityCoefficient = 1.0;          // Calculated later via Pitzer
   double concentration = 0.0;                // mol/L (can be initialized later)
   bool isPrimary = false;
};
// Represents a solid phase with a solubility product
struct SolidPhase {
   std::string name;                          // Solid name (e.g., UO2(C2O4)·2H2O)
   double logKsp;                             // log10 solubility product
   std::map<std::string, int> stoichiometry;  // {ion -> coefficient}
};
// Represents the solution system: composition, state variables, etc.
struct Solution {
   std::map<std::string, double> concentrations;  // Free species concentrations (mol/L)
   double ionicStrength = 0.0;                    // Calculated from species charges and concentrations
   double temperature = 298.15;                   // K (default: 25°C)
   double pH = 0.0;                               // Solution pH
};
// Loads aqueous species from a database file
std::vector<AqueousSpecies> loadSpecies(const std::string& filename);
// (optional future) Load solid phases
std::vector<SolidPhase> loadSolids(const std::string& filename);
#endif
