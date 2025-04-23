#include <iostream>
#include "species.h"
#include "pitzer.h"
#include "solver.h"
#include "speciation.h"
int main() {
   // Load species from database file
   std::vector<AqueousSpecies> species = loadSpecies("database.txt");
   // Mark primary species manually
   for (auto& sp : species) {
       if (sp.name == "UO2++" || sp.name == "F-" || sp.name == "H+") {
           sp.isPrimary = true;
       }
   }
   // Define solution and total input
   Solution solution;
   TotalInput input;
   input.totals["U"] = 1e-3;  // 1 mM total uranium
   input.totals["F"] = 2e-3;  // 2 mM total fluoride
   input.totals["H"] = 0.0;   // H is not mass-balanced; only used for charge
   // Initial guesses for free species
   solution.concentrations["UO2++"] = 5e-4;
   solution.concentrations["F-"] = 1e-3;
   solution.concentrations["H+"] = 1e-3;  // pH ≈ 3
   // Calculate activity coefficients
   calculateActivityCoefficients(species, solution);
   // Solve for speciation and pH
   solveSpeciation(species, solution, input);
   return 0;
}
