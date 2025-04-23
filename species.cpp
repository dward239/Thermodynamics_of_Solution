#include "species.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

// Load aqueous species from a simple text database
std::vector<AqueousSpecies> loadSpecies(const std::string& filename) {
   std::vector<AqueousSpecies> speciesList;
   std::ifstream file(filename);
   std::string line;
   if (!file.is_open()) {
       std::cerr << "Error: Could not open species database file: " << filename << std::endl;
       return speciesList;
   }
   while (getline(file, line)) {
       if (line.empty() || line[0] == '#') continue; // Skip comments
       std::istringstream iss(line);
       AqueousSpecies sp;
       int num_components;
       iss >> sp.name >> sp.charge >> sp.logK >> num_components;
       for (int i = 0; i < num_components; ++i) {
           std::string component;
           int coeff;
           iss >> component >> coeff;
           sp.stoichiometry[component] = coeff;
       }
       speciesList.push_back(sp);
   }
   return speciesList;
}


// Load solid phases from a similar file format
std::vector<SolidPhase> loadSolids(const std::string& filename) {
   std::vector<SolidPhase> solids;
   std::ifstream file(filename);
   std::string line;
   if (!file.is_open()) {
       std::cerr << "Error: Could not open solids database file: " << filename << std::endl;
       return solids;
   }
   while (getline(file, line)) {
       if (line.empty() || line[0] == '#') continue;
       std::istringstream iss(line);
       SolidPhase sp;
       int num_components;
       iss >> sp.name >> sp.logKsp >> num_components;
       for (int i = 0; i < num_components; ++i) {
           std::string component;
           int coeff;
           iss >> component >> coeff;
           sp.stoichiometry[component] = coeff;
       }
       solids.push_back(sp);
   }
   return solids;
}
