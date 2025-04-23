#include "speciation.h"
#include <iostream>
#include <cmath>
#include <limits>
#include <set>
#include <map>
#include <Eigen/Dense>
// === Main Speciation Solver Interface ===
void solveSpeciation(std::vector<AqueousSpecies>& species,
                     Solution& solution,
                     const TotalInput& input) {
    std::cout << "\n=== Speciation Solver ===\n";
    std::cout << "Total element inputs (mol/L):\n";
    for (const auto& pair : input.totals) {
        std::cout << "  " << pair.first << ": " << pair.second << "\n";
    }
    auto stoichMatrix = buildStoichiometricMatrix(species);
    newtonSolve(species, solution, input, stoichMatrix);
    std::cout << "\nFinal species concentrations:\n";
    for (const auto& sp : species) {
        std::cout << "  " << sp.name << ": " << solution.concentrations[sp.name] << "\n";
    }
}
// === Element Detection from Species ===
std::set<std::string> detectUniqueElements(const std::vector<AqueousSpecies>& species) {
    std::set<std::string> elements;
    for (const auto& sp : species) {
        for (const auto& comp : sp.stoichiometry) {
            const std::string& label = comp.first;
            bool isElement = true;
            for (char c : label) {
                if (isdigit(c) || c == '+' || c == '-' || c == '(' || c == ')') {
                    isElement = false;
                    break;
                }
            }
            if (isElement) {
                elements.insert(label);
            }
        }
    }
    return elements;
}
// === Build Stoichiometric Matrix ===
std::map<std::string, std::map<std::string, int>>
buildStoichiometricMatrix(const std::vector<AqueousSpecies>& species) {
    std::map<std::string, std::map<std::string, int>> matrix;
    for (const auto& sp : species) {
        for (const auto& comp : sp.stoichiometry) {
            matrix[comp.first][sp.name] += comp.second;
        }
    }
    return matrix;
}
// === Compute Mass Balance Residuals ===
std::map<std::string, double>
computeMassBalanceResiduals(const std::vector<AqueousSpecies>& species,
                            const Solution& solution,
                            const TotalInput& input,
                            const std::map<std::string, std::map<std::string, int>>& stoichMatrix) {
    std::map<std::string, double> residuals;
    for (const auto& elementPair : stoichMatrix) {
        const std::string& element = elementPair.first;
        double total = input.totals.count(element) ? input.totals.at(element) : 0.0;
        double sum = 0.0;
        for (const auto& speciesPair : elementPair.second) {
            const std::string& speciesName = speciesPair.first;
            int coeff = speciesPair.second;
            auto it = solution.concentrations.find(speciesName);
            if (it != solution.concentrations.end()) {
                sum += coeff * it->second;
            }
        }
        residuals[element] = sum - total;
    }
    return residuals;
}
// === Newton-Raphson Solver ===
void newtonSolve(std::vector<AqueousSpecies>& species,
                Solution& solution,
                const TotalInput& input,
                const std::map<std::string, std::map<std::string, int>>& stoichMatrix) {
   const double tol = 1e-8;
   const int max_iter = 50;
   const double perturb = 1e-6;
   std::vector<std::string> speciesNames;
   for (const auto& sp : species) {
       if (sp.isPrimary) {
           speciesNames.push_back(sp.name);
           if (solution.concentrations.count(sp.name) == 0) {
               solution.concentrations[sp.name] = 1e-10;
           }
       }
   }
   for (int iter = 0; iter < max_iter; ++iter) {
       evaluateMassAction(species, solution);
       auto residuals = computeMassBalanceResiduals(species, solution, input, stoichMatrix);
       double chargeResidual = computeChargeBalance(species, solution);
       std::vector<std::string> elements;
       for (const auto& r : residuals) {
           elements.push_back(r.first);
       }
       elements.push_back("Charge");
       Eigen::VectorXd f(elements.size());
       for (size_t i = 0; i < elements.size() - 1; ++i) {
           f(i) = residuals[elements[i]];
       }
       f(elements.size() - 1) = chargeResidual;
       if (f.norm() < tol) {
           std::cout << "Converged in " << iter << " iterations.\n";
           return;
       }
       Eigen::MatrixXd J(elements.size(), speciesNames.size());
       for (size_t j = 0; j < speciesNames.size(); ++j) {
           std::string sj = speciesNames[j];
           double original = solution.concentrations[sj];
           double delta = std::max(original * 0.01, perturb);
           solution.concentrations[sj] += delta;
           evaluateMassAction(species, solution);
           auto perturbed = computeMassBalanceResiduals(species, solution, input, stoichMatrix);
           double pertCharge = computeChargeBalance(species, solution);
           solution.concentrations[sj] = original;
           evaluateMassAction(species, solution);
           for (size_t i = 0; i < elements.size() - 1; ++i) {
               double df = perturbed[elements[i]] - residuals[elements[i]];
               J(i, j) = df / delta;
           }
           J(elements.size() - 1, j) = (pertCharge - chargeResidual) / delta;
       }
       Eigen::VectorXd delta = J.colPivHouseholderQr().solve(-f);
       for (size_t j = 0; j < speciesNames.size(); ++j) {
           std::string sj = speciesNames[j];
           solution.concentrations[sj] += delta(j);
           if (solution.concentrations[sj] < 0) {
               solution.concentrations[sj] = 1e-12;
           }
       }
   }
   std::cout << "Newton-Raphson did not converge after " << max_iter << " iterations.\n";
}
// === Mass Action Law Evaluation ===
void evaluateMassAction(std::vector<AqueousSpecies>& species, const Solution& solution) {
    std::map<std::string, double> activities;
    for (const auto& sp : species) {
        if (sp.isPrimary) {
            double conc = solution.concentrations.count(sp.name) ? solution.concentrations.at(sp.name) : 1e-12;
            double act = conc * sp.activityCoefficient;
            activities[sp.name] = act;
        }
    }
    for (auto& sp : species) {
        if (!sp.isPrimary) {
            double log_conc = sp.logK;
            for (const auto& pair : sp.stoichiometry) {
                const std::string& component = pair.first;
                int coeff = pair.second;
                if (activities.count(component)) {
                    log_conc += coeff * std::log10(activities.at(component));
                } else {
                    log_conc += coeff * std::log10(1e-12);
                }
            }
            sp.concentration = std::pow(10.0, log_conc);
        }
        // Write back into solution
        solution.concentrations[sp.name] = sp.concentration;
    }
}
