#include "pitzer.h"
#include "solver.h"
#include "species.h"
#include "speciation.h"
#include <cmath>
#include <exception>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    try {
        if (argc >= 3 && std::string(argv[1]) == "ideal") {
            const auto problem = loadIdealEquilibriumProblem(argv[2]);
            const auto result = solveIdealEquilibrium(problem);
            std::cout << "Ideal-solution mass-action equilibrium (illustrative model)\n";
            std::cout << "Newton iterations: " << result.iterations << "\n";
            for (const auto& [name, concentration] : result.concentrations)
                std::cout << "  " << name << ": " << concentration << " mol/L\n";
            for (const auto& [name, residual] : result.balanceResiduals)
                std::cout << "  Mass residual (" << name << "): " << residual << " mol/L\n";
            std::cout << "  Charge residual (diagnostic only): "
                      << result.chargeResidual << " mol charge/L\n";
            return 0;
        }
        if (argc < 2) {
            std::cout << "Thermodynamics of Solution — experimental engineering software\n"
                      << "  thermo ideal examples/association.txt\n"
                      << "  thermo <legacy-species-database> [pitzer-parameters]\n"
                      << "This is not a validated process or safety analysis code.\n";
            return 0;
        }
        auto species = loadSpecies(argv[1]);
        if (species.empty()) {
            std::cerr << "No species loaded. Supply a populated legacy database.\n";
            return 1;
        }
        std::cout << "Loaded " << species.size() << " legacy species.\n";
        if (argc >= 3) {
            auto parameters = loadPitzerParameters(argv[2]);
            std::cout << "Loaded " << parameters.size() << " illustrative interaction pairs.\n";
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
