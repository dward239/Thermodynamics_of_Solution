#include "pitzer.h"
#include "solver.h"
#include "species.h"
#include <exception>
#include <iostream>
#include <map>

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            std::cout << "Thermodynamics of Solution: experimental numerical components\n"
                      << "Usage: thermo <species-database> [pitzer-parameters]\n"
                      << "An equilibrium/speciation solver is not yet implemented.\n";
            return 0;
        }
        auto species = loadSpecies(argv[1]);
        if (species.empty()) {
            std::cerr << "No species loaded. Supply a populated database.\n";
            return 1;
        }
        std::cout << "Loaded " << species.size() << " species.\n";
        if (argc >= 3) {
            auto parameters = loadPitzerParameters(argv[2]);
            std::cout << "Loaded " << parameters.size() << " interaction pairs.\n";
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
