#include "species.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <set>

namespace {
std::ifstream openInput(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Unable to open data file: " + path);
    return input;
}
bool skipLine(const std::string& line) {
    const auto pos = line.find_first_not_of(" \t\r");
    return pos == std::string::npos || line[pos] == '#';
}
void requireEnd(std::istringstream& in, const std::string& filename, int lineNumber) {
    std::string trailing;
    if (in >> trailing && trailing[0] != '#') {
        throw std::runtime_error(filename + ":" + std::to_string(lineNumber) + ": unexpected trailing field");
    }
}
}

std::vector<AqueousSpecies> loadSpecies(const std::string& filename) {
    auto input = openInput(filename);
    std::vector<AqueousSpecies> result;
    std::set<std::string> names;
    std::string line;
    int lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (skipLine(line)) continue;
        std::istringstream stream(line);
        AqueousSpecies species{};
        int n = 0;
        if (!(stream >> species.name >> species.charge >> species.logK >> n) ||
            n < 0 || !std::isfinite(species.charge) || !std::isfinite(species.logK))
            throw std::runtime_error(filename + ":" + std::to_string(lineNumber) + ": invalid aqueous species header");
        for (int i = 0; i < n; ++i) {
            std::string component;
            int coefficient;
            if (!(stream >> component >> coefficient) || coefficient == 0 ||
                !species.stoichiometry.emplace(component, coefficient).second)
                throw std::runtime_error(filename + ":" + std::to_string(lineNumber) + ": invalid stoichiometry");
        }
        requireEnd(stream, filename, lineNumber);
        if (!names.insert(species.name).second)
            throw std::runtime_error(filename + ":" + std::to_string(lineNumber) + ": duplicate species");
        result.push_back(std::move(species));
    }
    return result;
}

std::vector<SolidPhase> loadSolids(const std::string& filename) {
    auto input = openInput(filename);
    std::vector<SolidPhase> result;
    std::set<std::string> names;
    std::string line;
    int lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        if (skipLine(line)) continue;
        std::istringstream stream(line);
        SolidPhase solid{};
        int n = 0;
        if (!(stream >> solid.name >> solid.logKsp >> n) ||
            n < 0 || !std::isfinite(solid.logKsp))
            throw std::runtime_error(filename + ":" + std::to_string(lineNumber) + ": invalid solid header");
        for (int i = 0; i < n; ++i) {
            std::string component;
            int coefficient;
            if (!(stream >> component >> coefficient) || coefficient == 0 ||
                !solid.stoichiometry.emplace(component, coefficient).second)
                throw std::runtime_error(filename + ":" + std::to_string(lineNumber) + ": invalid solid stoichiometry");
        }
        requireEnd(stream, filename, lineNumber);
        if (!names.insert(solid.name).second)
            throw std::runtime_error(filename + ":" + std::to_string(lineNumber) + ": duplicate solid");
        result.push_back(std::move(solid));
    }
    return result;
}
