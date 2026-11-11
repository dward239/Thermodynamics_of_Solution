#pragma once
#include <map>
#include <string>
#include <vector>

namespace thermo {
enum class ActivityModel { Ideal, Davies };
enum class Closure { Totals, Charge, FixedPH };

struct Basis {
    std::string name;
    int charge = 0;
    double totalMolal = 0.0;
    bool adjustable = false;
};
struct Formation {
    std::string name;
    int charge = 0;
    double logK25 = 0.0;
    double deltaHkJmol = 0.0;
    std::map<std::string, int> powers; // signed exponents of primary activities
};
struct Model {
    double temperatureC = 25.0;
    ActivityModel activity = ActivityModel::Ideal;
    Closure closure = Closure::Totals;
    std::string adjustedBasis;
    double fixedPH = 7.0;
    std::vector<Basis> basis;
    std::vector<Formation> formations;
};
struct Species {
    std::string name;
    int charge = 0;
    std::vector<int> powers;
    double logK = 0.0;
};
struct Prepared {
    Model model;
    std::vector<Species> species; // primary followed by secondary
    int adjustedIndex = -1;
    int hydrogenIndex = -1;
};
Model readModel(const std::string& filename);
Prepared prepare(Model model);
double temperatureLogK(double logK25, double deltaHkJmol, double temperatureC);
} // namespace thermo
