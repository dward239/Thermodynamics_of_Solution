#pragma once
#include <map>
#include <string>
#include <vector>

namespace aqueous {

enum class ActivityModel { Ideal, Davies };
struct Component {
    std::string name;                 // Basis aqueous species, e.g., H+, CO2*, Na+
    int charge = 0;
    bool fixedTotal = true;           // Exactly one free component activates electroneutrality
    double totalMolL = 0.0;
};
struct Complex {
    std::string name;
    int charge = 0;
    double logKReference = 0.0;       // dimensionless, at referenceTemperatureK
    double deltaHkJMol = 0.0;         // assumed constant enthalpy for formation
    std::map<std::string, int> basisStoichiometry; // signed coefficients allowed
};
struct Solid {
    std::string name;
    double logKspReference = 0.0; // dissolution ion activity product at reference T
    double deltaHkJMol = 0.0;    // constant dissolution enthalpy assumption
    std::map<std::string, int> dissolutionStoichiometry; // positive basis species coefficients
};
struct Model {
    double temperatureK = 298.15;
    double referenceTemperatureK = 298.15;
    ActivityModel activity = ActivityModel::Ideal;
    std::vector<Component> components;
    std::vector<Complex> complexes;
    std::vector<Solid> solids;       // supported: at most one ideal, neutral solid
};
struct Options {
    int maxIterations = 150;
    double relativeTolerance = 1e-9;
};
struct Result {
    std::map<std::string, double> concentrationMolL;
    std::map<std::string, double> activity;
    std::map<std::string, double> componentResidualMolL;
    std::map<std::string, double> solidPrecipitatedMolL; // moles solid / liter of original solution
    std::map<std::string, double> saturationIndex;      // log10(IAP/Ksp)
    double chargeResidualMolL = 0.0;
    double ionicStrengthMolL = 0.0;
    double pH = 0.0;
    bool hasPH = false;
    int iterations = 0;
    double normalizedResidual = 0.0;
    double temperatureK = 298.15;
    ActivityModel model = ActivityModel::Ideal;
};
Model readModel(const std::string& filename);
Result solve(const Model& model, const Options& options = {});
double logKAtTemperature(const Complex& reaction, double temperatureK, double referenceTemperatureK);
double daviesLogGamma(int charge, double ionicStrengthMolL, double temperatureK);
std::string toJson(const Result& result);
void writeCsv(const Result& result, const std::string& filename);
}
