#include "speciation.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {
constexpr double LN10 = 2.302585092994045684;
struct SpeciesEntry {
    std::string name;
    std::vector<int> nu;
    double log10K;
    int charge;
};
struct Evaluation {
    std::vector<double> values;
    std::vector<double> residuals;
    std::vector<std::vector<double>> jacobian;
    double norm = 0.0;
};
bool solveLinear(std::vector<std::vector<double>> matrix,
                 std::vector<double> rhs, std::vector<double>& step) {
    const std::size_t n = rhs.size();
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i)
            if (std::abs(matrix[i][k]) > std::abs(matrix[pivot][k])) pivot = i;
        if (!std::isfinite(matrix[pivot][k]) || std::abs(matrix[pivot][k]) < 1e-14)
            return false;
        std::swap(matrix[k], matrix[pivot]);
        std::swap(rhs[k], rhs[pivot]);
        const double diag = matrix[k][k];
        for (std::size_t i = k + 1; i < n; ++i) {
            const double factor = matrix[i][k] / diag;
            for (std::size_t j = k + 1; j < n; ++j)
                matrix[i][j] -= factor * matrix[k][j];
            rhs[i] -= factor * rhs[k];
        }
    }
    step.assign(n, 0.0);
    for (std::size_t k = n; k-- > 0;) {
        double sum = rhs[k];
        for (std::size_t j = k + 1; j < n; ++j) sum -= matrix[k][j] * step[j];
        step[k] = sum / matrix[k][k];
        if (!std::isfinite(step[k])) return false;
    }
    return true;
}
std::vector<SpeciesEntry> validate(const IdealEquilibriumProblem& problem) {
    if (problem.components.empty()) throw std::invalid_argument("At least one primary component is required");
    const std::size_t n = problem.components.size();
    std::map<std::string, std::size_t> indices;
    std::vector<SpeciesEntry> entries;
    for (std::size_t i = 0; i < n; ++i) {
        const auto& c = problem.components[i];
        if (c.name.empty() || !indices.emplace(c.name, i).second)
            throw std::invalid_argument("Duplicate or empty primary component name");
        if (!std::isfinite(c.total) || c.total <= 0.0)
            throw std::invalid_argument("Primary totals must be finite and positive: " + c.name);
        SpeciesEntry e{c.name, std::vector<int>(n, 0), 0.0, c.charge};
        e.nu[i] = 1;
        entries.push_back(std::move(e));
    }
    for (const auto& complex : problem.complexes) {
        if (complex.name.empty() || !indices.emplace(complex.name, n + entries.size()).second)
            throw std::invalid_argument("Duplicate or empty species name: " + complex.name);
        if (!std::isfinite(complex.log10K) || complex.stoichiometry.empty())
            throw std::invalid_argument("Invalid formation equilibrium: " + complex.name);
        SpeciesEntry e{complex.name, std::vector<int>(n, 0), complex.log10K, complex.charge};
        long long charge = 0;
        for (const auto& [name, coefficient] : complex.stoichiometry) {
            const auto it = indices.find(name);
            if (it == indices.end() || it->second >= n || coefficient <= 0)
                throw std::invalid_argument("Invalid primary/stoichiometry in complex " + complex.name);
            e.nu[it->second] = coefficient;
            charge += static_cast<long long>(coefficient) * problem.components[it->second].charge;
        }
        if (charge != complex.charge)
            throw std::invalid_argument("Complex charge inconsistent with primary components: " + complex.name);
        entries.push_back(std::move(e));
    }
    return entries;
}
bool evaluate(const IdealEquilibriumProblem& problem,
              const std::vector<SpeciesEntry>& entries,
              const std::vector<double>& logs, Evaluation& out) {
    const std::size_t n = logs.size();
    out.values.assign(entries.size(), 0.0);
    out.residuals.assign(n, -1.0);
    out.jacobian.assign(n, std::vector<double>(n, 0.0));
    for (std::size_t s = 0; s < entries.size(); ++s) {
        const auto& e = entries[s];
        double logValue = e.log10K;
        for (std::size_t j = 0; j < n; ++j) logValue += e.nu[j] * logs[j];
        if (!std::isfinite(logValue) || logValue > 250.0 || logValue < -300.0) return false;
        const double value = std::pow(10.0, logValue);
        if (!std::isfinite(value)) return false;
        out.values[s] = value;
        for (std::size_t i = 0; i < n; ++i) {
            if (e.nu[i] == 0) continue;
            const double ratio = value / problem.components[i].total;
            out.residuals[i] += e.nu[i] * ratio;
            for (std::size_t j = 0; j < n; ++j)
                out.jacobian[i][j] += LN10 * e.nu[i] * e.nu[j] * ratio;
        }
    }
    out.norm = 0.0;
    for (double residual : out.residuals) {
        if (!std::isfinite(residual)) return false;
        out.norm = std::max(out.norm, std::abs(residual));
    }
    return true;
}
void requireEnd(std::istringstream& in, int lineNo) {
    std::string extra;
    if (in >> extra && extra[0] != '#')
        throw std::runtime_error("Unexpected extra field at line " + std::to_string(lineNo));
}
}

IdealEquilibriumProblem loadIdealEquilibriumProblem(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open ideal model: " + filename);
    IdealEquilibriumProblem model;
    std::string line;
    int lineNo = 0;
    while (std::getline(file, line)) {
        ++lineNo;
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream in(line);
        std::string kind;
        in >> kind;
        if (kind == "component") {
            PrimaryComponent c;
            if (!(in >> c.name >> c.total >> c.charge))
                throw std::runtime_error("Invalid component at line " + std::to_string(lineNo));
            requireEnd(in, lineNo);
            model.components.push_back(std::move(c));
        } else if (kind == "complex") {
            SecondaryComplex s;
            int n = 0;
            if (!(in >> s.name >> s.charge >> s.log10K >> n) || n <= 0)
                throw std::runtime_error("Invalid complex at line " + std::to_string(lineNo));
            for (int j = 0; j < n; ++j) {
                std::string name;
                int coefficient = 0;
                if (!(in >> name >> coefficient) ||
                    !s.stoichiometry.emplace(name, coefficient).second)
                    throw std::runtime_error("Invalid complex component at line " + std::to_string(lineNo));
            }
            requireEnd(in, lineNo);
            model.complexes.push_back(std::move(s));
        } else {
            throw std::runtime_error("Unknown record type at line " + std::to_string(lineNo));
        }
    }
    validate(model);
    return model;
}

EquilibriumResult solveIdealEquilibrium(const IdealEquilibriumProblem& problem,
                                        const EquilibriumOptions& options) {
    if (options.maxIterations < 1 || !std::isfinite(options.relativeTolerance) ||
        options.relativeTolerance <= 0.0 || options.relativeTolerance >= 1.0)
        throw std::invalid_argument("Invalid equilibrium solver options");
    const auto entries = validate(problem);
    const std::size_t n = problem.components.size();
    std::vector<double> logs(n);
    for (std::size_t i = 0; i < n; ++i)
        logs[i] = std::log10(problem.components[i].total) - 0.5;
    Evaluation current;
    if (!evaluate(problem, entries, logs, current))
        throw std::runtime_error("Initial equilibrium evaluation is not finite");
    int iterations = 0;
    while (current.norm > options.relativeTolerance && iterations < options.maxIterations) {
        std::vector<double> rhs = current.residuals, step;
        for (double& r : rhs) r = -r;
        if (!solveLinear(current.jacobian, rhs, step))
            throw std::runtime_error("Equilibrium Jacobian is singular");
        bool accepted = false;
        for (int attempt = 0; attempt < 40; ++attempt) {
            const double alpha = std::ldexp(1.0, -attempt);
            std::vector<double> trial = logs;
            for (std::size_t j = 0; j < n; ++j) trial[j] += alpha * step[j];
            Evaluation candidate;
            if (evaluate(problem, entries, trial, candidate) && candidate.norm < current.norm) {
                logs = std::move(trial);
                current = std::move(candidate);
                accepted = true;
                break;
            }
        }
        if (!accepted) throw std::runtime_error("Equilibrium line search failed");
        ++iterations;
    }
    if (current.norm > options.relativeTolerance)
        throw std::runtime_error("Equilibrium solver did not converge in allowed iterations");
    EquilibriumResult result;
    result.iterations = iterations;
    for (std::size_t s = 0; s < entries.size(); ++s) {
        result.concentrations[entries[s].name] = current.values[s];
        result.chargeResidual += entries[s].charge * current.values[s];
    }
    for (std::size_t i = 0; i < n; ++i)
        result.balanceResiduals[problem.components[i].name] =
            current.residuals[i] * problem.components[i].total;
    return result;
}
