#pragma once
#include "aqueous_engine.h"
#include <string>
#include <vector>

namespace aqueous {
struct SweepPoint {
    double parameterValue = 0;
    Result result;
};

// Parameter = "temperature" or "total"; componentName required only for "total".
// On any invalid model/failed step throws, not partial success.
std::vector<SweepPoint> runSweep(
    const Model& model, const std::string& parameter,
    const std::string& componentName, double start, double stop,
    int steps, const Options& options = {});
std::string sweepCsv(const std::vector<SweepPoint>& points,
                     const std::string& parameter);
}
