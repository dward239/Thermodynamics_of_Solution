#pragma once
#include "species.h"
#include <map>
#include <string>

// Screening calculations only. No multiphase equilibrium solver is implemented.
double calculateSaturationIndex(const SolidPhase& solid,
    const std::map<std::string, double>& activities);
std::map<std::string, double> estimateIdealActivities(const Solution& solution);
