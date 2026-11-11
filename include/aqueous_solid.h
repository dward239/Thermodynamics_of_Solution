#pragma once
#include "aqueous_engine.h"
namespace aqueous {
// Restricted equilibrium for one neutral pure solid in the ideal model.
// No free/charge-balance basis components or competing solids; see README.
Result solveSingleIdealSolid(const Model& model, const Options& options);
}
