/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_COSTCALCULATOREDITOR_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_COSTCALCULATOREDITOR_H

#include "ReactionPathCostCalculator.h"
#include <Utils/UniversalSettings/GenericInstanceEditor.h>

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

/*!
 * Generic editor for cost calculators.
 */
using CostCalculatorEditor = Utils::UniversalSettings::GenericInstanceEditorWithDefaultConstructor<ReactionPathCostCalculator>;

template<typename T>
using CostCalculatorEditorImpl =
    Utils::UniversalSettings::GenericInstanceEditorWithDefaultConstructorImpl<CostCalculatorEditor, T>;

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_COSTCALCULATOREDITOR_H
