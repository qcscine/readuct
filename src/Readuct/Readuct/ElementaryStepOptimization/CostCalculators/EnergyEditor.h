/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_ENERGYEDITOR_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_ENERGYEDITOR_H

#include "CostCalculatorEditor.h"
#include "Energy.h"

// TODO: Improve documentation

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

/*!
 * Editor for the Energy cost calculator.
 */
class EnergyEditor : public CostCalculatorEditorImpl<Energy> {
 public:
  Utils::UniversalSettings::DescriptorCollection getSettingDescriptors() const override;

 private:
  void applyImpl(Energy& instance, const Utils::UniversalSettings::ValueCollection& values) const override;
  Utils::UniversalSettings::ValueCollection getAppliedSettingsImpl(const Energy& instance) const override;
};

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // ELEMENTARYSTEPOPTIMIZATION_COSTBASEDOPTIMIZATION_ENERGYEDITOR_H
