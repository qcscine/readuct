/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "EnergyEditor.h"
#include <Utils/UniversalSettings/DescriptorCollection.h>
#include <Utils/UniversalSettings/ValueCollection.h>

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

Utils::UniversalSettings::DescriptorCollection EnergyEditor::getSettingDescriptors() const {
  Utils::UniversalSettings::DescriptorCollection descriptorCollection{"Energy along path"};
  return descriptorCollection;
}

void EnergyEditor::applyImpl(Energy& /*instance*/, const Utils::UniversalSettings::ValueCollection& /*values*/) const {
}

Utils::UniversalSettings::ValueCollection EnergyEditor::getAppliedSettingsImpl(const Energy& /*instance*/) const {
  return {};
}

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine
