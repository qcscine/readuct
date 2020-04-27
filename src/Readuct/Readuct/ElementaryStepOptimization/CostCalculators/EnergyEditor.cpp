/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "EnergyEditor.h"
#include <Utils/UniversalSettings/DescriptorCollection.h>
#include <Utils/UniversalSettings/ValueCollection.h>

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

namespace {
constexpr const char* description = "Energy along path";
}

Utils::UniversalSettings::DescriptorCollection EnergyEditor::getSettingDescriptors() const {
  Utils::UniversalSettings::DescriptorCollection descriptorCollection{description};
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
