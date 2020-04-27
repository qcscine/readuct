/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ElasticPathEditor.h"
#include <Utils/UniversalSettings/DescriptorCollection.h>
#include <Utils/UniversalSettings/ValueCollection.h>

namespace Scine {
namespace Readuct {

using namespace Utils::UniversalSettings;

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

namespace {
constexpr const char* description = "Elastic path cost calculator settings";
}

DescriptorCollection ElasticPathEditor::getSettingDescriptors() const {
  DescriptorCollection descriptors(description);
  return descriptors;
}

void ElasticPathEditor::applyImpl(ElasticPath& /*instance*/, const ValueCollection& /*values*/) const {
}

ValueCollection ElasticPathEditor::getAppliedSettingsImpl(const ElasticPath& /*instance*/) const {
  ValueCollection values;
  return values;
}

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine
