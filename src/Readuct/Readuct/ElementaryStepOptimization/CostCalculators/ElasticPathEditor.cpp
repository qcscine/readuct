/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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

DescriptorCollection ElasticPathEditor::getSettingDescriptors() const {
  DescriptorCollection descriptors("Elastic path cost calculator settings");
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
