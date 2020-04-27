/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_ELASTICPATHEDITOR_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_ELASTICPATHEDITOR_H

#include "CostCalculatorEditor.h"
#include "ElasticPath.h"

// TODO: Improve documentation

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

/*!
 * Cost calculator editor for ElasticPath.
 */
class ElasticPathEditor : public CostCalculatorEditorImpl<ElasticPath> {
 public:
  Utils::UniversalSettings::DescriptorCollection getSettingDescriptors() const override;

 private:
  void applyImpl(ElasticPath& instance, const Utils::UniversalSettings::ValueCollection& values) const override;
  Utils::UniversalSettings::ValueCollection getAppliedSettingsImpl(const ElasticPath& instance) const override;
};

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // ELEMENTARYSTEPOPTIMIZATION_COSTBASEDOPTIMIZATION_ELASTICPATHEDITOR_H
