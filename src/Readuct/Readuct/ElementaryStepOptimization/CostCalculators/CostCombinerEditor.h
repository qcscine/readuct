/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_COSTCOMBINEREDITOR_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_COSTCOMBINEREDITOR_H

#include "CostCalculatorEditor.h"
#include "CostCombiner.h"
#include <Utils/UniversalSettings/DescriptorCollection.h>
#include <Utils/UniversalSettings/ValueCollection.h>

// TODO: Improve documentation

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

/* typedef for the combiner class from the editor types instead of from the classes themselves. */
template<typename E1, typename E2>
using Combiner = CostCombiner<typename E1::InstanceClass, typename E2::InstanceClass>;

/*!
 * Cost calculator editor for CostCombiner
 * \tparam E1 First cost calculator editor
 * \tparam E2 Second cost calculator editor
 */
template<typename E1, typename E2>
class CostCombinerEditor : public CostCalculatorEditorImpl<Combiner<E1, E2>> {
  static constexpr const char* description = "Combined cost calculator";
  static constexpr const char* factorKey = "first_weight";
  static constexpr const char* factorDescription = "Weight for first calculator (between 0 and 1)";

 public:
  CostCombinerEditor(std::string key1, std::string key2)
    : calculator1Key_(std::move(key1)), calculator2Key_(std::move(key2)) {
  }

  Utils::UniversalSettings::DescriptorCollection getSettingDescriptors() const override {
    Utils::UniversalSettings::DescriptorCollection descriptors(description);

    Utils::UniversalSettings::DoubleDescriptor weight(factorDescription);
    weight.setMinimum(0);
    weight.setMaximum(1);
    weight.setDefaultValue(Combiner<E1, E2>::defaultContributionFromFirstCalculator);

    auto descriptors1 = editor1_.getSettingDescriptors();
    auto descriptors2 = editor2_.getSettingDescriptors();

    descriptors.push_back(calculator1Key_, std::move(descriptors1));
    descriptors.push_back(calculator2Key_, std::move(descriptors2));
    descriptors.push_back(factorKey, std::move(weight));
    return descriptors;
  }

 private:
  void applyImpl(Combiner<E1, E2>& instance, const Utils::UniversalSettings::ValueCollection& values) const override {
    editor1_.apply(instance.firstCalculator(), values.getCollection(calculator1Key_));
    editor2_.apply(instance.secondCalculator(), values.getCollection(calculator2Key_));
    instance.setFirstCalculatorContribution(values.getDouble(factorKey));
  }

  Utils::UniversalSettings::ValueCollection getAppliedSettingsImpl(const Combiner<E1, E2>& instance) const override {
    Utils::UniversalSettings::ValueCollection values;

    auto firstCalculatorValues = editor1_.getAppliedSettings(instance.firstCalculator());
    auto secondCalculatorValues = editor2_.getAppliedSettings(instance.secondCalculator());

    values.addCollection(calculator1Key_, std::move(firstCalculatorValues));
    values.addCollection(calculator2Key_, std::move(secondCalculatorValues));
    values.addDouble(factorKey, instance.getFirstCalculatorContribution());
    return values;
  }

  std::string calculator1Key_;
  std::string calculator2Key_;
  E1 editor1_;
  E2 editor2_;
};

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_COSTCOMBINEREDITOR_H
