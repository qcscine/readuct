/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_COSTCOMBINER_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_COSTCOMBINER_H

#include "ReactionPathCostCalculator.h"
#include <iostream>

// TODO: Improve documentation

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

/*!
 * Cost calculator combining two other calculators
 * \tparam C1 First cost calculator
 * \tparam C2 Second cost calculator
 */
template<typename C1, typename C2>
class CostCombiner : public ReactionPathCostCalculator {
  static_assert(std::is_base_of<ReactionPathCostCalculator, C1>::value,
                "C1 must be a descendant of ReactionPathCostCalculator");
  static_assert(std::is_base_of<ReactionPathCostCalculator, C2>::value,
                "C2 must be a descendant of ReactionPathCostCalculator");

 public:
  /*! Give the contribution of the first calculator; must be between 0.0 and 1.0. */
  void setFirstCalculatorContribution(double f) {
    assert(0.0 <= f && f <= 1.0);
    contributionFromFirstCalculator_ = f;
  }

  double getFirstCalculatorContribution() const {
    return contributionFromFirstCalculator_;
  }

  const ReactionPathCostCalculator& firstCalculator() const {
    return calculator1_;
  }
  ReactionPathCostCalculator& firstCalculator() {
    return calculator1_;
  }

  const ReactionPathCostCalculator& secondCalculator() const {
    return calculator2_;
  }
  ReactionPathCostCalculator& secondCalculator() {
    return calculator2_;
  }

  static constexpr double defaultContributionFromFirstCalculator = 0.5;

 private:
  std::unique_ptr<ReactionPathCostCalculator> cloneImpl() const override {
    return std::make_unique<CostCombiner<C1, C2>>(*this);
  }

  bool energiesRequiredImpl() const override {
    return calculator1_.energiesRequired() || calculator2_.energiesRequired();
  }

  void calculateCostImpl(const Utils::BSplines::BSpline& spline, const EnergiesAndGradientsAlongSpline& energyValues) override {
    calculator1_.calculateCost(spline, energyValues);
    calculator2_.calculateCost(spline, energyValues);
  }

  double getCostImpl() const override {
    double f1 = contributionFromFirstCalculator_;
    double f2 = 1.0 - f1;

    return f1 * calculator1_.getCost() + f2 * calculator2_.getCost();
  }

  Eigen::MatrixXd getCostDerivativesImpl() const override {
    double f1 = contributionFromFirstCalculator_;
    double f2 = 1.0 - f1;

    return f1 * calculator1_.getCostDerivatives() + f2 * calculator2_.getCostDerivatives();
  }

  double contributionFromFirstCalculator_ = defaultContributionFromFirstCalculator;
  C1 calculator1_{};
  C2 calculator2_{};
};

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // ELEMENTARYSTEPOPTIMIZATION_COSTBASEDOPTIMIZATION_COSTCOMBINER_H