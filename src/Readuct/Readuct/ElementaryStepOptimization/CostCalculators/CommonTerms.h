/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_COMMONTERMS_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_COMMONTERMS_H

#include <Eigen/Core>
#include <vector>

// TODO: Improve documentation

namespace Scine {

namespace Utils {
namespace AutomaticDifferentiation {
class FirstND;
}
} // namespace Utils

namespace Utils {
namespace BSplines {
class BSpline;
}
} // namespace Utils

namespace Readuct {

namespace ElementaryStepOptimization {
class RecurringProfileCalculator;

namespace CostBasedOptimization {

/*!
 * Frequent terms in cost functions.
 * The functions return objects for the cost with derivatives of the inner control points.
 */
class CommonTerms {
 public:
  static Utils::AutomaticDifferentiation::FirstND squaredNormOfSplineDerivative(const Utils::BSplines::BSpline& spline,
                                                                                double u);
  static Utils::AutomaticDifferentiation::FirstND singlePointEnergy(const Utils::BSplines::BSpline& spline, double u,
                                                                    double energy, const Eigen::VectorXd& gradients);
  static Utils::AutomaticDifferentiation::FirstND squaredNormOfGradient(const Utils::BSplines::BSpline& spline,
                                                                        double u, const Eigen::VectorXd& gradients,
                                                                        const Eigen::MatrixXd& hessian);
  static Utils::AutomaticDifferentiation::FirstND
  dotProductOfGradientAndSplineDerivative(const Utils::BSplines::BSpline& spline, double u,
                                          const Eigen::VectorXd& gradients, const Eigen::MatrixXd& hessian);
  static Utils::AutomaticDifferentiation::FirstND
  dotProductOfPositionAndSplineDerivative(const Utils::BSplines::BSpline& spline, double uPosition, double uTangent);
  static Utils::AutomaticDifferentiation::FirstND distanceSquared(const Utils::BSplines::BSpline& spline, double uA, double uB);
  static Utils::AutomaticDifferentiation::FirstND distanceSquaredToPosition(const Utils::BSplines::BSpline& spline,
                                                                            double u, const Eigen::VectorXd& pos);

  static Utils::AutomaticDifferentiation::FirstND dotProductOfPositions(const Utils::BSplines::BSpline& spline,
                                                                        double uA, double uB);

  static Utils::AutomaticDifferentiation::FirstND
  energyOfHighestPointFromQuadraticInterpolation(RecurringProfileCalculator& profileCalculator,
                                                 const std::vector<double>& energies, const Utils::BSplines::BSpline& spline);

  static Utils::AutomaticDifferentiation::FirstND energyAlongSpline(RecurringProfileCalculator& profileCalculator,
                                                                    const Utils::BSplines::BSpline& spline, double u);
  /*! cost: d(||C'(u)||^2)/du. */
  static Utils::AutomaticDifferentiation::FirstND derivativeOfSquaredTangentNorm(const Utils::BSplines::BSpline& spline,
                                                                                 double u);

 private:
  static double interpolateExtremumUValue(double u0, double u1, double u2, double e0, double e1, double e2);
  static int getIndexForHighestEnergy(const std::vector<double>& energies);
  static double getUValueWithMaxEnergy(const std::vector<double>& coordinates, const std::vector<double>& energies, int maxIndex);
};

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATION_COMMONTERMS_H