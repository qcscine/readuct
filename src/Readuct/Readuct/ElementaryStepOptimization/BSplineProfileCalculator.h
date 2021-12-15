/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_BSPLINEPROFILECALCULATOR_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_BSPLINEPROFILECALCULATOR_H

#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <vector>

// TODO: Improve documentation

namespace Scine {

namespace Core {
class Calculator;
}

namespace Utils {
namespace BSplines {
class BSpline;
}
} // namespace Utils

namespace Readuct {

namespace ElementaryStepOptimization {

/*!
 * This class calculates the energies and forces for points along a bspline.
 */
class BSplineProfileCalculator {
 public:
  explicit BSplineProfileCalculator(Core::Calculator& calculator);

  void calculateAllEnergies(const Utils::BSplines::BSpline& spline, const std::vector<double>& uValues,
                            std::vector<double>& energies);
  void calculateAllEnergiesAndGradients(const Utils::BSplines::BSpline& spline, const std::vector<double>& uValues,
                                        std::vector<double>& energies, std::vector<Utils::GradientCollection>& gradients);
  void calculateAllUpToSecondDerivative(const Utils::BSplines::BSpline& spline, const std::vector<double>& uValues,
                                        std::vector<double>& energies, std::vector<Utils::GradientCollection>& gradients,
                                        std::vector<Eigen::MatrixXd>& hessians);

  void calculateEnergy(const Utils::BSplines::BSpline& spline, double u, double& energy);
  void calculateEnergyAndGradients(const Utils::BSplines::BSpline& spline, double u, double& energy,
                                   Utils::GradientCollection& gradients);
  void calculateUpToSecondDerivatives(const Utils::BSplines::BSpline& spline, double u, double& energy,
                                      Utils::GradientCollection& gradients, Eigen::MatrixXd& hessians);

  Utils::PositionCollection positionsForUValue(const Utils::BSplines::BSpline& spline, double u) const;

 private:
  Core::Calculator& calculator_;
};

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // ELEMENTARYSTEPOPTIMIZATION_BSPLINEPROFILECALCULATOR_H
