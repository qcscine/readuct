/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TypeConverter.h"
#include <Utils/Math/BSplines/BSpline.h>

namespace Scine {
namespace Readuct {
namespace ElementaryStepOptimization {

Eigen::MatrixXd TypeConverter::getInnerControlPointMatrix(const Utils::BSplines::BSpline& spline) {
  const auto& controlPoints = spline.getControlPointMatrix();
  auto numberInnerControlPoints = spline.controlPointCount() - 2;
  return controlPoints.middleRows(1, numberInnerControlPoints);
}

void TypeConverter::setInnerControlPoints(Utils::BSplines::BSpline& spline, const Eigen::MatrixXd& m) {
  Eigen::MatrixXd controlPoints = spline.getControlPointMatrix();
  auto numberInnerControlPoints = spline.controlPointCount() - 2;
  controlPoints.middleRows(1, numberInnerControlPoints) = m;

  auto newSpline = Utils::BSplines::BSpline{spline.getKnotVector(), controlPoints, spline.getDegree()};
  spline = std::move(newSpline);
}

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine
