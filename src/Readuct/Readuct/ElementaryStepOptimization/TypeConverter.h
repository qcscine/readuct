/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_TYPECONVERTER_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_TYPECONVERTER_H

#include "MSVCCompatibility.h"
#include <Eigen/Core>

// TODO: Improve documentation

namespace Scine {
namespace Utils {
namespace BSplines {
class BSpline;
}
} // namespace Utils
namespace Readuct {

namespace ElementaryStepOptimization {

/*!
 * This class converts inner control points and eigen types.
 */
class SCINE_DLLEXPORT TypeConverter {
 public:
  static Eigen::MatrixXd getInnerControlPointMatrix(const Utils::BSplines::BSpline& spline);
  static void setInnerControlPoints(Utils::BSplines::BSpline& spline, const Eigen::MatrixXd& m);
};

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // ELEMENTARYSTEPOPTIMIZATION_TYPECONVERTER_H
