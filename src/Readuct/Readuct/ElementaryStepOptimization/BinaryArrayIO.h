/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_BINARYARRAYIO_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_BINARYARRAYIO_H

#include <istream>
#include <ostream>
#include <vector>

// TODO: Improve documentation

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

/*!
 * Class for the input and output of arrays in binary format.
 */
class BinaryArrayIO {
 public:
  static std::vector<double> readDoubles(std::istream& in);

  static void writeDoubles(std::ostream& out, const std::vector<double>& v);

 private:
};

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine

#endif // CPPUTIL_BINARYARRAYIO_H