/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "BinaryArrayIO.h"

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

std::vector<double> BinaryArrayIO::readDoubles(std::istream& in) {
  int32_t size;
  in.read(reinterpret_cast<char*>(&size), sizeof(int32_t)); // NOLINT

  std::vector<double> array(size);
  in.read(reinterpret_cast<char*>(array.data()), size * sizeof(double)); // NOLINT
  return array;
}

void BinaryArrayIO::writeDoubles(std::ostream& out, const std::vector<double>& v) {
  auto size = static_cast<int32_t>(v.size());
  out.write(reinterpret_cast<char*>(&size), sizeof(int32_t));                // NOLINT
  out.write(reinterpret_cast<const char*>(v.data()), size * sizeof(double)); // NOLINT
}

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine
