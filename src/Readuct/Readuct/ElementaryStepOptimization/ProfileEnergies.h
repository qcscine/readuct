/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_PROFILEENERGIES_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_PROFILEENERGIES_H

#include <cassert>
#include <utility>
#include <vector>

// TODO: Improve documentation

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

/*!
 * Class for collection of (x, y) values, where x is a BSpline coordinate (between 0 and 1), and y is the
 * corresponding energy.
 */
class ProfileEnergies {
 public:
  using Array = std::vector<double>;

  ProfileEnergies() = default;
  ProfileEnergies(Array x, Array y);

  void setValues(Array x, Array y);

  const Array& getCoordinates() const;
  const Array& getEnergies() const;
  std::pair<double, double> getPair(int index) const;

  int size() const;
  bool empty() const;

 private:
  Array x_;
  Array y_;
};

inline ProfileEnergies::ProfileEnergies(ProfileEnergies::Array x, ProfileEnergies::Array y) {
  setValues(std::move(x), std::move(y));
}

inline void ProfileEnergies::setValues(ProfileEnergies::Array x, ProfileEnergies::Array y) {
  assert(x.size() == y.size());
  x_ = std::move(x);
  y_ = std::move(y);
}

inline int ProfileEnergies::size() const {
  return static_cast<int>(x_.size());
}

inline bool ProfileEnergies::empty() const {
  return size() == 0;
}

inline const ProfileEnergies::Array& ProfileEnergies::getCoordinates() const {
  return x_;
}

inline const ProfileEnergies::Array& ProfileEnergies::getEnergies() const {
  return y_;
}

inline std::pair<double, double> ProfileEnergies::getPair(int index) const {
  assert(0 <= index && index < size());
  return std::make_pair(x_[index], y_[index]);
}

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine

#endif // ELEMENTARYSTEPOPTIMIZATION_PROFILEENERGIES_H
