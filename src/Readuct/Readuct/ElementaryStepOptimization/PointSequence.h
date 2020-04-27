/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_POINTSEQUENCE_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_POINTSEQUENCE_H

#include <cassert>
#include <vector>

// TODO: Improve documentation

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

/*!
 * Class representing a collection of points that can be used, for instance, for an integration grid.
 * TODO: allow for non-uniform points.
 */
class PointSequence {
 public:
  PointSequence() = default;
  PointSequence(double min, double max, int numberPoints);

  int count() const;

  /*! whether all integration points are equidistant. */
  bool isUniform() const;

  /*! Distance between two integration points, provided the points are uniformly distributed. */
  double interval() const;

  double operator[](int index) const;
  double at(int index) const;

  double min() const;
  double max() const;

  const std::vector<double>& underlyingArray() const;

 private:
  using container = std::vector<double>;

  container points_;
  double interval_{0};
};

inline PointSequence::PointSequence(double min, double max, int numberPoints) {
  assert(max >= min);
  assert(numberPoints > 1);

  interval_ = (max - min) / (numberPoints - 1);

  points_.resize(numberPoints);

  for (int i = 0; i < numberPoints; ++i) {
    points_[i] = min + i * interval_;
  }
}

inline int PointSequence::count() const {
  return static_cast<int>(points_.size());
}

inline bool PointSequence::isUniform() const {
  return true;
}

inline double PointSequence::interval() const {
  assert(isUniform());
  return interval_;
}

inline double PointSequence::operator[](int index) const {
  return points_[static_cast<container::size_type>(index)];
}

inline double PointSequence::at(int index) const {
  return points_.at(static_cast<container::size_type>(index));
}

inline double PointSequence::min() const {
  return points_.front();
}

inline double PointSequence::max() const {
  return points_.back();
}

inline const std::vector<double>& PointSequence::underlyingArray() const {
  return points_;
}

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine

#endif // ELEMENTARYSTEPOPTIMIZATION_POINTSEQUENCE_H