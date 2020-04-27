/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_REACTIONPROFILE_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_REACTIONPROFILE_H

#include "ProfileEnergies.h"
#include <Utils/Math/BSplines/MolecularSpline.h>

// TODO: Improve documentation

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

/*!
 * This class describes a reaction path with the corresponding energies, as a MolecularSpline
 * and, optionally, the corresponding energies.
 */
class ReactionProfile {
 public:
  explicit ReactionProfile(Utils::BSplines::MolecularSpline spline = {}, ProfileEnergies energies = {});

  void setProfileEnergies(ProfileEnergies energies);

  bool hasEnergies() const;

  Utils::BSplines::MolecularSpline& getMolecularSpline();
  const Utils::BSplines::MolecularSpline& getMolecularSpline() const;

  ProfileEnergies& getProfileEnergies();
  const ProfileEnergies& getProfileEnergies() const;

 private:
  Utils::BSplines::MolecularSpline spline_;
  ProfileEnergies energies_;
};

inline ReactionProfile::ReactionProfile(Utils::BSplines::MolecularSpline spline, ProfileEnergies energies)
  : spline_(std::move(spline)), energies_(std::move(energies)) {
}

inline void ReactionProfile::setProfileEnergies(ProfileEnergies energies) {
  energies_ = std::move(energies);
}

inline bool ReactionProfile::hasEnergies() const {
  return !energies_.empty();
}

inline Utils::BSplines::MolecularSpline& ReactionProfile::getMolecularSpline() {
  return spline_;
}

inline const Utils::BSplines::MolecularSpline& ReactionProfile::getMolecularSpline() const {
  return spline_;
}

inline ProfileEnergies& ReactionProfile::getProfileEnergies() {
  return energies_;
}

inline const ProfileEnergies& ReactionProfile::getProfileEnergies() const {
  return energies_;
}

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine

#endif // ELEMENTARYSTEPOPTIMIZATION_REACTIONPROFILE_H