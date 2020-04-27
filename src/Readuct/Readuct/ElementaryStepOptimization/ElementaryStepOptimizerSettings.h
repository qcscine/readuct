/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_ELEMENTARYSTEPOPTIMIZERSETTINGS_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_ELEMENTARYSTEPOPTIMIZERSETTINGS_H

#include "Utils/Settings.h"

namespace Scine {
namespace Readuct {
namespace ElementaryStepOptimization {

/**
 * @brief Settings for an ElementaryStepOptimizer.
 *
 * Uses template arguments in order to automatically include the
 * settings of underlying objects into the given settings.
 *
 * @tparam OptimizerType The underlying Optimizer class.
 * @tparam ConvergenceCheckType The underlying ConvergenceCheck class.
 */
template<class OptimizerType, class ConvergenceCheckType>
class ElementaryStepOptimizerSettings : public Utils::Settings {
 public:
  /**
   * @brief Construct a new ElementaryStepOptimizerSettings object.
   *
   * Sets the default values of the settings to the current values set in the objects
   * given to the constructor.
   *
   * @param optimizer The optimizer.
   * @param check The convergence check criteria.
   */
  ElementaryStepOptimizerSettings(const OptimizerType& optimizer, const ConvergenceCheckType& check)
    : Settings("ElementaryStepOptimizerSettings") {
    optimizer.addSettingsDescriptors(this->_fields);
    check.addSettingsDescriptors(this->_fields);

    Utils::UniversalSettings::IntDescriptor num_integration_points(
        "The number of integration points used to optimize the spline.");
    num_integration_points.setDefaultValue(21);
    this->_fields.push_back("num_integration_points", num_integration_points);

    this->resetToDefaults();
  }
};

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine

#endif // READUCT_ELEMENTARYSTEPOPTIMIZATION_ELEMENTARYSTEPOPTIMIZERSETTINGS_H
