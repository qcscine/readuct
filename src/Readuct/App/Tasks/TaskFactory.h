/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_TASKFACTORY_H_
#define READUCT_TASKFACTORY_H_

#include "Tasks/AfirOptimizationTask.h"
#include "Tasks/BSplineInterpolationTask.h"
#include "Tasks/BondOrderTask.h"
#include "Tasks/GeometryOptimizationTask.h"
#include "Tasks/HessianTask.h"
#include "Tasks/IrcTask.h"
#include "Tasks/SinglePointTask.h"
#include "Tasks/Task.h"
#include "Tasks/TsOptimizationTask.h"
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace Scine {
namespace Readuct {

/**
 * @brief A factory generating Tasks by name (mainly).
 */
class TaskFactory {
 public:
  /// @brief Has only static functions.
  TaskFactory() = delete;
  /**
   * @brief Contstructs a Task with a given set of input and output systems.
   *
   * @param name   The name of the requested task.
   * @param input  The input system names for the task.
   * @param output The ouput system names for the task.
   * @return std::unique_ptr<Task> The requested task.
   */
  static std::unique_ptr<Task> produce(std::string name, const std::vector<std::string>& input,
                                       const std::vector<std::string>& output) {
    std::unique_ptr<Task> task;
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    if (name == "OPT" || name == "GEOOPT" || name == "GEOMETRYOPTIMIZATION" || name == "GEOMETRY_OPTIMIZATION") {
      task = std::make_unique<GeometryOptimizationTask>(input, output);
    }
    else if (name == "AFIR" || name == "AFIROPT" || name == "AFIROPTIMIZATION" || name == "AFIR_OPTIMIZATION") {
      task = std::make_unique<AfirOptimizationTask>(input, output);
    }
    else if (name == "TS" || name == "TSOPT" || name == "TRANSITIONSTATE_OPTIMIZATION" ||
             name == "TRANSITION_STATE_OPTIMIZATION") {
      task = std::make_unique<TsOptimizationTask>(input, output);
    }
    else if (name == "SP" || name == "SINGLEPOINT" || name == "ENERGY" || name == "SINGLE_POINT") {
      task = std::make_unique<SinglePointTask>(input, output);
    }
    else if (name == "IRC" || name == "IRCOPT") {
      task = std::make_unique<IrcTask>(input, output);
    }
    else if (name == "HESSIAN" || name == "FREQUENCY_ANALYSIS" || name == "FREQUENCYANALYSIS" || name == "FREQ" ||
             name == "FREQUENCY" || name == "FREQUENCIES") {
      task = std::make_unique<HessianTask>(input, output);
    }
    else if (name == "BONDS" || name == "BOND_ORDERS" || name == "BONDORDERS" || name == "BOS" || name == "BO") {
      task = std::make_unique<BondOrderTask>(input, output);
    }
    else if (name == "BSPLINE_INTERPOLATION" || name == "BSPLINEINTERPOLATION" || name == "BSPLINE") {
      task = std::make_unique<BSplineInterpolationTask>(input, output);
    }
    else {
      throw std::runtime_error("The requested task '" + name + "' is not available.\n");
    }
    return task;
  }
}; // namespace Readuct

} // namespace Readuct
} // namespace Scine

#endif // READUCT_TASKFACTORY_H_
