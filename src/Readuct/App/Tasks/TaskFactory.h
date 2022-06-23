/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
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
#include "Tasks/NtOptimization2Task.h"
#include "Tasks/NtOptimizationTask.h"
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
   * @param output The output system names for the task.
   * @param logger The logger to/through which all text output will be handled.
   * @return std::unique_ptr<Task> The requested task.
   */
  static std::unique_ptr<Task> produce(std::string name, const std::vector<std::string>& input,
                                       const std::vector<std::string>& output, std::shared_ptr<Core::Log> logger = nullptr) {
    std::unique_ptr<Task> task;
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    if (name == "OPT" || name == "GEOOPT" || name == "GEOMETRYOPTIMIZATION" || name == "GEOMETRY_OPTIMIZATION") {
      task = std::make_unique<GeometryOptimizationTask>(input, output, logger);
    }
    else if (name == "AFIR" || name == "AFIROPT" || name == "AFIROPTIMIZATION" || name == "AFIR_OPTIMIZATION") {
      task = std::make_unique<AfirOptimizationTask>(input, output, logger);
    }
    else if (name == "TS" || name == "TSOPT" || name == "TRANSITIONSTATE_OPTIMIZATION" ||
             name == "TRANSITION_STATE_OPTIMIZATION") {
      task = std::make_unique<TsOptimizationTask>(input, output, logger);
    }
    else if (name == "SP" || name == "SINGLEPOINT" || name == "ENERGY" || name == "SINGLE_POINT") {
      task = std::make_unique<SinglePointTask>(input, output, logger);
    }
    else if (name == "IRC" || name == "IRCOPT") {
      task = std::make_unique<IrcTask>(input, output, logger);
    }
    else if (name == "NT" || name == "NEWTONTRAJECTORY" || name == "NTOPTIMIZATION" ||
             name == "NEWTONTRAJECTORYOPTIMIZATION" || name == "NTOPT" || name == "NT1") {
      task = std::make_unique<NtOptimizationTask>(input, output);
    }
    else if (name == "NT2" || name == "NEWTONTRAJECTORY2" || name == "NTOPTIMIZATION2" ||
             name == "NEWTONTRAJECTORYOPTIMIZATION2" || name == "NTOPT2") {
      task = std::make_unique<NtOptimization2Task>(input, output);
    }
    else if (name == "HESSIAN" || name == "FREQUENCY_ANALYSIS" || name == "FREQUENCYANALYSIS" || name == "FREQ" ||
             name == "FREQUENCY" || name == "FREQUENCIES") {
      task = std::make_unique<HessianTask>(input, output, logger);
    }
    else if (name == "BONDS" || name == "BOND_ORDERS" || name == "BONDORDERS" || name == "BOS" || name == "BO") {
      task = std::make_unique<BondOrderTask>(input, output, logger);
    }
    else if (name == "BSPLINE_INTERPOLATION" || name == "BSPLINEINTERPOLATION" || name == "BSPLINE") {
      task = std::make_unique<BSplineInterpolationTask>(input, output, logger);
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
