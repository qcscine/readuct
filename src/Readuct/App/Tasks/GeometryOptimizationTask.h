/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_GEOMETRYOPTIMIZATIONTASK_H_
#define READUCT_GEOMETRYOPTIMIZATIONTASK_H_

/* Readuct */
#include "Tasks/Task.h"
/* Scine */
#include <Utils/GeometryOptimization/CoordinateSystem.h>
#include <Utils/GeometryOptimization/GeometryOptimization.h>
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Optimizer/GradientBased/Bfgs.h>
#include <Utils/Optimizer/GradientBased/Lbfgs.h>
#include <Utils/Optimizer/GradientBased/SteepestDescent.h>
#include <Utils/Optimizer/HessianBased/NewtonRaphson.h>
#include <Utils/UniversalSettings/SettingsNames.h>
/* External includes */
#include <boost/exception/diagnostic_information.hpp>
#include <boost/filesystem.hpp>
#include <cstdio>
#include <fstream>

namespace Scine {
namespace Readuct {

class GeometryOptimizationTask : public Task {
 public:
  /**
   * @brief Construct a new GeometryOptimizationTask.
   * @param input  The input system names for the task.
   * @param output The output system names for the task.
   * @param logger The logger to/through which all text output will be handled.
   */
  GeometryOptimizationTask(std::vector<std::string> input, std::vector<std::string> output,
                           std::shared_ptr<Core::Log> logger = nullptr)
    : Task(std::move(input), std::move(output), std::move(logger)) {
  }

  std::string name() const override {
    return "Geometry Optimization";
  }

  bool run(SystemsMap& systems, Utils::UniversalSettings::ValueCollection taskSettings, bool testRunOnly = false) const final {
    warningIfMultipleInputsGiven();
    warningIfMultipleOutputsGiven();
    // Get/Copy Calculator
    std::shared_ptr<Core::Calculator> calc;
    if (!testRunOnly) { // leave out in case of task chaining --> attention calc is NULL
      // Note: _input is guaranteed not to be empty by Task constructor
      calc = copyCalculator(systems, _input.front(), name());
    }

    // Generate optimizer
    auto optimizertype = taskSettings.extract("optimizer", std::string{"BFGS"});
    std::transform(optimizertype.begin(), optimizertype.end(), optimizertype.begin(), ::toupper);
    std::shared_ptr<Utils::GeometryOptimizerBase> optimizer;
    if (optimizertype == "LBFGS") {
      auto tmp = std::make_shared<Utils::GeometryOptimizer<Utils::Lbfgs>>(*calc);
      // Default convergence options
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "BFGS") {
      auto tmp = std::make_shared<Utils::GeometryOptimizer<Utils::Bfgs>>(*calc);
      // Default convergence options
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "SD" || optimizertype == "STEEPESTDESCENT") {
      auto tmp = std::make_shared<Utils::GeometryOptimizer<Utils::SteepestDescent>>(*calc);
      // Default convergence options
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "NR" || optimizertype == "NEWTONRAPHSON") {
      auto tmp = std::make_shared<Utils::GeometryOptimizer<Utils::NewtonRaphson>>(*calc);
      // Default convergence options
      optimizer = std::move(tmp);
    }
    else {
      throw std::runtime_error(
          "Unknown Optimizer requested for a geometry optimization, available are: SD, NR, BFGS and LBFGS!");
    }

    // Read and delete special settings
    bool stopOnError = stopOnErrorExtraction(taskSettings);
    // Apply settings
    auto settings = optimizer->getSettings();
    settings.merge(taskSettings);
    if (!testRunOnly && systems.at(_input[0])->getStructure()->size() == 1 &&
        Utils::CoordinateSystemInterpreter::getCoordinateSystemFromString(
            settings.getString("geoopt_coordinate_system")) != Utils::CoordinateSystem::Cartesian) {
      /* For monoatomic systems avoid transformation to internals resulting into
       * empty optimization parameters. Necessary because observer and
       * convergence check use transformed parameters
       */
      _logger->warning << "  Cannot treat a monoatomic system in internal coordinates. Switching to Cartesians."
                       << Core::Log::endl
                       << Core::Log::endl;
      settings.modifyString("geoopt_coordinate_system", Utils::CoordinateSystemInterpreter::getStringFromCoordinateSystem(
                                                            Utils::CoordinateSystem::Cartesian));
    }
    if (!settings.valid()) {
      settings.throwIncorrectSettings();
    }
    optimizer->setSettings(settings);
    if (!testRunOnly && !Utils::GeometryOptimization::settingsMakeSense(*optimizer)) {
      throw std::logic_error("The given calculator settings are too inaccurate for the given convergence criteria of "
                             "this optimization Task");
    }

    // If no errors encountered until here, the basic settings should be alright
    if (testRunOnly) {
      return true;
    }

    // Add observer
    // Trajectory stream
    const std::string& outputSystem = (!_output.empty() ? _output.front() : _input.front());
    using Writer = Utils::XyzStreamHandler;
    boost::filesystem::path dir(outputSystem);
    boost::filesystem::create_directory(dir);
    boost::filesystem::path trjfile(outputSystem + ".opt.trj.xyz");
    std::ofstream trajectory((dir / trjfile).string(), std::ofstream::out);
    double oldEnergy = 0.0;
    Eigen::VectorXd oldParams;
    auto cout = _logger->output;
    auto func = [&](const int& cycle, const double& energy, const Eigen::VectorXd& params) {
      if (oldParams.size() != params.size()) {
        oldParams = params;
      }
      if (cycle == 1) {
        cout.printf("%7s %16s %16s %16s %16s\n", "Cycle", "Energy", "Energy Diff.", "Step RMS", "Max. Step");
      }
      auto diff = (params - oldParams).eval();
      cout.printf("%7d %+16.9f %+16.9f %+16.9f %+16.9f\n", cycle, energy, energy - oldEnergy,
                  sqrt(diff.squaredNorm() / diff.size()), diff.cwiseAbs().maxCoeff());
      oldEnergy = energy;
      oldParams = params;
      auto structure = calc->getStructure();
      Writer::write(trajectory, *structure, std::to_string(energy));
    };
    optimizer->addObserver(func);

    // Run optimization
    auto structure = calc->getStructure();
    int cycles = 0;
    try {
      cycles = optimizer->optimize(*structure, *_logger);
    }
    catch (...) {
      Writer::write(trajectory, *calc->getStructure());
      trajectory.close();
      _logger->error << "Optimization failed with error!" << Core::Log::endl;
      if (stopOnError) {
        throw;
      }
      _logger->error << boost::current_exception_diagnostic_information() << Core::Log::endl;
      return false;
    }
    trajectory.close();

    int maxiter = settings.getInt("convergence_max_iterations");
    if (cycles < maxiter) {
      cout << Core::Log::endl
           << "    Converged after " << cycles << " iterations." << Core::Log::endl
           << Core::Log::endl;
    }
    else {
      cout << Core::Log::endl
           << "    Stopped after " << maxiter << " iterations." << Core::Log::endl
           << Core::Log::endl;
      if (stopOnError) {
        throw std::runtime_error("Problem: Structure optimization did not converge.");
      }
    }

    // Print/Store result
    systems[outputSystem] = calc;
    boost::filesystem::path xyzfile(outputSystem + ".xyz");
    std::ofstream xyz((dir / xyzfile).string(), std::ofstream::out);
    Writer::write(xyz, *(calc->getStructure()));
    xyz.close();

    return cycles < maxiter;
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_GEOMETRYOPTIMIZATIONTASK_H_
