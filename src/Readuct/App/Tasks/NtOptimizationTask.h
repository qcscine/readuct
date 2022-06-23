/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_NTOPTIMIZATIONTASK_H_
#define READUCT_NTOPTIMIZATIONTASK_H_

/* Readuct */
#include "Tasks/Task.h"
/* Scine */
#include <Utils/GeometryOptimization/NtOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/Yaml.h>
/* External */
#include <boost/exception/diagnostic_information.hpp>
#include <boost/filesystem.hpp>
/* std c++ */
#include <cstdio>
#include <fstream>
#include <iostream>

namespace Scine {
namespace Readuct {

class NtOptimizationTask : public Task {
 public:
  /**
   * @brief Construct a new NtOptimizationTask.
   * @param input  The input system names for the task.
   * @param output The output system names for the task.
   * @param logger The logger to/through which all text output will be handled.
   */
  NtOptimizationTask(std::vector<std::string> input, std::vector<std::string> output, std::shared_ptr<Core::Log> logger = nullptr)
    : Task(std::move(input), std::move(output), std::move(logger)) {
  }

  std::string name() const override {
    return "NT1 Optimization";
  }

  bool run(SystemsMap& systems, Utils::UniversalSettings::ValueCollection taskSettings, bool testRunOnly = false) const final {
    warningIfMultipleInputsGiven();
    warningIfMultipleOutputsGiven();

    bool silentCalculator = taskSettings.extract("silent_stdout_calculator", true);
    std::shared_ptr<Core::Calculator> calc;
    if (!testRunOnly) { // leave out in case of task chaining --> attention calc is NULL
      // Note: _input is guaranteed not to be empty by Task constructor
      calc = copyCalculator(systems, _input.front(), name());
      Utils::CalculationRoutines::setLog(*calc, true, true, !silentCalculator);

      // Check system size
      if (calc->getStructure()->size() == 1) {
        throw std::runtime_error("Cannot calculate NT1 optimization for monoatomic systems.");
      }
    }

    // Generate optimizer
    auto optimizer = std::make_unique<Utils::NtOptimizer>(*calc);
    // Read and delete special settings
    bool stopOnError = stopOnErrorExtraction(taskSettings);
    // Apply user settings
    auto settings = optimizer->getSettings();
    settings.merge(taskSettings);
    if (!settings.valid()) {
      settings.throwIncorrectSettings();
    }
    optimizer->setSettings(settings);

    // If no errors encountered until here, the basic settings should be alright
    if (testRunOnly) {
      return true;
    }

    // Add observer
    // Trajectory stream
    using Writer = Utils::XyzStreamHandler;
    const std::string& outputSystem = ((!_output.empty()) ? _output[0] : _input[0]);
    boost::filesystem::path dir(outputSystem);
    boost::filesystem::create_directory(dir);
    boost::filesystem::path trjfile(outputSystem + ".nt.trj.xyz");
    std::ofstream trajectory((dir / trjfile).string(), std::ofstream::out);
    double oldEnergy = 0.0;
    Eigen::VectorXd oldParams;
    auto func = [&](const int& cycle, const double& energy, const Eigen::VectorXd& params) {
      if (oldParams.size() != params.size()) {
        oldParams = params;
      }
      if (cycle == 1) {
        printf("%7s %16s %16s %16s %16s\n", "Cycle", "Energy", "Energy Diff.", "Step RMS", "Max. Step");
      }
      auto diff = (params - oldParams).eval();
      printf("%7d %+16.9f %+16.9f %+16.9f %+16.9f\n", cycle, energy, energy - oldEnergy,
             sqrt(diff.squaredNorm() / diff.size()), diff.cwiseAbs().maxCoeff());
      oldEnergy = energy;
      oldParams = params;
      auto structure = calc->getStructure();
      Writer::write(trajectory, *structure, std::to_string(energy));
    };
    optimizer->addObserver(func);
    auto cout = _logger->output;

    // Run optimization
    auto structure = calc->getStructure();
    int cycles = 0;
    try {
      cycles = optimizer->optimize(*structure, *_logger);
    }
    catch (...) {
      trajectory.close();
      if (stopOnError) {
        throw;
      }
      // check if thrown exception corresponds to no actual calculation failures but simply no guess found
      // --> no thrown exception intended in this case if stopOnError equals false
      std::string noGuess = "No transition state guess was found in Newton Trajectory scan";
      size_t found = boost::current_exception_diagnostic_information().find(noGuess);
      if (found == std::string::npos) { // did not find harmless exception -> throw
        throw;
      }
      cout << Core::Log::endl << "    No TS guess found in NT1 scan." << Core::Log::endl << Core::Log::endl;
      return false;
    }

    trajectory.close();

    cout << Core::Log::endl
         << "    Found TS guess after " << cycles << " iterations." << Core::Log::endl
         << Core::Log::endl;

    // Print/Store results
    systems[outputSystem] = calc;
    boost::filesystem::path xyzfile(outputSystem + ".xyz");
    std::ofstream xyz((dir / xyzfile).string(), std::ofstream::out);
    Writer::write(xyz, *(calc->getStructure()));
    xyz.close();

    return true;
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_NTOPTIMIZATIONTASK_H_
