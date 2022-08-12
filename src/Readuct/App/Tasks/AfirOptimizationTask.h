/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_AFIROPTIMIZATIONTASK_H_
#define READUCT_AFIROPTIMIZATIONTASK_H_

/* Readuct */
#include "Tasks/Task.h"
/* Scine */
#include <Utils/GeometryOptimization/AfirOptimizer.h>
#include <Utils/GeometryOptimization/GeometryOptimization.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Optimizer/GradientBased/Bfgs.h>
#include <Utils/Optimizer/GradientBased/Lbfgs.h>
#include <Utils/Optimizer/GradientBased/SteepestDescent.h>
#include <Utils/Optimizer/HessianBased/NewtonRaphson.h>
/* External */
#include <boost/exception/diagnostic_information.hpp>
#include <boost/filesystem.hpp>
/* std c++ */
#include <cstdio>
#include <fstream>
#include <iostream>

namespace Scine {
namespace Readuct {

class AfirOptimizationTask : public Task {
 public:
  /**
   * @brief Construct a new AfirOptimizationTask.
   * @param input  The input system names for the task.
   * @param output The output system names for the task.
   * @param logger The logger to/through which all text output will be handled.
   */
  AfirOptimizationTask(std::vector<std::string> input, std::vector<std::string> output, std::shared_ptr<Core::Log> logger = nullptr)
    : Task(std::move(input), std::move(output), std::move(logger)) {
  }

  std::string name() const override {
    return "AFIR Optimization";
  }

  bool run(SystemsMap& systems, Utils::UniversalSettings::ValueCollection taskSettings, bool testRunOnly = false,
           std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>
               observers = {}) const final {
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
        throw std::runtime_error("Cannot calculate AFIR for monoatomic systems.");
      }
    }

    // Generate optimizer
    auto optimizertype = taskSettings.extract("optimizer", std::string{"BFGS"});
    std::transform(optimizertype.begin(), optimizertype.end(), optimizertype.begin(), ::toupper);
    std::shared_ptr<Utils::AfirOptimizerBase> optimizer;
    if (optimizertype == "LBFGS") {
      auto tmp = std::make_shared<Utils::AfirOptimizer<Utils::Lbfgs>>(*calc);
      tmp->optimizer.useTrustRadius = true;
      tmp->optimizer.trustRadius = 0.1;
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "BFGS") {
      auto tmp = std::make_shared<Utils::AfirOptimizer<Utils::Bfgs>>(*calc);
      tmp->optimizer.useTrustRadius = true;
      tmp->optimizer.trustRadius = 0.1;
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "SD" || optimizertype == "STEEPESTDESCENT") {
      auto tmp = std::make_shared<Utils::AfirOptimizer<Utils::SteepestDescent>>(*calc);
      optimizer = std::move(tmp);
    }
    else {
      throw std::runtime_error("Unknown Optimizer requested for AFIR, available are: SD, BFGS and LBFGS!");
    }

    // Read and delete special settings
    bool stopOnError = stopOnErrorExtraction(taskSettings);

    // Apply user settings
    auto settings = optimizer->getSettings();
    settings.merge(taskSettings);
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
    using Writer = Utils::XyzStreamHandler;
    const std::string& outputSystem = (!_output.empty() ? _output.front() : _input.front());
    boost::filesystem::path dir(outputSystem);
    boost::filesystem::create_directory(dir);
    boost::filesystem::path trjfile(outputSystem + ".afir.trj.xyz");
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

    // Add custom observers
    auto customObservers = [&calc, &observers](const int& cycle, const double& /*energy*/, const Eigen::VectorXd& /*params*/) {
      for (auto& observer : observers) {
        auto atoms = calc->getStructure();
        Utils::Results& results = calc->results();
        observer(cycle, *atoms, results, "afir_scan");
      }
    };
    optimizer->addObserver(customObservers);

    // Run optimization
    auto structure = calc->getStructure();
    int cycles = 0;
    try {
      cycles = optimizer->optimize(*structure, *_logger);
    }
    catch (...) {
      Writer::write(trajectory, *calc->getStructure());
      trajectory.close();
      _logger->error << "AFIR Optimization failed with error!" << Core::Log::endl;
      if (stopOnError) {
        throw;
      }
      _logger->error << boost::current_exception_diagnostic_information() << Core::Log::endl;
      return false;
    }
    trajectory.close();

    int maxiter = settings.getInt("convergence_max_iterations");
    if (maxiter > cycles) {
      cout << Core::Log::endl
           << "    Converged after " << cycles << " iterations." << Core::Log::endl
           << Core::Log::endl;
    }
    else {
      cout << Core::Log::endl
           << "    Stopped after " << maxiter << " iterations." << Core::Log::endl
           << Core::Log::endl;
      if (stopOnError) {
        throw std::runtime_error("Problem: AFIR optimization did not converge.");
      }
    }

    // Print/Store results
    boost::filesystem::path xyzfile(outputSystem + ".xyz");
    std::ofstream xyz((dir / xyzfile).string(), std::ofstream::out);
    Writer::write(xyz, *(calc->getStructure()));
    xyz.close();
    systems[outputSystem] = std::move(calc);

    return cycles < maxiter;
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_AFIROPTIMIZATIONTASK_H_
