/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_IRCTASK_H_
#define READUCT_IRCTASK_H_

/* Readuct */
#include "Tasks/Task.h"
/* Scine */
#include <Utils/GeometricDerivatives/NormalModeAnalysis.h>
#include <Utils/GeometryOptimization/GeometryOptimization.h>
#include <Utils/GeometryOptimization/IrcOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Optimizer/GradientBased/Lbfgs.h>
#include <Utils/Optimizer/GradientBased/SteepestDescent.h>
#include <Utils/Optimizer/HessianBased/NewtonRaphson.h>
/* std c++ */
#include <boost/exception/diagnostic_information.hpp>
#include <boost/filesystem.hpp>
#include <cstdio>
#include <fstream>

namespace Scine {
namespace Readuct {

/**
 * @brief A task to run an IRC optimization using a given normal mode.
 */
class IrcTask : public Task {
 public:
  /**
   * @brief Construct a new IrcTask.
   * @param input  The input system names for the task.
   * @param output The output system names for the task.
   * @param logger The logger to/through which all text output will be handled.
   */
  IrcTask(std::vector<std::string> input, std::vector<std::string> output, std::shared_ptr<Core::Log> logger = nullptr)
    : Task(std::move(input), std::move(output), std::move(logger)) {
  }

  std::string name() const override {
    return "IRC Optimizations";
  }

  bool run(SystemsMap& systems, Utils::UniversalSettings::ValueCollection taskSettings, bool testRunOnly = false) const final {
    warningIfMultipleInputsGiven();
    // Warn if only one output system was specified
    if (_output.size() == 1) {
      _logger->warning
          << "  Warning: To store IRC results in a new location two output systems have to be specified.\n";
    }

    bool silentCalculator = taskSettings.extract("silent_stdout_calculator", true);
    std::shared_ptr<Core::Calculator> calc;
    if (!testRunOnly) { // leave out in case of task chaining --> attention calc is NULL
      // Note: _input is guaranteed not to be empty by Task constructor
      calc = copyCalculator(systems, _input.front(), name());
      Utils::CalculationRoutines::setLog(*calc, true, true, !silentCalculator);

      // Check system size
      if (calc->getStructure()->size() == 1) {
        throw std::runtime_error("Cannot perform IRC task for monoatomic systems.");
      }
    }

    // Generate optimizer
    std::string optimizertype = taskSettings.extract("optimizer", std::string{"SD"});
    int mode = taskSettings.extract("irc_mode", 0);
    // Read and delete special settings
    bool stopOnError = stopOnErrorExtraction(taskSettings);
    std::transform(optimizertype.begin(), optimizertype.end(), optimizertype.begin(), ::toupper);
    std::shared_ptr<Utils::IrcOptimizerBase> optimizer;
    if (optimizertype == "LBFGS") {
      auto tmp = std::make_shared<Utils::IrcOptimizer<Utils::Lbfgs>>(*calc);
      // Default convergence options
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "BFGS") {
      auto tmp = std::make_shared<Utils::IrcOptimizer<Utils::Bfgs>>(*calc);
      // Default convergence options
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "SD" || optimizertype == "STEEPESTDESCENT") {
      auto tmp = std::make_shared<Utils::IrcOptimizer<Utils::SteepestDescent>>(*calc);
      // Default convergence options
      optimizer = std::move(tmp);
    }
    else {
      throw std::runtime_error(
          "Unknown Optimizer requested for an IRC optimization, available are: SD, BFGS and LBFGS!");
    }
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

    // Get/find mode
    if (!calc->results().has<Utils::Property::Hessian>()) {
      calc->setRequiredProperties(Utils::Property::Hessian);
      calc->calculate("Hessian Calculation");
    }
    auto hessian = calc->results().get<Utils::Property::Hessian>();
    auto system = calc->getStructure();
    auto modes = Utils::NormalModeAnalysis::calculateNormalModes(hessian, system->getElements(), system->getPositions());
    if (mode < 0 || mode > modes.size() - 1) {
      throw std::runtime_error(
          "The chosen normal mode number is smaller than 0 or larger than the number of modes present!");
    }
    auto ircMode = modes.getMode(mode);
    auto ircModeVector = Eigen::Map<const Eigen::VectorXd>(ircMode.data(), ircMode.cols() * ircMode.rows());

    /*========================*
     *  Forward optimization
     *========================*/
    auto cout = _logger->output;
    cout << "\n    IRC Optimization: Forward\n\n";
    std::shared_ptr<Scine::Core::Calculator> forwardCalc = std::move(calc->clone());
    // Add observer
    double oldEnergy = 0.0;
    Eigen::VectorXd oldParams;
    // Trajectory stream
    using Writer = Utils::XyzStreamHandler;

    const std::string& partialOutput = ((_output.size() > 1) ? _output[0] : _input[0]);
    boost::filesystem::path dirF(partialOutput);
    boost::filesystem::create_directory(dirF);
    boost::filesystem::path trjfileF(partialOutput + ".irc.forward.trj.xyz");
    std::ofstream trajectoryF((dirF / trjfileF).string(), std::ofstream::out);
    auto forward = [&](const int& cycle, const double& energy, const Eigen::VectorXd& params) {
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
      Writer::write(trajectoryF, *structure, std::to_string(energy));
    };
    optimizer->addObserver(forward);

    // Run optimization
    auto structure = systems.at(_input[0])->getStructure();
    int cycles = 0;
    bool forwardAborted = false;
    try {
      cycles = optimizer->optimize(*structure, *_logger, ircModeVector, true);
    }
    catch (...) {
      Writer::write(trajectoryF, *calc->getStructure());
      _logger->error << "Forward IRC failed with error!" << Core::Log::endl;
      if (stopOnError) {
        trajectoryF.close();
        throw;
      }
      _logger->error << boost::current_exception_diagnostic_information() << Core::Log::endl;
      forwardAborted = true;
    }
    trajectoryF.close();

    int maxiter = settings.getInt("convergence_max_iterations");
    const bool forwardConverged = cycles < maxiter && !forwardAborted;
    if (forwardConverged) {
      cout << Core::Log::endl
           << "    Converged after " << cycles << " iterations." << Core::Log::endl
           << Core::Log::endl;
    }
    else if (!forwardAborted) {
      cout << Core::Log::endl
           << "    Stopped after " << maxiter << " iterations." << Core::Log::endl
           << Core::Log::endl;
      if (stopOnError) {
        throw std::runtime_error("Problem: IRC optimization did not converge.");
      }
    }

    // only write separate file with single last frame if no error occurred
    if (!forwardAborted) {
      // Print/Store results
      if (_output.size() > 1) {
        systems[_output[0]] = std::shared_ptr<Core::Calculator>(calc->clone().release());
        boost::filesystem::path xyzfile(_output[0] + ".xyz");
        std::ofstream xyz((dirF / xyzfile).string(), std::ofstream::out);
        Writer::write(xyz, *(calc->getStructure()));
        xyz.close();
      }
      else {
        boost::filesystem::path xyzfile(_input[0] + ".irc.forward.xyz");
        std::ofstream xyz((dirF / xyzfile).string(), std::ofstream::out);
        Writer::write(xyz, *(calc->getStructure()));
        xyz.close();
      }
    }

    /*=========================*
     *  Backward optimization
     *=========================*/
    cout << "\n    IRC Optimization: Backward\n\n";
    // Add observer
    oldEnergy = 0.0;
    oldParams.resize(0);
    // Reset optimizer
    optimizer->setSettings(settings);
    optimizer->reset();
    // Trajectory stream
    boost::filesystem::path dirB(((_output.size() > 1) ? _output[1] : _input[0]));
    boost::filesystem::create_directory(dirB);
    boost::filesystem::path trjfileB(((_output.size() > 1) ? _output[1] : _input[0]) + ".irc.backward.trj.xyz");
    std::ofstream trajectoryB((dirB / trjfileB).string(), std::ofstream::out);
    auto backward = [&](const int& cycle, const double& energy, const Eigen::VectorXd& params) {
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
      Writer::write(trajectoryB, *structure, std::to_string(energy));
    };
    optimizer->clearObservers();
    optimizer->addObserver(backward);

    // Run optimization
    structure = systems.at(_input[0])->getStructure();
    cycles = 0;
    try {
      cycles = optimizer->optimize(*structure, *_logger, ircModeVector, false);
    }
    catch (...) {
      Writer::write(trajectoryB, *calc->getStructure());
      trajectoryB.close();
      _logger->error << "Calculation in IRC backward optimization failed!" << Core::Log::endl;
      if (stopOnError) {
        throw;
      }
      _logger->error << boost::current_exception_diagnostic_information() << Core::Log::endl;
      return false;
    }
    trajectoryB.close();

    maxiter = settings.getInt("convergence_max_iterations");
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
        throw std::runtime_error("Problem: IRC optimization did not converge.");
      }
    }

    // Print/Store results
    if (_output.size() > 1) {
      systems[_output[1]] = std::shared_ptr<Core::Calculator>(calc->clone().release());
      boost::filesystem::path xyzfile(_output[1] + ".xyz");
      std::ofstream xyz((dirB / xyzfile).string(), std::ofstream::out);
      Writer::write(xyz, *(calc->getStructure()));
      xyz.close();
    }
    else {
      boost::filesystem::path xyzfile(_input[0] + ".irc.backward.xyz");
      std::ofstream xyz((dirB / xyzfile).string(), std::ofstream::out);
      Writer::write(xyz, *(calc->getStructure()));
      xyz.close();
    }

    return cycles < maxiter && forwardConverged;
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_IRCTASK_H_
