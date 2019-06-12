/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_IRCTASK_H_
#define READUCT_IRCTASK_H_

/* Readuct */
#include "Tasks/Task.h"
/* Scine */
#include <Utils/GeometryOptimization/IRCOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XYZStreamHandler.h>
#include <Utils/IO/Yaml.h>
#include <Utils/Optimizer/GradientBased/LBFGS.h>
#include <Utils/Optimizer/GradientBased/SteepestDescent.h>
#include <Utils/Optimizer/HessianBased/NewtonRaphson.h>
/* std c++ */
#include <cstdio>
#include <fstream>

namespace Scine {
namespace Readuct {

/**
 * @brief A task to run an IRC optimization using a given normal mode.
 */
class IRCTask : public Task {
 public:
  IRCTask(std::vector<std::string> input, std::vector<std::string> output) : Task(input, output) {
  }

  std::string name() const override {
    return "IRC Optimizations";
  }

  virtual void run(std::map<std::string, std::shared_ptr<Core::Calculator>>& systems, const YAML::Node& taskSettings) const final {
    // Get/Copy Calculator
    std::shared_ptr<Core::Calculator> calc;
    if (systems.find(_input[0]) != systems.end()) {
      calc = std::shared_ptr<Core::Calculator>(systems.at(_input[0])->clone().release());
    }
    else {
      throw std::runtime_error("Missing system '" + _input[0] + "' in IRC optimization.");
    }

    // Generate optimizer
    std::string optimizertype = "LBFGS";
    int mode = 0;
    auto settingsCopy = taskSettings;
    if (auto type = settingsCopy["optimizer"]) {
      optimizertype = type.as<std::string>();
      settingsCopy.remove("optimizer");
    }
    if (auto m = settingsCopy["irc_mode"]) {
      mode = m.as<int>();
      settingsCopy.remove("irc_mode");
    }
    std::transform(optimizertype.begin(), optimizertype.end(), optimizertype.begin(), ::toupper);
    std::unique_ptr<Utils::IRCOptimizerBase> optimizer;
    if (optimizertype == "LBFGS") {
      auto tmp = std::make_unique<Utils::IRCOptimizer<Utils::LBFGS>>(*calc);
      // Default convergence options
      tmp->check.maxIter = 150;
      tmp->check.stepMaxCoeff = 5.0e-3;
      tmp->check.stepRMS = 1.0e-3;
      tmp->check.gradMaxCoeff = 5.0e-4;
      tmp->check.gradRMS = 1.0e-4;
      tmp->check.deltaValue = 1.0e-6;
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "SD" || optimizertype == "STEEPESTDESCENT") {
      auto tmp = std::make_unique<Utils::IRCOptimizer<Utils::SteepestDescent>>(*calc);
      // Default convergence options
      tmp->check.maxIter = 150;
      tmp->check.stepMaxCoeff = 5.0e-3;
      tmp->check.stepRMS = 1.0e-3;
      tmp->check.gradMaxCoeff = 5.0e-4;
      tmp->check.gradRMS = 1.0e-4;
      tmp->check.deltaValue = 1.0e-6;
      optimizer = std::move(tmp);
    }
    else {
      throw std::runtime_error("Unknown Optimizer requested for an IRC optimization, available are: SD and LBFGS!");
    }
    // Apply user settings
    auto settings = optimizer->getSettings();
    nodeToSettings(settings, settingsCopy);
    optimizer->setSettings(settings);

    // Get/find mode
    if (!calc->results().hasHessian()) {
      calc->setRequiredProperties(Utils::Property::Hessian);
      calc->calculate("Hessian Calculation");
    }
    auto hessian = calc->results().getHessian();
    auto system = calc->getStructure();
    Utils::NormalModeAnalyzer nma(hessian, system->getElements(), system->getPositions());
    auto modes = nma.calculateNormalModes();
    if (mode < 0 || mode > modes.size() - 1) {
      throw std::runtime_error(
          "The chosen normal mode number is smaller than 0 or larger than the number of modes present!");
    }
    auto ircMode = modes.getMode(mode);
    auto ircModeVector = Eigen::Map<const Eigen::VectorXd>(ircMode.data(), ircMode.cols() * ircMode.rows());

    /*========================*
     *  Forward optimization
     *========================*/
    printf("\n    IRC Optimization: Forward\n\n");
    // Add observer
    int counter = 0;
    double oldEnergy = 0.0;
    Eigen::VectorXd oldParams;
    // Trajectory stream
    Utils::XYZStreamHandler writer;
    boost::filesystem::path dirF(((_output.size() > 1) ? _output[0] : _input[0]));
    boost::filesystem::create_directory(dirF);
    boost::filesystem::path trjfileF(((_output.size() > 1) ? _output[0] : _input[0]) + ".irc.forward.trj");
    std::ofstream trajectoryF((dirF / trjfileF).string(), std::ofstream::out);
    auto forward = [&](const int& cycle, const double& energy, const Eigen::VectorXd& params) {
      if (oldParams.size() == 0) {
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
      writer.write(trajectoryF, *structure);
    };
    optimizer->addObserver(forward);

    // Run optimization
    auto structure = systems.at(_input[0])->getStructure();
    int cycles = optimizer->optimize(*structure, ircModeVector, true);
    int maxiter = settings.getInt("convergence_max_iterations");
    if (maxiter > cycles) {
      std::cout << std::endl << "    Converged after " << cycles << " iterations." << std::endl << std::endl;
    }
    else {
      std::cout << std::endl << "    Stopped after " << maxiter << " iterations." << std::endl << std::endl;
    }

    // Print/Store results
    trajectoryF.close();
    if (_output.size() > 1) {
      systems[_output[0]] = calc;
      boost::filesystem::path xyzfile(_output[0] + ".xyz");
      std::ofstream xyz((dirF / xyzfile).string(), std::ofstream::out);
      writer.write(xyz, *structure);
      xyz.close();
    }
    else {
      boost::filesystem::path xyzfile(_input[0] + ".irc.forward.xyz");
      std::ofstream xyz((dirF / xyzfile).string(), std::ofstream::out);
      writer.write(xyz, *structure);
      xyz.close();
    }

    /*=========================*
     *  Backward optimization
     *=========================*/
    printf("\n    IRC Optimization: Backward\n\n");
    // Add observer
    counter = 0;
    oldEnergy = 0.0;
    oldParams.resize(0);
    // Trajectory stream
    boost::filesystem::path dirB(((_output.size() > 1) ? _output[1] : _input[0]));
    boost::filesystem::create_directory(dirB);
    boost::filesystem::path trjfileB(((_output.size() > 1) ? _output[1] : _input[0]) + ".irc.backward.trj");
    std::ofstream trajectoryB((dirB / trjfileB).string(), std::ofstream::out);
    auto backward = [&](const int& cycle, const double& energy, const Eigen::VectorXd& params) {
      if (oldParams.size() == 0) {
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
      writer.write(trajectoryB, *structure);
    };
    optimizer->clearObservers();
    optimizer->addObserver(backward);

    // Run optimization
    structure = systems.at(_input[0])->getStructure();
    cycles = optimizer->optimize(*structure, ircModeVector, false);
    maxiter = settings.getInt("convergence_max_iterations");
    if (maxiter > cycles) {
      std::cout << std::endl << "    Converged after " << cycles << " iterations." << std::endl << std::endl;
    }
    else {
      std::cout << std::endl << "    Stopped after " << maxiter << " iterations." << std::endl << std::endl;
    }

    // Print/Store results
    trajectoryB.close();
    if (_output.size() > 1) {
      systems[_output[1]] = calc;
      boost::filesystem::path xyzfile(_output[1] + ".xyz");
      std::ofstream xyz((dirB / xyzfile).string(), std::ofstream::out);
      writer.write(xyz, *structure);
      xyz.close();
    }
    else {
      boost::filesystem::path xyzfile(_input[0] + ".irc.backward.xyz");
      std::ofstream xyz((dirB / xyzfile).string(), std::ofstream::out);
      writer.write(xyz, *structure);
      xyz.close();
    }
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_IRCTASK_H_
