/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_TSOPTIMIZATIONTASK_H_
#define READUCT_TSOPTIMIZATIONTASK_H_

/* Readuct */
#include "Tasks/Task.h"
/* Scine */
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XYZStreamHandler.h>
#include <Utils/IO/Yaml.h>
#include <Utils/Optimizer/GradientBased/Bofill.h>
#include <Utils/Optimizer/HessianBased/EigenvectorFollowing.h>
/* External */
#include <boost/filesystem.hpp>
/* std c++ */
#include <cstdio>
#include <fstream>
#include <iostream>

namespace Scine {
namespace Readuct {

class TSOptimizationTask : public Task {
 public:
  TSOptimizationTask(std::vector<std::string> input, std::vector<std::string> output) : Task(input, output) {
  }

  std::string name() const override {
    return "TS Optimization";
  }

  virtual void run(std::map<std::string, std::shared_ptr<Core::Calculator>>& systems, const YAML::Node& taskSettings) const final {
    // Get/Copy Calculator
    std::shared_ptr<Core::Calculator> calc;
    if (systems.find(_input[0]) != systems.end()) {
      calc = std::shared_ptr<Core::Calculator>(systems.at(_input[0])->clone().release());
    }
    else {
      throw std::runtime_error("Missing system '" + _input[0] + "' in TS optimization.");
    }

    // Generate optimizer
    std::string optimizertype = "BOFILL";
    auto settingsCopy = taskSettings;
    if (auto type = settingsCopy["optimizer"]) {
      optimizertype = type.as<std::string>();
      settingsCopy.remove("optimizer");
    }
    std::transform(optimizertype.begin(), optimizertype.end(), optimizertype.begin(), ::toupper);
    std::unique_ptr<Utils::GeometryOptimizerBase> optimizer;
    if (optimizertype == "BOFILL") {
      auto tmp = std::make_unique<Utils::GeometryOptimizer<Utils::Bofill>>(*calc);
      tmp->transformCoordinates = true;
      // Default convergence options
      tmp->check.maxIter = 150;
      tmp->check.stepMaxCoeff = 2.0e-3;
      tmp->check.stepRMS = 1.0e-3;
      tmp->check.gradMaxCoeff = 2.0e-4;
      tmp->check.gradRMS = 1.0e-4;
      tmp->check.deltaValue = 1.0e-6;
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "EV" || optimizertype == "EVF" || optimizertype == "EIGENVECTORFOLLOWING" ||
             optimizertype == "EIGENVECTOR_FOLLOWING") {
      auto tmp = std::make_unique<Utils::GeometryOptimizer<Utils::EigenvectorFollowing>>(*calc);
      tmp->transformCoordinates = true;
      // Default convergence options
      tmp->check.maxIter = 150;
      tmp->check.stepMaxCoeff = 2.0e-3;
      tmp->check.stepRMS = 1.0e-3;
      tmp->check.gradMaxCoeff = 2.0e-4;
      tmp->check.gradRMS = 1.0e-4;
      tmp->check.deltaValue = 1.0e-6;
      optimizer = std::move(tmp);
    }
    else {
      throw std::runtime_error("Unknown Optimizer requested for TS optimization, available are: Bofill and EVF!");
    }
    // Apply settings
    auto settings = optimizer->getSettings();
    nodeToSettings(settings, settingsCopy);
    optimizer->setSettings(settings);

    // Add observer
    int counter = 0;
    // Trajectory stream
    Utils::XYZStreamHandler writer;
    boost::filesystem::path dir(((_output.size() > 0) ? _output[0] : _input[0]));
    boost::filesystem::create_directory(dir);
    boost::filesystem::path trjfile(((_output.size() > 0) ? _output[0] : _input[0]) + ".tsopt.trj");
    std::ofstream trajectory((dir / trjfile).string(), std::ofstream::out);
    printf("%7s %16s %16s %16s %16s\n", "Cycle", "Energy", "Energy Diff.", "Step RMS", "Max. Step");
    double oldEnergy = 0.0;
    Eigen::VectorXd oldParams;
    auto func = [&](const int& cycle, const double& energy, const Eigen::VectorXd& params) {
      if (oldParams.size() == 0) {
        oldParams = params;
      }
      auto diff = (params - oldParams).eval();
      printf("%7d %+16.9f %+16.9f %+16.9f %+16.9f\n", cycle, energy, energy - oldEnergy,
             sqrt(diff.squaredNorm() / diff.size()), diff.cwiseAbs().maxCoeff());
      oldEnergy = energy;
      oldParams = params;
      auto structure = calc->getStructure();
      writer.write(trajectory, *structure);
    };
    optimizer->addObserver(func);

    // Run optimization
    auto structure = calc->getStructure();
    int cycles = optimizer->optimize(*structure);
    int maxiter = settings.getInt("convergence_max_iterations");
    if (maxiter > cycles) {
      std::cout << std::endl << "    Converged after " << cycles << " iterations." << std::endl << std::endl;
    }
    else {
      std::cout << std::endl << "    Stopped after " << maxiter << " iterations." << std::endl << std::endl;
    }

    // Print/Store results
    trajectory.close();
    if (_output.size() > 0) {
      systems[_output[0]] = calc;
      boost::filesystem::path xyzfile(_output[0] + ".xyz");
      std::ofstream xyz((dir / xyzfile).string(), std::ofstream::out);
      writer.write(xyz, *structure);
      xyz.close();
    }
    else {
      systems[_input[0]] = calc;
      boost::filesystem::path xyzfile(_input[0] + ".xyz");
      std::ofstream xyz((dir / xyzfile).string(), std::ofstream::out);
      writer.write(xyz, *structure);
      xyz.close();
    }
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_TSOPTIMIZATIONTASK_H_
