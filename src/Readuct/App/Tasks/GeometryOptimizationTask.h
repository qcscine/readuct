/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_GEOMETRYOPTIMIZATIONTASK_H_
#define READUCT_GEOMETRYOPTIMIZATIONTASK_H_

/* Readuct */
#include "Tasks/Task.h"
/* Scine */
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/Yaml.h>
#include <Utils/Optimizer/GradientBased/Bfgs.h>
#include <Utils/Optimizer/GradientBased/Lbfgs.h>
#include <Utils/Optimizer/GradientBased/SteepestDescent.h>
#include <Utils/Optimizer/HessianBased/NewtonRaphson.h>
/* std c++ */
#include <cstdio>
#include <fstream>

namespace Scine {
namespace Readuct {

class GeometryOptimizationTask : public Task {
 public:
  GeometryOptimizationTask(std::vector<std::string> input, std::vector<std::string> output) : Task(input, output) {
  }

  std::string name() const override {
    return "Geometry Optimization";
  }

  virtual bool run(std::map<std::string, std::shared_ptr<Core::Calculator>>& systems, const YAML::Node& taskSettings) const final {
    // Get/Copy Calculator
    std::shared_ptr<Core::Calculator> calc;
    if (systems.find(_input[0]) != systems.end()) {
      calc = std::shared_ptr<Core::Calculator>(systems.at(_input[0])->clone().release());
    }
    else {
      throw std::runtime_error("Missing system '" + _input[0] + "' in geometry optimization.");
    }

    // Generate optimizer
    std::string optimizertype = "BFGS";
    auto settingsCopy = taskSettings;
    if (auto type = settingsCopy["optimizer"]) {
      optimizertype = type.as<std::string>();
      settingsCopy.remove("optimizer");
    }
    std::transform(optimizertype.begin(), optimizertype.end(), optimizertype.begin(), ::toupper);
    std::unique_ptr<Utils::GeometryOptimizerBase> optimizer;
    if (optimizertype == "LBFGS") {
      auto tmp = std::make_unique<Utils::GeometryOptimizer<Utils::Lbfgs>>(*calc);
      // Default convergence options
      tmp->check.maxIter = 150;
      tmp->check.stepMaxCoeff = 2.0e-3;
      tmp->check.stepRMS = 1.0e-3;
      tmp->check.gradMaxCoeff = 2.0e-4;
      tmp->check.gradRMS = 1.0e-4;
      tmp->check.deltaValue = 1.0e-6;
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "BFGS") {
      auto tmp = std::make_unique<Utils::GeometryOptimizer<Utils::Bfgs>>(*calc);
      // Default convergence options
      tmp->check.maxIter = 150;
      tmp->check.stepMaxCoeff = 2.0e-3;
      tmp->check.stepRMS = 1.0e-3;
      tmp->check.gradMaxCoeff = 2.0e-4;
      tmp->check.gradRMS = 1.0e-4;
      tmp->check.deltaValue = 1.0e-6;
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "SD" || optimizertype == "STEEPESTDESCENT") {
      auto tmp = std::make_unique<Utils::GeometryOptimizer<Utils::SteepestDescent>>(*calc);
      // Default convergence options
      tmp->check.maxIter = 150;
      tmp->check.stepMaxCoeff = 2.0e-3;
      tmp->check.stepRMS = 1.0e-3;
      tmp->check.gradMaxCoeff = 2.0e-4;
      tmp->check.gradRMS = 1.0e-4;
      tmp->check.deltaValue = 1.0e-6;
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "NR" || optimizertype == "NEWTONRAPHSON") {
      auto tmp = std::make_unique<Utils::GeometryOptimizer<Utils::NewtonRaphson>>(*calc);
      // Default convergence options
      tmp->check.maxIter = 150;
      tmp->check.stepMaxCoeff = 2.0e-3;
      tmp->check.stepRMS = 1.0e-3;
      tmp->check.gradMaxCoeff = 2.0e-4;
      tmp->check.gradRMS = 1.0e-4;
      tmp->check.deltaValue = 1.0e-6;
      tmp->transformCoordinates = true;
      optimizer = std::move(tmp);
    }
    else {
      throw std::runtime_error(
          "Unknown Optimizer requested for a geometry optimization, available are: SD, NR, BFGS and LBFGS!");
    }

    // Read and delete special setttings
    bool allowUnconverged = false;
    if (auto m = settingsCopy["allow_unconverged"]) {
      allowUnconverged = m.as<bool>();
      settingsCopy.remove("allow_unconverged");
    }
    // Apply settings
    auto settings = optimizer->getSettings();
    nodeToSettings(settings, settingsCopy);
    optimizer->setSettings(settings);
    // Add observer
    int counter = 0;
    // Trajectory stream
    Utils::XyzStreamHandler writer;
    boost::filesystem::path dir(((_output.size() > 0) ? _output[0] : _input[0]));
    boost::filesystem::create_directory(dir);
    boost::filesystem::path trjfile(((_output.size() > 0) ? _output[0] : _input[0]) + ".opt.trj.xyz");
    std::ofstream trajectory((dir / trjfile).string(), std::ofstream::out);
    double oldEnergy = 0.0;
    Eigen::VectorXd oldParams;
    auto func = [&](const int& cycle, const double& energy, const Eigen::VectorXd& params) {
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
    else if (allowUnconverged) {
      std::cout << std::endl << "    Stopped after " << maxiter << " iterations." << std::endl << std::endl;
    }
    else {
      std::cout << std::endl << "    Stopped after " << maxiter << " iterations." << std::endl << std::endl;
      throw std::runtime_error("Problem: Structure optimization did not converge.");
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
    if (cycles >= maxiter)
      return false;
    return true;
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_GEOMETRYOPTIMIZATIONTASK_H_
