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
#include "Utils/Geometry/GeometryUtilities.h"
#include "Utils/Geometry/InternalCoordinates.h"
#include <Utils/GeometricDerivatives/HessianUtilities.h>
#include <Utils/GeometricDerivatives/NormalModeAnalyzer.h>
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/Yaml.h>
#include <Utils/Optimizer/GradientBased/Dimer.h>
#include <Utils/Optimizer/HessianBased/Bofill.h>
#include <Utils/Optimizer/HessianBased/EigenvectorFollowing.h>
/* External */
#include <boost/filesystem.hpp>
/* std c++ */
#include <cstdio>
#include <fstream>
#include <iostream>

namespace Scine {
namespace Readuct {

class TsOptimizationTask : public Task {
 public:
  TsOptimizationTask(std::vector<std::string> input, std::vector<std::string> output) : Task(input, output) {
  }

  std::string name() const override {
    return "TS Optimization";
  }

  virtual bool run(std::map<std::string, std::shared_ptr<Core::Calculator>>& systems, const YAML::Node& taskSettings) const final {
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

    if (optimizertype == "DIMER") {
      if (auto type = settingsCopy["dimer_calculate_hessian_once"]) {
        useEigenvectorForFirstStep_ = type.as<bool>();
        settingsCopy.remove("dimer_calculate_hessian_once");
      }
      if (auto type = settingsCopy["dimer_guess_vector_file"]) {
        tangentFileName_ = type.as<std::string>();
        useTangent_ = true;
        settingsCopy.remove("dimer_guess_vector_file");
      }
      if (auto type = settingsCopy["dimer_discrete_guesses"]) {
        useDiscreteGuess_ = true;
        auto fileNames = type.as<std::vector<std::string>>();
        if (fileNames.size() != 2)
          throw std::runtime_error("Problem with discrete guesses. You may only give list with two filepaths");
        filepath1_ = fileNames.at(0);
        filepath2_ = fileNames.at(1);
        settingsCopy.remove("dimer_discrete_guesses");
      }
    }

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
      tmp->check.requirement = 3;
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "EF" || optimizertype == "EV" || optimizertype == "EVF" ||
             optimizertype == "EIGENVECTORFOLLOWING" || optimizertype == "EIGENVECTOR_FOLLOWING") {
      auto tmp = std::make_unique<Utils::GeometryOptimizer<Utils::EigenvectorFollowing>>(*calc);
      tmp->transformCoordinates = true;
      // Default convergence options
      tmp->check.maxIter = 150;
      tmp->check.stepMaxCoeff = 2.0e-3;
      tmp->check.stepRMS = 1.0e-3;
      tmp->check.gradMaxCoeff = 2.0e-4;
      tmp->check.gradRMS = 1.0e-4;
      tmp->check.deltaValue = 1.0e-6;
      tmp->optimizer.trustRadius = 1.0e-1;
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "DIMER") {
      auto tmp = std::make_unique<Utils::GeometryOptimizer<Utils::Dimer>>(*calc);
      // Default convergence options
      tmp->check.maxIter = 150;
      tmp->check.stepMaxCoeff = 2.0e-3;
      tmp->check.stepRMS = 1.0e-3;
      tmp->check.gradMaxCoeff = 2.0e-4;
      tmp->check.gradRMS = 1.0e-4;
      tmp->check.deltaValue = 1.0e-6;
      tmp->check.requirement = 3;

      bool transform = tmp->transformCoordinates;
      if (auto type = settingsCopy["geoopt_transform_coordinates"])
        transform = type.as<bool>();
      /* Generate and pass initial eigenvector and tell optimizer to not perform the first rotation */
      if (useEigenvectorForFirstStep_) {
        tmp->optimizer.skipFirstRotation = true;
        if (!calc->results().has<Utils::Property::Hessian>()) {
          calc->setRequiredProperties(Utils::Property::Hessian);
          calc->calculate("Hessian Calculation");
        }
        auto hessian = calc->results().get<Utils::Property::Hessian>();
        auto system = calc->getStructure();
        auto elements = system->getElements();
        auto positions = system->getPositions();
        bool massWeighted = false;
        Utils::NormalModeAnalyzer nma(hessian, system->getElements(), system->getPositions());
        auto modes = nma.calculateNormalModes();
        auto mode = modes.getMode(0);
        if (transform) {
          Utils::InternalCoordinates internalCoordinates(*system);
          tmp->optimizer.guessVector = std::make_unique<Eigen::VectorXd>(internalCoordinates.coordinatesToInternal(mode));
        }
        else {
          tmp->optimizer.guessVector =
              std::make_unique<Eigen::VectorXd>(Eigen::Map<const Eigen::VectorXd>(mode.data(), mode.cols() * mode.rows()));
        }
      }
      /* Read tangent of spline from file to give as guess vector */
      else if (useTangent_) {
        std::vector<double> container;
        double num = 0.0;
        std::ifstream file(tangentFileName_, std::ios::in);
        if (!file.is_open())
          throw std::runtime_error("Problem when opening file " + tangentFileName_);
        while (file >> num) {
          container.push_back(num);
        }
        file.close();
        Eigen::Map<Eigen::VectorXd> tangentVector(container.data(), container.size());
        if (transform) {
          auto system = calc->getStructure();
          Utils::InternalCoordinates internalCoordinates(*system);
          Utils::PositionCollection tangentPosition =
              Eigen::Map<Utils::PositionCollection>(tangentVector.data(), tangentVector.size() / 3, 3);
          tmp->optimizer.guessVector =
              std::make_unique<Eigen::VectorXd>(internalCoordinates.coordinatesToInternal(tangentPosition));
        }
        else
          tmp->optimizer.guessVector = std::make_unique<Eigen::VectorXd>(tangentVector);
      }
      /* Read 2 structures from files to form guess vector from discrete difference of the two */
      else if (useDiscreteGuess_) {
        Utils::XyzStreamHandler handler;
        Utils::PositionCollection position1;
        Utils::PositionCollection position2;
        Utils::AtomCollection system;
        std::filebuf fb;
        if (fb.open(filepath1_, std::istream::in)) {
          std::istream is(&fb);
          system = handler.read(is);
          position1 = system.getPositions();
          fb.close();
        }
        if (fb.open(filepath2_, std::istream::in)) {
          std::istream is(&fb);
          system = handler.read(is);
          position2 = system.getPositions();
          fb.close();
        }
        Utils::PositionCollection positionDiff = position2 - position1;
        Eigen::Map<Eigen::VectorXd> mode(positionDiff.data(), positionDiff.size());
        if (transform) {
          Utils::InternalCoordinates internalCoordinates(system);
          tmp->optimizer.guessVector =
              std::make_unique<Eigen::VectorXd>(internalCoordinates.coordinatesToInternal(positionDiff));
        }
        else
          tmp->optimizer.guessVector = std::make_unique<Eigen::VectorXd>(mode);
      }
      optimizer = std::move(tmp);
    }
    else {
      throw std::runtime_error(
          "Unknown Optimizer requested for TS optimization, available are: Bofill, EVF and Dimer!");
    }
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
    boost::filesystem::path trjfile(((_output.size() > 0) ? _output[0] : _input[0]) + ".tsopt.trj.xyz");
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
    else if (allowUnconverged) {
      std::cout << std::endl << "    Stopped after " << maxiter << " iterations." << std::endl << std::endl;
    }
    else {
      std::cout << std::endl << "    Stopped after " << maxiter << " iterations." << std::endl << std::endl;
      throw std::runtime_error("Problem: TS optimization did not converge.");
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

  mutable bool useEigenvectorForFirstStep_ = false;
  mutable bool useTangent_ = false;
  mutable bool useDiscreteGuess_ = false;
  mutable std::string filepath1_ = "";
  mutable std::string filepath2_ = "";
  mutable std::string tangentFileName_ = "";
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_TSOPTIMIZATIONTASK_H_
