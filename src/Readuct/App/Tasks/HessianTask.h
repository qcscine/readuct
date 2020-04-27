/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_HESSIANTASK_H_
#define READUCT_HESSIANTASK_H_

/* Readuct */
#include "Tasks/Task.h"
/* Scine */
#include <Utils/GeometricDerivatives/NormalModeAnalyzer.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/IO/FormattedString.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/IO/Yaml.h>
#include <Utils/Properties/Thermochemistry/ThermochemistryCalculator.h>
/* External */
#include <boost/filesystem.hpp>
/* std c++ */
#include <cstdio>
#include <fstream>
#include <iostream>

namespace Scine {
namespace Readuct {

class HessianTask : public Task {
 public:
  HessianTask(std::vector<std::string> input, std::vector<std::string> output) : Task(input, output) {
  }

  std::string name() const override {
    return "Hessian Calculation";
  }

  virtual bool run(std::map<std::string, std::shared_ptr<Core::Calculator>>& systems, const YAML::Node& taskSettings) const final {
    // Get Calculator
    std::shared_ptr<Core::Calculator> calc;
    if (systems.find(_input[0]) != systems.end()) {
      calc = systems.at(_input[0]);
    }
    else {
      throw std::runtime_error("Missing system '" + _input[0] + "' in Hessian Task.");
    }

    // Setup output directory
    boost::filesystem::path dir(_input[0]);
    boost::filesystem::create_directory(dir);

    // Get/calculate Hessian
    if (!calc->results().has<Utils::Property::Hessian>() or !calc->results().has<Utils::Property::Thermochemistry>()) {
      calc->setRequiredProperties(Utils::Property::Hessian | Utils::Property::Thermochemistry);
      calc->calculate("Hessian Calculation");
    }
    auto hessian = calc->results().get<Utils::Property::Hessian>();

    // Print frequencies
    auto system = calc->getStructure();
    Utils::NormalModeAnalyzer nma(hessian, system->getElements(), system->getPositions());
    auto modes = nma.calculateNormalModes();
    auto wavenumbers = modes.getWaveNumbers();
    printf("  Vib. Frequencies:\n");
    printf("  [Rot. and trans. vib. removed, imaginary freq. shown as negative.]\n\n");
    printf("  %4s %8s\n", "#", "cm^-1");
    for (unsigned int i = 0; i < wavenumbers.size(); i++) {
      printf("  %4d %+8.1f\n", i + 1, wavenumbers[i]);
    }
    std::cout << std::endl;

    // Print Imaginary modes
    for (unsigned int i = 0; i < wavenumbers.size(); i++) {
      if (!(wavenumbers[i] < 0.0))
        break;
      auto trj = modes.getModeAsMolecularTrajectory(i, *system, 1.0);
      std::string filename = Utils::format("%s.vibmode-%05d.trj.xyz", _input[0].c_str(), i + 1);
      std::ofstream xyz((dir / filename).string(), std::ofstream::out);
      Utils::MolecularTrajectoryIO::write(Utils::MolecularTrajectoryIO::format::xyz, xyz, trj);
      xyz.close();
    }

    // Print thermochemistry results
    auto thermo = calc->results().get<Utils::Property::Thermochemistry>();
    printf("  Thermochemical Data:\n");
    printf("  %35s %+14.9f K\n", "Temperature (T):", thermo.overall.temperature);
    printf("  %35s %+14.9f E_h\n", "Zero Point Vibrational Energy:", thermo.overall.zeroPointVibrationalEnergy);
    printf("  %35s %+14.9f E_h/K\n", "Entropy Correction (dS):", thermo.overall.entropy);
    printf("  %35s %+14.9f E_h\n", "Entropy Correction (-TdS):", -thermo.overall.temperature * thermo.overall.entropy);
    printf("  %35s %+14.9f E_h\n", "Enthalpy Correction (dH):", thermo.overall.enthalpy);
    printf("  %35s %+14.9f E_h\n", "Gibbs Energy Correction (dG):", thermo.overall.gibbsFreeEnergy);
    if (calc->results().has<Utils::Property::Energy>()) {
      printf("  %35s %+14.9f E_h\n", "Electronic Energy (E):", calc->results().get<Utils::Property::Energy>());
      printf("  %35s %+14.9f E_h\n",
             "Gibbs Energy (E + dG):", thermo.overall.gibbsFreeEnergy + calc->results().get<Utils::Property::Energy>());
    }
    std::cout << std::endl << std::endl;
    return true;
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_HESSIANTASK_H_
