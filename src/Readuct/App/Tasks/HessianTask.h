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
#include <Utils/IO/Yaml.h>
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

  virtual void run(std::map<std::string, std::shared_ptr<Core::Calculator>>& systems, const YAML::Node& taskSettings) const final {
    // Get Calculator
    std::shared_ptr<Core::Calculator> calc;
    if (systems.find(_input[0]) != systems.end()) {
      calc = systems.at(_input[0]);
    }
    else {
      throw std::runtime_error("Missing system '" + _input[0] + "' in Hessian Task.");
    }

    // Get/calculate Hessian
    if (!calc->results().hasHessian()) {
      calc->setRequiredProperties(Utils::Property::Hessian);
      calc->calculate("Hessian Calculation");
    }
    auto hessian = calc->results().getHessian();

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
    std::cout << std::endl << std::endl;
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_HESSIANTASK_H_
