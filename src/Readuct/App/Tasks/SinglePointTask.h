/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_SINGLEPOINTTASK_H_
#define READUCT_SINGLEPOINTTASK_H_

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

class SinglePointTask : public Task {
 public:
  SinglePointTask(std::vector<std::string> input, std::vector<std::string> output) : Task(input, output) {
  }

  std::string name() const override {
    return "Singe Point Calculation";
  }

  virtual void run(std::map<std::string, std::shared_ptr<Core::Calculator>>& systems, const YAML::Node& taskSettings) const final {
    // Get/Copy Calculator
    std::shared_ptr<Core::Calculator> calc;
    if (systems.find(_input[0]) != systems.end()) {
      calc = std::shared_ptr<Core::Calculator>(systems.at(_input[0])->clone().release());
    }
    else {
      throw std::runtime_error("Missing system '" + _input[0] + "' in Single Point Task.");
    }

    // Get/calculate Energy
    if (!calc->results().hasEnergy()) {
      calc->setRequiredProperties(Utils::Property::Energy);
      calc->calculate("Energy Calculation");
    }
    auto energy = calc->results().getEnergy();

    // Print
    printf("  The (electronic) energy is: %+16.9f hartree\n\n\n", energy);
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_SINGLEPOINTTASK_H_
