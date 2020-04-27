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

  virtual bool run(std::map<std::string, std::shared_ptr<Core::Calculator>>& systems, const YAML::Node& taskSettings) const final {
    // Get/Copy Calculator
    std::shared_ptr<Core::Calculator> calc;
    if (systems.find(_input[0]) != systems.end()) {
      calc = systems.at(_input[0]);
    }
    else {
      throw std::runtime_error("Missing system '" + _input[0] + "' in Single Point Task.");
    }

    // Check for charges
    bool chargesAvailable = calc->possibleProperties().containsSubSet(Utils::Property::AtomicCharges);
    if (auto m = taskSettings["require_charges"]) {
      bool requireCharges = m.as<bool>();
      if (requireCharges and !chargesAvailable)
        throw std::runtime_error("Charges required, but chosen calculator does not provide them.");
    }

    // Calculate Energy, and possibly charges
    if (chargesAvailable) {
      calc->setRequiredProperties(Utils::Property::Energy | Utils::Property::AtomicCharges);
    }
    else {
      calc->setRequiredProperties(Utils::Property::Energy);
    }
    try {
      calc->calculate("Energy Calculation");
    }
    catch (...) {
      return false;
    }
    auto energy = calc->results().get<Utils::Property::Energy>();

    // Print Energy
    printf("  The (electronic) energy is: %+16.9f hartree\n\n", energy);
    // Print Charges
    if (chargesAvailable) {
      auto charges = calc->results().get<Utils::Property::AtomicCharges>();
      auto atomColl = calc->getStructure();
      auto elements = atomColl->getElements();
      printf("  Atomic Partial Charges:\n\n", energy);
      for (size_t i = 0; i < charges.size(); i++) {
        printf("  %5d %-6s: %+11.6f\n", i, Utils::ElementInfo::symbol(elements[i]).c_str(), charges[i]);
      }
    }
    printf("\n\n");
    return true;
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_SINGLEPOINTTASK_H_
