/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_BONDORDERTASK_H_
#define READUCT_BONDORDERTASK_H_

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

class BondOrderTask : public Task {
 public:
  BondOrderTask(std::vector<std::string> input, std::vector<std::string> output) : Task(input, output) {
  }

  std::string name() const override {
    return "Bond Order Calculation";
  }

  virtual bool run(std::map<std::string, std::shared_ptr<Core::Calculator>>& systems, const YAML::Node& taskSettings) const final {
    // Get/Copy Calculator
    std::shared_ptr<Core::Calculator> calc;
    if (systems.find(_input[0]) != systems.end()) {
      calc = systems.at(_input[0]);
    }
    else {
      throw std::runtime_error("Missing system '" + _input[0] + "' in Bond Order Task.");
    }

    // Get/calculate Energy
    if (!calc->results().has<Utils::Property::BondOrderMatrix>()) {
      calc->setRequiredProperties(Utils::Property::BondOrderMatrix | Utils::Property::Energy);
      try {
        calc->calculate("Bond Order Calculation");
      }
      catch (...) {
        return false;
      }
    }

    auto energy = calc->results().get<Utils::Property::Energy>();
    auto bos = calc->results().get<Utils::Property::BondOrderMatrix>();

    // Print
    printf("  The (electronic) energy is: %+16.9f hartree\n\n\n", energy);

    printf("  Showing all bond orders >0.3 a.u.:\n\n");
    printf("  Atom#1 Atom#2 Bond Order\n\n");
    auto& mat = bos.getMatrix();
    for (int i = 0; i < mat.outerSize(); i++) {
      for (typename Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it) {
        if (it.value() > 0.3 && it.row() < it.col())
          printf("  %6d %6d %+16.9f\n", it.row(), it.col(), it.value());
      }
    }
    printf("\n\n");
    return true;
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_BONDORDERTASK_H_
