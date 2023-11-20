/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "io.h"
#include <Core/Interfaces/Calculator.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <map>
#include <memory>
#include <string>
#include <tuple>

using namespace Scine;

void init_io(pybind11::module& m) {
  m.def(
      "load_yaml",
      [&](const std::string& filename)
          -> std::tuple<std::map<std::string, std::tuple<std::string, std::string, std::shared_ptr<Core::Calculator>>>,
                        std::vector<std::tuple<std::string, std::vector<std::string>, std::vector<std::string>>>,
                        std::vector<Utils::UniversalSettings::ValueCollection>> {
        const auto& [systems, tasks, tasksettings] = Readuct::loadYamlFile(filename);
        std::vector<std::tuple<std::string, std::vector<std::string>, std::vector<std::string>>> ioPairs;
        for (const auto& t : tasks) {
          ioPairs.push_back(std::make_tuple(t->name(), t->input(), t->output()));
        }
        return std::make_tuple(systems, ioPairs, tasksettings);
      },
      "Load the system map, the task, input, and output names and the task settings from a ReaDuct yaml input file");
}
