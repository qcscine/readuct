/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Tasks/AFIROptimizationTask.h"
#include "Tasks/GeometryOptimizationTask.h"
#include "Tasks/HessianTask.h"
#include "Tasks/IRCTask.h"
#include "Tasks/SinglePointTask.h"
#include "Tasks/TSOptimizationTask.h"
#include "Tasks/Task.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/IO/Yaml.h>
#include <pybind11/pybind11.h>
#include <iostream>
// #include <pybind11/stl.h>

using namespace Scine;

template<class TaskName>
pybind11::dict run(pybind11::dict inputSystems, pybind11::list inputNames, pybind11::kwargs kwargs) {
  // Read all systems
  std::map<std::string, std::shared_ptr<Core::Calculator>> systems;
  for (auto item : inputSystems) {
    std::string key = item.first.cast<std::string>();
    std::shared_ptr<Core::Calculator> value = item.second.cast<std::shared_ptr<Core::Calculator>>();
    systems[key] = value;
  }
  // Read all input system names
  std::vector<std::string> input;
  for (auto item : inputNames) {
    std::string name = item.cast<std::string>();
    input.push_back(name);
  }
  // Read all settings
  YAML::Node settingsnode;
  std::vector<std::string> output;
  for (auto item : kwargs) {
    std::string key = item.first.cast<std::string>();
    if (!key.compare("output")) {
      pybind11::list outputNames = item.second.cast<pybind11::list>();
      for (auto item : outputNames) {
        std::string name = item.cast<std::string>();
        output.push_back(name);
      }
    }
    else {
      std::string value = pybind11::str(item.second).cast<std::string>();
      settingsnode[key] = value;
    }
  }
  // Run
  auto task = std::make_unique<TaskName>(input, output);
  task->run(systems, settingsnode);
  // Update systems
  for (auto item : systems) {
    inputSystems[pybind11::cast(item.first)] = item.second;
  }

  return inputSystems;
}

void init_tasks(pybind11::module& m) {
  m.def("run_single_point_task", &run<Readuct::SinglePointTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"));
  m.def("run_sp_task", &run<Readuct::SinglePointTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"));
  m.def("run_optimization_task", &run<Readuct::GeometryOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"));
  m.def("run_opt_task", &run<Readuct::GeometryOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"));
  m.def("run_hessian_task", &run<Readuct::HessianTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"));
  m.def("run_freq_task", &run<Readuct::HessianTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"));
  m.def("run_transition_state_optimization_task", &run<Readuct::TSOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"));
  m.def("run_tsopt_task", &run<Readuct::TSOptimizationTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"));
  m.def("run_irc_task", &run<Readuct::IRCTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"));
  m.def("run_afir_task", &run<Readuct::AFIROptimizationTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"));
}
