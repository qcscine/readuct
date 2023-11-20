/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Tasks/AfirOptimizationTask.h"
#include "Tasks/BSplineInterpolationTask.h"
#include "Tasks/BondOrderTask.h"
#include "Tasks/GeometryOptimizationTask.h"
#include "Tasks/HessianTask.h"
#include "Tasks/IrcTask.h"
#include "Tasks/NtOptimization2Task.h"
#include "Tasks/NtOptimizationTask.h"
#include "Tasks/SinglePointTask.h"
#include "Tasks/Task.h"
#include "Tasks/TsOptimizationTask.h"
#include <Core/Interfaces/Calculator.h>
#include <Python.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine;

namespace {

using SystemsMap = std::map<std::string, std::shared_ptr<Core::Calculator>>;

template<class TaskType>
std::pair<SystemsMap, bool>
run(pybind11::dict& systems, std::vector<std::string> inputNames, bool testRunOnly,
    std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>> observers,
    const pybind11::kwargs& kwargs) {
  Utils::UniversalSettings::ValueCollection settingValues;

  std::vector<std::string> output;
  std::shared_ptr<Core::Log> logger = nullptr;
  for (auto item : kwargs) {
    std::string key = item.first.cast<std::string>();
    if (key == "output") {
      pybind11::list outputNames = item.second.cast<pybind11::list>();
      for (auto item : outputNames) {
        std::string name = item.cast<std::string>();
        output.push_back(name);
      }
    }
    else if (key == "logger") {
      logger = item.second.cast<std::shared_ptr<Core::Log>>();
    }
    else if (pybind11::isinstance<pybind11::bool_>(item.second)) {
      bool value = item.second.cast<bool>();
      settingValues.addBool(key, value);
    }
    else if (pybind11::isinstance<pybind11::str>(item.second)) {
      std::string value = item.second.cast<std::string>();
      settingValues.addString(key, value);
    }
    else if (pybind11::isinstance<pybind11::int_>(item.second)) {
      int value = item.second.cast<int>();
      settingValues.addInt(key, value);
    }
    else if (pybind11::isinstance<pybind11::float_>(item.second)) {
      double value = item.second.cast<double>();
      settingValues.addDouble(key, value);
    }
    else if (pybind11::isinstance<pybind11::list>(item.second)) {
      pybind11::list list = item.second.cast<pybind11::list>();
      if (!list.empty()) {
        // Figure out type from first item
        const auto& firstItem = *list.begin();
        if (pybind11::isinstance<pybind11::int_>(firstItem)) {
          std::vector<int> ints;
          for (const auto& item : list) {
            ints.push_back(item.cast<int>());
          }
          settingValues.addIntList(key, ints);
        }
        else if (pybind11::isinstance<pybind11::str>(firstItem)) {
          std::vector<std::string> strings;
          for (const auto& item : list) {
            strings.push_back(item.cast<std::string>());
          }
          settingValues.addStringList(key, strings);
        }
        else {
          throw std::runtime_error("CollectionLists are unimplemented here!");
        }
      }
    }
    else {
      throw std::runtime_error(key + " could not be converted from Python into a C++ variable, check its type!");
    }
  }

  // Construct task
  TaskType task(inputNames, std::move(output), std::move(logger));

  // we need to make sure that we do not cast a pybind11::none to a shared_ptr<Core::Calculator>,
  // which would segfault, but we want to return the systems map with still the None values.
  auto pybindNone = pybind11::none();
  std::vector<std::string> noneKeys;

  SystemsMap systemsMap;
  for (auto& item : systems) {
    std::string key = item.first.cast<std::string>();
    if (item.second.is(pybindNone)) {
      if (std::find(inputNames.begin(), inputNames.end(), key) != inputNames.end()) {
        throw std::runtime_error("The system '" + key + "' is 'None' in the systems map!");
      }
      noneKeys.push_back(key);
      continue;
    }
    auto cls = item.second.cast<std::shared_ptr<Core::Calculator>>();
    auto sptr = cls->shared_from_this();
    systemsMap.emplace(key, sptr);
  }
  // check if all input systems are present
  for (const auto& name : inputNames) {
    if (systemsMap.find(name) == systemsMap.end()) {
      throw std::runtime_error("The system '" + name + "' is not present in the systems map!");
    }
  }

  // Run the task without GIL if possible
  bool canReleaseGIL = std::all_of(inputNames.begin(), inputNames.end(),
                                   [&systemsMap](const auto& name) { return systemsMap[name]->allowsPythonGILRelease(); }) &&
                       observers.empty();
  if (canReleaseGIL) {
    pybind11::gil_scoped_release release;
  }
  bool success = false;
  try {
    // run the task
    success = task.run(systemsMap, settingValues, testRunOnly, observers);
    if (canReleaseGIL) {
      pybind11::gil_scoped_acquire acquire;
    }
  }
  catch (...) {
    if (canReleaseGIL) {
      pybind11::gil_scoped_acquire acquire;
    }
    throw;
  }

  // add the None values back to the systems map
  for (const auto& key : noneKeys) {
    systemsMap.emplace(key, nullptr);
  }

  return {systemsMap, success};
}

} // namespace

void init_tasks(pybind11::module& m) {
  m.def("run_single_point_task", &run<Readuct::SinglePointTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());
  m.def("run_sp_task", &run<Readuct::SinglePointTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"),
        pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());

  m.def("run_optimization_task", &run<Readuct::GeometryOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());
  m.def("run_opt_task", &run<Readuct::GeometryOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());

  m.def("run_bond_order_task", &run<Readuct::BondOrderTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());

  m.def("run_hessian_task", &run<Readuct::HessianTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());
  m.def("run_freq_task", &run<Readuct::HessianTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"),
        pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());

  m.def("run_transition_state_optimization_task", &run<Readuct::TsOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());
  m.def("run_tsopt_task", &run<Readuct::TsOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());

  m.def("run_irc_task", &run<Readuct::IrcTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"),
        pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());

  m.def("run_nt_task", &run<Readuct::NtOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());

  m.def("run_nt2_task", &run<Readuct::NtOptimization2Task>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());

  m.def("run_afir_task", &run<Readuct::AfirOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());

  m.def("run_bspline_task", &run<Readuct::BSplineInterpolationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false,
        pybind11::arg("observers") =
            std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>());
}
