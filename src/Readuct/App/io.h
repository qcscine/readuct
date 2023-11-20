/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_IO_H_
#define READUCT_IO_H_

/* Readuct */
#include "Tasks/TaskFactory.h"
/* Scine */
#include <Core/Interfaces/Calculator.h>
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Utils/CalculatorBasics/CalculationRoutines.h>
#include <Utils/Geometry.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/Yaml.h>

/* Boost program arguments */
#include "boost/exception/diagnostic_information.hpp"

/* C++ std */
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace Scine {
namespace Readuct {

std::tuple<std::map<std::string, std::tuple<std::string, std::string, std::shared_ptr<Core::Calculator>>>,
           std::vector<std::shared_ptr<Task>>, std::vector<Utils::UniversalSettings::ValueCollection>>
loadYamlFile(const std::string& filename) {
  // Initialize map for all systems (name - their calculators)
  std::map<std::string, std::tuple<std::string, std::string, std::shared_ptr<Core::Calculator>>> systems;

  // Initialize Logger
  auto logger = std::make_shared<Core::Log>();

  auto input = YAML::LoadFile(filename);
  // Check for invalid top level input sections
  std::vector<std::string> keywords{"systems", "tasks"};
  Utils::checkYamlKeyRecognition(input, keywords);

  // Load systems
  auto systemsInput = input["systems"];
  for (size_t i = 0; i < systemsInput.size(); i++) {
    // Check for invalid system input sections
    std::vector<std::string> keywords{"name", "method_family", "path", "program", "settings"};
    Utils::checkYamlKeyRecognition(systemsInput[i], keywords);

    auto current = systemsInput[i];
    if (!current["name"]) {
      throw std::logic_error("System no. " + std::to_string(i + 1) + " is missing a name.\n");
    }
    std::string name = current["name"].as<std::string>();
    if (!current["method_family"]) {
      throw std::logic_error("A method_family is missing for the system: '" + name + "'.\n");
    }
    std::string methodFamily = current["method_family"].as<std::string>();
    if (!current["path"]) {
      throw std::logic_error("An input path is missing for the system: '" + name + "'.\n");
    }
    std::string path = current["path"].as<std::string>();
    std::string program;
    if (auto node = current["program"]) {
      program = current["program"].as<std::string>();
    }
    // Load molecule
    auto readResults = Utils::ChemicalFileHandler::read(path);

    // Generate Calculator
    auto calc = Utils::CalculationRoutines::getCalculator(methodFamily, program);
    // Apply settings to Calculator
    if (auto settingsnode = current["settings"]) {
      nodeToSettings(calc->settings(), settingsnode);
    }
    // Set initial structure
    calc->setStructure(readResults.first);
    if (!calc->settings().valid()) {
      calc->settings().throwIncorrectSettings();
    }
    if (program.empty()) {
      program = "Any";
    }
    Utils::CalculationRoutines::inputPreparation(methodFamily, program);
    systems[name] = std::make_tuple(methodFamily, program, calc);
  }

  // Generate Tasks
  std::vector<std::shared_ptr<Task>> tasks;
  std::vector<Utils::UniversalSettings::ValueCollection> tasksettings;
  auto tasksInput = input["tasks"];
  for (auto current : tasksInput) {
    // Check for invalid task input sections
    std::vector<std::string> keywords{"type", "input", "output", "settings"};
    Utils::checkYamlKeyRecognition(current, keywords);
    std::string type = current["type"].as<std::string>();
    // Parse inputs
    std::vector<std::string> inputs;
    auto inputnames = current["input"];
    if (inputnames.size() == 0) {
      throw std::logic_error("Missing system (input name) in task " + type + "\n");
    }
    for (auto inputname : inputnames) {
      inputs.push_back(inputname.as<std::string>());
    }
    // Parse outputs if present
    std::vector<std::string> outputs;
    if (auto outputsnames = current["output"]) {
      for (auto outputname : outputsnames) {
        outputs.push_back(outputname.as<std::string>());
      }
    }
    // Generate task
    tasks.emplace_back(TaskFactory::produce(type, inputs, outputs, logger));
    // Get task settings
    if (current["settings"]) {
      tasksettings.push_back(Utils::deserializeValueCollection(current["settings"]));
    }
    else {
      tasksettings.emplace_back();
    }
  }

  // Test task dependencies
  std::vector<std::string> existing;
  existing.reserve(systems.size());
  for (const auto& s : systems) {
    existing.push_back(s.first);
  }
  for (size_t i = 0; i < tasks.size(); i++) {
    for (const auto& required : tasks[i]->input()) {
      if (std::find(existing.begin(), existing.end(), required) == existing.end()) {
        throw std::logic_error("Task No. " + std::to_string(i + 1) + " requires a system named '" + required +
                               "' which won't be present at that stage.\n");
      }
    }
    for (const auto& generated : tasks[i]->output()) {
      existing.push_back(generated);
    }
  }

  // Run tasks in testing mode to check all settings without calculations and writing files
  // if this goes through all basic settings are correct, otherwise we hopefully have a meaningful error in task
  std::map<std::string, std::shared_ptr<Core::Calculator>> calcs;
  for (const auto& s : systems) {
    calcs[s.first] = std::get<2>(s.second);
  }
  for (size_t i = 0; i < tasks.size(); i++) {
    auto name = tasks[i]->name();
    try {
      tasks[i]->run(calcs, tasksettings[i], true);
    }
    catch (...) {
      throw std::logic_error("Encountered the following error in task " + name + ":\n" +
                             boost::current_exception_diagnostic_information());
    }
  }
  return std::make_tuple(systems, tasks, tasksettings);
}

} // namespace Readuct
} // namespace Scine

#endif // READUCT_IO_H_
