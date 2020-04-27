/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Readuct */
#include "Tasks/TaskFactory.h"
/* Scine */
#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Utils/Geometry.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/Logger.h>
#include <Utils/IO/Yaml.h>

/* Boost program arguments */
#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

/* C++ std */
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace Scine;
using namespace Scine::Readuct;

int main(int argc, char* argv[]) {
  // Initialize global map for all systems (name - their calculators)
  std::map<std::string, std::shared_ptr<Core::Calculator>> systems;

  namespace po = boost::program_options;

  // Arguments
  po::options_description optionsDescription("Recognized options");
  optionsDescription.add_options()("help,h", "Produce this help message")(
      "log,l", po::value<std::string>()->default_value("info"),
      "Sets the logging verbosity level")("config,c", po::value<std::string>(), "YAML input file to read");

  po::variables_map optionsMap;
  po::positional_options_description positionalDescription;
  positionalDescription.add("config", 1);
  positionalDescription.add("log", 1);
  po::store(po::command_line_parser(argc, argv)
                .options(optionsDescription)
                .positional(positionalDescription)
                .style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise)
                .run(),
            optionsMap);
  po::notify(optionsMap);

  // Handle help
  if (optionsMap.count("help") > 0 || optionsMap.count("config") == 0) {
    std::cout << optionsDescription << std::endl;
    return 1;
  }

  // Load input file
  const std::string filename = optionsMap["config"].as<std::string>();
  if (!boost::filesystem::exists(filename)) {
    std::cout << "Specified config file does not exist!\n";
    return 1;
  }

  auto input = YAML::LoadFile(filename);

  // Handle logging
  const std::string loggingVerbosity = optionsMap["log"].as<std::string>();
  Utils::Log::startConsoleLogging(loggingVerbosity);

  // Header
  std::cout << R"(#=============================================================================#)" << std::endl;
  std::cout << R"(|   ______   _______     ___       _______   __    __    ______  __________   |)" << std::endl;
  std::cout << R"(|  |   _  \ |   ____|   /   \     |       \ |  |  |  |  /      ||          |  |)" << std::endl;
  std::cout << R"(|  |  |_)  ||  |__     /  ^  \    |  .--.  ||  |  |  | |  ,----'`---|  |---`  |)" << std::endl;
  std::cout << R"(|  |      / |   __|   /  /_\  \   |  |  |  ||  |  |  | |  |         |  |      |)" << std::endl;
  std::cout << R"(|  |  |\  \ |  |____ /  _____  \  |  '--'  ||  `--'  | |  `----.    |  |      |)" << std::endl;
  std::cout << R"(|  | _| `._||_______/__/     \__\ |_______/  \______/   \______|    |__|      |)" << std::endl;
  std::cout << R"(|                                                                             |)" << std::endl;
  std::cout << R"(#=============================================================================#)" << std::endl;
  std::cout << std::endl << std::endl;

  // Load module manager
  auto& manager = Core::ModuleManager::getInstance();

  // Load systems
  auto systemsInput = input["systems"];
  for (size_t i = 0; i < systemsInput.size(); i++) {
    auto current = systemsInput[i];
    if (!current["name"]) {
      std::cout << "System no. " << i + 1 << " is missing a name.\n";
      return 1;
    }
    std::string name = current["name"].as<std::string>();
    if (!current["method_family"]) {
      std::cout << "A method_family is missing for the system: '" << name << "'.\n";
      return 1;
    }
    std::string method_family = current["method_family"].as<std::string>();
    if (!current["path"]) {
      std::cout << "An input path is missing for the system: '" << name << "'.\n";
      return 1;
    }
    std::string path = current["path"].as<std::string>();
    std::string program;
    if (auto node = current["program"]) {
      program = current["program"].as<std::string>();
    }
    // Load molecule
    auto readResults = Utils::ChemicalFileHandler::read(path);
    // Generate Calculator
    std::shared_ptr<Core::Calculator> calc;
    try {
      for (auto& x : program)
        x = std::tolower(x);
      program[0] = std::toupper(program[0]);
      for (auto& x : method_family)
        x = std::toupper(x);
      calc = manager.get<Core::Calculator>(Core::Calculator::supports(method_family), program);
    }
    catch (...) {
      if (program.empty()) {
        std::cout << "No SCINE module providing '" << method_family << "' is currently loaded.\n";
        std::cout << "Please add the module to the SCINE_MODULE_PATH in order for it to be accessible.\n";
        return 1;
      }
      else {
        std::cout << "No SCINE module named '" << program << "' providing '" << method_family << "' is currently loaded.\n";
        std::cout << "Please add the module to the SCINE_MODULE_PATH in order for it to be accessible.\n";
        return 1;
      }
    }
    // Apply settings to Calculator
    if (auto settingsnode = current["settings"]) {
      nodeToSettings(calc->settings(), settingsnode);
    }
    // Set initial structure
    calc->setStructure(readResults.first);
    systems[name] = calc;
  }

  // Generate Tasks
  std::vector<std::unique_ptr<Task>> tasks;
  std::vector<YAML::Node> tasksettings;
  auto tasksInput = input["tasks"];
  for (size_t i = 0; i < tasksInput.size(); i++) {
    auto current = tasksInput[i];
    std::string type = current["type"].as<std::string>();
    // Parse inputs
    std::vector<std::string> inputs;
    auto inputnames = current["input"];
    if (inputnames.size() == 0) {
      std::cout << "Missing system (input name) in task " + type + "\n";
      return 1;
    }
    for (size_t i = 0; i < inputnames.size(); i++) {
      inputs.push_back(inputnames[i].as<std::string>());
    }
    // Parse outputs if present
    std::vector<std::string> outputs;
    if (auto outputsnames = current["output"]) {
      for (size_t i = 0; i < outputsnames.size(); i++) {
        outputs.push_back(outputsnames[i].as<std::string>());
      }
    }
    // Generate task
    tasks.emplace_back(TaskFactory::produce(type, inputs, outputs));
    // Get task settings
    if (current["settings"]) {
      tasksettings.push_back(current["settings"]);
    }
    else {
      tasksettings.push_back(YAML::Node());
    }
  }

  // Test task dependencies
  std::vector<std::string> existing;
  for (auto const& s : systems)
    existing.push_back(s.first);
  for (size_t i = 0; i < tasks.size(); i++) {
    for (auto& required : tasks[i]->input()) {
      if (std::find(existing.begin(), existing.end(), required) == existing.end()) {
        std::cout << "Task No. " << i + 1 << " requires a system named '" << required
                  << "' which won't be present at that stage." << std::endl;
        return 1;
      }
    }
    for (auto& generated : tasks[i]->output()) {
      existing.push_back(generated);
    }
  }

  // Run tasks
  for (size_t i = 0; i < tasks.size(); i++) {
    auto name = tasks[i]->name();
    printf("  #============#    #==%s==#\n", std::string(name.size(), '=').c_str());
    printf("  |  Task %3d  |    |  %s  |\n", static_cast<int>(i + 1), name.c_str());
    printf("  #============#    #==%s==#\n", std::string(name.size(), '=').c_str());
    std::cout << std::endl << std::endl;
    const auto& inpt = tasks[i]->input();
    const auto& outpt = tasks[i]->output();

    std::cout << "  Input System(s):  " + tasks[i]->input()[0] << std::endl;
    for (unsigned int i = 1; i < inpt.size(); i++) {
      std::cout << "                    " + inpt[i] << std::endl;
    }
    std::cout << std::endl;
    if (outpt.size() > 0) {
      std::cout << "  Output System(s): " + outpt[0] << std::endl;
      for (unsigned int i = 1; i < outpt.size(); i++) {
        std::cout << "                    " + outpt[i] << std::endl;
      }
      std::cout << std::endl;
    }
    tasks[i]->run(systems, tasksettings[i]);
  }

  std::cout << "  Completed all tasks, exiting!" << std::endl;
  return 0;
}
