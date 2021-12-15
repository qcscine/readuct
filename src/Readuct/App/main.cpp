/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Readuct */
#include "Tasks/TaskFactory.h"
/* Scine */
#include <Core/Interfaces/Calculator.h>
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Utils/Geometry.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/Yaml.h>

/* Boost program arguments */
#include "boost/exception/diagnostic_information.hpp"
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
  optionsDescription.add_options()("help,h", "Produce this help message")("config,c", po::value<std::string>(),
                                                                          "YAML input file to read");

  po::variables_map optionsMap;
  po::positional_options_description positionalDescription;
  positionalDescription.add("config", 1);
  po::store(po::command_line_parser(argc, argv)
                .options(optionsDescription)
                .positional(positionalDescription)
                .style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise)
                .run(),
            optionsMap);
  po::notify(optionsMap);

  // Initialize Logger
  auto logger = std::make_shared<Core::Log>();
  auto cout = logger->output;

  // Handle help
  if (optionsMap.count("help") > 0 || optionsMap.count("config") == 0) {
    cout << optionsDescription << Core::Log::endl;
    return 1;
  }

  // Load input file
  const std::string filename = optionsMap["config"].as<std::string>();
  if (!boost::filesystem::exists(filename)) {
    cout << "Specified config file does not exist!\n";
    return 1;
  }

  auto input = YAML::LoadFile(filename);

  // Check for invalid top level input sections
  std::vector<std::string> keywords{"systems", "tasks"};
  Scine::Utils::checkYamlKeyRecognition(input, keywords);

  // Header
  cout << R"(#=============================================================================#)" << Core::Log::endl;
  cout << R"(|   ______   _______     ___       _______   __    __    ______  __________   |)" << Core::Log::endl;
  cout << R"(|  |   _  \ |   ____|   /   \     |       \ |  |  |  |  /      ||          |  |)" << Core::Log::endl;
  cout << R"(|  |  |_)  ||  |__     /  ^  \    |  .--.  ||  |  |  | |  ,----'`---|  |---`  |)" << Core::Log::endl;
  cout << R"(|  |      / |   __|   /  /_\  \   |  |  |  ||  |  |  | |  |         |  |      |)" << Core::Log::endl;
  cout << R"(|  |  |\  \ |  |____ /  _____  \  |  '--'  ||  `--'  | |  `----.    |  |      |)" << Core::Log::endl;
  cout << R"(|  | _| `._||_______/__/     \__\ |_______/  \______/   \______|    |__|      |)" << Core::Log::endl;
  cout << R"(|                                                                             |)" << Core::Log::endl;
  cout << R"(#=============================================================================#)" << Core::Log::endl;
  cout << Core::Log::endl << Core::Log::endl;

  // Load module manager
  auto& manager = Core::ModuleManager::getInstance();

  // Load systems
  auto systemsInput = input["systems"];
  for (size_t i = 0; i < systemsInput.size(); i++) {
    // Check for invalid system input sections
    std::vector<std::string> keywords{"name", "method_family", "path", "program", "settings"};
    Scine::Utils::checkYamlKeyRecognition(systemsInput[i], keywords);

    auto current = systemsInput[i];
    if (!current["name"]) {
      cout << "System no. " << i + 1 << " is missing a name.\n";
      return 1;
    }
    std::string name = current["name"].as<std::string>();
    if (!current["method_family"]) {
      cout << "A method_family is missing for the system: '" << name << "'.\n";
      return 1;
    }
    std::string method_family = current["method_family"].as<std::string>();
    if (!current["path"]) {
      cout << "An input path is missing for the system: '" << name << "'.\n";
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
      for (auto& x : program) {
        x = std::tolower(x);
      }
      program[0] = std::toupper(program[0]);
      for (auto& x : method_family) {
        x = std::toupper(x);
      }
      calc = manager.get<Core::Calculator>(Core::Calculator::supports(method_family), program);
    }
    catch (...) {
      if (program.empty()) {
        cout << "No SCINE module providing '" << method_family << "' is currently loaded.\n";
        cout << "Please add the module to the SCINE_MODULE_PATH in order for it to be accessible.\n";
        return 1;
      }
      cout << "No SCINE module named '" << program << "' providing '" << method_family << "' is currently loaded.\n";
      cout << "Please add the module to the SCINE_MODULE_PATH in order for it to be accessible.\n";
      return 1;
    }
    // Apply settings to Calculator
    if (auto settingsnode = current["settings"]) {
      nodeToSettings(calc->settings(), settingsnode);
    }
    // Set initial structure
    calc->setStructure(readResults.first);
    if (!calc->settings().valid()) {
      calc->settings().throwIncorrectSettings();
    }
    systems[name] = calc;
  }

  // Generate Tasks
  std::vector<std::unique_ptr<Task>> tasks;
  std::vector<Utils::UniversalSettings::ValueCollection> tasksettings;
  auto tasksInput = input["tasks"];
  for (auto current : tasksInput) {
    // Check for invalid task input sections
    std::vector<std::string> keywords{"type", "input", "output", "settings"};
    Scine::Utils::checkYamlKeyRecognition(current, keywords);
    std::string type = current["type"].as<std::string>();
    // Parse inputs
    std::vector<std::string> inputs;
    auto inputnames = current["input"];
    if (inputnames.size() == 0) {
      cout << "Missing system (input name) in task " + type + "\n";
      return 1;
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
        cout << "Task No. " << i + 1 << " requires a system named '" << required
             << "' which won't be present at that stage." << Core::Log::endl;
        return 1;
      }
    }
    for (const auto& generated : tasks[i]->output()) {
      existing.push_back(generated);
    }
  }

  // Run tasks in testing mode to check all settings without calculations and writing files
  // if this goes through all basic settings are correct, otherwise we hopefully have a meaningful error in task
  for (size_t i = 0; i < tasks.size(); i++) {
    auto name = tasks[i]->name();
    try {
      tasks[i]->run(systems, tasksettings[i], true);
    }
    catch (...) {
      throw std::logic_error("Encountered the following error in task " + name + ":\n" +
                             boost::current_exception_diagnostic_information());
    }
  }

  // Run tasks
  for (size_t i = 0; i < tasks.size(); i++) {
    auto name = tasks[i]->name();
    cout.printf("  #============#    #==%s==#\n", std::string(name.size(), '=').c_str());
    cout.printf("  |  Task %3d  |    |  %s  |\n", static_cast<int>(i + 1), name.c_str());
    cout.printf("  #============#    #==%s==#\n", std::string(name.size(), '=').c_str());
    cout << Core::Log::endl << Core::Log::endl;
    const auto& inpt = tasks[i]->input();
    const auto& outpt = tasks[i]->output();

    cout << "  Input System(s):  " + tasks[i]->input()[0] << Core::Log::endl;
    for (unsigned int i = 1; i < inpt.size(); i++) {
      cout << "                    " + inpt[i] << Core::Log::endl;
    }
    cout << Core::Log::endl;
    if (!outpt.empty()) {
      cout << "  Output System(s): " + outpt[0] << Core::Log::endl;
      for (unsigned int i = 1; i < outpt.size(); i++) {
        cout << "                    " + outpt[i] << Core::Log::endl;
      }
      cout << Core::Log::endl;
    }
    tasks[i]->run(systems, tasksettings[i]);
  }

  cout << "  Completed all tasks, exiting!" << Core::Log::endl;
  return 0;
}
