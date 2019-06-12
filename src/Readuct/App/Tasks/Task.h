/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_TASK_H_
#define READUCT_TASK_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace YAML {
class Node;
}
namespace Scine {
namespace Core {
class Calculator;
}
namespace Readuct {

/**
 * @brief The base class for all tasks in Readuct.
 */
class Task {
 public:
  Task(std::vector<std::string> input, std::vector<std::string> output) : _input(input), _output(output) {
  }

  virtual std::string name() const = 0;
  virtual void run(std::map<std::string, std::shared_ptr<Core::Calculator>>& systems, const YAML::Node& taskSettings) const = 0;

  const std::vector<std::string>& input() const {
    return _input;
  };
  const std::vector<std::string>& output() const {
    return _output;
  };

 protected:
  std::vector<std::string> _input;
  std::vector<std::string> _output;
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_TASK_H_
