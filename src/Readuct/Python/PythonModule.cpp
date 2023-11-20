/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* Includes */
#include <pybind11/pybind11.h>

void init_tasks(pybind11::module& m);
void init_io(pybind11::module& m);

PYBIND11_MODULE(scine_readuct, m) {
  m.doc() = "Pybind11 Bindings for SCINE ReaDuct";

  // Requires other modules to function properly
  auto utils = pybind11::module::import("scine_utilities");

  init_tasks(m);
  init_io(m);
}
