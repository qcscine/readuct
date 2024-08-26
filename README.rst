SCINE - ReaDuct
===============

Introduction
------------

SCINE ReaDuct is a command-line tool that allows you to carry out

- single point calculations,
- bond order calculations,
- Hessian calculations,
- structure optimizations,
- single-ended transition state searches,
- double-ended B-Spline transition state searches,
- intrinsic reaction coordinate (IRC) calculations,
- artificial force induced reaction (AFIR) calculations, and
- Newton trajectory scans searching for transition state guesses.

For these calculations, it relies on a backend program to provide the necessary
quantum chemical properties (such as nuclear gradients). Currently, SCINE Sparrow,
XTB, CP2K, Gaussian, ORCA, Serenity, and Turbomole are supported as backend programs.

License and Copyright Information
---------------------------------

ReaDuct is distributed under the BSD 3-clause "New" or "Revised" License.
For more license and copyright information, see the file ``LICENSE.txt`` in the
repository.

Installation and Usage
----------------------

For instructions on how to install and use ReaDuct as well as for a detailed
documentation of the entire functionality of ReaDuct, please consult the user
manual found in the ``manual`` directory in in the repository.
Alternatively the manual can also be found on the official GitHub website
and on the SCINE website.

How to Cite
-----------

When publishing results obtained with ReaDuct, please cite the corresponding
release as archived on `Zenodo <https://zenodo.org/record/3244107>`_ (DOI
10.5281/zenodo.3244107; please use the DOI of the respective release).

In addition, we kindly request you to cite the following article when using ReaDuct:

A. C. Vaucher, M. Reiher, "Minimum Energy Paths and Transition States by Curve
Optimization", *J. Chem. Theory Comput.*, **2018**, *16*, 3091.

Furthermore, when publishing results obtained with any SCINE module, please cite the following paper:

T. Weymuth, J. P. Unsleber, P. L. Türtscher, M. Steiner, J.-G. Sobez, C. H. Müller, M. Mörchen,
V. Klasovita, S. A. Grimmel, M. Eckhoff, K.-S. Csizi, F. Bosia, M. Bensberg, M. Reiher,
"SCINE—Software for chemical interaction networks", *J. Chem. Phys.*, **2024**, *160*, 222501
(DOI `10.1063/5.0206974 <https://doi.org/10.1063/5.0206974>`_).

Support and Contact
-------------------

In case you should encounter problems or bugs, please write a short message
to scine@phys.chem.ethz.ch.

Third-Party Libraries Used
--------------------------

SCINE ReaDuct makes use of the following third-party libraries:

- `Boost <https://www.boost.org/>`_
- `Cereal <https://uscilab.github.io/cereal/>`_
- `Eigen <http://eigen.tuxfamily.org>`_
- `IRC <https://github.com/rmeli/irc>`_
- `Google Test <https://github.com/google/googletest>`_
- `pybind11 <https://github.com/pybind/pybind11>`_
- `yaml-cpp <https://github.com/jbeder/yaml-cpp>`_
