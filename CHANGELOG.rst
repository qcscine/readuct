Changelog
=========

Release 6.0.0
-------------

- Add test for QM/MM transition state optimization
- Improve automatic mode selection by introducing a mode score considering a weighted sum of
  the contributions of the relevant atoms and the wavenumber of the mode.
- Add the option to export the selected mode of a transition state optimization.
- Enable thermochemistry calculations for single atoms
- The one-electron integrals may be written to file through a dedicated integral evaluation task
- Separated energy and gradient contributions in QM/MM can be requested with the keyword "require_partial_energies",
  and "require_partial_gradients" in the single-point task.

Release 5.1.0
-------------

- Enable QM/MM methods for all tasks
- Improve support for compilation on Windows (MSVC)
- Add Python binding for reading YAML input files
- Update address in license

Release 5.0.0
-------------

- Avoid spin propensity calculations for calculators without spin multiplicity setting
- Write the trajectory with the trajectory with the structures and energies optimized during the B-spline
  optimization to file (BSpline Task)

Release 4.1.0
-------------

- Allow for custom observers in optimization tasks
- Various improvements

Release 4.0.0
-------------

- Add 2nd Newton trajectory scan algorithm (NT2)
- Update automatic TS mode picking to respect frequencies to some degree
- Deprecate BondOrderTask and allow bond order calculation in SinglePointTask
- Add option for spin propensity check
- Add option to optimize periodic boundaries in geometry optimization

Release 3.0.0
-------------

- Add Newton trajectory scans, searching for transition state guesses
- Add option to automatically select mode for transition state search based on important atoms
- Improved sanity checks before calculations
- Removed 'allow_unconverged' option and replaced with 'stop_on_error' option, which is now available in all tasks
- Allow Conan distribution
- Allow PyPI distribution

Release 2.0.0
-------------

- Interface to Gaussian
- Add B-Spline interpolation task; this allows for a
  linear interpolation between two structures, and
  for the optimization of a reaction path as well as for
  the extraction of a guess transition state
- Add dimer algorithm for transition state optimization
- Print various thermochemical quantities after a
  Hessian matrix has been calculated
- Add BFGS optimization algorithm
- Add G-DIIS convergence acceleration
- Improve the Bofill algorithm for transition state optimization
- Various bugfixes and improvements

Release 1.0.0
-------------

Initial release with the following features:

- Support the following tasks:

  - Single point calculations
  - Hessian calculations
  - Structure optimizaitons
  - Transition state searches
  - Intrinsic reaction coordinate (IRC) calculations
  - Artifial force induced reaction (AFIR) calculations

- Support for task chaining
- Interfaces to SCINE Sparrow and ORCA
- Python bindings

