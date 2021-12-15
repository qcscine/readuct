Changelog
=========

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

