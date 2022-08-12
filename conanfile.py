from dev.conan.base import ScineConan


class ScineReaductConan(ScineConan):
    name = "scine_readuct"
    version = "4.1.0"
    url = "https://github.com/qcscine/readuct"
    description = """
SCINE ReaDuct is a command-line tool that allows to carry out:
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
quantum chemical properties (such as nuclear gradients). Currently, SCINE Sparrow
XTB, CP2K, Gaussian, ORCA, Serenity, and Turbomole are supported as backend programs."""
    options = {
        "shared": [True, False],
        "python": [True, False],
        "tests": [True, False],
        "coverage": [True, False],
        "microarch": ["detect", "none"],
    }
    default_options = {
        "shared": True,
        "python": False,
        "tests": False,
        "coverage": False,
        "microarch": "none"
    }
    exports = "dev/conan/*.py"
    exports_sources = [
        "dev/cmake/*", "src/*", "CMakeLists.txt", "README.rst",
        "LICENSE.txt", "dev/conan/hook.cmake", "dev/conan/glue/*"
    ]
    requires = ["scine_utilities/6.0.0",
                "boost/[>1.65.0]",
                "yaml-cpp/0.6.3"]
    cmake_name = "Readuct"
    cmake_definitions = {
        "BUILD_SPARROW": lambda self: self.options.tests
    }

    def configure(self):
        if self.options.tests:
            setattr(self.options["scine_sparrow"],
                    "python", self.options.python)

        super().configure()

    def build_requirements(self):
        if self.options.tests:
            self.build_requires("scine_sparrow/3.1.0")

        super().build_requirements()
