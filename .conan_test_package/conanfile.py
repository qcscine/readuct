__copyright__ = """This file is part of SCINE ReaDuct.
This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

from conans import ConanFile


class TestPackageConan(ConanFile):
    def build(self):
        pass

    def test(self):
        if self.options["scine_readuct"].python:
            self.output.info("Trying to import 'scine_readuct'")
            import scine_readuct
            self.output.info("Import worked")
