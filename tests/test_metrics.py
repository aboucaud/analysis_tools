# This file is part of analysis_tools.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import annotations

from unittest import TestCase, main

import lsst.utils.tests
import numpy as np
from lsst.analysis.tools.atools import (
    MagnitudeTool,
    MatchedRefCoaddDiffMagMetric,
    MatchedRefCoaddDiffMetric,
    MatchedRefCoaddDiffPositionMetric,
)


class TestDiffMatched(TestCase):
    def setUp(self) -> None:
        super().setUp()
        self.band_default = "analysisTools"

    def _testMatchedRefCoaddMetricDerived(self, type_metric: type[MatchedRefCoaddDiffMetric], **kwargs):
        for compute_chi in (False, True):
            tester = type_metric(**kwargs, compute_chi=compute_chi)
            # tester.getInputSchema won't work properly before finalizing
            tester.finalize()

            keys = set(k[0] for k in tester.getInputSchema())
            self.assertGreater(len(keys), 0)
            self.assertGreater(len(list(tester.configureMetrics())), 0)
            data = {key.format(band=self.band_default): np.arange(5) for key in keys}
            output = tester(data)
            self.assertGreater(len(output), 0)

    def testMatchedRefCoaddMetric(self):
        tester = MatchedRefCoaddDiffMetric(unit="", name_prefix="")
        tester.finalize()
        self.assertGreater(len(list(tester.getInputSchema())), 0)
        self.assertGreater(len(list(tester.configureMetrics())), 0)

    def testMatchedRefCoaddDiffMagMetric(self):
        self._testMatchedRefCoaddMetricDerived(
            MatchedRefCoaddDiffMagMetric,
            fluxes={"cmodel": MagnitudeTool.fluxes_default.cmodel_err},
            mag_y="cmodel",
            name_prefix="",
            unit="",
        )

    def testMatchedRefCoaddDiffPositionMetric(self):
        for variable in ("x", "y"):
            self._testMatchedRefCoaddMetricDerived(MatchedRefCoaddDiffPositionMetric, variable=variable)


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    main()
