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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
from __future__ import annotations

__all__ = ("SkySourceSkyPlot", "SkySourceHistPlot")

from ..actions.plot.histPlot import HistPanel, HistPlot
from ..actions.plot.skyPlot import SkyPlot
from ..actions.vector import CalcSn, LoadVector
from ..actions.vector.selectors import SkySourceSelector, SnSelector
from ..interfaces import AnalysisTool


class SkySourceSkyPlot(AnalysisTool):
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.skySourceSelector = SkySourceSelector()

        # TODO: Can we make these defaults somewhere?
        self.process.buildActions.x = LoadVector()
        self.process.buildActions.x.vectorKey = "coord_ra"
        self.process.buildActions.y = LoadVector()
        self.process.buildActions.y.vectorKey = "coord_dec"
        self.process.buildActions.z = LoadVector()
        self.process.buildActions.z.vectorKey = "ap09Flux"
        self.process.buildActions.statMask = SnSelector()
        self.process.buildActions.statMask.threshold = -1e12
        self.process.buildActions.statMask.fluxType = "psfFlux"

        self.produce.plot = SkyPlot()
        self.produce.plot.plotTypes = ["any"]
        self.produce.plot.plotName = "skySource"
        self.produce.plot.xAxisLabel = "R.A. (degrees)"
        self.produce.plot.yAxisLabel = "Dec. (degrees)"
        self.produce.plot.zAxisLabel = "Sky Source ap09Flux (nJy)"
        self.produce.plot.plotOutlines = False
        self.produce.plot.fixAroundZero = True


class SkySourceHistPlot(AnalysisTool):
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.skySourceSelector = SkySourceSelector()

        self.process.buildActions.hist_psf_flux = LoadVector(vectorKey="psfFlux")
        self.process.buildActions.hist_ap09_flux = LoadVector(vectorKey="ap09Flux")
        self.process.buildActions.hist_psf_sn = CalcSn(fluxType="psfFlux")
        self.process.buildActions.hist_ap09_sn = CalcSn(fluxType="ap09Flux")

        self.produce.plot = HistPlot()

        self.produce.plot.panels["panel_flux"] = HistPanel()
        self.produce.plot.panels["panel_flux"].label = "Flux (nJy)"
        self.produce.plot.panels["panel_flux"].hists = dict(
            hist_psf_flux="psfFlux",
            hist_ap09_flux="ap09Flux",
        )
        self.produce.plot.panels["panel_flux"].rangeType = "sigmaMad"
        self.produce.plot.panels["panel_flux"].lowerRange = 3.5
        self.produce.plot.panels["panel_flux"].upperRange = 3.5
        self.produce.plot.panels["panel_flux"].referenceValue = 0.0
        self.produce.plot.panels["panel_flux"].validate()

        self.produce.plot.panels["panel_sn"] = HistPanel()
        self.produce.plot.panels["panel_sn"].label = "S/N"
        self.produce.plot.panels["panel_sn"].hists = dict(
            hist_psf_sn="psf S/N",
            hist_ap09_sn="ap09 S/N",
        )
        self.produce.plot.panels["panel_sn"].rangeType = "sigmaMad"
        self.produce.plot.panels["panel_sn"].lowerRange = 3.5
        self.produce.plot.panels["panel_sn"].upperRange = 3.5
        self.produce.plot.panels["panel_sn"].referenceValue = 0.0
        self.produce.plot.panels["panel_sn"].histDensity = True
        self.produce.plot.panels["panel_sn"].validate()
