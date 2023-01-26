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

__all__ = (
    "TargetRefCatDeltaRAScatterPlot",
    "TargetRefCatDeltaDecScatterPlot",
    "TargetRefCatDeltaRASkyPlot",
    "TargetRefCatDeltaDecSkyPlot",
    "TargetRefCatDeltaPhotom",
    "TargetRefCatDeltaSkyPlotPhotom",
)

from lsst.pex.config import Field

from ..actions.plot.barPlots import BarPanel, BarPlot
from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction, ScatterPlotWithTwoHists
from ..actions.plot.skyPlot import SkyPlot
from ..actions.vector import (
    AstromDiff,
    CoaddPlotFlagSelector,
    DownselectVector,
    ExtinctionCorrectedMagDiff,
    LoadVector,
    MagColumnNanoJansky,
    SnSelector,
    StarSelector,
    VectorSelector,
)
from ..analysisParts.baseFluxRatio import BasePsfApRatio
from ..analysisParts.baseSources import BaseSources
from ..analysisParts.genericPrep import CoaddPrep, VisitPrep
from ..interfaces import AnalysisPlot


class TargetRefCatDelta(AnalysisPlot):
    """Plot the difference in milliseconds between a target catalog and a
    reference catalog for the coordinate set in `setDefaults`.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )

    def coaddContext(self) -> None:
        self.prep = CoaddPrep()
        self.process.buildActions.starSelector.vectorKey = "{band}_extendedness"
        self.process.buildActions.mags = MagColumnNanoJansky(vectorKey="{band}_psfFlux")
        self.process.filterActions.psfFlux = DownselectVector(
            vectorKey="{band}_psfFlux", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.filterActions.psfFluxErr = DownselectVector(
            vectorKey="{band}_psfFluxErr", selector=VectorSelector(vectorKey="starSelector")
        )

    def visitContext(self) -> None:
        self.parameterizedBand = False
        self.prep = VisitPrep()
        self.process.buildActions.starSelector.vectorKey = "extendedness"
        self.process.buildActions.mags = MagColumnNanoJansky(vectorKey="psfFlux")
        self.process.filterActions.psfFlux = DownselectVector(
            vectorKey="psfFlux", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.filterActions.psfFluxErr = DownselectVector(
            vectorKey="psfFluxErr", selector=VectorSelector(vectorKey="starSelector")
        )

    def setDefaults(self, vectorKey):
        super().setDefaults()

        self.process.buildActions.starSelector = StarSelector()

        self.process.filterActions.xStars = DownselectVector(
            vectorKey="mags", selector=VectorSelector(vectorKey="starSelector")
        )

        self.process.calculateActions.stars = ScatterPlotStatsAction(vectorKey="yStars")
        self.process.calculateActions.stars.lowSNSelector.fluxType = "psfFlux"
        self.process.calculateActions.stars.highSNSelector.fluxType = "psfFlux"
        self.process.calculateActions.stars.fluxType = "psfFlux"

        self.produce = ScatterPlotWithTwoHists()

        self.produce.plotTypes = ["stars"]
        self.produce.xAxisLabel = "PSF Magnitude (mag)"
        self.produce.magLabel = "PSF Magnitude (mag)"


class TargetRefCatDeltaAstrom(TargetRefCatDelta):
    """Plot the difference in milliseconds between a target catalog and a
    reference catalog for the coordinate set in `setDefaults`.
    """

    def setDefaults(self, vectorKey):
        super().setDefaults(vectorKey)

        coordStr = vectorKey.lower()
        self.process.buildActions.astromDiff = AstromDiff(
            col1=f"coord_{coordStr}_target", col2=f"coord_{coordStr}_ref"
        )

        self.process.filterActions.yStars = DownselectVector(
            vectorKey="astromDiff", selector=VectorSelector(vectorKey="starSelector")
        )

        self.produce.yAxisLabel = f"${vectorKey}_{{target}} - {vectorKey}_{{ref}}$ (marcsec)"


class TargetRefCatDeltaPhotom(TargetRefCatDelta):
    """Plot the difference in milliseconds between a target catalog and a
    reference catalog for the coordinate set in `setDefaults`.
    """

    def setDefaults(self, vectorKey):
        super().setDefaults(vectorKey)

        self.process.buildActions.magDiff = ExtinctionCorrectedMagDiff()
        self.process.buildActions.magDiff.magDiff.col1 = "{band}" + f"_{vectorKey}_target"
        self.process.buildActions.magDiff.magDiff.col2 ="{band}_mag_ref"
        self.process.buildActions.magDiff.magDiff.fluxUnits2 = "ABmag"

        self.process.filterActions.yStars = DownselectVector(
            vectorKey="magDiff", selector=VectorSelector(vectorKey="starSelector")
        )

        self.produce.yAxisLabel = f"Output Mag - Ref Mag (mmag)"


class TargetRefCatDeltaPsfScatterPlot(TargetRefCatDeltaPhotom):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults(vectorKey="psfFlux")


class TargetRefCatDeltaRAScatterPlot(TargetRefCatDeltaAstrom):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults(vectorKey="RA")


class TargetRefCatDeltaDecScatterPlot(TargetRefCatDeltaAstrom):
    """Plot the difference in milliseconds between the Dec of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults(vectorKey="Dec")


class TargetRefCatDeltaSkyPlot(AnalysisPlot):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference between the RA or Dec of the target and reference catalog as
    the color.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )

    def coaddContext(self) -> None:
        self.prep = CoaddPrep()

        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux"

    def visitContext(self) -> None:
        self.parameterizedBand = False
        self.prep = VisitPrep()

        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "psfFlux"

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.xStars = LoadVector()
        self.process.buildActions.xStars.vectorKey = "coord_ra_target"
        self.process.buildActions.yStars = LoadVector()
        self.process.buildActions.yStars.vectorKey = "coord_dec_target"

        self.produce = SkyPlot()
        self.produce.plotTypes = ["stars"]
        self.produce.xAxisLabel = "R.A. (degrees)"
        self.produce.yAxisLabel = "Dec. (degrees)"
        self.produce.plotOutlines = False

class TargetRefCatDeltaSkyPlotAstrom(TargetRefCatDeltaSkyPlot):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference between the RA or Dec of the target and reference catalog as
    the color.
    """

    def setDefaults(self, vectorKey):
        super().setDefaults()

        coordStr = vectorKey.lower()
        self.process.buildActions.zStars = AstromDiff(
            col1=f"coord_{coordStr}_target", col2=f"coord_{coordStr}_ref"
        )

        self.produce.plotName = f"astromDiffSky_{vectorKey}"
        self.produce.zAxisLabel = f"${vectorKey}_{{target}} - {vectorKey}_{{ref}}$ (marcsec)"


class TargetRefCatDeltaSkyPlotPhotom(TargetRefCatDeltaSkyPlot):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference between the RA or Dec of the target and reference catalog as
    the color.
    """

    def setDefaults(self, vectorKey):
        super().setDefaults()

        self.process.buildActions.zStars = ExtinctionCorrectedMagDiff()
        self.process.buildActions.zStars.magDiff.col1 = "{band}" + f"_{vectorKey}_target"
        self.process.buildActions.zStars.magDiff.col2 = "{band}_mag_ref"
        self.process.buildActions.zStars.magDiff.fluxUnits2="ABmag"

        self.produce.plotName = "photomDiffSky_{band}"
        self.produce.zAxisLabel = "Output Mag - Ref Mag (mmag)"


class TargetRefCatDeltaPsfSkyPlot(TargetRefCatDeltaSkyPlotPhotom):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog as a function of RA and Dec.
    """

    def setDefaults(self):
        super().setDefaults(vectorKey="psfFlux")


class TargetRefCatDeltaRASkyPlot(TargetRefCatDeltaSkyPlotAstrom):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog as a function of RA and Dec.
    """

    def setDefaults(self):
        super().setDefaults(vectorKey="RA")


class TargetRefCatDeltaDecSkyPlot(TargetRefCatDeltaSkyPlotAstrom):
    """Plot the difference in milliseconds between the Dec of a target catalog
    and a reference catalog as a function of RA and Dec.
    """

    def setDefaults(self):
        super().setDefaults(vectorKey="Dec")
