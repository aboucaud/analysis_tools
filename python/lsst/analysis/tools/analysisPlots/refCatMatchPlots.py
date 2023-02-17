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
    "TargetRefCatDeltaScatterPhotom",
    "TargetRefCatDeltaScatterAstrom",
    "TargetRefCatDeltaSkyPlotPhotom",
    "TargetRefCatDeltaPsfScatterPlot",
    "TargetRefCatDeltaCModelScatterPlot",
    "TargetRefCatDeltaPsfSkyPlot",
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
from ..contexts import CoaddContext, RefMatchContext


class TargetRefCatDelta(AnalysisPlot):
    """Plot the difference in milliseconds between a target catalog and a
    reference catalog for the coordinate set in `setDefaults`.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )

    def coaddContext(self) -> None:
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux_target"
        self.prep.selectors.snSelector.threshold = 200
        self.prep.selectors.starSelector.vectorKey = "{band}_extendedness_target"
        self.process.buildActions.xStars = MagColumnNanoJansky()
        self.process.buildActions.xStars.vectorKey = "{band}_psfFlux_target"
        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux_target"

        self.process.calculateActions.stars = ScatterPlotStatsAction(vectorKey="yStars")
        self.process.calculateActions.stars.lowSNSelector.fluxType = "{band}_psfFlux_target"
        self.process.calculateActions.stars.highSNSelector.fluxType = "{band}_psfFlux_target"
        self.process.calculateActions.stars.fluxType = "{band}_psfFlux_target"

    def visitContext(self) -> None:
        self.parameterizedBand = False
        self.prep = VisitPrep()
        self.process.buildActions.starSelector.vectorKey = "extendedness_target"
        self.process.buildActions.mags = MagColumnNanoJansky(vectorKey="psfFlux_target")
        self.process.filterActions.psfFlux = DownselectVector(
            vectorKey="psfFlux_target", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.filterActions.psfFluxErr = DownselectVector(
            vectorKey="psfFluxErr_target", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "psfFlux_target"


    def setDefaults(self, vectorKey):
        super().setDefaults()

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.starSelector = StarSelector()


class TargetRefCatDeltaScatterAstrom(TargetRefCatDelta):
    """Plot the difference in milliseconds between a target catalog and a
    reference catalog for the coordinate set in `setDefaults`.
    """

    def setDefaults(self, vectorKey):
        super().setDefaults(vectorKey)

        coordStr = vectorKey.lower()
        self.process.buildActions.astromDiff = AstromDiff(
            col1=f"coord_{coordStr}_target", col2=f"coord_{coordStr}_ref"
        )

        #self.process.buildActions.yStars = LoadVector()
        #self.process.buildActions.yStars.vectorKey = "coord_ra_target"

        self.produce = ScatterPlotWithTwoHists()
        self.produce.plotTypes = ["stars"]
        self.produce.magLabel = "PSF Magnitude (mag)"
        self.produce.yAxisLabel = f"${vectorKey}_{{target}} - {vectorKey}_{{ref}}$ (marcsec)"
        self.produce.xAxisLabel = "PSF Magnitude (mag)"
        self.applyContext(CoaddContext)
        self.applyContext(RefMatchContext)


class TargetRefCatDeltaScatterPhotom(TargetRefCatDelta):
    """Plot the difference in milliseconds between a target catalog and a
    reference catalog for the coordinate set in `setDefaults`.
    """

    def setDefaults(self, vectorKey):
        super().setDefaults(vectorKey)

        self.process.buildActions.yStars = ExtinctionCorrectedMagDiff()
        self.process.buildActions.yStars.magDiff.col1 = "{band}" + f"_{vectorKey}"
        self.process.buildActions.yStars.magDiff.col2 ="{band}_mag_ref"
        self.process.buildActions.yStars.magDiff.fluxUnits2 = "mag(AB)"
        self.process.buildActions.yStars.ebvCol = "ebv_target"

        self.produce = ScatterPlotWithTwoHists()
        self.produce.plotTypes = ["stars"]
        self.produce.magLabel = "PSF Magnitude (mag)"
        self.produce.xAxisLabel = "PSF Magnitude (mag)"
        self.produce.yAxisLabel = f"Output Mag - Ref Mag (mmag)"
        self.applyContext(CoaddContext)
        self.applyContext(RefMatchContext)


class TargetRefCatDeltaPsfScatterPlot(TargetRefCatDeltaScatterPhotom):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults(vectorKey="psfFlux_target")

class TargetRefCatDeltaCModelScatterPlot(TargetRefCatDeltaScatterPhotom):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults(vectorKey="cModelFlux_target")


class TargetRefCatDeltaRAScatterPlot(TargetRefCatDeltaScatterAstrom):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults(vectorKey="RA")

class TargetRefCatDeltaDecScatterPlot(TargetRefCatDeltaScatterAstrom):
    """Plot the difference in milliseconds between the Dec of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults(vectorKey="Dec")


class TargetRefCatDeltaSkyPlot(TargetRefCatDelta):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference between the RA or Dec of the target and reference catalog as
    the color.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )


    def setDefaults(self, vectorKey):
        super().setDefaults(vectorKey)

        self.process.buildActions.xStars = LoadVector()
        self.process.buildActions.xStars.vectorKey = "coord_ra_target"
        self.process.buildActions.yStars = LoadVector()
        self.process.buildActions.yStars.vectorKey = "coord_dec_target"

        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux_target"

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
        super().setDefaults(vectorKey)

        coordStr = vectorKey.lower()
        self.process.buildActions.zStars = AstromDiff(
            col1=f"coord_{coordStr}_target", col2=f"coord_{coordStr}_ref"
        )

        self.produce.plotName = f"astromDiffSky_{vectorKey}"
        self.produce.zAxisLabel = f"${vectorKey}_{{target}} - {vectorKey}_{{ref}}$ (marcsec)"
        self.applyContext(CoaddContext)
        self.applyContext(RefMatchContext)


class TargetRefCatDeltaSkyPlotPhotom(TargetRefCatDeltaSkyPlot):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference between the RA or Dec of the target and reference catalog as
    the color.
    """

    def setDefaults(self, vectorKey):
        super().setDefaults(vectorKey)

        self.applyContext(RefMatchContext)

        self.process.buildActions.zStars = ExtinctionCorrectedMagDiff()
        self.process.buildActions.zStars.magDiff.col1 = "{band}" + f"_{vectorKey}"
        self.process.buildActions.zStars.magDiff.col2 = "{band}_mag_ref"
        self.process.buildActions.zStars.magDiff.fluxUnits2= "mag(AB)"
        self.process.buildActions.zStars.ebvCol = "ebv_target"

        self.produce.plotName = "photomDiffSky_{band}"
        self.produce.zAxisLabel = "Output Mag - Ref Mag (mmag)"
        self.applyContext(CoaddContext)
        self.applyContext(RefMatchContext)


class TargetRefCatDeltaPsfSkyPlot(TargetRefCatDeltaSkyPlotPhotom):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog as a function of RA and Dec.
    """

    def setDefaults(self):
        super().setDefaults(vectorKey="psfFlux_target")

class TargetRefCatDeltaCModelSkyPlot(TargetRefCatDeltaSkyPlotPhotom):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog as a function of RA and Dec.
    """

    def setDefaults(self):
        super().setDefaults(vectorKey="cModelFlux_target")


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
