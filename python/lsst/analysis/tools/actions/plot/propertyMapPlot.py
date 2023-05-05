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

__all__ = ("PropertyMapPlot",)

from typing import Mapping, Optional

import matplotlib.patheffects as pathEffects
import matplotlib.pyplot as plt
import numpy as np
import skyproj
from lsst.pex.config import Field, ListField
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from ...interfaces import KeyedData, PlotAction, Scalar, Vector
from ...statistics import nansigmaMad
from .plotUtils import (addPlotInfo, mkColormap, plotProjectionWithBinning,
                        sortAllArrays)


class PropertyMapPlot(PlotAction):

    plotName = Field[str](doc="The name for the plot.", optional=False)

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        # self._validateInput(data, **kwargs)
        return self.makePlot(data, **kwargs)

    def makePlot(
        self,
        data: KeyedData,
        tract: int,
        **kwargs,
    ) -> Figure:
        """Make the survey property map plot.

        Parameters
        ----------
        data : `KeyedData`
            The HealSparseMap to plot the points from.

        tract: `int`
            The tract that the data comes from

        **kwargs :
            Additional keyword arguments to pass to the plot

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The resulting figure.
        """

        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(111)

        sp = skyproj.GnomonicSkyproj(ax=ax, lon_0=tract.ctr_coord.getRa().asDegrees(),
                                     lat_0=tract.ctr_coord.getDec().asDegrees())

        _ = sp.draw_hspmap(data, zoom=True)

        sp.draw_colorbar(label="Exposure Time (s)")  # TODO: make it for multiple bands

        return fig
