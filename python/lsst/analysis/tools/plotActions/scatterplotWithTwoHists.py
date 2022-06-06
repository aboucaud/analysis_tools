from __future__ import annotations

import numpy as np
import scipy.stats as sps

from functools import partial
from itertools import chain
from typing import Mapping, Optional, NamedTuple, cast

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.path import Path
from matplotlib.collections import PolyCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import PatchCollection


from .plotUtils import mkColormap

from lsst.pex.config import Field
from lsst.pex.config.listField import ListField

from ..interfaces import PlotAction, KeyedDataSchema, KeyedData, Scalar, Vector

cmapPatch = plt.cm.coolwarm.copy()
cmapPatch.set_bad(color="none")

sigmaMad = partial(sps.median_abs_deviation, scale="normal")  # type: ignore


def _validatePlotTypes(value):
    return value in ("stars", "galaxies", "unknown", "any", "mag")


class _StatsContainer(NamedTuple):
    median: Scalar
    sigmaMad: Scalar
    count: Scalar
    aproxMag: Scalar


class ScatterPlotWithTwoHists(PlotAction):
    yLims = ListField(
        doc="ylimits of the plot, if not specified determined from data",
        dtype=float,
        length=2,
        optional=True,
    )

    xLims = ListField(
        doc="xlimits of the plot, if not specified determined from data", dtype=float, length=2, optional=True
    )
    xAxisLabel = Field(doc="Label to use for the x axis", dtype=str, optional=False)
    yAxisLabel = Field(doc="Label to use for the y axis", dtype=str, optional=False)
    magLabel = Field(doc="Label to use for the magnitudes used for SNR", dtype=str, optional=False)
    nBins = Field(doc="Number of bins on x axis", dtype=float, default=40.0)
    plot2DHist = Field(
        doc="Plot a 2D histogram in dense areas of points on the scatter plot."
        "Doesn't look great if plotting multiple datasets on top of each other.",
        default=True,
        dtype=bool,
    )
    plotTypes = ListField(
        doc="Selection of types of objects to plot. Can take any combination of"
        " stars, galaxies, unknown, mag, any",
        dtype=str,
        optional=False,
        itemCheck=_validatePlotTypes,
    )
    highThreshold = Field(doc="The Value used as high threshold in statistics", dtype=float, optional=False)
    lowThreshold = Field(doc="The Value used as low threshold in statistics", dtype=float, optional=False)

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.stats = ("median", "sigmaMad", "count", "approxMag")

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        base = []
        if "stars" in self.plotTypes:  # type: ignore
            base.append(("xsStars", Vector))
            base.append(("ysStars", Vector))
            base.append(("starHighSNMask", Vector))
            base.append(("starLowSNMask", Vector))
            # statistics
            for name in self.stats:
                base.append((f"{{band}}_highSNStars_{name}".format(**kwargs), Scalar))
                base.append((f"{{band}}_lowSNStars_{name}".format(**kwargs), Scalar))
        if "galaxies" in self.plotTypes:  # type: ignore
            base.append(("xsGalaxies", Vector))
            base.append(("ysGalaxies", Vector))
            base.append(("galaxyHighSNMask", Vector))
            base.append(("galaxyLowSNMask", Vector))
            # statistics
            for name in self.stats:
                base.append((f"{{band}}_highSNGalaxies_{name}".format(**kwargs), Scalar))
                base.append((f"{{band}}_lowSNGalaxies_{name}".format(**kwargs), Scalar))
        if "unknown" in self.plotTypes:  # type: ignore
            base.append(("xsUnknown", Scalar))
            base.append(("ysUnknown", Scalar))
            base.append(("unknownHighSNMask", Vector))
            base.append(("unknownLowSNMask", Vector))
            # statistics
            for name in self.stats:
                base.append(f"{{band}}_highSNUnknown_{name}".format(**kwargs))
                base.append(f"{{band}}_lowSNUnknown_{name}".format(**kwargs))
        if "any" in self.plotTypes:  # type: ignore
            base.extend((("x", Vector), ("y", Vector)))
            base.append(("anyHighSNMask", Vector))
            base.append(("anySNMask", Vector))
            # statistics
            for name in self.stats:
                base.append((f"{{band}}_highSNAny_{name}".format(**kwargs), Scalar))
                base.append((f"{{band}}_lowSNAny_{name}".format(**kwargs), Scalar))
        return base

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        self._validateInput(data)
        return self.makePlot(data, **kwargs)
        # table is a dict that needs: x, y, run, skymap, filter, tract,

    def _validateInput(self, data: KeyedData, **kwargs) -> None:
        needed = self.getInputSchema()
        if remainder := {key.format(**kwargs) for key, _ in needed} - {
            key.format(**kwargs) for key in data.keys()
        }:
            raise ValueError(f"Task needs keys {remainder} but they were not found in input")
        # Can't do this until Vector is converted to a proper protocol
        # for name, typ in needed:
        #     if not isinstance(colType := type(data[name.format(**kwargs)]), typ):
        #         raise ValueError(f"Data keyed by {name} has type {colType} but action requires type {typ}")

    def _scatterPlot(
        self, data: KeyedData, fig: Figure, gs: gridspec.GridSpec, **kwargs
    ) -> tuple[Axes, Optional[PolyCollection]]:
        # Main scatter plot
        ax = fig.add_subplot(gs[1:, :-1])

        newBlues = mkColormap(["paleturquoise", "midnightBlue"])
        newReds = mkColormap(["lemonchiffon", "firebrick"])

        binThresh = 5

        yBinsOut = []
        linesForLegend = []

        toPlotList = []
        histIm = None
        if "stars" in self.plotTypes:  # type: ignore
            highArgs = {}
            lowArgs = {}
            for name in self.stats:
                highArgs[name] = f"{{band}}_highSNStars_{name}".format(**kwargs)
                lowArgs[name] = f"{{band}}_lowSNStars_{name}".format(**kwargs)
            highStats = _StatsContainer(**highArgs)
            lowStats = _StatsContainer(**lowArgs)

            toPlotList.append(
                (
                    data["xsStars"],
                    data["ysStars"],
                    data["starHighSNMask"],
                    data["starLowSNMask"],
                    "midnightblue",
                    newBlues,
                    highStats,
                    lowStats,
                )
            )
        if "galaxies" in self.plotTypes:  # type: ignore
            highArgs = {}
            lowArgs = {}
            for name in self.stats:
                highArgs[name] = f"{{band}}_highSNGalaxies_{name}".format(**kwargs)
                lowArgs[name] = f"{{band}}_lowSNGalaxies_{name}".format(**kwargs)
            highStats = _StatsContainer(**highArgs)
            lowStats = _StatsContainer(**lowArgs)

            toPlotList.append(
                (
                    data["xsGalaxies"],
                    data["ysGalaxies"],
                    data["galaxyHighSNMask"],
                    data["galaxyLowSNMask"],
                    "firebrick",
                    newReds,
                    highStats,
                    lowStats,
                )
            )
        if "unknown" in self.plotTypes:  # type: ignore
            highArgs = {}
            lowArgs = {}
            for name in self.stats:
                highArgs[name] = f"{{band}}_highSNUnknown_{name}".format(**kwargs)
                lowArgs[name] = f"{{band}}_lowSNUnknown_{name}".format(**kwargs)
            highStats = _StatsContainer(**highArgs)
            lowStats = _StatsContainer(**lowArgs)

            toPlotList.append(
                (
                    data["xsUnknown"],
                    data["ysUnknown"],
                    data["unknownHighSNMask"],
                    data["unknownLowSNMask"],
                    "green",
                    None,
                    highStats,
                    lowStats,
                )
            )
        if "any" in self.plotTypes:  # type: ignore
            highArgs = {}
            lowArgs = {}
            for name in self.stats:
                highArgs[name] = f"{{band}}_highSNUnknown_{name}".format(**kwargs)
                lowArgs[name] = f"{{band}}_lowSNUnknown_{name}".format(**kwargs)
            highStats = _StatsContainer(**highArgs)
            lowStats = _StatsContainer(**lowArgs)

            toPlotList.append(
                (
                    data["x"],
                    data["y"],
                    data["anyHighSNMask"],
                    data["anyLowSNMask"],
                    "purple",
                    None,
                    highStats,
                    lowStats,
                )
            )

        highStats: _StatsContainer
        lowStats: _StatsContainer
        xMin = None
        for (j, (xs, ys, highSn, lowSn, color, cmap, highStats, lowStats)) in enumerate(toPlotList):
            sigMadYs = sigmaMad(ys, nan_policy="omit")
            if len(xs) < 2:
                (medLine,) = ax.plot(
                    xs, np.nanmedian(ys), color, label=f"Median: {np.nanmedian(ys):0.3g}", lw=0.8
                )
                linesForLegend.append(medLine)
                sigMads = np.array([sigmaMad(ys, nan_policy="omit")] * len(xs))
                (sigMadLine,) = ax.plot(
                    xs,
                    np.nanmedian(ys) + 1.0 * sigMads,
                    color,
                    alpha=0.8,
                    lw=0.8,
                    label=r"$\sigma_{MAD}$: " + f"{sigMads[0]:0.3g}",
                )
                ax.plot(xs, np.nanmedian(ys) - 1.0 * sigMads, color, alpha=0.8)
                linesForLegend.append(sigMadLine)
                histIm = None
                continue

            [xs1, xs25, xs50, xs75, xs95, xs97] = np.nanpercentile(xs, [1, 25, 50, 75, 95, 97])
            xScale = (xs97 - xs1) / 20.0  # This is ~5% of the data range

            # 40 was used as the default number of bins because it looked good
            xEdges = np.arange(
                np.nanmin(xs) - xScale,
                np.nanmax(xs) + xScale,
                (np.nanmax(xs) + xScale - (np.nanmin(xs) - xScale)) / self.nBins,
            )
            medYs = np.nanmedian(ys)
            fiveSigmaHigh = medYs + 5.0 * sigMadYs
            fiveSigmaLow = medYs - 5.0 * sigMadYs
            binSize = (fiveSigmaHigh - fiveSigmaLow) / 101.0
            yEdges = np.arange(fiveSigmaLow, fiveSigmaHigh, binSize)

            counts, xBins, yBins = np.histogram2d(xs, ys, bins=(xEdges, yEdges))
            yBinsOut.append(yBins)
            countsYs = np.sum(counts, axis=1)

            ids = np.where((countsYs > binThresh))[0]
            xEdgesPlot = xEdges[ids][1:]
            xEdges = xEdges[ids]

            if len(ids) > 1:
                # Create the codes needed to turn the sigmaMad lines
                # into a path to speed up checking which points are
                # inside the area.
                codes = np.ones(len(xEdgesPlot) * 2) * Path.LINETO
                codes[0] = Path.MOVETO
                codes[-1] = Path.CLOSEPOLY

                meds = np.zeros(len(xEdgesPlot))
                threeSigMadVerts = np.zeros((len(xEdgesPlot) * 2, 2))
                sigMads = np.zeros(len(xEdgesPlot))

                for (i, xEdge) in enumerate(xEdgesPlot):
                    ids = np.where((xs < xEdge) & (xs > xEdges[i]) & (np.isfinite(ys)))[0]
                    med = np.median(ys[ids])
                    sigMad = sigmaMad(ys[ids])
                    meds[i] = med
                    sigMads[i] = sigMad
                    threeSigMadVerts[i, :] = [xEdge, med + 3 * sigMad]
                    threeSigMadVerts[-(i + 1), :] = [xEdge, med - 3 * sigMad]

                (medLine,) = ax.plot(xEdgesPlot, meds, color, label="Running Median")
                linesForLegend.append(medLine)

                # Make path to check which points lie within one sigma mad
                threeSigMadPath = Path(threeSigMadVerts, codes)

                # Add lines for the median +/- 3 * sigma MAD
                (threeSigMadLine,) = ax.plot(
                    xEdgesPlot,
                    threeSigMadVerts[: len(xEdgesPlot), 1],
                    color,
                    alpha=0.4,
                    label=r"3$\sigma_{MAD}$",
                )
                ax.plot(xEdgesPlot[::-1], threeSigMadVerts[len(xEdgesPlot) :, 1], color, alpha=0.4)

                # Add lines for the median +/- 1 * sigma MAD
                (sigMadLine,) = ax.plot(
                    xEdgesPlot, meds + 1.0 * sigMads, color, alpha=0.8, label=r"$\sigma_{MAD}$"
                )
                linesForLegend.append(sigMadLine)
                ax.plot(xEdgesPlot, meds - 1.0 * sigMads, color, alpha=0.8)

                # Add lines for the median +/- 2 * sigma MAD
                (twoSigMadLine,) = ax.plot(
                    xEdgesPlot, meds + 2.0 * sigMads, color, alpha=0.6, label=r"2$\sigma_{MAD}$"
                )
                linesForLegend.append(twoSigMadLine)
                linesForLegend.append(threeSigMadLine)
                ax.plot(xEdgesPlot, meds - 2.0 * sigMads, color, alpha=0.6)

                # Check which points are outside 3 sigma MAD of the median
                # and plot these as points.
                inside = threeSigMadPath.contains_points(np.array([xs, ys]).T)
                ax.plot(xs[~inside], ys[~inside], ".", ms=3, alpha=0.3, mfc=color, mec=color, zorder=-1)

                # Add some stats text
                xPos = 0.65 - 0.4 * j
                bbox = dict(edgecolor=color, linestyle="--", facecolor="none")
                highThresh = self.highThreshold
                statText = f"S/N > {highThresh} Stats ({self.magLabel} < {highStats.aproxMag:0.3g})\n"
                highStatsStr = (
                    f"Median: {highStats.median:0.3g}    "
                    + r"$\sigma_{MAD}$: "
                    + f"{highStats.sigmaMad:0.3g}    "
                    + r"N$_{points}$: "
                    + f"{highStats.count}"
                )
                statText += highStatsStr
                fig.text(xPos, 0.090, statText, bbox=bbox, transform=fig.transFigure, fontsize=6)

                bbox = dict(edgecolor=color, linestyle=":", facecolor="none")
                lowThresh = self.lowThreshold
                statText = f"S/N > {lowThresh} Stats ({self.magLabel} < {lowStats.aproxMag:0.3g})\n"
                lowStatsStr = (
                    f"Median: {lowStats.median:0.3g}    "
                    + r"$\sigma_{MAD}$: "
                    + f"{lowStats.sigmaMad:0.3g}    "
                    + r"N$_{points}$: "
                    + f"{lowStats.count}"
                )
                statText += lowStatsStr
                fig.text(xPos, 0.020, statText, bbox=bbox, transform=fig.transFigure, fontsize=6)

                if self.plot2DHist:
                    histIm = ax.hexbin(xs[inside], ys[inside], gridsize=75, cmap=cmap, mincnt=1, zorder=-2)

                # If there are not many sources being used for the
                # statistics then plot them individually as just
                # plotting a line makes the statistics look wrong
                # as the magnitude estimation is iffy for low
                # numbers of sources.
                if np.sum(highSn) < 100 and np.sum(highSn) > 0:
                    ax.plot(xs[highSn], ys[highSn], marker="x", ms=4, mec="w", mew=2, ls="none")
                    (highSnLine,) = ax.plot(
                        xs[highSn], ys[highSn], color=color, marker="x", ms=4, ls="none", label="High SN"
                    )
                    linesForLegend.append(highSnLine)
                    xMin = np.min(xs[highSn])
                else:
                    ax.axvline(highStats.aproxMag, color=color, ls="--")

                if np.sum(lowSn) < 100 and np.sum(lowSn) > 0:
                    ax.plot(xs[lowSn], ys[lowSn], marker="+", ms=4, mec="w", mew=2, ls="none")
                    (lowSnLine,) = ax.plot(
                        xs[lowSn], ys[lowSn], color=color, marker="+", ms=4, ls="none", label="Low SN"
                    )
                    linesForLegend.append(lowSnLine)
                    if xMin is None or xMin > np.min(xs[lowSn]):
                        xMin = np.min(xs[lowSn])
                else:
                    ax.axvline(lowStats.aproxMag, color=color, ls=":")

            else:
                ax.plot(xs, ys, ".", ms=5, alpha=0.3, mfc=color, mec=color, zorder=-1)
                meds = np.array([np.nanmedian(ys)] * len(xs))
                (medLine,) = ax.plot(xs, meds, color, label=f"Median: {np.nanmedian(ys):0.3g}", lw=0.8)
                linesForLegend.append(medLine)
                sigMads = np.array([sigmaMad(ys, nan_policy="omit")] * len(xs))
                (sigMadLine,) = ax.plot(
                    xs,
                    meds + 1.0 * sigMads,
                    color,
                    alpha=0.8,
                    lw=0.8,
                    label=r"$\sigma_{MAD}$: " + f"{sigMads[0]:0.3g}",
                )
                ax.plot(xs, meds - 1.0 * sigMads, color, alpha=0.8)
                linesForLegend.append(sigMadLine)
                histIm = None

        # Set the scatter plot limits
        if len(cast(Vector, data["ysStars"])) > 0:
            plotMed = np.nanmedian(data["ysStars"])
        else:
            plotMed = np.nanmedian(data["ysGalaxies"])
        if len(xs) < 2:
            meds = [np.median(ys)]
        if self.yLims:
            ax.set_ylim(self.yLims[0], self.yLims[1])  # type: ignore
        else:
            numSig = 4
            yLimMin = plotMed - numSig * sigMadYs
            yLimMax = plotMed + numSig * sigMadYs
            while (yLimMax < np.max(meds) or yLimMin > np.min(meds)) and numSig < 10:
                numSig += 1

            numSig += 1
            yLimMin = plotMed - numSig * sigMadYs
            yLimMax = plotMed + numSig * sigMadYs
            ax.set_ylim(yLimMin, yLimMax)

        if self.xLims:
            ax.set_xlim(self.xLims[0], self.xLims[1])  # type: ignore
        elif len(xs) > 2:
            if xMin is None:
                xMin = xs1 - 2 * xScale
            ax.set_xlim(xMin, xs97 + 2 * xScale)

        # Add a line legend
        ax.legend(
            handles=linesForLegend,
            ncol=4,
            fontsize=6,
            loc="upper left",
            framealpha=0.9,
            edgecolor="k",
            borderpad=0.4,
            handlelength=1,
        )

        # Add axes labels
        ax.set_ylabel(self.yAxisLabel, fontsize=10, labelpad=10)
        ax.set_xlabel(self.xAxisLabel, fontsize=10, labelpad=2)

        return ax, histIm

    def makePlot(self, data: KeyedData, **kwargs):
        """Makes a generic plot with a 2D histogram and collapsed histograms of
        each axis.
        Parameters
        ----------
        catPlot : `pandas.core.frame.DataFrame`
            The catalog to plot the points from.
        plotInfo : `dict`
            A dictionary of information about the data being plotted with keys:
                ``"run"``
                    The output run for the plots (`str`).
                ``"skymap"``
                    The type of skymap used for the data (`str`).
                ``"filter"``
                    The filter used for this data (`str`).
                ``"tract"``
                    The tract that the data comes from (`str`).
        sumStats : `dict`
            A dictionary where the patchIds are the keys which store the R.A.
            and dec of the corners of the patch, along with a summary
            statistic for each patch.
        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The resulting figure.
        Notes
        -----
        Uses the axisLabels config options `x` and `y` and the axisAction
        config options `xAction` and `yAction` to plot a scatter
        plot of the values against each other. A histogram of the points
        collapsed onto each axis is also plotted. A summary panel showing the
        median of the y value in each patch is shown in the upper right corner
        of the resultant plot. The code uses the selectorActions to decide
        which points to plot and the statisticSelector actions to determine
        which points to use for the printed statistics.
        """
        if not self.plotTypes:
            noDataFig = plt.Figure()
            noDataFig.text(0.3, 0.5, "No data to plot after selectors applied")
            noDataFig = addPlotInfo(noDataFig, plotInfo)
            return noDataFig

        fig = plt.figure(dpi=300)
        gs = gridspec.GridSpec(4, 4)

        # add the various plot elements
        ax, imhist = self._scatterPlot(data, fig, gs, **kwargs)
        self._makeTopHistogram(data, fig, gs, ax)

    def _makeTopHistogram(
        self, data: KeyedData, figure: Figure, gs: gridspec.GridSpec, ax: Axes, **kwargs
    ) -> None:
        # Top histogram
        totalX = []
        if "stars" in self.plotTypes:  # type: ignore
            totalX.append(data["xsStars"])
        if "galaxies" in self.plotTypes:  # type: ignore
            totalX.append(data["xsGalaxies"])
        if "unknown" in self.plotTypes:  # type: ignore
            totalX.append(data["xsUknown"])
        if "any" in self.plotTypes:  # type: ignore
            totalX.append(data["x"])

        # cheat to get the total count while iterating once
        totalXChained = list(chain(totalX))

        topHist = figure.add_subplot(gs[0, :-1], sharex=ax)
        topHist.hist(
            totalXChained, bins=100, color="grey", alpha=0.3, log=True, label=f"All ({len(totalXChained)})"
        )
        if "galaxies" in self.plotTypes:  # type: ignore
            topHist.hist(
                data["xsGalaxies"],
                bins=100,
                color="firebrick",
                histtype="step",
                log=True,
                label=f"Galaxies ({len(cast(Vector, data['xsGalaxies']))})",
            )
        if "stars" in self.plotTypes:  # type: ignore
            topHist.hist(
                data["xsStars"],
                bins=100,
                color="midnightblue",
                histtype="step",
                log=True,
                label=f"Stars ({len(cast(Vector, data['xsStars']))})",
            )
        topHist.axes.get_xaxis().set_visible(False)
        topHist.set_ylabel("Number", fontsize=8)
        topHist.legend(fontsize=6, framealpha=0.9, borderpad=0.4, loc="lower left", ncol=3, edgecolor="k")

        # Side histogram

    def _makeSideHistogram(
        self,
        data: KeyedData,
        figure: Figure,
        gs: gridspec.Gridspec,
        ax: Axes,
        histIm: Optional[PolyCollection],
        **kwargs,
    ) -> None:
        sideHist = figure.add_subplot(gs[1:, -1], sharey=ax)

        totalY = []
        if "stars" in self.plotTypes:  # type: ignore
            totalY.append(data["ysStars"])
        if "galaxies" in self.plotTypes:  # type: ignore
            totalY.append(data["ysGalaxies"])
        if "unknown" in self.plotTypes:  # type: ignore
            totalY.append(data["ysUknown"])
        if "any" in self.plotTypes:  # type: ignore
            totalY.append(data["y"])
        totalYChained = [y for y in chain(totalY) if y == y]

        # cheat to get the total count while iterating once
        yLimMin, yLimMax = ax.get_ylim()
        bins = np.linspace(yLimMin, yLimMax)
        sideHist.hist(
            totalYChained,
            bins=bins,
            color="grey",
            alpha=0.3,
            orientation="horizontal",
            log=True,
        )
        if "galaxies" in self.plotTypes:  # type: ignore
            sideHist.hist(
                [g for g in cast(Vector, data["ysGalaxies"]) if g == g],
                bins=bins,
                color="firebrick",
                histtype="step",
                orientation="horizontal",
                log=True,
            )
            sideHist.hist(
                cast(Vector, data["ysGalaxies"])[data["galaxyHighSNMask"]],
                bins=bins,
                color="firebrick",
                histtype="step",
                orientation="horizontal",
                log=True,
                ls="--",
            )
            sideHist.hist(
                cast(Vector, data["ysGalaxies"])[data["galaxyLowSNMask"]],
                bins=bins,
                color="firebrick",
                histtype="step",
                orientation="horizontal",
                log=True,
                ls=":",
            )

        if "stars" in self.plotTypes:  # type: ignore
            sideHist.hist(
                [s for s in cast(Vector, data["ysStars"]) if s == s],
                bins=bins,
                color="midnightblue",
                histtype="step",
                orientation="horizontal",
                log=True,
            )
            sideHist.hist(
                cast(Vector, data["ysStars"])[data["starHighSNMask"]],
                bins=bins,
                color="midnightblue",
                histtype="step",
                orientation="horizontal",
                log=True,
                ls="--",
            )
            sideHist.hist(
                cast(Vector, data["ysStars"])[data["starLowSNMask"]],
                bins=bins,
                color="midnightblue",
                histtype="step",
                orientation="horizontal",
                log=True,
                ls=":",
            )

        sideHist.axes.get_yaxis().set_visible(False)
        sideHist.set_xlabel("Number", fontsize=8)
        if self.plot2DHist and histIm is not None:
            divider = make_axes_locatable(sideHist)
            cax = divider.append_axes("right", size="8%", pad=0)
            plt.colorbar(histIm, cax=cax, orientation="vertical", label="Number of Points Per Bin")

        # Corner plot of patches showing summary stat in each
        axCorner = plt.gcf().add_subplot(gs[0, -1])
        axCorner.yaxis.tick_right()
        axCorner.yaxis.set_label_position("right")
        axCorner.xaxis.tick_top()
        axCorner.xaxis.set_label_position("top")
        axCorner.set_aspect("equal")

        patches = []
        colors = []
        for dataId in sumStats.keys():
            (corners, stat) = sumStats[dataId]
            ra = corners[0][0].asDegrees()
            dec = corners[0][1].asDegrees()
            xy = (ra, dec)
            width = corners[2][0].asDegrees() - ra
            height = corners[2][1].asDegrees() - dec
            patches.append(Rectangle(xy, width, height))
            colors.append(stat)
            ras = [ra.asDegrees() for (ra, dec) in corners]
            decs = [dec.asDegrees() for (ra, dec) in corners]
            axCorner.plot(ras + [ras[0]], decs + [decs[0]], "k", lw=0.5)
            cenX = ra + width / 2
            cenY = dec + height / 2
            if dataId != "tract":
                axCorner.annotate(dataId, (cenX, cenY), color="k", fontsize=4, ha="center", va="center")

        # Set the bad color to transparent and make a masked array
        colors = np.ma.array(colors, mask=np.isnan(colors))
        collection = PatchCollection(patches, cmap=cmapPatch)
        collection.set_array(colors)
        axCorner.add_collection(collection)

        axCorner.set_xlabel("R.A. (deg)", fontsize=7)
        axCorner.set_ylabel("Dec. (deg)", fontsize=7)
        axCorner.tick_params(axis="both", labelsize=6, length=0, pad=1.5)
        axCorner.invert_xaxis()

        # Add a colorbar
        pos = axCorner.get_position()
        cax = fig.add_axes([pos.x0, pos.y0 + 0.23, pos.x1 - pos.x0, 0.025])
        plt.colorbar(collection, cax=cax, orientation="horizontal")
        cax.text(
            0.5,
            0.5,
            "Median Value",
            color="k",
            transform=cax.transAxes,
            rotation="horizontal",
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=6,
        )
        cax.tick_params(
            axis="x", labelsize=6, labeltop=True, labelbottom=False, bottom=False, top=True, pad=0.5, length=2
        )

        plt.draw()
        plt.subplots_adjust(wspace=0.0, hspace=0.0, bottom=0.22, left=0.21)
        fig = plt.gcf()
        fig = addPlotInfo(fig, plotInfo)

        return fig