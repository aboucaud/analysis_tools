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

__all__ = ("AstrometricCatalogMatchConfig", "AstrometricCatalogMatchTask")

import astropy.units as units
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
from lsst.pipe.tasks.loadReferenceCatalog import LoadReferenceCatalogTask
from lsst.meas.algorithms import ReferenceObjectLoader
from lsst.pipe.tasks.configurableActions import ConfigurableActionStructField
from lsst.skymap import BaseSkyMap

from ..actions.vector import (
    CoaddPlotFlagSelector,
    GalaxySelector,
    SnSelector,
    StarSelector,
    VisitPlotFlagSelector,
)

from ..tasks.catalogMatch import CatalogMatchConfig, CatalogMatchTask, CatalogMatchConnections

class AstrometricCatalogMatchConfig(CatalogMatchConfig, pipelineConnections=CatalogMatchConnections):

    bands = pexConfig.ListField[str](
        doc="The bands to persist downstream",
        default=["g", "r", "i", "z", "y"],
    )

    def setDefaults(self):
        super().setDefaults()
        self.referenceCatalogLoader.doApplyColorTerms = False
        self.referenceCatalogLoader.refObjLoader.requireProperMotion = False
        self.referenceCatalogLoader.refObjLoader.anyFilterMapsToThis = "phot_g_mean"

class AstrometricCatalogMatchTask(CatalogMatchTask):
    """Match a tract-level catalog to a reference catalog"""

    ConfigClass = AstrometricCatalogMatchConfig
    _DefaultName = "analysisToolsAstrometricCatalogMatch"


    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        # Docs inherited from base class

        inputs = butlerQC.get(inputRefs)
        columns = self.prepColumns(self.config.bands)
        table = inputs["catalog"].get(parameters={"columns": columns})
        inputs["catalog"] = table

        tract = butlerQC.quantum.dataId["tract"]

        loaderTask = LoadReferenceCatalogTask(
            config=self.config.referenceCatalogLoader,
            dataIds=[ref.dataId for ref in inputRefs.refCat],
            name=inputs["refCat"][0].ref.datasetType.name,
            refCats=inputs["refCat"],
        )

        skymap = inputs.pop("skymap")
        loadedRefCat = self._loadRefCat(loaderTask, skymap[tract])
        outputs = self.run(catalog=inputs["catalog"], loadedRefCat=loadedRefCat, bands=self.config.bands)

        butlerQC.put(outputs, outputRefs)


class CatalogMatchVisitConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("visit",),
    defaultTemplates={"targetCatalog": "sourceTable_visit", "refCatalog": "gaia_dr2_20200414"},
):

    catalog = pipeBase.connectionTypes.Input(
        doc="The visit-wide catalog to make plots from.",
        storageClass="DataFrame",
        name="sourceTable_visit",
        dimensions=("visit",),
        deferLoad=True,
    )

    refCat = pipeBase.connectionTypes.PrerequisiteInput(
        doc="The astrometry reference catalog to match to loaded input catalog sources.",
        name="gaia_dr2_20200414",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True,
    )

    visitSummaryTable = pipeBase.connectionTypes.Input(
        doc="A summary table of the ccds in the visit",
        storageClass="ExposureCatalog",
        name="finalVisitSummary",
        dimensions=("visit",),
    )

    matchedCatalog = pipeBase.connectionTypes.Output(
        doc="Catalog with matched target and reference objects with separations",
        name="{targetCatalog}_{refCatalog}_match",
        storageClass="DataFrame",
        dimensions=("visit",),
    )


class CatalogMatchVisitConfig(CatalogMatchConfig, pipelineConnections=CatalogMatchVisitConnections):
    selectorActions = ConfigurableActionStructField(
        doc="Which selectors to use to narrow down the data for QA plotting.",
        default={"flagSelector": VisitPlotFlagSelector},
    )

    extraColumns = pexConfig.ListField[str](
        doc="Other catalog columns to persist to downstream tasks",
        default=["psfFlux", "psfFluxErr"],
    )

    def setDefaults(self):
        # sourceSelectorActions.sourceSelector is StarSelector
        self.sourceSelectorActions.sourceSelector.vectorKey = "extendedness"
        # extraColumnSelectors.selector1 is SnSelector
        self.extraColumnSelectors.selector1.fluxType = "psfFlux"
        # extraColumnSelectors.selector2 is GalaxySelector
        self.extraColumnSelectors.selector2.vectorKey = "extendedness"
        self.referenceCatalogLoader.doApplyColorTerms = False
        self.referenceCatalogLoader.refObjLoader.requireProperMotion = False


class CatalogMatchVisitTask(CatalogMatchTask):
    """Match a visit-level catalog to a reference catalog"""

    ConfigClass = CatalogMatchVisitConfig
    _DefaultName = "analysisToolsCatalogMatchVisit"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        # Docs inherited from base class

        inputs = butlerQC.get(inputRefs)

        columns = ["coord_ra", "coord_dec", "detector"] + self.config.extraColumns.list()
        for selectorAction in [
            self.config.selectorActions,
            self.config.sourceSelectorActions,
            self.config.extraColumnSelectors,
        ]:
            for selector in selectorAction:
                selectorSchema = selector.getFormattedInputSchema()
                columns += [s[0] for s in selectorSchema]

        dataFrame = inputs["catalog"].get(parameters={"columns": columns})
        inputs["catalog"] = dataFrame

        self.refObjLoader = ReferenceObjectLoader(
            dataIds=[ref.datasetRef.dataId for ref in inputRefs.refCat],
            refCats=inputs.pop("refCat"),
            name=self.config.connections.refCat,
            log=self.log,
        )

        self.setRefCat(inputs.pop("visitSummaryTable"))

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    def setRefCat(self, visitSummaryTable):
        """Make a reference catalog with coordinates in degrees

        Parameters
        ----------
        visitSummaryTable : `lsst.afw.table.ExposureCatalog`
            The table of visit information
        """
        # Get convex hull around the detectors, then get its center and radius
        corners = []
        for visSum in visitSummaryTable:
            for (ra, dec) in zip(visSum["raCorners"], visSum["decCorners"]):
                corners.append(lsst.geom.SpherePoint(ra, dec, units=lsst.geom.degrees).getVector())
        visitBoundingCircle = lsst.sphgeom.ConvexPolygon.convexHull(corners).getBoundingCircle()
        center = lsst.geom.SpherePoint(visitBoundingCircle.getCenter())
        radius = visitBoundingCircle.getOpeningAngle()

        # Get the observation date of the visit
        obsDate = visSum.getVisitInfo().getDate()
        epoch = Time(obsDate.toPython())

        # Load the reference catalog in the skyCircle of the detectors, then
        # convert the coordinates to degrees and convert the catalog to a
        # dataframe
        skyCircle = self.refObjLoader.loadSkyCircle(center, radius, "i", epoch=epoch)
        refCat = skyCircle.refCat

        refCat["coord_ra"] = (refCat["coord_ra"] * units.radian).to(units.degree).to_value()
        refCat["coord_dec"] = (refCat["coord_dec"] * units.radian).to(units.degree).to_value()
        self.refCat = refCat.asAstropy().to_pandas()
