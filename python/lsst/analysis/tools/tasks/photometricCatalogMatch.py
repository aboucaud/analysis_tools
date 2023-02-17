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

__all__ = ("PhotometricCatalogMatchConfig", "PhotometricCatalogMatchTask")


import numpy as np
from astropy.time import Time
from astropy.table import Table, hstack
from smatch import Matcher

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.tasks.loadReferenceCatalog import LoadReferenceCatalogTask
from lsst.skymap import BaseSkyMap
from lsst.pipe.tasks.configurableActions import ConfigurableActionStructField
import lsst.geom
from ..actions.vector import (
    CoaddPlotFlagSelector,
    GalaxySelector,
    SnSelector,
    StarSelector,
)


class PhotometricCatalogMatchConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("tract", "skymap"),
    defaultTemplates={"targetCatalog": "objectTable_tract", "refCatalog": "ps1_pv3_3pi_20170110"},
):
    catalog = pipeBase.connectionTypes.Input(
        doc="The tract-wide catalog to make plots from.",
        storageClass="ArrowAstropy",
        name="{targetCatalog}",
        dimensions=("tract", "skymap"),
        deferLoad=True,
    )

    refCat = pipeBase.connectionTypes.PrerequisiteInput(
        doc="The reference catalog to match to loaded input catalog sources.",
        name="ps1_pv3_3pi_20170110",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True,
    )

    skymap = pipeBase.connectionTypes.Input(
        doc="The skymap for the tract",
        storageClass="SkyMap",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        dimensions=("skymap",),
    )

    matchedCatalog = pipeBase.connectionTypes.Output(
        doc="Catalog with matched target and reference objects with separations",
        name="{targetCatalog}_{refCatalog}_match",
        storageClass="ArrowAstropy",
        dimensions=("tract", "skymap"),
    )


class PhotometricCatalogMatchConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=PhotometricCatalogMatchConnections
):

    referenceCatalogLoader = pexConfig.ConfigurableField(
        target=LoadReferenceCatalogTask,
        doc="Reference catalog loader",
    )

    epoch = pexConfig.Field[float](
        doc="Epoch to which reference objects are shifted.",
        default=2015.0,
    )

    filterNames = pexConfig.ListField[str](
        doc="Physical filter names to persist downstream.",
        default=["HSC-G", "HSC-R", "HSC-I", "HSC-Z", "HSC-Y"],
    )

    selectorBand = pexConfig.Field[str](
        doc="Band to use when selecting objects, primarily for extendedness.",
        default="i",
    )

    selectorActions = ConfigurableActionStructField(
        doc="Which selectors to use to narrow down the data for QA plotting.",
        default={"flagSelector": CoaddPlotFlagSelector},
    )

    sourceSelectorActions = ConfigurableActionStructField(
        doc="What types of sources to use.",
        default={"sourceSelector": StarSelector},
    )

    extraColumnSelectors = ConfigurableActionStructField(
        doc="Other selectors that are not used in this task, but whose columns; may be needed downstream.",
        default={"selector1": SnSelector, "selector2": GalaxySelector},
    )

    extraColumns = pexConfig.ListField[str](
        doc="Other catalog columns to persist to downstream tasks.",
        default=["x", "y", "ebv"],
    )

    extraPerBandColumns = pexConfig.ListField[str](
        doc="Other columns to load that should be loaded for each band individually.",
        default=["cModelFlux"],
    )

    raColumn = pexConfig.Field[str](doc="RA column.", default="coord_ra")
    decColumn = pexConfig.Field[str](doc="Dec column.", default="coord_dec")
    patchColumn = pexConfig.Field[str](doc="Patch column.", default="patch")

    # in the validation, we need to confirm that all of the filter names are in the filter map.

    # do a setDefaults and set however you like.
    def setDefaults(self):
        super().setDefaults()
        self.referenceCatalogLoader.doReferenceSelection = False


class PhotometricCatalogMatchTask(pipeBase.PipelineTask):
    ConfigClass = PhotometricCatalogMatchConfig
    _DefaultName = "analysisToolsPhotometricCatalogMatch"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Make any subtasks?  matching?

    def runQuantum(self, butlerQC, inputRefs, outputRefs):

        inputs = butlerQC.get(inputRefs)

        bands = []
        for filterName in self.config.filterNames:
            bands.append(self.config.referenceCatalogLoader.refObjLoader.filterMap[filterName])

        bandColumns = []
        for band in bands:
            for col in self.config.extraPerBandColumns:
                bandColumns.append(band + "_" + col)

        columns = [
            self.config.raColumn,
            self.config.decColumn,
            self.config.patchColumn,
        ] + self.config.extraColumns.list() + bandColumns

        for selectorAction in [
            self.config.selectorActions,
            self.config.sourceSelectorActions,
            self.config.extraColumnSelectors,
        ]:
            for selector in selectorAction:
                for band in bands:
                    selectorSchema = selector.getFormattedInputSchema(band=band)
                    columns += [s[0] for s in selectorSchema]

        table = inputs["catalog"].get(parameters={"columns": columns})
        inputs["catalog"] = table

        tract = butlerQC.quantum.dataId["tract"]

        loaderTask = LoadReferenceCatalogTask(
            config=self.config.referenceCatalogLoader,
            dataIds=[ref.dataId for ref in inputRefs.refCat],
            refCats=inputs.pop("refCat"),
            name=self.config.connections.refCat,
        )

        skymap = inputs.pop("skymap")
        loadedRefCat = self._loadRefCat(loaderTask, skymap[tract])

        outputs = self.run(catalog=inputs["catalog"], loadedRefCat=loadedRefCat, bands=bands)

        butlerQC.put(outputs, outputRefs)

    def _loadRefCat(self, loaderTask, tractInfo):
        """docstring"""
        boundingCircle = tractInfo.getOuterSkyPolygon().getBoundingCircle()
        center = lsst.geom.SpherePoint(boundingCircle.getCenter())
        radius = boundingCircle.getOpeningAngle()

        epoch = Time(self.config.epoch, format="decimalyear")

        # This is always going to return degrees.
        loadedRefCat = loaderTask.getSkyCircleCatalog(center, radius, self.config.filterNames, epoch=epoch)

        return Table(loadedRefCat)

    def run(self, *, catalog, loadedRefCat, bands):
        """docstring"""
        # Apply the selectors to the catalog
        mask = np.ones(len(catalog), dtype=bool)
        for selector in self.config.selectorActions:
            mask &= selector(catalog, bands=bands)

        for selector in self.config.sourceSelectorActions:
            mask &= selector(catalog, band=self.config.selectorBand).astype(bool)

        targetCatalog = catalog[mask]

        if (len(targetCatalog) == 0) or (len(loadedRefCat) == 0):
            refMatchIndices = np.array([], dtype=np.int64)
            targetMatchIndices = np.array([], dtype=np.int64)
            dists = np.array([], dtype=np.float64)
        else:
            # Run the matcher.

            # This all assumes that everything is in degrees.
            # Which I think is okay, but the current task allows different units.
            # Need to configure match radius, either in this task or a subtask.
            with Matcher(loadedRefCat["ra"], loadedRefCat["dec"]) as m:
                idx, refMatchIndices, targetMatchIndices, dists = m.query_radius(
                    targetCatalog[self.config.raColumn],
                    targetCatalog[self.config.decColumn],
                    1.0 / 3600.0,
                    return_indices=True,
                )

            # Convert degrees to arcseconds.
            dists *= 3600.0

        targetCatalogMatched = targetCatalog[targetMatchIndices]
        loadedRefCatMatched = loadedRefCat[refMatchIndices]

        targetCols = targetCatalogMatched.columns.copy()
        for col in targetCols:
            targetCatalogMatched.rename_column(col, col + "_target")
        refCols = loadedRefCatMatched.columns.copy()
        for col in refCols:
            loadedRefCatMatched.rename_column(col, col + "_ref")

        for (i, band) in enumerate(bands):
            loadedRefCatMatched[band + "_mag_ref"] = loadedRefCatMatched["refMag_ref"][:, i]
            loadedRefCatMatched[band + "_magErr_ref"] = loadedRefCatMatched["refMagErr_ref"][:, i]
        loadedRefCatMatched.remove_column("refMag_ref")
        loadedRefCatMatched.remove_column("refMagErr_ref")
        tMatched = hstack([targetCatalogMatched, loadedRefCatMatched], join_type="exact")
        tMatched["matchDistance"] = dists

        return pipeBase.Struct(matchedCatalog=tMatched)
