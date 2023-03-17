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

from ..tasks.catalogMatch import CatalogMatchConfig, CatalogMatchTask, CatalogMatchConnections

class PhotometricCatalogMatchConnections(CatalogMatchConnections):
    pass

class PhotometricCatalogMatchConfig(
    CatalogMatchConfig, pipelineConnections=PhotometricCatalogMatchConnections
):

    filterNames = pexConfig.ListField[str](
        doc="Physical filter names to persist downstream.",
        default=["HSC-G", "HSC-R", "HSC-I", "HSC-Z", "HSC-Y"],
    )

    # in the validation, we need to confirm that all of the filter names are in the filter map.

    # do a setDefaults and set however you like.
    def setDefaults(self):
        super().setDefaults()
        self.referenceCatalogLoader.doReferenceSelection = False
        self.referenceCatalogLoader.doApplyColorTerms = True


class PhotometricCatalogMatchTask(CatalogMatchTask):
    """A wrapper task to provide the information that
    is specific to the photometric reference catalog.
    """
    ConfigClass = PhotometricCatalogMatchConfig
    _DefaultName = "analysisToolsPhotometricCatalogMatch"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        """Run the matching to the photometric reference
        catalog.

        Parameters
        ----------
        `butlerQC` : lsst.pipe.base.butlerQuantumContext.ButlerQuantumContext
        `inputRefs` : lsst.pipe.base.connections.InputQuantizedConnection
        `outputRefs` : lsst.pipe.base.connections.OutputQuantizedConnection

        """

        inputs = butlerQC.get(inputRefs)
        bands = []
        for filterName in self.config.filterNames:
            bands.append(self.config.referenceCatalogLoader.refObjLoader.filterMap[filterName])

        columns = self.prepColumns(bands)
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
        outputs = self.run(catalog=inputs["catalog"], loadedRefCat=loadedRefCat, bands=bands)

        butlerQC.put(outputs, outputRefs)
