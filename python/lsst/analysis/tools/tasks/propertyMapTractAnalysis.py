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

__all__ = [
    "PropertyMapTractAnalysisConfig",
    "PropertyMapTractAnalysisTask",
]


from lsst.daf.butler import DataCoordinate
from lsst.pipe.base import (ButlerQuantumContext, InputQuantizedConnection,
                            OutputQuantizedConnection)
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap

from ..interfaces import (AnalysisBaseConfig, AnalysisBaseConnections,
                          AnalysisPipelineTask)


class PropertyMapTractAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract", "band"),
    defaultTemplates={"outputName": "propertyMapTract", "propertyMapName": "deepCoadd_exposure_time_map_sum"},
):

    # Either:
    # " band is in the dimensions and multiple=False "
    # OR
    # " band is not in the dimensions and multiple=True "

    skymap = ct.Input(
        doc="The skymap that covers the tract that the data is from.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )

    data = ct.Input(
        doc="Sum-value map of exposure time",
        name="{propertyMapName}",
        storageClass="HealSparseMap",
        dimensions=("tract", "skymap", "band"),
        multiple=False,
        deferLoad=False,  # True if your input are > 4GB (we are limited to 4GB per quantum)
    )


class PropertyMapTractAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=PropertyMapTractAnalysisConnections
):
    pass


class PropertyMapTractAnalysisTask(AnalysisPipelineTask):
    ConfigClass = PropertyMapTractAnalysisConfig
    _DefaultName = "propertyMapTractAnalysis"

    def runQuantum(
        self,
        butlerQC: ButlerQuantumContext,
        inputRefs: InputQuantizedConnection,
        outputRefs: OutputQuantizedConnection,
    ) -> None:
        inputs = butlerQC.get(inputRefs)
        if "skymap" in inputs.keys():
            skymap = inputs["skymap"]
        else:
            skymap = None
        outputs = self.run(data=inputs["data"], skymap=skymap)
        butlerQC.put(outputs, outputRefs)
