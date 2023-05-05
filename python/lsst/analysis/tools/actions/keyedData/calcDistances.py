__all__ = ("CalcDistances",)

import astropy.units as u
import numpy as np
import pandas as pd
from lsst.pex.config import Field
from lsst.pex.config.configurableActions import ConfigurableActionField

from ...interfaces import KeyedData, KeyedDataAction, KeyedDataSchema, Vector, VectorAction
from ..vector import LoadVector


class CalcDistances(KeyedDataAction):
    """Calculate distances in a matched catalog"""

    groupKey = Field[str](doc="Column key to use for forming groups", default="obj_index")
    annulus = Field[float](doc="Radial distance of the annulus in arcmin", default=5.0)
    width = Field[float](doc="Width of annulus in arcmin", default=2.0)
    mags = ConfigurableActionField[VectorAction](doc="Action which supplies magnitudes")
    minMag = Field[float](doc="Minimum magnitude", default=17)
    maxMag = Field[float](doc="Maximum magnitude", default=21.5)

    def getInputSchema(self) -> KeyedDataSchema:
        return tuple(self.mags.getInputSchema()) + ((self.groupKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        # This part comes from AMxTask/astromResiduals in Faro
        D = self.annulus * u.arcmin
        width = self.width * u.arcmin
        annulus = ((D + (width / 2) * np.array([-1, +1])) * u.arcmin).to(u.radian)

        # This part copies calcRmsDistances/calcSepOutliers in Faro
        inMagRange = (self.mags(data) > self.minMag) & (self.mags(data) < self.maxMag)

        # TODO: this line needs to add ra, dec,
        df = pd.DataFrame({"groupKey": data[self.groupKey], "value": LoadVector(data, **kwargs)})
        # result = df.groupby("groupKey")["value"].aggregate(self.func)
        meanRa = df.groupby("groupKey")["coord_ra"].aggregate("mean")
        meanDec = df.groupby("groupKey")["coord_dec"].aggregate("mean")
        print(inMagRange.sum(), meanRa, meanDec, annulus)
        import ipdb

        ipdb.set_trace()
