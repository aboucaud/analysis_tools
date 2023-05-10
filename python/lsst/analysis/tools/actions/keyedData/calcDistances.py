__all__ = ("CalcDistances",)

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from lsst.faro.utils.separations import matchVisitComputeDistance
from lsst.pex.config import Field

from ...interfaces import KeyedData, KeyedDataAction, KeyedDataSchema, Vector


class CalcDistances(KeyedDataAction):
    """Calculate distances in a matched catalog"""

    groupKey = Field[str](doc="Column key to use for forming groups", default="obj_index")
    visitKey = Field[str](doc="Column key to use for matching visits", default="visit")
    raKey = Field[str](doc="RA column key", default="coord_ra")
    decKey = Field[str](doc="Dec column key", default="coord_dec")
    annulus = Field[float](doc="Radial distance of the annulus in arcmin", default=5.0)
    width = Field[float](doc="Width of annulus in arcmin", default=2.0)
    threshAD = Field[float](doc="Threshold in mas for AFx calculation.", default=20.0)
    threshAF = Field[float](
        doc="Percentile of differences that can vary by more than threshAD.", default=10.0
    )

    def getInputSchema(self) -> KeyedDataSchema:
        return (  # tuple(self.mags.getInputSchema()) + (
            (self.groupKey, Vector),
            (self.raKey, Vector),
            (self.decKey, Vector),
            (self.visitKey, Vector),
        )

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        # This part comes from AMxTask/astromResiduals in Faro
        print("Annulus, data:", self.annulus, len(data["coord_ra"]))
        D = self.annulus * u.arcmin
        width = self.width * u.arcmin
        annulus = (D + (width / 2) * np.array([-1, +1])).to(u.radian)

        # TODO: this line needs to add ra, dec,
        df = pd.DataFrame(
            {
                "groupKey": data[self.groupKey],
                "coord_ra": data[self.raKey],
                "coord_dec": data[self.decKey],
                "visit": data[self.visitKey],
            }
        )

        meanRa = df.groupby("groupKey")["coord_ra"].aggregate("mean")
        meanDec = df.groupby("groupKey")["coord_dec"].aggregate("mean")

        catalog = SkyCoord(meanRa.to_numpy() * u.degree, meanDec.to_numpy() * u.degree)
        idx, idxCatalog, d2d, d3d = catalog.search_around_sky(catalog, annulus[1])
        inAnnulus = d2d > annulus[0]
        idx = idx[inAnnulus]
        idxCatalog = idxCatalog[inAnnulus]

        rmsDistances = []
        sepResiduals = []
        sepResiduals2 = []
        tt1 = 0
        tt2 = 0
        tt3 = 0
        tt4 = 0
        import time

        for id in range(len(meanRa[:100])):
            match_inds = idx == id
            match_ids = idxCatalog[match_inds & (idxCatalog != id)]
            if match_ids.sum() == 0:
                continue
            if id % 100 == 0:
                print(id, len(match_ids))
            object_srcs = df.loc[df["groupKey"] == meanRa.index[id]]

            object_visits = object_srcs["visit"].to_numpy()
            object_ras = (object_srcs["coord_ra"].to_numpy() * u.degree).to(u.radian).value
            object_decs = (object_srcs["coord_dec"].to_numpy() * u.degree).to(u.radian).value
            if len(object_srcs) <= 1:
                continue
            object_srcs = object_srcs.set_index("visit")
            object_srcs.sort_index(inplace=True)

            for id2 in match_ids:
                match_srcs = df.loc[df["groupKey"] == meanRa.index[id2]]
                match_visits = match_srcs["visit"].to_numpy()
                match_ras = (match_srcs["coord_ra"].to_numpy() * u.degree).to(u.radian).value
                match_decs = (match_srcs["coord_dec"].to_numpy() * u.degree).to(u.radian).value

                t2 = time.time()
                separations = matchVisitComputeDistance(
                    object_visits, object_ras, object_decs, match_visits, match_ras, match_decs
                )

                t3 = time.time()

                if len(separations) > 1:
                    rmsDist = np.std(separations, ddof=1)
                    rmsDistances.append(rmsDist)
                if len(separations) > 2:
                    sepResids = np.abs(separations - np.median(separations))
                    sepResiduals.append(sepResids)
                    sepResiduals2.append(separations - np.median(separations))
                t4 = time.time()

                # tt1 += t1 - t0
                # tt2 += t2 - t1
                tt3 += t3 - t2
                tt4 += t4 - t3
            if id % 100 == 0:
                # print(t1 - t0, t2- t1, t3- t2, t4 - t3)
                print(t3 - t2, t4 - t3)
        print(tt1, tt2, tt3, tt4)
        if len(rmsDistances) == 0:
            AMx = np.nan * u.marcsec
        else:
            AMx = (np.median(rmsDistances) * u.radian).to(u.marcsec)

        if len(sepResiduals) <= 1:
            AFx = np.nan * u.percent
            ADx = np.nan * u.marcsec
            absDiffSeparations = np.array([]) * u.marcsec
            absDiffSeparations2 = np.array([]) * u.marcsec

        else:
            sepResiduals = np.concatenate(sepResiduals)
            sepResiduals2 = np.concatenate(sepResiduals2)
            absDiffSeparations = ((sepResiduals - np.median(sepResiduals)) * u.radian).to(u.marcsec)
            absDiffSeparations2 = (abs(sepResiduals2 - np.median(sepResiduals2)) * u.radian).to(u.marcsec)
            afThreshhold = 100.0 - self.threshAF
            ADx = np.percentile(absDiffSeparations, afThreshhold)
            AFx = 100 * np.mean(np.abs(absDiffSeparations) > self.threshAD * u.marcsec) * u.percent

        distanceParams = {
            "rmsDistances": (np.array(rmsDistances) * u.radian).to(u.marcsec).value,
            "separationResiduals": absDiffSeparations2.value,
            "AMx": AMx.value,
            "ADx": ADx.value,
            "AFx": AFx.value,
        }
        print(distanceParams)
        return distanceParams
