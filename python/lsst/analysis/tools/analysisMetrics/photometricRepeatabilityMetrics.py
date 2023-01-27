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

__all__ = ("StellarPhotometricRepeatabilityMetric",)

from ..analysisParts.photometricRepeatability import StellarPhotometricRepeatabilityMixin
from ..interfaces import AnalysisMetric


class StellarPhotometricRepeatabilityMetric(AnalysisMetric, StellarPhotometricRepeatabilityMixin):
    """Compute photometric repeatability from multiple measurements of a set of
    stars. First, a set of per-source quality criteria are applied. Second,
    the individual source measurements are grouped together by object index
    and per-group quantities are computed (e.g., a representative S/N for the
    group based on the median of associated per-source measurements). Third,
    additional per-group criteria are applied. Fourth, summary statistics are
    computed for the filtered groups.
    """

    def setDefaults(self):
        super().setDefaults()

        self.produce.units = {  # type: ignore
            "photRepeatStdev": "mmag",
            "photRepeatOutlier": "percent",
            "photRepeatNsources": "ct",
        }
        self.produce.newNames = {
            "photRepeatStdev": "{band}_stellarPhotRepeatStdev",
            "photRepeatOutlier": "{band}_stellarPhotRepeatOutlierFraction",
            "photRepeatNsources": "{band}_ct",
        }
