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
    "Scalar",
    "ScalarType",
    "KeyedData",
    "KeyedDataTypes",
    "KeyedDataSchema",
    "Vector",
    "PlotTypes",
    "KeyedResults",
)

from abc import ABCMeta
from numbers import Number
from typing import Any, Iterable, Mapping, MutableMapping, TypeAlias

import numpy as np
from healsparse import HealSparseMap
from lsst.verify import Measurement
from matplotlib.figure import Figure
from numpy.typing import NDArray


class ScalarMeta(ABCMeta):
    def __instancecheck__(cls: ABCMeta, instance: Any) -> Any:
        return isinstance(instance, tuple(cls.mro()[1:]))


class Scalar(Number, np.number, metaclass=ScalarMeta):  # type: ignore
    """This is an interface only class, and is intended to abstract around all
    the various types of numbers used in Python.

    This has been tried many times with various levels of success in python,
    and this is another attempt. However, as this class is only intended as an
    interface, and not something concrete to use it works.

    Users should not directly instantiate from this class, instead they should
    use a built in python number type, or a numpy number.
    """

    def __init__(self) -> None:
        raise NotImplementedError("Scalar is only an interface and should not be directly instantiated")


ScalarType = type[Scalar]
"""A type alias for the Scalar interface."""

Vector = NDArray
"""A Vector is an abstraction around the NDArray interface, things that 'quack'
like an NDArray should be considered a Vector.
"""

KeyedData = MutableMapping[str, Vector | Scalar | HealSparseMap]
"""KeyedData is an interface where either a `Vector` or `Scalar` can be
retrieved using a key which is of str type.
"""

KeyedDataTypes = MutableMapping[str, type[Vector] | ScalarType | type[HealSparseMap]]
r"""A mapping of str keys to the Types which are valid in `KeyedData` objects.
This is useful in conjunction with `AnalysisAction`\ 's ``getInputSchema`` and
``getOutputSchema`` methods.
"""

KeyedDataSchema = Iterable[tuple[str, type[Vector] | ScalarType | type[HealSparseMap]]]
r"""An interface that represents a type returned by `AnalysisAction`\ 's
``getInputSchema`` and ``getOutputSchema`` methods.
"""

PlotTypes = Figure
"""An interface that represents the various plot types analysis tools supports.
"""

KeyedResults = Mapping[str, PlotTypes | Measurement]
"""A mapping of the return types for an analysisTool."""

MetricResultType: TypeAlias = Mapping[str, Measurement] | Measurement
"""A type alias for the return type of a MetricAction."""

PlotResultType: TypeAlias = Mapping[str, PlotTypes] | PlotTypes
"""A type alias for the return type of a PlotAction."""
