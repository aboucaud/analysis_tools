# This file is part of analysis_tools.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import annotations

__all__ = ("SasquatchDatastore",)

"""Sasquatch datastore"""
import logging
from collections.abc import Iterable, Mapping, Sequence
from typing import TYPE_CHECKING, Any, ClassVar

from lsst.daf.butler import (
    DatasetRef,
    DatasetRefURIs,
    DatasetTypeNotSupportedError,
    DatastoreRecordData,
    StorageClass,
)
from lsst.daf.butler.datastores.genericDatastore import GenericBaseDatastore
from lsst.daf.butler.registry.interfaces import DatastoreRegistryBridge
from lsst.resources import ResourcePath

from . import SasquatchDispatcher

if TYPE_CHECKING:
    from lsst.daf.butler import Config, DatasetType, LookupKey
    from lsst.daf.butler.registry.interfaces import DatasetIdRef, DatastoreRegistryBridgeManager


log = logging.getLogger(__name__)


class SasquatchDatastore(GenericBaseDatastore):
    """Basic Datastore for writing to an in Sasquatch instance.

    This Datastore is currently write only, meaning that it can dispatch data
    to a Sasquatch instance, but at the present can not be used to retrieve
    values.


    Parameters
    ----------
    config : `DatastoreConfig` or `str`
        Configuration.
    bridgeManager : `DatastoreRegistryBridgeManager`
        Object that manages the interface between `Registry` and datastores.
    butlerRoot : `str`, optional
        Unused parameter.
    """

    defaultConfigFile: ClassVar[str | None] = "sasquatchDatastore.yaml"
    """Path to configuration defaults. Accessed within the ``configs`` resource
    or relative to a search path. Can be None if no defaults specified.
    """

    restProxyUrl: str
    """Url which points to the http rest proxy. This is where datasets will be
    dispatched to.
    """

    accessToken: str
    """Access token which is used to authenticate to the restProxy.
    """

    namespace: str
    """The namespace in Sasquatch where the uploaded metrics will be
    dispatched.
    """

    def __init__(
        self,
        config: Config | str,
        bridgeManager: DatastoreRegistryBridgeManager,
        butlerRoot: str | None = None,
    ):
        super().__init__(config, bridgeManager)

        # Name ourselves either using an explicit name or a name
        # derived from the (unexpanded) root.
        self.name = self.config.get("name", "{}@{}".format(type(self).__name__, self.config["restProxyUrl"]))
        log.debug("Creating datastore %s", self.name)

        self._bridge = bridgeManager.register(self.name, ephemeral=False)

        self.restProxyUrl = self.config["restProxyUrl"]

        self.accessToken = self.config.get("accessToken", "na")

        self.namespace = self.config.get("namespace", "lsst.dm")

        self._dispatcher = SasquatchDispatcher(self.restProxyUrl, self.accessToken, self.namespace)

    @property
    def bridge(self) -> DatastoreRegistryBridge:
        return self._bridge

    def put(self, inMemoryDataset: Any, ref: DatasetRef) -> None:
        if self.constraints.isAcceptable(ref):
            self._dispatcher.dispatchRef(inMemoryDataset, ref)
        else:
            log.debug("Could not put dataset type %s with Sasquatch datastore", ref.datasetType)
            raise DatasetTypeNotSupportedError(
                f"Could not put dataset type {ref.datasetType} with Sasquatch datastore"
            )

    def addStoredItemInfo(self, refs: Iterable[DatasetRef], infos: Iterable[Any]) -> None:
        raise NotImplementedError()

    def getStoredItemsInfo(self, ref: DatasetRef) -> Sequence[Any]:
        raise NotImplementedError()

    def removeStoredItemInfo(self, ref: DatasetRef) -> None:
        raise NotImplementedError()

    def trash(self, ref: DatasetRef | Iterable[DatasetRef], ignore_errors: bool = True) -> None:
        log.debug("Sasquatch datastore does not support trashing skipping %s", ref)
        raise FileNotFoundError()

    def emptyTrash(self, ignore_errors: bool = True) -> None:
        log.debug("Sasquatch datastore does not support trash, nothing to empty")

    def forget(self, ref: Iterable[DatasetRef]) -> None:
        pass

    def exists(self, datasetRef: DatasetRef) -> bool:
        # sasquatch is not currently searchable
        return False

    def knows(self, ref: DatasetRef) -> bool:
        return False

    def get(
        self,
        datasetRef: DatasetRef,
        parameters: Mapping[str, Any] | None = None,
        storageClass: StorageClass | str | None = None,
    ) -> Any:
        raise FileNotFoundError()

    def validateConfiguration(
        self, entities: Iterable[DatasetRef | DatasetType | StorageClass], logFailures: bool = False
    ) -> None:
        """Validate some of the configuration for this datastore.

        Parameters
        ----------
        entities : iterable of `DatasetRef`, `DatasetType`, or `StorageClass`
            Entities to test against this configuration.  Can be differing
            types.
        logFailures : `bool`, optional
            If `True`, output a log message for every validation error
            detected.

        Raises
        ------
        DatastoreValidationError
            Raised if there is a validation problem with a configuration.
            All the problems are reported in a single exception.

        Notes
        -----
        This method is a no-op.
        """
        return

    def validateKey(self, lookupKey: LookupKey, entity: DatasetRef | DatasetType | StorageClass) -> None:
        # Docstring is inherited from base class.
        return

    def getLookupKeys(self) -> set[LookupKey]:
        # Docstring is inherited from base class.
        return self.constraints.getLookupKeys()

    def needs_expanded_data_ids(
        self,
        transfer: str | None,
        entity: DatasetRef | DatasetType | StorageClass | None = None,
    ) -> bool:
        # Docstring inherited.
        return False

    def import_records(self, data: Mapping[str, DatastoreRecordData]) -> None:
        # Docstring inherited from the base class.
        return

    def export_records(self, refs: Iterable[DatasetIdRef]) -> Mapping[str, DatastoreRecordData]:
        # Docstring inherited from the base class.

        # Sasquatch Datastore records cannot be exported or imported.
        return {}

    def getURI(self, datasetRef: DatasetRef, predict: bool = False) -> ResourcePath:
        raise NotImplementedError()

    def getURIs(self, datasetRef: DatasetRef, predict: bool = False) -> DatasetRefURIs:
        raise NotImplementedError()

    def retrieveArtifacts(
        self,
        refs: Iterable[DatasetRef],
        destination: ResourcePath,
        transfer: str = "auto",
        preserve_path: bool = True,
        overwrite: bool = False,
    ) -> list[ResourcePath]:
        raise NotImplementedError()

    @classmethod
    def setConfigRoot(cls, root: str, config: Config, full: Config, overwrite: bool = True) -> None:
        pass
