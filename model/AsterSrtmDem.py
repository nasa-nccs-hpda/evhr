import logging
import os
import pathlib

from core.model.Envelope import Envelope

from evhr.model.InputDem import InputDem


# ----------------------------------------------------------------------------
# class AsterSrtmDem
#
# EvhrToa -> AsterSrtmDem
# ----------------------------------------------------------------------------
class AsterSrtmDem(InputDem):
    """
    Class which handles verification and mosaicking of adapt-only DEM given
    a footprints shapefile of the DEM.
    """

    SRTM_FOOTPRINTS_PATH: str = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 'SRTM/srtm.shp'
    )
    ASTER_FOOTPRINTS_PATH: str = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 'ASTERGDEM/astergdem.shp'
    )

    # ------------------------------------------------------------------------
    # __init__
    # ------------------------------------------------------------------------
    def __init__(self, demDir: pathlib.Path,
                 logger: logging.Logger = None) -> None:
        super().__init__(self.SRTM_FOOTPRINTS_PATH, demDir, logger)

    # ------------------------------------------------------------------------
    # mosaicAndClipDemTiles
    # ------------------------------------------------------------------------
    def mosaicAndClipDemTiles(self,
                              outDemName: pathlib.Path,
                              envelope: Envelope) -> None:
        # ---
        # SRTM was collected between -54 and 60 degrees of latitude. Use
        # ASTERGDEM where SRTM is unavailable
        # ---
        if envelope.uly() < -54.0 or envelope.uly() > 60.0 or \
                envelope.lry() < -54.0 or envelope.lry() > 60.0:

            if self._logger:
                self._logger.info(f'{envelope} is outside bounds or SRTM ' +
                                  'DEM, switching to ASTERGDEM')

            self._demFootprintsPath = pathlib.Path(self.ASTER_FOOTPRINTS_PATH)

        super().mosaicAndClipDemTiles(outDemName, envelope)
