import logging
import os
import pathlib

from core.model.Envelope import Envelope
from core.model.SystemCommand import SystemCommand

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

        self._logger = logger

        # Default to SRTM DEM.
        self._demFootprintsPath = pathlib.Path(self.SRTM_FOOTPRINTS_PATH)

        self._demDir = pathlib.Path(demDir)

        if not self._demFootprintsPath.exists():

            raise FileNotFoundError(self._demFootprintsPath)

        if self._logger:

            self._logger.info('DEM footprints shapefile: ' +
                              str(self._demFootprintsPath))

        if not self._demDir.exists():

            errorMessage = f'{self._demDir} is expected to be created as'
            ' part of the EvhrToA process but does not exist.'

            raise FileNotFoundError(errorMessage)

        self._validateDem()

    # ------------------------------------------------------------------------
    # mosaicAndClipDemTiles
    # ------------------------------------------------------------------------
    def mosaicAndClipDemTiles(self, outDemName: pathlib.Path,
                              envelope: Envelope) -> None:

        if self._logger:
            self._logger.info(f'Creating DEM {str(outDemName)}')

        outDemNameTemp = outDemName.replace('.tif', '-temp.tif')

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

            self._validateDem()

        features = self._clipShp(envelope)

        if not features or len(features) == 0:

            msg = 'Clipping rectangles to supplied DEM footprints did not' + \
                f' return any features. Corners: ({str(envelope.ulx())},' + \
                f' {str(envelope.uly())}), ({str(envelope.lrx())},' + \
                f' {str(envelope.lry())})'

            return RuntimeError(msg)

        tiles = []

        for feature in features:

            tileFile = str(feature.
                           getElementsByTagName('ogr:location')[0].
                           firstChild.
                           data)

            tiles.append(tileFile)

        # Mosaic the tiles.
        cmd = 'gdal_merge.py' + \
              ' -o ' + outDemNameTemp + \
              ' -ul_lr' + \
              ' ' + str(envelope.ulx()) + \
              ' ' + str(envelope.uly()) + \
              ' ' + str(envelope.lrx()) + \
              ' ' + str(envelope.lry()) + \
              ' ' + ' '.join(tiles)

        SystemCommand(cmd, self._logger, True)

        # Run mosaicked DEM through geoid correction
        cmd = self.BASE_SP_CMD + \
            'dem_geoid ' + \
            outDemNameTemp + ' --geoid EGM96 -o ' + \
            outDemName.strip('-adj.tif') + \
            ' --reverse-adjustment'

        SystemCommand(cmd, self._logger, True)

        for log in self._demDir.glob('*log*.txt'):
            log.unlink()  # remove dem_geoid log file(s)
