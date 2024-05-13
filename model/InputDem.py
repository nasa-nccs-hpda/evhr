import logging
import pathlib
import tempfile
from xml.dom import minidom

from osgeo import ogr
from osgeo import gdal

from core.model.Envelope import Envelope
from core.model.SystemCommand import SystemCommand


# ----------------------------------------------------------------------------
# class InputDem
#
# EvhrToa -> InputDem
# ----------------------------------------------------------------------------
class InputDem(object):
    """
    Class which handles verification and mosaicking of user-supplied DEM given
    a footprints shapefile of the DEM.
    """

    BASE_SP_CMD = '/opt/StereoPipeline/bin/'

    # ------------------------------------------------------------------------
    # __init__
    # ------------------------------------------------------------------------
    def __init__(self, demFootprintsPath: pathlib.Path, demDir: pathlib.Path,
                 logger: logging.Logger = None) -> None:

        self._logger = logger

        self._demFootprintsPath = pathlib.Path(demFootprintsPath)

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

        features = self._clipShp(envelope)

        if not features or len(features) == 0:

            msg = 'Clipping rectangles to supplied DEM footprints did not'
            f' return any features. Corners: ({str(envelope.ulx())},'
            f' {str(envelope.uly())}), ({str(envelope.lrx())},'
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

    # ------------------------------------------------------------------------
    # _clipShp
    # ------------------------------------------------------------------------
    def _clipShp(self, envelope: Envelope):

        if self._logger:
            self._logger.info(f'Clipping shapefile {self._demFootprintsPath}')

        # Create a temporary file for the clip output.
        tempClipFile = tempfile.mkstemp()[1]

        # ---
        # To filter scenes that only overlap the AoI slightly, decrease both
        # corners of the query AoI.
        # ---
        MIN_OVERLAP_IN_DEGREES = 0.02
        ulx = float(envelope.ulx()) + MIN_OVERLAP_IN_DEGREES
        uly = float(envelope.uly()) - MIN_OVERLAP_IN_DEGREES
        lrx = float(envelope.lrx()) - MIN_OVERLAP_IN_DEGREES
        lry = float(envelope.lry()) + MIN_OVERLAP_IN_DEGREES

        # Clip.  The debug option somehow prevents an occasional seg. fault!
        cmd = 'ogr2ogr' + \
              ' -f "GML"' + \
              ' -spat' + \
              ' ' + str(ulx) + \
              ' ' + str(lry) + \
              ' ' + str(lrx) + \
              ' ' + str(uly) + \
              ' -spat_srs' + \
              ' "' + envelope.GetSpatialReference().ExportToProj4() + '"' + \
              ' --debug on'

        # ---
        # This was in the EVHR project, but I do not remember why it is
        # important.
        # ---
        MAXIMUM_SCENES = 5
        cmd += ' -limit ' + str(MAXIMUM_SCENES)

        cmd += ' "' + tempClipFile + '"' + \
               ' "' + str(self._demFootprintsPath) + '"'

        SystemCommand(cmd, self._logger, True)

        xml = minidom.parse(tempClipFile)
        features = xml.getElementsByTagName('gml:featureMember')

        return features

    # -------------------------------------------------------------------------
    # _validateDem
    # -------------------------------------------------------------------------
    def _validateDem(self) -> None:

        demShapeFilePath = self._demFootprintsPath

        # Test that ogr can open the footprints file and access features.
        try:

            demShapeSource = ogr.Open(str(demShapeFilePath))

            demShapeLayer = demShapeSource.GetLayer()

            demShapeSampleFeature = demShapeLayer.GetNextFeature()

            demShapeSampleLocation = demShapeSampleFeature.GetField("location")

        except Exception as e:

            errorMessage = 'Error opening or accessing layers and' + \
                f' fields of {str(demShapeFilePath)}, Error: {e}'

            raise RuntimeError(errorMessage)

        # Get a path to the DEM, check if it exists
        demShapeSamplePath = pathlib.Path(demShapeSampleLocation)

        if not demShapeSamplePath.exists():

            errorMessage = 'DEM footprints shapefile must point to a file' + \
                f' that exist. Cannot find: {demShapeSamplePath}'

            raise FileNotFoundError(errorMessage)

        # Attempt to open the sampled DEM with GDAL.
        try:

            dataset = gdal.Open(str(demShapeSamplePath), gdal.GA_ReadOnly)

            if dataset is None:

                errorMessage = 'Dataset is none. Unable to open' + \
                    f' {str(demShapeSamplePath)}'

                raise RuntimeError(errorMessage)

            dataset = None

        except Exception as e:

            errorMessage = f'Error opening {str(demShapeSamplePath)}' + \
                ' using GDAL. Please pass a valid footprints file which' + \
                ' contains a field "location" that is a path to a valid' + \
                f' DEM readable by GDAL. Error: {e}'

            raise RuntimeError(errorMessage)
