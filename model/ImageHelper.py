from pathlib import Path
import sys

import numpy

from osgeo import gdal


# ----------------------------------------------------------------------------
# Class ImageHelper
#
# TODO: add accessors using @property.
# ----------------------------------------------------------------------------
class ImageHelper(object):

    # ------------------------------------------------------------------------
    # init
    # ------------------------------------------------------------------------
    def __init__(self):

        self._dataset: gdal.Dataset = None
        self._inputFile: Path = None
        self._redBandId: int = 1
        self._greenBandId: int = 1
        self._blueBandId: int = 1
        self._redBand: numpy.ndarray = None
        self._greenBand: numpy.ndarray = None
        self._blueBand: numpy.ndarray = None
        self._noDataValue: float = -9999.0
        self._minValue: float = sys.float_info.max
        self._maxValue: float = sys.float_info.min

    # ------------------------------------------------------------------------
    # initFromDataset
    # ------------------------------------------------------------------------
    def initFromDataset(self,
                        dataset: gdal.Dataset,
                        noDataValue: float,
                        redBandId: int = 1,
                        greenBandId: int = 1,
                        blueBandId: int = 1,
                        ) -> None:

        self._dataset = dataset

        self._completeInitialization(noDataValue,
                                     redBandId,
                                     greenBandId,
                                     blueBandId)

    # ------------------------------------------------------------------------
    # initFromFile
    # ------------------------------------------------------------------------
    def initFromFile(self,
                     inputFile: Path,
                     noDataValue: float,
                     redBandId: int = 1,
                     blueBandId: int = 1,
                     greenBandId: int = 1,
                     ) -> None:

        # TODO: check validity and comment.
        self._inputFile = inputFile
        self._dataset: gdal.Dataset = gdal.Open(str(self._inputFile))

        self._completeInitialization(noDataValue,
                                     redBandId,
                                     greenBandId,
                                     blueBandId)

    # ------------------------------------------------------------------------
    # completeInitialization
    # ------------------------------------------------------------------------
    def _completeInitialization(self,
                                noDataValue: float,
                                redBandId: int = 1,
                                greenBandId: int = 1,
                                blueBandId: int = 1,
                                ) -> None:

        self._redBandId = redBandId
        self._greenBandId = greenBandId
        self._blueBandId = blueBandId

        # Generate overviews, if they do not exist.
        if self._dataset.GetRasterBand(1).GetOverviewCount() == 0:
            dummy = self. _dataset.BuildOverviews()

        # ---
        # Read the bands.
        # ---
        self._redBand: numpy.ndarray = \
            self._dataset.GetRasterBand(self._redBandId).ReadAsArray()

        self._greenBand: numpy.ndarray = \
            self._dataset.GetRasterBand(self._greenBandId).ReadAsArray()

        self._blueBand: numpy.ndarray = \
            self._dataset.GetRasterBand(self._blueBandId).ReadAsArray()

        # ---
        # Initialize the no-data value.
        # ---
        self._noDataValue = self._dataset.GetRasterBand(self._redBandId). \
            GetNoDataValue() or noDataValue

        # ---
        # Compute the minimum and maximum pixels values to help the renderer.
        # ---
        forExtremes = self._redBand[self._redBand != self._noDataValue]
        self._minValue: float = forExtremes.min()
        self._maxValue: float = forExtremes.max()

    # ------------------------------------------------------------------------
    # getCorners
    #
    # Why is this not built into gdal?
    # ------------------------------------------------------------------------
    def getCorners(self):

        minx, xres, xskew, maxy, yskew, yres = \
            self._dataset.GetGeoTransform()

        maxx = minx + (self._dataset.RasterXSize * xres)
        miny = maxy + (self._dataset.RasterYSize * yres)

        return (minx, miny, maxx, maxy)

    # ------------------------------------------------------------------------
    # getRgbIndexList
    # ------------------------------------------------------------------------
    def getRgbIndexList(self) -> list:

        return [self._redBandId, self._greenBandId, self._blueBandId]

    # ------------------------------------------------------------------------
    # getRgbBands
    # ------------------------------------------------------------------------
    def getRgbBands(self) -> list:

        return [self._redBand, self._greenBand, self._blueBand]

    # ------------------------------------------------------------------------
    # __str__
    # ------------------------------------------------------------------------
    def __str__(self):

        return ('Input file: ' + str(self._inputFile) +
                '\nMin. pixel: ' + str(self._minValue) +
                '\nMax. pixel: ' + str(self._maxValue) +
                '\nNo-data value: ' + str(self._noDataValue) +
                '\nRed band index: ' + str(self._redBandId) +
                '\nGreen band index: ' + str(self._greenBandId) +
                '\nBlue band index: ' + str(self._blueBandId) +
                '\nCorners: ' + str(self.getCorners())
                )
