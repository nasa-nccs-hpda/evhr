
import filecmp
import logging
import os
from pathlib import Path
import shutil
import tempfile

from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from osgeo.osr import CoordinateTransformation
from osgeo.osr import SpatialReference

from core.model.BaseFile import BaseFile
from core.model.DgFile import DgFile
from core.model.Envelope import Envelope
from core.model.GeospatialImageFile import GeospatialImageFile
from core.model.SystemCommand import SystemCommand

from evhr.model.ToaCalculation import ToaCalculation
from evhr.model.InputDem import InputDem
from evhr.model.AsterSrtmDem import AsterSrtmDem
from evhr.model.EvhrUtils import EvhrUtils


# -----------------------------------------------------------------------------
# class EvhrToA
#
# run
#   _collectImagesByStrip
#   _computeEnvelope
#   _runOneStrip
#       _createStrip
#           _scenesToStripFromBandList
#           _stripToToa
#               _orthoOne
#                   _createDemForOrthos
#       _stripToToa
#           _ToaCalculation.run
#           _mergeBands
# -----------------------------------------------------------------------------
class EvhrToA(object):

    BASE_SP_CMD = '/opt/StereoPipeline/bin/'
    MAXIMUM_SCENES = 100
    NO_DATA_VALUE = -10001

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self, outDir, inputDemPath=None, panResolution=1,
                 panSharpen=False, logger=None):

        self._logger = logger

        # Output directory
        if not os.path.exists(outDir):
            os.mkdir(outDir)

        self._outDir = BaseFile(outDir).fileName()  # BaseFile tests validity.

        if not os.path.isdir(self._outDir):
            raise RuntimeError(self._outDir + ' must be a directory.')

        self._outSrsProj4 = None
        self._panResolution = panResolution
        self._panSharpen = panSharpen

        # Ensure the ortho and toa directories exist.
        self._bandDir = os.path.join(self._outDir, '1-bands')
        self._stripDir = os.path.join(self._outDir, '2-strips')
        self._demDir = os.path.join(self._outDir, '3-dems')
        self._orthoDir = os.path.join(self._outDir, '4-orthos')
        self._toaDir = os.path.join(self._outDir, '5-toas')

        for d in [self._bandDir, self._stripDir, self._demDir,
                  self._orthoDir, self._toaDir]:

            if not os.path.exists(d):
                os.mkdir(d)

        if self._logger:

            self._logger.info('Output directory: ' + self._outDir)
            self._logger.info('Band directory: ' + self._bandDir)
            self._logger.info('Strip directory: ' + self._stripDir)
            self._logger.info('DEM directory: ' + self._demDir)
            self._logger.info('Ortho image directory: ' + self._orthoDir)
            self._logger.info('ToA directory: ' + self._toaDir)

        if inputDemPath:

            self._inputDem = InputDem(inputDemPath,
                                      self._demDir,
                                      self._logger)
        else:

            if self._logger:
                self._logger.info('No user-supplied DEM, defaulting to' +
                                  ' ADAPT SRTM/ASTERGDEM')

            self._inputDem = AsterSrtmDem(self._demDir, self._logger)

    # -------------------------------------------------------------------------
    # collectImagesByStrip
    #
    # run
    #    _collectImagesByStrip()
    #
    # This returns a map of strip names to images in the strip.
    #
    # {strip1: [image, image, ...], strip2: [image, image, ...], ...}
    # -------------------------------------------------------------------------
    def _collectImagesByStrip(self, dgScenes: list) -> dict:

        if self._logger:

            self._logger.info('In collectImagesByStrip')
            self._logger.info('Collating ' + str(len(dgScenes)) + ' scenes')

        stripsWithScenes = {}

        for dgf in dgScenes:

            stripID = dgf.getStripName()

            if stripID:

                if stripID not in stripsWithScenes:

                    stripsWithScenes[stripID] = []

                if dgf not in stripsWithScenes[stripID]:

                    stripsWithScenes[stripID].append(dgf)

                else:
                    if self._logger:
                        self._logger.info(str(dgf) + ' already added')

            else:

                if self._logger:
                    self._logger.warn('Unable to get strip name for: ' +
                                      stripID)

        return stripsWithScenes

    # -------------------------------------------------------------------------
    # _computeEnvelope
    # -------------------------------------------------------------------------
    def _computeEnvelope(self, sceneList: list) -> Envelope:

        if self._logger:
            self._logger.info('In _computeEnvelope')

        env = sceneList[0].envelope()

        for dgf in sceneList[1:]:

            env.addPoint(dgf._ulx, dgf._uly, 0.0, dgf.srs())
            env.addPoint(dgf._lrx, dgf._lry, 0.0, dgf.srs())

        return env

    # -------------------------------------------------------------------------
    # _createDemForOrthos
    #
    # run -> runOneStrip -> stripToToa -> orthoOne -> createDemForOrthos
    # -------------------------------------------------------------------------
    @staticmethod
    def _createDemForOrthos(inputDem, envelope, demDir, logger):

        if logger:
            logger.info('Creating DEM for orthorectification.')

        # If there is already a clipped DEM for this bounding box, use it.
        demName = 'dem-' + \
            str(envelope.ulx()) + '-' + \
            str(envelope.uly()) + '-' + \
            str(envelope.lrx()) + '-' + \
            str(envelope.lry()) + '-' + \
            str(envelope.GetSpatialReference().GetAuthorityCode(None)) + \
            '-adj.tif'

        demName = os.path.join(demDir, demName)

        if os.path.exists(demName):
            return demName

        # Expand the bounding box before clipping the DEM.
        tempEnv = Envelope()

        tempEnv.addPoint(envelope.ulx(),
                         envelope.uly(),
                         0.0,
                         envelope.GetSpatialReference())

        tempEnv.addPoint(envelope.lrx(),
                         envelope.lry(),
                         0.0,
                         envelope.GetSpatialReference())

        xUlx, xUly, xLrx, xLry = tempEnv.expandByPercentage(10)

        # Mosaic SRTM tiles to cover this AoI.
        xEnv = Envelope()
        xEnv.addPoint(xUlx, xUly, 0.0, tempEnv.GetSpatialReference())
        xEnv.addPoint(xLrx, xLry, 0.0, tempEnv.GetSpatialReference())

        inputDem.mosaicAndClipDemTiles(demName, xEnv)

        return demName

    # -------------------------------------------------------------------------
    # _createStrip
    # _runOneStrip -> _createStrip
    # -------------------------------------------------------------------------
    @staticmethod
    def _createStrip(stripName,
                     stripScenes,
                     bandDir,
                     stripDir,
                     logger) -> list:

        if logger:

            logger.info('Extracting bands and mosaicking to strips' +
                        ' for {} ({} input scenes)'.
                        format(stripName, len(stripScenes)))

        bands = ['BAND_P'] if 'P1BS' in stripName else \
            stripScenes[0].bandNameList

        stripBandList = EvhrToA._scenesToStripFromBandList(stripName,
                                                           stripScenes,
                                                           bands,
                                                           bandDir,
                                                           stripDir,
                                                           logger)

        return stripBandList

    # -------------------------------------------------------------------------
    # getUtmSrs
    #
    # This method returns the UTM SRS definition in PROJ4 format.
    # -------------------------------------------------------------------------
    def _getUtmSrs(self, envelope: Envelope) -> str:

        # If it is already in UTM, just return the PROJ4 string.
        srs = envelope.GetSpatialReference()

        if srs.IsProjected() and 'UTM' in srs.GetAttrValue('PROJCS'):
            return envelope.GetSpatialReference().ExportToProj4()

        # If SRS is not 4326, convert coordinates
        targetSRS = SpatialReference()
        targetSRS.ImportFromEPSG(4326)

        if not srs.IsSame(targetSRS):

            targetSRS.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            xform = CoordinateTransformation(srs, targetSRS)
            envelope.Transform(xform)

        # Check if AOI is within UTM boundary
        if envelope.uly() >= 84.0 or envelope.lry() <= -80.0:

            msg = 'Cannot process request with AoI outside of (-80, 84) ' + \
                  'degrees latitude.  uly = ' + \
                  str(envelope.uly()) + ' lry = ' + \
                  str(envelope.lry())

            raise RuntimeError(msg)

        # Clip the UTM Shapefile for this bounding box.
        clipFile = tempfile.mkdtemp()

        UTM_FILE = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'UTM_Zone_Boundaries/UTM_Zone_Boundaries.shp')

        cmd = 'ogr2ogr' + \
              ' -clipsrc' + \
              ' ' + str(envelope.ulx()) + \
              ' ' + str(envelope.lry()) + \
              ' ' + str(envelope.lrx()) + \
              ' ' + str(envelope.uly()) + \
              ' -f "ESRI Shapefile"' + \
              ' -select "Zone_Hemi"' + \
              ' "' + clipFile + '"' + \
              ' "' + UTM_FILE + '"'

        SystemCommand(cmd, logger=self._logger, raiseException=True)

        # Read clipped shapefile
        driver = ogr.GetDriverByName('ESRI Shapefile')
        ds = driver.Open(clipFile, 0)
        layer = ds.GetLayer()

        # Is the envelope too large?
        if len(layer) > 3:

            raise RuntimeError('The area requested spans more than 3 UTM ' +
                               ' zones.')

        maxArea = 0

        for feature in layer:

            area = feature.GetGeometryRef().GetArea()

            if area > maxArea:

                maxArea = area
                zone, hemi = feature.GetField('Zone_Hemi').split(',')

        # Configure proj.4 string
        proj4 = '+proj=utm +zone={} +ellps=WGS84 ' \
                '+datum=WGS84 +units=m +no_defs'.format(zone)

        if hemi.upper() == 'S':
            proj4 += ' +south'

        # Remove temporary clipFile and its auxiliary files
        driver.DeleteDataSource(clipFile)

        return proj4

    # -------------------------------------------------------------------------
    # mergeBands
    #
    # stripToToA -> mergeBands
    # -------------------------------------------------------------------------
    @staticmethod
    def _mergeBands(bandFiles, outFileName, logger) -> None:

        if logger:
            logger.info('Merging bands into ' + str(outFileName))

        if len(bandFiles) == 0:

            logger.warn('No band files to merge.')
            return

        outFileNameNonCogTemp = outFileName.replace('.tif', '-noncog.tif')

        bandDs = gdal.Open(bandFiles[0])
        driver = gdal.GetDriverByName('GTiff')
        options = ['COMPRESS=LZW',
                   'BIGTIFF=YES']

        ds = driver.Create(outFileNameNonCogTemp,
                           bandDs.RasterXSize,
                           bandDs.RasterYSize,
                           len(bandFiles),
                           gdal.GDT_Int16,
                           options=options)

        try:

            ds.SetSpatialRef(bandDs.GetSpatialRef())
            ds.SetProjection(bandDs.GetProjectionRef())
            ds.SetGeoTransform(bandDs.GetGeoTransform())

            # Add the bands.
            outBandNum = 0

            for bandFile in bandFiles:

                # Create the output band.
                outBandNum += 1
                outBand = ds.GetRasterBand(outBandNum)
                outBand.SetNoDataValue(EvhrToA.NO_DATA_VALUE)

                baseName = os.path.basename(bandFile)
                nameNoExt = baseName.split('.')[0]
                bandName = '-'.join(nameNoExt.split('_')[-2:])
                outBand.SetDescription(bandName)

                # Write the input raster to the output band.
                outBand.WriteRaster(0,
                                    0,
                                    bandDs.RasterXSize,
                                    bandDs.RasterYSize,
                                    gdal.Open(bandFile).ReadRaster())

            # Flush data to disk
            ds = None

            # Translate geotiff to cloud-optimized-geotiff format
            EvhrUtils.createCloudOptimizedGeotiff(outFileName,
                                                  outFileNameNonCogTemp,
                                                  logger)

            # Delete the non-cog TOA file
            os.remove(outFileNameNonCogTemp)

        except Exception:
            os.remove(outFileNameNonCogTemp)

        # Delete the band files.
        for bandFile in bandFiles:
            os.remove(bandFile)

    # --------------------------------------------------------------------------
    # _orthoOne
    #
    # run -> processStrips -> runOneStrip -> stripToToa -> orthoOne
    # --------------------------------------------------------------------------
    @staticmethod
    def _orthoOne(bandFile, inputDem, orthoDir, demDir, outSrsProj4,
                  mapproject_threads, panResolution, logger):

        baseName = os.path.splitext(os.path.basename(bandFile.fileName()))[0]
        orthoFile = os.path.join(orthoDir, baseName + '-ortho.tif')
        orthoDg = None

        if not os.path.exists(orthoFile):

            if logger:

                logger.info('Orthorectifying ' + str(bandFile) + ' to ' +
                            orthoFile)

            try:

                clippedDEM = EvhrToA._createDemForOrthos(inputDem,
                                                         bandFile.envelope(),
                                                         demDir,
                                                         logger)

            except RuntimeError as e:

                raise RuntimeError(str(e) + ' Band file: ' + str(bandFile)) \
                    from e

            # Orthorectify.
            orthoFileTemp = orthoFile.replace('.tif', '-temp.tif')
            outRes = panResolution if bandFile.isPanchromatic() else 2

            try:
                if logger:
                    msg = 'Using: {} threads for mapproject'.format(
                        mapproject_threads)
                    logger.info(msg)

                cmd = EvhrToA.BASE_SP_CMD + \
                    'mapproject --nodata-value 0' + \
                    ' --threads={}'.format(mapproject_threads) + \
                    ' --num-processes=1 -t rpc' + \
                    ' --mpp={}'.format(outRes) + \
                    ' --t_srs "{}"'.format(outSrsProj4) + \
                    ' ' + clippedDEM + \
                    ' ' + bandFile.fileName() + \
                    ' ' + bandFile.xmlFileName + \
                    ' ' + orthoFileTemp

                SystemCommand(cmd, logger, True)

            except Exception as e:

                msg = 'Encountered exception executing mapproject: ' + \
                    '{}'.format(e)

                if logger:
                    logger.error(msg)

                if os.path.exists(orthoFile):
                    os.remove(orthoFile)

                raise RuntimeError(msg) from e

            try:
                # Convert NoData to settings value, set output type to Int16
                cmd = EvhrToA.BASE_SP_CMD + \
                    'image_calc -c "var_0" {} -d int16   \
                    --output-nodata-value {} -o {}'. \
                    format(orthoFileTemp, EvhrToA.NO_DATA_VALUE, orthoFile)

                SystemCommand(cmd, logger, True)

            except Exception as e:

                msg = 'Encountered exception executing image_calc: ' + \
                    '{}'.format(e)

                if logger:
                    logger.error(msg)

                if os.path.exists(orthoFile):
                    os.remove(orthoFile)

                raise RuntimeError(msg) from e

            if os.path.exists(orthoFileTemp):
                os.remove(orthoFileTemp)

            # Copy xml to accompany ortho file (needed for TOA)
            shutil.copy(bandFile.xmlFileName,
                        orthoFile.replace('.tif', '.xml'))

            orthoDg = DgFile(orthoFile)
            orthoDg.setBandName(bandFile.getBandName())

        else:

            try:
                orthoDg = DgFile(orthoFile)

            except Exception as e:

                msg = 'Encountered exception creating DgFile ' + \
                    'from {}: {}'.format(orthoFile, e)

                if logger:
                    logger.error(msg)

                os.remove(orthoFile)

                raise RuntimeError(msg) from e

        return orthoDg

    # -------------------------------------------------------------------------
    # processStrips
    #
    # run -> processStrips
    # -> runOneStrip
    # -------------------------------------------------------------------------
    def processStrips(self,
                      stripsWithScenes: dict,
                      bandDir,
                      stripDir,
                      orthoDir,
                      demDir,
                      toaDir,
                      outSrsProj4,
                      panResolution,
                      panSharpen,
                      inputDem,
                      logger) -> list:

        toaNames = []

        for key in iter(stripsWithScenes):

            toaName = EvhrToA._runOneStrip(key,
                                           stripsWithScenes[key],
                                           bandDir,
                                           stripDir,
                                           orthoDir,
                                           demDir,
                                           toaDir,
                                           outSrsProj4,
                                           panResolution,
                                           panSharpen,
                                           inputDem,
                                           logger)

            toaNames.append(toaName)

        return toaNames

    # -------------------------------------------------------------------------
    # removeDuplicates
    # -------------------------------------------------------------------------
    def _removeDuplicates(self, dgScenes: list) -> list:

        if len(dgScenes) == 1:
            return dgScenes

        # ---
        # Footprints might erroneously not be normalized.  Remove duplicates
        # now.  Dg_mosaic would find them later and throw an exception.
        # ---
        numBefore = len(dgScenes)
        dgScenes = list(set(dgScenes))
        numAfter = len(dgScenes)

        if self._logger:
            self._logger.info('Number of input scenes ' + str(numBefore))

            self._logger.info('Set operation removed ' +
                              str(numBefore - numAfter) +
                              ' scenes.')

        # ---
        # There are also cases where NTF files of the same name exist in two
        # different directories, the .../X1BS/... and the .../M1BS/...
        # directories.  Their NTF files differ, yet they have identical XML
        # files.  Remove the X1BS version.  First, sort the list to make
        # searching easier.
        # ---
        sortedScenes = dgScenes.copy()
        sortedScenes.sort(reverse=True)  # Put X1BS first, so M1BS is added.
        dgScenes = []
        numScenes = len(sortedScenes)
        secondToLastMatched = False

        for i in range(numScenes-1):

            # ---
            # If the base name of this item matches the base name of the
            # next item, check the XML.
            # ---
            ntf1 = sortedScenes[i].fileName()
            ntf2 = sortedScenes[i+1].fileName()

            if os.path.basename(ntf1) == os.path.basename(ntf2):

                xml1 = os.path.splitext(ntf1)[0] + '.xml'
                xml2 = os.path.splitext(ntf2)[0] + '.xml'

                # If the contents of the XML differs, keep them.
                if not filecmp.cmp(xml1, xml2, False):

                    dgScenes.append(sortedScenes[i])

                else:
                    if self._logger:
                        self._logger.info('XML match: ' + ntf1 + ' ' + ntf2)

                    if i == numScenes - 2:
                        secondToLastMatched = True
            else:

                dgScenes.append(sortedScenes[i])

        # ---
        # The loop above compares the current scene to the one after, so the
        # last scene in the list is never checked--or added to the filtered
        # list.  If the last item checked, which is the second-to-last item,
        # was a match, then the last item should not be added.  Otherwise,
        # add it.
        # ---
        if not secondToLastMatched:
            dgScenes.append(sortedScenes[-1])

        numAfter = len(dgScenes)

        if self._logger:
            self._logger.info(str(numAfter) +
                              ' after removal of matching XMLs.')

        return dgScenes

    # -------------------------------------------------------------------------
    # run
    # -> collectImagesByStrip
    # -> computeEnvelope
    # -> getUtmSrs
    # -> processStrips
    # -------------------------------------------------------------------------
    def run(self, inDgScenes: list = None) -> list:

        if self._logger:
            self._logger.info('In run')

        if not inDgScenes:
            raise RuntimeError('No valid scenes found.')

        dgScenes: list = self._removeDuplicates(inDgScenes)
        envelope: Envelope = self._computeEnvelope(list(dgScenes))
        self._outSrsProj4 = self._getUtmSrs(envelope)
        stripsWithDgScenes: dict = self._collectImagesByStrip(dgScenes)

        # ---
        # Process the strips.
        # ---
        toaNames = self.processStrips(stripsWithDgScenes,
                                      self._bandDir,
                                      self._stripDir,
                                      self._orthoDir,
                                      self._demDir,
                                      self._toaDir,
                                      self._outSrsProj4,
                                      self._panResolution,
                                      self._panSharpen,
                                      self._inputDem,
                                      self._logger)

        toaPaths = [Path(f) for f in toaNames]

        return toaPaths

    # -------------------------------------------------------------------------
    # runOneStrip
    #
    # run -> processStrips -> runOneStrip
    # -> createStrip
    # -> stripToToA
    # -------------------------------------------------------------------------
    @staticmethod
    def _runOneStrip(stripID, scenes, bandDir, stripDir, orthoDir, demDir,
                     toaDir, outSrsProj4, panResolution, panSharpen, inputDem,
                     logger, thisToaIsForPanSharpening=False):

        if logger:
            logger.info('In runOneStrip')

        toaName = os.path.join(toaDir, stripID + '-toa.tif')

        if thisToaIsForPanSharpening:
            toaName = toaName.replace('_M1BS_', '_P1BS_')

        if not os.path.exists(toaName):

            imageForEachBandInStrip: list = EvhrToA._createStrip(stripID,
                                                                 scenes,
                                                                 bandDir,
                                                                 stripDir,
                                                                 logger)

            mapproject_threads = 4

            # StripToToa detects an empty imageForEachBandInStrip, and stops.
            EvhrToA._stripToToa(imageForEachBandInStrip,
                                toaName,
                                inputDem,
                                orthoDir,
                                demDir,
                                toaDir,
                                outSrsProj4,
                                mapproject_threads,
                                panResolution,
                                logger)

        if panSharpen and \
            scenes[0].isMultispectral() and \
                os.path.exists(toaName):

            EvhrToA._runPanSharpening(toaName, stripID, scenes, bandDir,
                                      stripDir, orthoDir,
                                      demDir, toaDir, outSrsProj4,
                                      panResolution, inputDem, logger)

        return toaName

    # -------------------------------------------------------------------------
    # runPanSharpening
    # -------------------------------------------------------------------------
    @staticmethod
    def _runPanSharpening(toaName: str,
                          stripID: str,
                          scenes: list,
                          bandDir: str,
                          stripDir: str,
                          orthoDir: str,
                          demDir: str,
                          toaDir: str,
                          outSrsProj4: str,
                          panResolution: float,
                          inputDem: InputDem,
                          logger: logging.RootLogger):

        if logger:
            logger.info('In _runPanSharpening')

        # ---
        # Get the panchromatic mates of the ToA.  By definition, panchromatic
        # mates should have the same name as their multispectral counterparts
        # except the substring "M1BS" should be "P1BS".
        # ---
        panNames = [m.fileName().replace('-M1BS-', '-P1BS-') for m in scenes]
        panDgMates = []

        for panName in panNames:

            if not os.path.exists(panName):

                logger.warning('Panchromatic mate, ' +
                               str(panName) +
                               ', does not exist.  ' +
                               'Pansharpening will not be applied.')

            else:
                panDgMates.append(DgFile(panName))

        if not panDgMates:

            logger.warning(
                'There are not panchromatic scenes for catalog ID ' +
                str(stripID))

            return

        # ---
        # Make a ToA from the panchromatic mates.  If any of these scenes were
        # passed into the application from the user, they will not be
        # redundantly processed.  Each step checks for the existence of the
        # file it is supposed to create before running its code.
        #
        # This should look a lot like the regular EVHR process.
        #
        # Note, if we could predict the name of the pan. sharpened ToA, we
        # could skip its creation when it already exists.
        # ---
        if logger:
            logger.info('Creating pan. ToA for sharpening')

        panSharpeningToA = EvhrToA._runOneStrip(stripID,
                                                panDgMates,
                                                bandDir,
                                                stripDir,
                                                orthoDir,
                                                demDir,
                                                toaDir,
                                                outSrsProj4,
                                                panResolution,
                                                False,
                                                inputDem,
                                                logger,
                                                thisToaIsForPanSharpening=True)

        # ---
        # Warp the pan. ToA to match the multispectral ToA it will sharpen.
        # ---
        warpedPanSharpeningToA = \
            panSharpeningToA.replace('-toaForPanSharp.tif', '-warpedPsToa.tif')

        if logger:
            logger.info('Creating ' + str(warpedPanSharpeningToA))

        if not os.path.exists(warpedPanSharpeningToA):

            inputEnv = \
                GeospatialImageFile(toaName, None, None, logger).envelope()

            cmd = 'gdalwarp' + \
                  ' -dstnodata -10001 -co COMPRESS=LZW -co BIGTIFF=YES' + \
                  ' -te ' + str(inputEnv.ulx()) + ' ' + str(inputEnv.lry()) + \
                  ' ' + str(inputEnv.lrx()) + ' ' + str(inputEnv.uly()) + \
                  ' ' + panSharpeningToA + \
                  ' ' + warpedPanSharpeningToA

            SystemCommand(cmd, logger, True)

        # ---
        # Pan sharpen.
        # ---
        psTempName = toaName.replace('-toa.tif', '-toa-sharpened.tif')

        if logger:
            logger.info('Creating ' + str(psTempName))

        if not os.path.exists(psTempName):

            cmd = 'gdal_pansharpen.py' + \
                  ' -nodata -10001 -co COMPRESS=LZW -co BIGTIFF=YES' + \
                  ' ' + warpedPanSharpeningToA + \
                  ' ' + toaName + \
                  ' ' + psTempName

            SystemCommand(cmd, logger, True)

            # Copy an XML for it.
            toaXmlName = toaName.replace('.tif', '.xml')
            psTempXmlName = psTempName.replace('.tif', '.xml')
            shutil.copy(toaXmlName, psTempXmlName)

    # ------------------------------------------------------------------------
    # scenesToStripFromBandList
    #
    # Input: a list of scenes all belonging to the same strip and list of
    # bands
    #
    # Output:  a mosaic of all the scenes for each band:  a mosaic containing
    # band1 from every scene, a mosaic containing band2 from every scene ...
    #
    # run -> processStrips -> runOneStrip -> createStrip
    # ------------------------------------------------------------------------
    @staticmethod
    def _scenesToStripFromBandList(stripName, stripScenes, bands, bandDir,
                                   stripDir, logger) -> list:

        stripBandList = []  # Length of list = number of bands

        for bandName in bands:

            if logger:

                logger.info('Reading band ' + bandName +
                            ' for all scenes in strip ' + stripName)

            bandScenes = [scene.getBand(bandDir, bandName)
                          for scene in stripScenes]

            bandScenesStr = ' '.join(bandScenes)

            stripBandFile = \
                os.path.join(stripDir,
                             '{}_{}.r100.tif'.format(stripName, bandName))

            if not os.path.exists(stripBandFile):

                if logger:

                    logger.info('Mosaicking all of band ' + bandName +
                                ' scenes in strip ' + stripName)

                cmd = EvhrToA.BASE_SP_CMD + \
                    'dg_mosaic ' + \
                    '--output-nodata-value 0' + \
                    ' --ignore-inconsistencies --output-prefix {} {}'. \
                    format(stripBandFile.replace('.r100.tif', ''),
                           bandScenesStr)

                try:
                    SystemCommand(cmd, logger, raiseException=True)

                except Exception as e:

                    if 'Input images have same input location which ' + \
                       'could be caused by repeated XML file or invalid ' + \
                       'TLC information.' in str(e):

                        if logger:

                            msg = 'Skipping ' + stripName + ' due to ' + \
                                  str(e)

                            logger.warn(msg)

                        break

                else:

                    stripBandDg = DgFile(stripBandFile)
                    stripBandDg.setBandName(bandName)
                    stripBandList.append(stripBandDg)

            else:

                stripBandDg = DgFile(stripBandFile)
                stripBandDg.setBandName(bandName)
                stripBandList.append(stripBandDg)

        # Return the list of band strips.
        return stripBandList

    # -------------------------------------------------------------------------
    # _stripToToa
    #
    # run -> runOneStrip -> stripToToa
    # -------------------------------------------------------------------------
    @staticmethod
    def _stripToToa(imageForEachBandInStrip, toaName, inputDem, orthoDir,
                    demDir, toaDir, outSrsProj4, mapproject_threads,
                    panResolution, logger):

        if logger:
            logger.info('In _stripToToa, processing ' + toaName)

        if os.path.exists(toaName):
            return

        if len(imageForEachBandInStrip) == 0:

            logger.warn('No strip images.')
            return

        toaBands = []

        for stripBand in imageForEachBandInStrip:

            orthoBandDg = EvhrToA._orthoOne(stripBand,
                                            inputDem,
                                            orthoDir,
                                            demDir,
                                            outSrsProj4,
                                            mapproject_threads,
                                            panResolution,
                                            logger)

            toaBands.append(ToaCalculation.run(orthoBandDg,
                                               toaDir,
                                               logger))

        EvhrToA._mergeBands(toaBands, toaName, logger)
        shutil.copy(orthoBandDg.xmlFileName, toaName.replace('.tif', '.xml'))

    # -------------------------------------------------------------------------
    # toaDir
    # -------------------------------------------------------------------------
    @property
    def toaDir(self) -> Path:
        return Path(self._toaDir)
