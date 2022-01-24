
import glob
import os
import shutil
import tempfile
from xml.dom import minidom

from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from osgeo.osr import CoordinateTransformation
from osgeo.osr import SpatialReference

from core.model.BaseFile import BaseFile
from core.model.DgFile import DgFile
from core.model.Envelope import Envelope
from core.model.FootprintsQuery import FootprintsQuery
from core.model.SystemCommand import SystemCommand

from evhr.model.ToaCalculation import ToaCalculation


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
#                       _mosaicAndClipDemTiles
#                           _clipShp
#       _stripToToa
#           _ToaCalculation.run
#           _mergeBands
# -----------------------------------------------------------------------------
class EvhrToA(object):

    BASE_SP_CMD = '/opt/StereoPipeline/bin/'

    # BASE_SP_CMD = '/opt/StereoPipeline/' + \
    #               'StereoPipeline-2.7.0-2020-07-29-x86_64-Linux/bin/'

    NO_DATA_VALUE = -9999

    MAPPROJECT_THREADS = 4

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self, outDir, logger=None):

        self._logger = logger

        # Output directory
        self._outDir = BaseFile(outDir).fileName()  # BaseFile tests validity.

        if not os.path.isdir(self._outDir):
            raise RuntimeError(self._outDir + ' must be a directory.')

        self._outSrsProj4 = None

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

    # --------------------------------------------------------------------------
    # clipShp
    #
    # run -> runOneStrip -> stripToToa -> orthoOne -> createDemForOrthos
    #     -> mosaicAndClipDemTiles -> clipShp
    # --------------------------------------------------------------------------
    @staticmethod
    def _clipShp(shpFile, envelope, logger):

        if logger:
            logger.info('Clipping Shapefile.')

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
               ' "' + shpFile + '"'

        SystemCommand(cmd, logger, True)

        xml = minidom.parse(tempClipFile)
        features = xml.getElementsByTagName('gml:featureMember')

        return features

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
    def _collectImagesByStrip(self, sceneFiles):

        if self._logger:
            self._logger.info('In collectImagesByStrip')

        # Aggregate the scenes into strips.
        stripsWithScenes = {}

        for sceneFile in sceneFiles:

            try:
                dgf = DgFile(sceneFile, self._logger)

            except RuntimeError as e:

                if self._logger:
                    self._logger.warn(e)

                continue

            stripID = dgf.getStripName()

            if stripID:

                if stripID not in stripsWithScenes:
                    stripsWithScenes[stripID] = []

                if sceneFile not in stripsWithScenes[stripID]:
                    stripsWithScenes[stripID].append(dgf)

            else:

                if self._logger:

                    self._logger.warn('Unable to get strip name for: ' +
                                      stripID)

        return stripsWithScenes

    # -------------------------------------------------------------------------
    # _computeEnvelope
    # -------------------------------------------------------------------------
    def _computeEnvelope(self, sceneList):

        if self._logger:
            self._logger.info('In _computeEnvelope')

        env = sceneList[0].envelope()

        for dgf in sceneList[1:]:

            env.addPoint(dgf._ulx, dgf._uly, 0.0, dgf.srs())
            env.addPoint(dgf._lrx, dgf._lry, 0.0, dgf.srs())

        return env

    # --------------------------------------------------------------------------
    # _createDemForOrthos
    #
    # run -> runOneStrip -> stripToToa -> orthoOne -> createDemForOrthos
    # --------------------------------------------------------------------------
    @staticmethod
    def _createDemForOrthos(envelope, demDir, logger):

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
        EvhrToA._mosaicAndClipDemTiles(demName, xEnv, demDir, logger)

        return demName

    # -------------------------------------------------------------------------
    # _createStrip
    # _runOneStrip -> _createStrip
    # -------------------------------------------------------------------------
    @staticmethod
    def _createStrip(stripName, stripScenes, bandDir, stripDir, logger):

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
    def _getUtmSrs(self, envelope):

        # If it is already in UTM, just return the PROJ4 string.
        srs = envelope.GetSpatialReference()

        if srs.IsProjected() and 'UTM' in srs.GetAttrValue('PROJCS'):
            return self.GetSpatialRef().ExportToProj4()

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

    # --------------------------------------------------------------------------
    # mergeBands
    #
    # stripToToA -> mergeBands
    # --------------------------------------------------------------------------
    @staticmethod
    def _mergeBands(bandFiles, outFileName, logger):

        if logger:
            logger.info('Merging bands into ' + str(outFileName))

        # cmd = 'gdal_merge.py -co COMPRESS=LZW -co BIGTIFF=YES -ot Int16 \
        #       -separate -init {} -a_nodata {} -o {} {}'. \
        #       format(EvhrToA.NO_DATA_VALUE,
        #              EvhrToA.NO_DATA_VALUE,
        #              outFileName,
        #              ' '.join(bandFiles))
        #
        # SystemCommand(cmd, logger, True)

        bandDs = gdal.Open(bandFiles[0])
        driver = gdal.GetDriverByName('GTiff')

        ds = driver.Create(outFileName,
                           bandDs.RasterXSize,
                           bandDs.RasterYSize,
                           len(bandFiles),
                           gdal.GDT_Int16,
                           options=['COMPRESS=LZW',
                                    'BIGTIFF=YES'])

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

        except Exception:
            os.remove(outFileName)

        # Delete the band files.
        for bandFile in bandFiles:
            os.remove(bandFile)

    # --------------------------------------------------------------------------
    # mosaicAndClipDemTiles
    #
    # To build the SRTM index file:
    # gdaltindex -t_srs "EPSG:4326" -src_srs_name SRS srtm.shp /att/pubrepo/DEM/SRTM/1-ArcSec/*.hgt
    #
    # To build the ASTERGDEM index file:
    # gdaltindex -t_srs "EPSG:4326" -src_srs_name SRS astergdem.shp /att/pubrepo/DEM/ASTERGDEM/v2/*dem.tif
    #
    # run -> runOneStrip -> stripToToa -> orthoOne -> createDemForOrthos
    #     -> mosaicAndClipDemTiles
    # --------------------------------------------------------------------------
    @staticmethod
    def _mosaicAndClipDemTiles(outDemName, envelope, demDir, logger):

        if logger:
            logger.info('Creating DEM ' + str(outDemName))

        outDemNameTemp = outDemName.replace('.tif', '-temp.tif')

        # ---
        # SRTM was collected between -54 and 60 degrees of latitude.  Use
        # ASTERGDEM where SRTM is unavailable.
        # ---
        SHP_INDEX = None

        if envelope.uly() >= -54.0 and envelope.uly() <= 60.0 and \
           envelope.lry() >= -54.0 and envelope.lry() <= 60.0:

            SHP_INDEX = \
                os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             'SRTM/srtm.shp')

        else:

            SHP_INDEX = \
                os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             'ASTERGDEM/astergdem.shp')

        # Get the tile Shapefile and intersect it with the AoI.
        features = EvhrToA._clipShp(SHP_INDEX, envelope, logger)

        if not features or len(features) == 0:

            msg = 'Clipping rectangle to SRTM did not return any ' + \
                  'features.  Corners: (' + str(envelope.ulx()) + ', ' + \
                  str(envelope.uly()) + '), (' + \
                  str(envelope.lrx()) + ', ' + str(envelope.lry()) + ')'

            raise RuntimeError(msg)

        # Get the list of tiles.
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

        SystemCommand(cmd, logger, True)

        # Run mosaicked DEM through geoid correction
        cmd = EvhrToA.BASE_SP_CMD + \
            'dem_geoid ' + \
            outDemNameTemp + ' --geoid EGM96 -o ' + \
            outDemName.strip('-adj.tif') + \
            ' --reverse-adjustment'

        SystemCommand(cmd, logger, True)

        for log in glob.glob(os.path.join(demDir, '*log*.txt')):
            os.remove(log)  # remove dem_geoid log file

    # --------------------------------------------------------------------------
    # _orthoOne
    #
    # run -> processStrips -> runOneStrip -> stripToToa -> orthoOne
    # --------------------------------------------------------------------------
    @staticmethod
    def _orthoOne(bandFile, orthoDir, demDir, outSrsProj4, mapproject_threads,
                  logger):

        baseName = os.path.splitext(os.path.basename(bandFile.fileName()))[0]
        orthoFile = os.path.join(orthoDir, baseName + '-ortho.tif')
        orthoDg = None

        if not os.path.exists(orthoFile):

            if logger:

                logger.info('Orthorectifying ' + str(bandFile) + ' to ' +
                            orthoFile)

            try:

                clippedDEM = EvhrToA._createDemForOrthos(bandFile.envelope(),
                                                         demDir,
                                                         logger)

            except RuntimeError as e:

                raise RuntimeError(str(e) + ' Band file: ' + str(bandFile))

            # Orthorectify.
            orthoFileTemp = orthoFile.replace('.tif', '-temp.tif')
            outRes = 2

            if bandFile.isPanchromatic():
                outRes = 1

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
                raise RuntimeError(msg)

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
                raise RuntimeError(msg)

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
                raise RuntimeError(msg)

        return orthoDg

    # -------------------------------------------------------------------------
    # processStrips
    #
    # run -> processStrips
    # -> runOneStrip
    # -------------------------------------------------------------------------
    def processStrips(self, stripsWithScenes, bandDir, stripDir, orthoDir,
                      demDir, toaDir, outSrsProj4, logger):

        for key in iter(stripsWithScenes):

            EvhrToA._runOneStrip(key,
                                 stripsWithScenes[key],
                                 bandDir,
                                 stripDir,
                                 orthoDir,
                                 demDir,
                                 toaDir,
                                 outSrsProj4,
                                 logger)

    # -------------------------------------------------------------------------
    # _queryScenes
    # -------------------------------------------------------------------------
    def _queryScenes(self, envelope):

        fpq = FootprintsQuery(logger=self._logger)
        fpq.addAoI(envelope)
        fpq.setMinimumOverlapInDegrees()

        MAXIMUM_SCENES = 100
        fpq.setMaximumScenes(MAXIMUM_SCENES)

        sceneFiles = fpq.getScenes()

        if not sceneFiles and self._logger:
            self._logger.error('There were no level 1B scenes.')

        return sceneFiles

    # -------------------------------------------------------------------------
    # run
    # -> queryScenes
    # -> collectImagesByStrip
    # -> computeEnvelope
    # -> getUtmSrs
    # -> processStrips
    # -------------------------------------------------------------------------
    def run(self, envelope=None, sceneList=None):

        if self._logger:
            self._logger.info('In run')

        # ---
        # We need both an envelope and a scene list.  If there is no envelope,
        # expect a scene list and compute its envelope.  If there is no scene
        # list, expect an envelope and query for scenes within it.
        # ---
        if not envelope and not sceneList:

            raise RuntimeError('Either an envelope or a scene list ' +
                               'must be provided.')

        if not sceneList:
            sceneList = self._queryScenes(envelope)

        else:

            # ---
            # Convert Path objects to strings.  The Path class is new in
            # Python 3.2.  We should use them extensively.  There isn't time
            # to do that now, so cast them to strings.
            # ---
            sceneList = [str(scene) for scene in sceneList]

        # ---
        # Footprints might erroneously not be normalized.  Remove duplicates
        # now.  Dg_mosaic would find them later and throw an exception.
        # ---
        sceneList = list(set(sceneList))
        sceneList.sort()

        # ---
        # _collectImagesByStrip and _computeEnvelope both instantiate each
        # scene as a DgFile.  It will be efficient to reuse the DgFiles from
        # _collectImagesByStrip, although this method will be a little more
        # difficult to read.  Relative to the strip processing below, this
        # only takes a fraction of the overall time.
        # ---
        stripsWithScenes = self._collectImagesByStrip(sceneList)

        if not envelope:

            sceneSet = set()

            for key in stripsWithScenes.keys():
                sceneSet |= set(stripsWithScenes[key])

            envelope = self._computeEnvelope(list(sceneSet))

        # The output SRS must be UTM.
        self._outSrsProj4 = self._getUtmSrs(envelope)

        self.processStrips(stripsWithScenes,
                           self._bandDir,
                           self._stripDir,
                           self._orthoDir,
                           self._demDir,
                           self._toaDir,
                           self._outSrsProj4,
                           self._logger)

    # -------------------------------------------------------------------------
    # runOneStrip
    #
    # run -> processStrips -> runOneStrip
    # -> createStrip
    # -> stripToToA
    # -------------------------------------------------------------------------
    @staticmethod
    def _runOneStrip(stripID, scenes, bandDir, stripDir, orthoDir, demDir,
                     toaDir, outSrsProj4, logger):

        if logger:
            logger.info('In runOneStrip')

        imageForEachBandInStrip = EvhrToA._createStrip(stripID,
                                                       scenes,
                                                       bandDir,
                                                       stripDir,
                                                       logger)

        toaName = os.path.join(toaDir, stripID + '-toa.tif')

        EvhrToA._stripToToa(imageForEachBandInStrip,
                            toaName,
                            orthoDir,
                            demDir,
                            toaDir,
                            outSrsProj4,
                            EvhrToA.MAPPROJECT_THREADS,
                            logger)

    # --------------------------------------------------------------------------
    # scenesToStripFromBandList
    #
    # Input: a list of scenes all belonging to the same strip and list of bands
    #
    # Output:  a mosaic of all the scenes for each band:  a mosaic containing
    # band1 from every scene, a mosaic containing band2 from every scene ...
    #
    # run -> processStrips -> runOneStrip -> createStrip
    # --------------------------------------------------------------------------
    @staticmethod
    def _scenesToStripFromBandList(stripName, stripScenes, bands, bandDir,
                                   stripDir, logger):

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

                SystemCommand(cmd, logger, True)

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
    def _stripToToa(imageForEachBandInStrip, toaName, orthoDir, demDir, toaDir,
                    outSrsProj4, mapproject_threads, logger):

        if logger:
            logger.info('In _stripToToa, processing ' + toaName)

        if os.path.exists(toaName):
            return

        toaBands = []

        for stripBand in imageForEachBandInStrip:

            orthoBandDg = EvhrToA._orthoOne(stripBand,
                                            orthoDir,
                                            demDir,
                                            outSrsProj4,
                                            mapproject_threads,
                                            logger)

            toaBands.append(ToaCalculation.run(orthoBandDg,
                                               toaDir,
                                               logger))

        EvhrToA._mergeBands(toaBands, toaName, logger)
        shutil.copy(orthoBandDg.xmlFileName, toaName.replace('.tif', '.xml'))
