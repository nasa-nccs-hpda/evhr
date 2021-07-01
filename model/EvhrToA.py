
import glob
import os
import shutil
import tempfile
from xml.dom import minidom

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

    # BASE_SP_CMD = '/usr/local/bin/singularity run -B ' + \
    #               '/usr/local/bin/singularity ' + \
    #               '/att/nobackup/iluser/containers/' + \
    #               'ilab-stereo-pipeline-3.0.0-sandbox/ '

    BASE_SP_CMD = '/opt/StereoPipeline/bin/'
    NO_DATA_VALUE = -9999

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self, envelope, outDir, logger=None):

        self._logger = logger
        self._envelope = envelope

        # Output directory
        self._outDir = BaseFile(outDir).fileName()  # BaseFile tests validity.

        if not os.path.isdir(self._outDir):
            raise RuntimeError(self._outDir + ' must be a directory.')

        # The output SRS must be UTM.
        self._outSrsProj4 = self._getUtmSrs()

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

            self._logger.info('Envelope: ' + str(self._envelope))
            self._logger.info('Output directory: ' + self._outDir)
            self._logger.info('Output SRS:' + self._outSrsProj4)
            self._logger.info('Band directory: ' + self._bandDir)
            self._logger.info('Strip directory: ' + self._stripDir)
            self._logger.info('DEM directory: ' + self._demDir)
            self._logger.info('Ortho image directory: ' + self._orthoDir)
            self._logger.info('ToA directory: ' + self._toaDir)

    # --------------------------------------------------------------------------
    # clipShp
    # --------------------------------------------------------------------------
    def _clipShp(self, shpFile, envelope):

        if self._logger:
            self._logger.info('Clipping Shapefile.')

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

        sCmd = SystemCommand(cmd, self._logger, True)

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
    def _collectImagesByStrip(self):

        if self._logger:
            self._logger.info('In collectImagesByStrip')

        fpq = FootprintsQuery(logger=self._logger)
        fpq.addAoI(self._envelope)
        fpq.setMinimumOverlapInDegrees()

        MAXIMUM_SCENES = 100
        fpq.setMaximumScenes(MAXIMUM_SCENES)

        sceneFiles = fpq.getScenes()
        sceneFiles.sort()

        if not sceneFiles and self._logger:
            self._logger.error('There were no level 1B scenes.')

        # Aggregate the scenes into strips.
        stripsWithScenes = {}

        for sceneFile in sceneFiles:

            dgf = DgFile(sceneFile, self._logger)
            stripID = dgf.getStripName()

            if stripID:

                if stripID not in stripsWithScenes:
                    stripsWithScenes[stripID] = []

                if sceneFile not in stripsWithScenes[stripID]:
                    stripsWithScenes[stripID].append(sceneFile)

            else:

                if self._logger:

                    self._logger.warn('Unable to get strip name for: ' +
                                      stripID)

        return stripsWithScenes

    # --------------------------------------------------------------------------
    # _createDemForOrthos
    #
    # --------------------------------------------------------------------------
    def _createDemForOrthos(self, envelope):

        if self._logger:
            self._logger.info('Creating DEM for orthorectification.')

        # If there is already a clipped DEM for this bounding box, use it.
        demName = 'dem-' + \
            str(envelope.ulx()) + '-' + \
            str(envelope.uly()) + '-' + \
            str(envelope.lrx()) + '-' + \
            str(envelope.lry()) + '-' + \
            str(envelope.GetSpatialReference().GetAuthorityCode(None)) + \
            '-adj.tif'

        demName = os.path.join(self._demDir, demName)

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
        self._mosaicAndClipDemTiles(demName, xEnv)

        return demName

    # -------------------------------------------------------------------------
    # _createStrip
    # -------------------------------------------------------------------------
    def _createStrip(self, stripName, stripScenes):

        if self._logger:

            self._logger.info('Extracting bands and mosaicking to strips' +
                              ' for {} ({} input scenes)'.
                              format(stripName, len(stripScenes)))

        bands = ['BAND_P'] if 'P1BS' in stripName else \
                ['BAND_B', 'BAND_G', 'BAND_R', 'BAND_N']

        stripBandList = \
            self._scenesToStripFromBandList(stripName, stripScenes, bands)

        return stripBandList

    # -------------------------------------------------------------------------
    # getUtmSrs
    #
    # This method returns the UTM SRS definition in PROJ4 format.
    # -------------------------------------------------------------------------
    def _getUtmSrs(self):

        # If it is already in UTM, just return the PROJ4 string.
        srs = self._envelope.GetSpatialReference()

        if srs.IsProjected() and 'UTM' in srs.GetAttrValue('PROJCS'):
            return self.GetSpatialReference().ExportToProj4()

        # If SRS is not 4326, convert coordinates
        targetSRS = SpatialReference()
        targetSRS.ImportFromEPSG(4326)

        if not srs.IsSame(targetSRS):

            targetSRS.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            xform = CoordinateTransformation(srs, targetSRS)
            self._envelope.Transform(xform)

        # Check if AOI is within UTM boundary
        if self._envelope.uly() >= 84.0 or self._envelope.lry() <= -80.0:

            msg = 'Cannot process request with AoI outside of (-80, 84) ' + \
                  'degrees latitude.  uly = ' + \
                  str(self._envelope.uly()) + ' lry = ' + \
                  str(self._envelope.lry())

            raise RuntimeError(msg)

        # Clip the UTM Shapefile for this bounding box.
        clipFile = tempfile.mkdtemp()

        UTM_FILE = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'UTM_Zone_Boundaries/UTM_Zone_Boundaries.shp')

        cmd = 'ogr2ogr' + \
              ' -clipsrc' + \
              ' ' + str(self._envelope.ulx()) + \
              ' ' + str(self._envelope.lry()) + \
              ' ' + str(self._envelope.lrx()) + \
              ' ' + str(self._envelope.uly()) + \
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
        proj4 = '+proj=utm +zone={} +ellps=WGS84 ' + \
                '+datum=WGS84 +units=m +no_defs'.format(zone)

        if hemi.upper() == 'S':
            proj4 += ' +south'

        # Remove temporary clipFile and its auxiliary files
        driver.DeleteDataSource(clipFile)

        return proj4

    # --------------------------------------------------------------------------
    # mergeBands
    # --------------------------------------------------------------------------
    def _mergeBands(self, bandFiles, outFileName):

        if self._logger:
            self._logger.info('Merging bands into ' + str(outFileName))

        cmd = 'gdal_merge.py -co COMPRESS=LZW -co BIGTIFF=YES -ot Int16 \
              -separate -init {} -a_nodata {} -o {} {}'. \
              format(EvhrToA.NO_DATA_VALUE,
                     EvhrToA.NO_DATA_VALUE,
                     outFileName,
                     ' '.join(bandFiles))

        sCmd = SystemCommand(cmd, self._logger, True)

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
    # --------------------------------------------------------------------------
    def _mosaicAndClipDemTiles(self, outDemName, envelope):

        if self._logger:
            self._logger.info('Creating DEM ' + str(outDemName))

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
        features = self._clipShp(SHP_INDEX, envelope)

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

        sCmd = SystemCommand(cmd, self._logger, True)

        # Run mosaicked DEM through geoid correction
        cmd = EvhrToA.BASE_SP_CMD + \
            'dem_geoid ' + \
            outDemNameTemp + ' --geoid EGM96 -o ' + \
            outDemName.strip('-adj.tif') + \
            ' --reverse-adjustment'

        sCmd = SystemCommand(cmd, self._logger, True)

        for log in glob.glob(os.path.join(self._demDir, '*log*.txt')):
            os.remove(log)  # remove dem_geoid log file

    # --------------------------------------------------------------------------
    # _orthoOne
    # --------------------------------------------------------------------------
    def _orthoOne(self, bandFile, origDgFile):

        baseName = os.path.splitext(os.path.basename(bandFile))[0]
        orthoFile = os.path.join(self._orthoDir, baseName + '-ortho.tif')

        if not os.path.exists(orthoFile):

            if self._logger:

                self._logger.info('Orthorectifying ' +
                                  str(bandFile) +
                                  ' to ' +
                                  orthoFile)

            try:

                clippedDEM = self._createDemForOrthos(origDgFile.envelope())

            except RuntimeError as e:

                msg = str(e) + ' Band file: ' + str(bandFile) + \
                      ' DgFile: ' + str(origDgFile.fileName)

                raise RuntimeError(msg)

            # Orthorectify.
            orthoFileTemp = orthoFile.replace('.tif', '-temp.tif')
            bandName = DgFile(bandFile).getBandName()
            outRes = 2

            if origDgFile.isPanchromatic():
                outRes = 1

            cmd = EvhrToA.BASE_SP_CMD + \
                'mapproject --nodata-value 0' + \
                ' --threads=2 -t rpc' + \
                ' --mpp={}'.format(outRes) + \
                ' --t_srs "{}"'.format(self._outSrsProj4) + \
                ' ' + clippedDEM + \
                ' ' + bandFile + \
                ' ' + origDgFile.xmlFileName + \
                ' ' + orthoFileTemp

            sCmd = SystemCommand(cmd, self._logger, True)

            # Convert NoData to settings value, set output type to Int16
            cmd = EvhrToA.BASE_SP_CMD + \
                'image_calc -c "var_0" {} -d int16   \
                --output-nodata-value {} -o {}'. \
                format(orthoFileTemp, EvhrToA.NO_DATA_VALUE, orthoFile)

            sCmd = SystemCommand(cmd, self._logger, True)

            # Copy xml to accompany ortho file (needed for TOA)
            shutil.copy(origDgFile.xmlFileName,
                        orthoFile.replace('.tif', '.xml'))

            DgFile(orthoFile).setBandName(bandName)

        return orthoFile

    # -------------------------------------------------------------------------
    # run
    # -------------------------------------------------------------------------
    def run(self):

        if self._logger:
            self._logger.info('In run')

        stripsWithScenes = self._collectImagesByStrip()

        for key in iter(stripsWithScenes):
            self._runOneStrip(key, stripsWithScenes[key])

    # -------------------------------------------------------------------------
    # runOne
    # -------------------------------------------------------------------------
    def _runOneStrip(self, stripID, scenes):

        if self._logger:
            self._logger.info('In runOneStrip')

        imageForEachBandInStrip = self._createStrip(stripID, scenes)
        toaName = os.path.join(self._toaDir, stripID + '-toa.tif')
        self._stripToToa(imageForEachBandInStrip, toaName)

    # --------------------------------------------------------------------------
    # scenesToStripFromBandList
    #
    # Input: a list of scenes all belonging to the same strip and list of bands
    #
    # Output:  a mosaic of all the scenes for each band:  a mosaic containing
    # band1 from every scene, a mosaic containing band2 from every scene ...
    # --------------------------------------------------------------------------
    def _scenesToStripFromBandList(self, stripName, stripScenes, bands):

        stripBandList = []  # Length of list = number of bands

        for bandName in bands:

            bandScenes = [DgFile(scene).getBand(self._bandDir, bandName)
                          for scene in stripScenes]

            bandScenesStr = ' '.join(bandScenes)

            stripBandFile = \
                os.path.join(self._stripDir,
                             '{}_{}.r100.tif'.format(stripName, bandName))

            if not os.path.exists(stripBandFile):

                # /opt/StereoPipeline/bin/dg_mosaic
                cmd = EvhrToA.BASE_SP_CMD + \
                      'dg_mosaic ' + \
                      '--output-nodata-value 0' + \
                      ' --ignore-inconsistencies --output-prefix {} {}'. \
                      format(stripBandFile.replace('.r100.tif', ''),
                             bandScenesStr)

                sCmd = SystemCommand(cmd, self._logger, True)

            DgFile(stripBandFile).setBandName(bandName)
            stripBandList.append(stripBandFile)

        # Return the list of band strips.
        return stripBandList

    # -------------------------------------------------------------------------
    # _stripToToa
    # -------------------------------------------------------------------------
    def _stripToToa(self, imageForEachBandInStrip, toaName):

        if self._logger:
            self._logger.info('In _stripToToa, processing ' + toaName)

        if os.path.exists(toaName):
            return

        toaBands = []

        for stripBand in imageForEachBandInStrip:

            dgStrip = DgFile(stripBand)
            orthoBand = self._orthoOne(stripBand, dgStrip)

            toaBands.append(ToaCalculation.run(orthoBand,
                                               self._toaDir,
                                               self._logger))

        self._mergeBands(toaBands, toaName)

        shutil.copy(DgFile(orthoBand).xmlFileName,
                    toaName.replace('.tif', '.xml'))
