
import os

from core.model.BaseFile import BaseFile
from core.model.SystemCommand import SystemCommand

from evhr.model.EvhrUtils import EvhrUtils


# -----------------------------------------------------------------------------
# Class DemCreator
# -----------------------------------------------------------------------------
class DemCreator(object):

    MAXIMUM_SCENES = 100

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self, outDir, logger, testMode=False):

        self._logger = logger
        self._testMode = testMode

        # Output directory
        self._outDir = BaseFile(outDir).fileName()  # BaseFile tests validity.

        if not os.path.isdir(self._outDir):
            raise RuntimeError(self._outDir + ' must be a directory.')

    # -------------------------------------------------------------------------
    # createCloudOptimizedGeotiffs
    # -------------------------------------------------------------------------
    @staticmethod
    def _createCloudOptimizedGeotiffs(workDir, logger):

        if logger:
            logger.info('In createCloudOptimizedGeotiffs')

        files = DemCreator._getFileList(workDir)
        cogs = []

        for f in files:

            name, _ = os.path.splitext(f)
            inFile = os.path.join(workDir, f)
            inFileNonCog = os.path.join(workDir, name + 'noncog.tif')

            # Rename inFile to reflect it is not the final CoG format raster.
            os.rename(inFile, inFileNonCog)

            # This is the final CoG file name
            cogName = os.path.join(workDir, f)

            if os.path.exists(cogName):
                cogs.append(cogName)
                os.remove(inFileNonCog)
                logger.info(f'{cogName} already exists.')
                continue

            logger.info(f'Translating {inFileNonCog} to CoG {cogName}')

            EvhrUtils.createCloudOptimizedGeotiff(cogName,
                                                  inFileNonCog,
                                                  logger)

            if not os.path.exists(cogName):
                errorMessage = f'CoG {cogName} was not created'
                logger.error(errorMessage)
                raise RuntimeError(errorMessage)

            cogs.append(cogName)
            os.remove(inFileNonCog)

        return cogs

    # -------------------------------------------------------------------------
    # demComplete
    # -------------------------------------------------------------------------
    @staticmethod
    def demComplete(workDir, logger=None):

        files = DemCreator._getFileList(workDir)

        print(files)
        for f in files:

            testFile = os.path.join(workDir, f)

            if not os.path.exists(testFile):

                if logger:
                    logger.warn('Output file does not exist: ' + testFile)

                return False

        return True

    # -------------------------------------------------------------------------
    # getFileList
    # -------------------------------------------------------------------------
    @staticmethod
    def _getFileList(workDir):

        pairName = os.path.basename(workDir)
        prefix = 'out'

        files = [prefix + '-DEM_1m.tif',
                 prefix + '-DEM_24m_hs_az315.tif',
                 prefix + '-DEM_24m.tif',
                 prefix + '-DEM_4m_hs_az315.tif',
                 prefix + '-DEM_4m.tif',
                 pairName + '_ortho_4m.tif',
                 pairName + '_ortho.tif']

        return files

    # -------------------------------------------------------------------------
    # processPairs
    #
    # This decompsition of the run method is in anticipation of a Celery
    # version of DemCreator.
    #
    # runPairs -> getPairs -> processPairs
    # -------------------------------------------------------------------------
    def processPairs(self, pairs):

        if self._logger:
            self._logger.info('In processPairs')

        for key in pairs:
            self._processPair(key, pairs[key],
                              self._outDir,
                              self._testMode,
                              self._logger)

    # -------------------------------------------------------------------------
    # processPair
    # -------------------------------------------------------------------------
    @staticmethod
    def _processPair(pairName, scenes, outDir, testMode, logger):

        if logger:
            logger.info('In _processPair')

        # Create the working directory, if necessary.
        workDir = os.path.join(outDir, pairName)

        if logger:
            logger.debug(f'Working directory: {workDir}')

        if not os.path.exists(workDir):
            os.mkdir(workDir)

        # If the DEM exists, do not proceed to make a new DEM.
        if not DemCreator.demComplete(workDir):

            # Copy the scenes to the working directory using sym links
            for scene in scenes:

                ext = os.path.splitext(scene)[1]  # could be .tif or .ntf
                dst = os.path.join(workDir, os.path.basename(scene))

                if not os.path.exists(dst):
                    os.symlink(scene, dst)

                dstXml = dst.replace(ext, '.xml')

                if not os.path.exists(dstXml):
                    os.symlink(scene.replace(ext, '.xml'), dstXml)

            try:
                DemCreator._runDgStereo(pairName,
                                        scenes,
                                        outDir,
                                        logger,
                                        testMode)

            except Exception as e:

                if not testMode:
                    raise RuntimeError(e)

            # ---
            # dg_stereo.sh leaves many files in its wake.  Presumably, it knows
            # when it has finished without a problem.  However, there have been
            # cases to the contrary.  At least until this is solved, use this
            # definitive, independent verification of success.
            # ---
            if not DemCreator.demComplete(workDir, logger):

                msg = 'DEM did not complete: ' + str(workDir)

                if testMode:

                    if logger:
                        logger.warn(msg)

                else:
                    raise RuntimeError(msg)

        # ---
        # Regardless of previous existing DEM, continue with COG generation
        # COG generation is now default behaviour
        # ---
        DemCreator._createCloudOptimizedGeotiffs(workDir, logger)

        if logger:
            logger.info('DEM completed in: ' + workDir)

    # -------------------------------------------------------------------------
    # runDgStereo
    # -------------------------------------------------------------------------
    @staticmethod
    def _runDgStereo(pairName, scenes, outDir, logger, testMode):

        DG_STEREO_DIR = '/opt/DgStereo'
        PAIR_NAME = pairName
        TEST = 'true'  # Keeps intermediate files
        ADAPT = 'true'
        MAP = 'false'
        RUN_PSTEREO = 'true'
        BATCH_NAME = '"' + pairName + '"'
        SGM = 'false'
        SUB_PIX_KNL = '15'
        ERODE_MAX = '24'
        COR_KNL_SIZE = '21'
        COR_TIME = '300'
        OUT_DIR = outDir
        QUERY = 'false'
        USE_NODE_LIST = 'false'
        NODES = os.path.join(DG_STEREO_DIR, 'nodeList.txt')
        # CROP_WINDOW = '"2000 2000 5000 5000"'  # xoff yoff xsize ysize
        CROP_WINDOW = '"0 15000 5000 5000"'  # xoff yoff xsize ysize

        DEM_APPLICATION = os.path.join(DG_STEREO_DIR,
                                       'evhr',
                                       'dg_stereo.sh')

        # Create the DEM.
        cmd = DEM_APPLICATION + \
            ' ' + PAIR_NAME + \
            ' ' + TEST + \
            ' ' + ADAPT + \
            ' ' + MAP + \
            ' ' + RUN_PSTEREO + \
            ' ' + BATCH_NAME + \
            ' _placeholder_for_rpcdem_' + \
            ' ' + USE_NODE_LIST + \
            ' ' + NODES + \
            ' ' + SGM + \
            ' ' + SUB_PIX_KNL + \
            ' ' + ERODE_MAX + \
            ' ' + COR_KNL_SIZE + \
            ' ' + COR_TIME + \
            ' ' + OUT_DIR + \
            ' ' + QUERY

        if testMode:
            cmd += ' ' + CROP_WINDOW

        SystemCommand(cmd, logger, True)

    # -------------------------------------------------------------------------
    # runPairs
    # -------------------------------------------------------------------------
    def runPairs(self, pairs):

        if self._logger:
            self._logger.info('In runPairs')
            self._logger.debug(f'Pairs: {pairs}')

        self.processPairs(pairs)
