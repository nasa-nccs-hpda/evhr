
import os

from osgeo.osr import SpatialReference

from core.model.BaseFile import BaseFile
from core.model.DgFile import DgFile
from core.model.Envelope import Envelope
from core.model.FootprintsQuery import FootprintsQuery
from core.model.SystemCommand import SystemCommand


# -----------------------------------------------------------------------------
# Class DemCreator
# -----------------------------------------------------------------------------
class DemCreator(object):

    MAXIMUM_SCENES = 100

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self, outDir, logger):

        self._logger = logger

        # Output directory
        self._outDir = BaseFile(outDir).fileName()  # BaseFile tests validity.

        if not os.path.isdir(self._outDir):
            raise RuntimeError(self._outDir + ' must be a directory.')

    # -------------------------------------------------------------------------
    # collatePairs
    #
    # run{env|scenes} -> getPairs -> collatePairs
    # -------------------------------------------------------------------------
    def _collatePairs(self, scenes):

        if self._logger:
            self._logger.info('In collatePairs')

        # ---
        # Now that dg_stereo.sh does not query redundantly, EDR must copy each
        # pair's files to the request directory for dg_stereo.sh to find them.
        # The first step is to associate the pair name with its files.
        # ---
        pairs = {}

        for fps in scenes:

            pairName = fps.pairName()

            if not pairName:

                raise RuntimeError('Scene ' +
                                   str(fps) +
                                   ' is not a member of a pair.')

            if pairName not in pairs:
                pairs[pairName] = []

            pairs[pairName].append(fps.fileName())

        return pairs

    # -------------------------------------------------------------------------
    # createCloudOptimizedGeotiffs
    # -------------------------------------------------------------------------
    def _createCloudOptimizedGeotiffs(self, workDir):

        files = self._getFileList(workDir)

        for f in files:

            name, ext = os.path.splitext(f)

            cmd = 'gdal_translate ' + \
                  f + \
                  ' ' + name + '-CoG.tif' + \
                  ' -co TILED=YES' + \
                  ' -co COPY_SRC_OVERVIEWS=YES' + \
                  ' -co COMPRESS=LZW'

            SystemCommand(cmd, self._logger, True)

    # -------------------------------------------------------------------------
    # demComplete
    #
    # out-DEM_1m.tif, out-DEM_24m_hs_az315.tif, out-DEM_24m.tif
    # out-DEM_4m_hs_az315.tif, out-DEM_4m.tif,
    # WV02_20161214_103001005F18FD00_103001005FC39D00_ortho_4m.tif,
    # WV02_20161214_103001005F18FD00_103001005FC39D00_ortho.tif
    # -------------------------------------------------------------------------
    @staticmethod
    def demComplete(workDir, logger=None):

        # pairName = os.path.basename(workDir)
        #
        # # files = ['out-DEM_1m.tif', 'out-DEM_24m_hs_az315.tif',
        # #          'out-DEM_24m.tif', 'out-DEM_4m_hs_az315.tif',
        # #          'out-DEM_4m.tif',
        # #          pairName + '_ortho_4m.tif',
        # #          pairName + '_ortho.tif']
        #
        # files = [pairName + '-DEM_1m.tif',
        #          pairName + '-DEM_24m_hs_az315.tif',
        #          pairName + '-DEM_24m.tif',
        #          pairName + '-DEM_4m_hs_az315.tif',
        #          pairName + '-DEM_4m.tif',
        #          pairName + '_ortho_4m.tif',
        #          pairName + '_ortho.tif']
        #

        files = DemCreator._getFileList(workDir)

        for f in files:

            testFile = os.path.join(workDir, f)

            if not os.path.exists(testFile):

                if logger:
                    logger.warn('Output file does not exist: ' + testFile)

                return False

        return True

    # -------------------------------------------------------------------------
    # findMissing
    #
    # run{env|scenes} -> getPairs -> findMates -> findMissing
    # -------------------------------------------------------------------------
    def _findMissing(self, pairName, pairScenes):

        dgScenes = [DgFile(s) for s in pairScenes]
        hasMate = False
        missing = []
        pairCats = set(pairName.split('_')[2:])

        while len(dgScenes) > 0:

            dgScene = dgScenes.pop()

            for potentialMate in dgScenes:

                if dgScene.isMate(pairName, potentialMate):

                    hasMate = True
                    break

            if not hasMate:
                missing += list(pairCats - set([dgScene.getCatalogId()]))

        return list(set(missing))

    # -------------------------------------------------------------------------
    # findMates
    #
    # run{env|scenes} -> getPairs -> findMates
    # -------------------------------------------------------------------------
    def _findMates(self, pairs):

        if self._logger:
            self._logger.info('In findMates')

        # ---
        # Ensure that each pair has its mates.  Querying Footprints takes a
        # long time, so identify all missing mates and query once.  This means
        # we need to check for missing scenes a second time, to determine
        # what Footprints didn't find.
        # ---
        pairsWithMissingScenes = []

        for pairName in pairs.keys():

            missing = self._findMissing(pairName, pairs[pairName])
            fpq = self._getBaseQuery()
            fpq.addCatalogID(missing)
            scenes = fpq.getScenes()

            # scenes = \
            #     fpq.getScenesFromResultsFile('/att/nobackup/rlgill/query2')

            # ---
            # These results can contain references to strip scenes that
            # weren't in the user's set of requested scenes.  Only keep ones
            # that are counterparts of input scenes.
            # ---
            keepers = []

            for pairScene in pairs[pairName]:

                for scene in scenes:

                    try:

                        if DgFile(scene.fileName()).isMate(pairName,
                                                           DgFile(pairScene)):

                            keepers.append(scene.fileName())
                            break

                    except RuntimeError:

                        # ---
                        # Some FP images are missing their XML files.  Skip
                        # these and keep going.
                        # ---
                        pass

                    except ValueError:

                        # ---
                        # Some FP images have incorrect stereo information.
                        # Skip these and keep going.
                        # ---
                        pass

            pairs[pairName] = list(set(pairs[pairName] + keepers))

            # Are any missing after searching?
            missing = self._findMissing(pairName, pairs[pairName])

            if len(missing) > 0:
                pairsWithMissingScenes.append(pairName)

        return pairs, pairsWithMissingScenes

    # -------------------------------------------------------------------------
    # getBaseQuery
    # -------------------------------------------------------------------------
    def _getBaseQuery(self):

        fpq = FootprintsQuery(logger=self._logger)
        fpq.setMinimumOverlapInDegrees()
        fpq.setMultispectralOff()
        fpq.setMaximumScenes(DemCreator.MAXIMUM_SCENES)
        fpq.setPairsOnly()
        return fpq

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
    # getPairs
    #
    # runScenes -> getPairs
    # -------------------------------------------------------------------------
    def _getPairs(self, scenes):

        if self._logger:
            self._logger.info('In getPairs')

        # Collate scenes into pairs.
        pairs = self._collatePairs(scenes)

        # Search for missing mates.
        pairs, pairsWithMissingScenes = self._findMates(pairs)

        # Remove unpaired scenes.
        numUnpairedScenes = 0

        for pair in pairsWithMissingScenes:

            numUnpairedScenes += len(pairs[pair])
            del pairs[pair]

        # Reconcile all this bookkeeping.
        self._reconcilePairing(pairs, scenes, numUnpairedScenes)

        return pairs

    # -------------------------------------------------------------------------
    # processPairs
    #
    # This decompsition of the run method is in anticipation of a Celery
    # version of DemCreator.
    #
    # run{env|scenes} -> getPairs -> processPairs
    # -------------------------------------------------------------------------
    def processPairs(self, pairs):

        if self._logger:
            self._logger.info('In processPairs')

        for key in pairs:
            self._processPair(key, pairs[key], self._outDir, self._logger)

    # -------------------------------------------------------------------------
    # processPair
    # -------------------------------------------------------------------------
    @staticmethod
    def _processPair(pairName, dgScenes, outDir, logger):

        if logger:
            logger.info('In _processPair')

        # Create the working directory, if necessary.
        workDir = os.path.join(outDir, pairName)

        if not os.path.exists(workDir):
            os.mkdir(workDir)

        # If the DEM exists, do not proceed.
        if DemCreator.demComplete(workDir):
            return

        # Copy the scenes to the working directory using sym links
        for scene in dgScenes:

            ext = os.path.splitext(scene)[1]  # could be .tif or .ntf
            dst = os.path.join(workDir, os.path.basename(scene))

            if not os.path.exists(dst):
                os.symlink(scene, dst)

            dstXml = dst.replace(ext, '.xml')

            if not os.path.exists(dstXml):
                os.symlink(scene.replace(ext, '.xml'), dstXml)

        # DEM application settings.
        DG_STEREO_DIR = '/opt/DgStereo'
        DEM_APPLICATION = os.path.join(DG_STEREO_DIR, 'evhr', 'dg_stereo.sh')
        PAIR_NAME = pairName
        TEST = 'false'  # 'true'
        ADAPT = 'true'
        MAP = 'false'
        RUN_PSTEREO = 'true'
        BATCH_NAME = '"' + pairName + '"'
        SGM = 'false'
        SUB_PIX_KNL = '15'
        ERODE_MAX = '24'
        COR_KNL_SIZE = '21'
        MYSTERY1 = '300'
        OUT_DIR = outDir
        QUERY = 'false'
        CROP_WINDOW = '"0 15000 5000 5000"'
        USE_NODE_LIST = 'false'
        NODES = os.path.join(DG_STEREO_DIR, 'nodeList.txt')

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
            ' ' + MYSTERY1 + \
            ' ' + OUT_DIR + \
            ' ' + QUERY + \
            ' ' + CROP_WINDOW

        SystemCommand(cmd, logger, True)

        # ---
        # dg_stereo.sh leaves many files in its wake.  Presumably, it knows
        # when it has finished without a problem.  However, there have been
        # cases to the contrary.  At least until this is solved, use this
        # definitive, independent verification of success.
        # ---
        if not DemCreator.demComplete(workDir, logger):

            if logger:
                logger.warn('DEM did not complete: ' + workDir)

        self._createCloudOptimizedGeotiffs(workDir)

    # -------------------------------------------------------------------------
    # reconcilePairing
    #
    # run{env|scenes} -> getPairs -> reconcilePairing
    # -------------------------------------------------------------------------
    def _reconcilePairing(self, pairs, scenes, numUnpairedScenes):

        if self._logger:
            self._logger.info('In reconcilePairing')

        numPairedScenes = 0

        for pair in pairs:
            numPairedScenes += len(pairs[pair])

        if self._logger:

            numQueriedScenes = len(scenes)

            self._logger.info('Queried scenes: ' +
                              str(numQueriedScenes) + '\n' +
                              'Unpaired scenes: ' +
                              str(numUnpairedScenes) + '\n' +
                              'Paired scenes: ' +
                              str(numPairedScenes) + '\n' +
                              'Pairs: ' + str(len(pairs)))

            for pair in sorted(pairs.keys()):
                self._logger.info(pair)

    # -------------------------------------------------------------------------
    # runEnv
    # -------------------------------------------------------------------------
    def runEnv(self, envelope):

        if self._logger:
            self._logger.info('In runEnv')

        # Envelope
        if not isinstance(envelope, Envelope):

            raise TypeError('The envelope argument must be of type ' +
                            'core.model.Envelope.')

        if not envelope.IsValid():

            raise RuntimeError('Envelope is invalid.')

        # The envelope must be in the geographic projection.
        srs = SpatialReference()
        srs.ImportFromEPSG(4326)
        envelope.TransformTo(srs)

        # Query for pairs of scenes.
        fpq = self._getBaseQuery()
        fpq.addAoI(envelope)
        fpq.setPairsOnly()
        fpScenes = fpq.getScenes()
        pairs = self._getPairs(fpScenes)
        self.processPairs(pairs)

        # import pandas as pd
        # savedPairs = \
        #     pd.read_pickle('/adapt/nobackup/people/rlgill/pairs.pickle')

        self.processPairs(savedPairs)

    # -------------------------------------------------------------------------
    # runScenes
    # -------------------------------------------------------------------------
    def runScenes(self, scenes):

        if self._logger:
            self._logger.info('In runScenes')

        # Convert Posix paths to strings.
        scenes = [str(s) for s in scenes]

        # Do not accept multispectral scenes.
        for scene in scenes:
            if DgFile(scene).isMultispectral():
                raise RuntimeError('Scenes must be panchromatic.')

        # Get FootprintsScene objects because they have pair information.
        fpq = self._getBaseQuery()
        fpq.addScenesFromNtf(scenes)
        fpScenes = fpq.getScenes()

        # fpScenes = \
        #     fpq.getScenesFromResultsFile('/att/nobackup/rlgill/query1')

        pairs = self._getPairs(fpScenes)
        self.processPairs(pairs)
