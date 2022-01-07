
import logging
import os
import re

from core.model.BaseFile import BaseFile
from core.model.DgFile import DgFile
from core.model.Envelope import Envelope
from core.model.FootprintsQuery import FootprintsQuery

# -----------------------------------------------------------------------------
# Class DemCreator
# -----------------------------------------------------------------------------
class DemCreator(object):
    
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
    # -------------------------------------------------------------------------
    def _collatePairs(self, scenes):
        
        # ---
        # Now that dg_stereo.sh does not query redundantly, EDR must copy each
        # pair's files to the request directory for dg_stereo.sh to find them.
        # The first step is to associate the pair name with its files.
        # ---
        pairs = {}

        for fps in scenes:
            
            pairName = fps.pairName()
            
            if not pairs.has_key(pairName):
                pairs[pairName] = []
                
            pairs[pairName].append(fps.fileName())
            
        return pairs
            
    # -------------------------------------------------------------------------
    # envelopeToSceneList
    #
    # If you do not have a list of scenes, only an Envelope, use this method
    # to get scenes from the Envelope.  Pass these scenes to run().
    # -------------------------------------------------------------------------
    def envelopeToSceneList(self, envelope):
        
        # Envelope
        if not isinstance(envelope, Envelope):
            
            raise TypeError('The envelope argument must be of type ' +
                            'core.model.Envelope.')
        
        if not envelope.IsValid():
            
            raise RuntimeError('Envelope is invalid.')
                             
        # The envelope must be in the geographic projection.
        srs = SpatialReference()
        srs.ImportFromEPSG('EPSG:4326')
        outEnv = envelope.TransformTo(srs)

        # Query for pairs of scenes.
        fpq = FootprintsQuery(logger=self.logger)
        
        fpq.addAoI(outEnv.ulx(),
                   outEnv.uly(),
                   outEnv.lrx(),
                   outEnv.lry(),
                   outEnv.GetSpatialReference())
        
        fpq.setMultispectralOff()
        fpq.setPairsOnly()
        fpScenes = fpq.getScenes()
        
        return fpScenes
        
    # -------------------------------------------------------------------------
    # findMates
    # -------------------------------------------------------------------------
    def _findMates(self, pairs):
        
        # Ensure that each pair has its mates.
        pairsWithMissingScenes = []
    
        for pairName in pairs.iterkeys():
        
            mates = pairName.split('_')[2:]
        
            for mate in mates:
            
                r = re.compile('.*' + mate + '.*')
            
                if not (filter(r.match, pairs[pairName])):

                    fpq = FootprintsQuery(logger=self.logger)
                    fpq.addCatalogID([mate])
                    mates = fpq.getScenes()
                
                    if mates:
                    
                        existingScenes = pairs[pairName]
                    
                        for fps in mates:
                        
                            if fps not in existingScenes:

                                pairs[pairName].append(fps.fileName())
                                fpScenes.append(fps)

                    else:
                    
                        if self._logger:

                            self._logger.warning('Pair ' +
                                                 pairName +
                                                 ' does not contain any' +
                                                 ' scenes for ' +
                                                 mate)

                        pairsWithMissingScenes.append(pairName)
                    
        return pairs, pairsWithMissingScenes
        
    # -------------------------------------------------------------------------
    # getPairs
    # -------------------------------------------------------------------------
    def _getPairs(self, dgScenes):
        
        # Collate scenes into pairs.
        pairs = self._collatePairs(dgScenes)
        
        # Search for missing mates.
        pairs, pairsWithMissingScenes = self._findMates(pairs)
        
        # Remove unpaired scenes.
        numUnpairedScenes = 0
        
        for pair in pairsWithMissingScenes:
            
            numUnpairedScenes += len(pairs[pair])
            del pairs[pair]
            
        # Reconcile all this bookkeeping.
        self._reconcilePairing(pairs, dgScenes)
        
    # -------------------------------------------------------------------------
    # reconcilePairing
    # -------------------------------------------------------------------------
    def _reconcilePairing(self, pairs, dgScenes):
        
        numPairedScenes = 0
        
        for pair in pairs:
            numPairedScenes += len(pairs[pair])
            
        if self._logger:
            
            numQueriedScenes = len(dgScenes)
            
            unaccountedScenes = \
                numQueriedScenes - numPairedScenes - numUnpairedScenes
            
            self._logger.info('Queried scenes: ' + \
                              str(numQueriedScenes) + '\n' + \
                              'Unpaired scenes: ' + \
                              str(numUnpairedScenes) + '\n' + \
                              'Paired scenes: ' + \
                              str(numPairedScenes) + '\n' + \
                              'Unaccounted scenes: ' + \
                              str(unaccountedScenes) + '\n' + \
                              'Pairs: ' + str(len(pairs)))
                             
            for pair in sorted(pairs.keys()):
                self._logger.info(pair)

    # -------------------------------------------------------------------------
    # run
    # -------------------------------------------------------------------------
    def run(self, dgScenes):
        
        pairs = self._getPairs(dgScenes)
        self.runPairs(pairs)

    # -------------------------------------------------------------------------
    # runPairs
    #
    # This decompsition of the run method is in anticipation of a Celery
    # version of DemCreator.
    # -------------------------------------------------------------------------
    def runPairs(self, pairs):
        
        for key in pairs:
            self._runPair(key, pairs[key], self._outDir)

    # -------------------------------------------------------------------------
    # runPair
    # -------------------------------------------------------------------------
    @staticmethod
    def _runPair(pairName, dgScenes, outDir):
        
        # If the DEM exists, do not proceed.
        demName = os.path.join(outDir, pairName + '.tif')
        
        if os.path.exists(demName):
            return  

        # Create the working directory, if necessary.
        workDir = os.path.join(outDir, pairName)
        
        if not os.path.exists(workDir):
            os.mkdir(workDir)

        # Copy the scenes to the working directory using sym links
        for scene in dgScenes:

            ext = os.path.splitext(scene)[1] # could be .tif or .ntf            
            dst = os.path.join(workDir, os.path.basename(scene))
            
            if not os.path.exists(dst):
                os.symlink(scene, dst)
                
            dstXml = dst.replace(ext, '.xml')
            
            if not os.path.exists(dstXml):
                os.symlink(scene.replace(ext, '.xml'), dstXml)

        # DEM application settings.
        PAIR_NAME     = pairName
        TEST          = 'false' #'true'
        ADAPT         = 'true'
        MAP           = 'false'
        RUN_PSTEREO   = 'true' 
        BATCH_NAME    = '"' + pairName + '"'
        SGM           = 'false'
        SUB_PIX_KNL   = '15'
        ERODE_MAX     = '24'
        COR_KNL_SIZE  = '21'
        MYSTERY1      = '300'
        OUT_DIR       = outDir
        QUERY         = 'false'
        CROP_WINDOW   = '"0 15000 5000 5000"'
        USE_NODE_LIST = 'true'
        NODES         = '/att/nobackup/mwooten3/EVHR_API/DgStereo/nodeList.txt'
        
        # Create the DEM.
        BASE_SP_CMD = '/opt/???/bin/'
        DEM_APPLICATION = 'dg_stereo.sh'
        