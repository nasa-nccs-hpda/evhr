
import os

from celery import group
from celery.utils.log import get_task_logger

from core.model.CeleryConfiguration import app
from evhr.model.EvhrToA import EvhrToA

logger = get_task_logger(__name__)


# -----------------------------------------------------------------------------
# class EvhrToaCelery
# -----------------------------------------------------------------------------
class EvhrToaCelery(EvhrToA):

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self, outDir, panResolution=1, unusedLogger=None):

        logger.info('In EvhrToaCelery.__init__')

        # Initialize the base class.
        super(EvhrToaCelery, self).__init__(outDir, 
                                            panResolution, 
                                            unusedLogger)

    # -------------------------------------------------------------------------
    # processStrips
    #
    # run -> processStrips
    # -------------------------------------------------------------------------
    def processStrips(self, stripsWithScenes, bandDir, stripDir, orthoDir,
                      demDir, toaDir, outSrsProj4, panResolution, 
                      unusedLogger):

        logger.info('In EvhrToaCelery.processStrips')

        wpi = group(EvhrToaCelery._runOneStrip.s(
            key,
            stripsWithScenes[key],
            bandDir,
            stripDir,
            orthoDir,
            demDir,
            toaDir,
            outSrsProj4,
            panResolution,
            logger) for key in iter(stripsWithScenes))

        result = wpi.apply_async()
        result.get()    # Waits for wpi to finish.

        return result

    # -------------------------------------------------------------------------
    # runOneStrip
    #
    # run -> processStrips -> runOneStrip
    # -------------------------------------------------------------------------
    @staticmethod
    @app.task(serializer='pickle')
    def _runOneStrip(stripID, scenes, bandDir, stripDir, orthoDir, demDir,
                     toaDir, outSrsProj4, panResolution, unusedLogger):

        logger.info('In EvhrToaCelery._runOneStrip')

        imageForEachBandInStrip = EvhrToA._createStrip(stripID,
                                                       scenes,
                                                       bandDir,
                                                       stripDir,
                                                       logger)

        toaName = os.path.join(toaDir, stripID + '-toa.tif')
        mapproject_threads = 1

        EvhrToA._stripToToa(imageForEachBandInStrip,
                            toaName,
                            orthoDir,
                            demDir,
                            toaDir,
                            outSrsProj4,
                            mapproject_threads,
                            panResolution,
                            logger)
