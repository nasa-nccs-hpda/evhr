
import os

from celery import group
from celery.contrib import rdb

from core.model.CeleryConfiguration import app
from evhr.model.EvhrToA import EvhrToA


# -----------------------------------------------------------------------------
# class EvhrToaCelery
# -----------------------------------------------------------------------------
class EvhrToaCelery(EvhrToA):

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self, outDir, logger=None):

        # Initialize the base class.
        super(EvhrToaCelery, self).__init__(outDir, logger)

    # -------------------------------------------------------------------------
    # processStrips
    #
    # run -> processStrips
    # -------------------------------------------------------------------------
    def processStrips(self, stripsWithScenes, bandDir, stripDir, orthoDir,
                      demDir, toaDir, outSrsProj4, logger):

        wpi = group(EvhrToaCelery._runOneStrip.s(
                            key,
                            stripsWithScenes[key],
                            bandDir,
                            stripDir,
                            orthoDir,
                            demDir,
                            toaDir,
                            outSrsProj4,
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
                     toaDir, outSrsProj4, logger):

        if logger:
            logger.info('In runOneStrip')

        imageForEachBandInStrip = EvhrToA._createStrip(stripID,
                                                       scenes,
                                                       bandDir,
                                                       stripDir,
                                                       logger)

        toaName = os.path.join(toaDir, stripID + '-toa.tif')

        EvhrToA._stripToToa(imageForEachBandInStrip, toaName, orthoDir,
                            demDir, toaDir, outSrsProj4, logger)
