
from celery import group
from celery.utils.log import get_task_logger

from evhr.model.CeleryConfiguration import app
from evhr.model.DemCreator import DemCreator

logger = get_task_logger(__name__)


# -----------------------------------------------------------------------------
# class DemCreatorCelery
# -----------------------------------------------------------------------------
class DemCreatorCelery(DemCreator):

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self, outDir, unusedLogger=None, testMode=False):

        if logger:
            logger.info('In DemCreatorCelery.__init__')

        # Initialize the base class.
        super(DemCreatorCelery, self).__init__(outDir, unusedLogger, testMode)

    # -------------------------------------------------------------------------
    # processPairs
    #
    # runPairs -> getPairs -> processPairs
    # -------------------------------------------------------------------------
    def processPairs(self, pairs):

        if logger:
            logger.info('In DemCreatorCelery.processPairs')

        wpi = group(DemCreatorCelery._processPair.s(
                        key,
                        pairs[key],
                        self._outDir,
                        self._testMode,
                        logger) for key in pairs)

        result = wpi.apply_async()
        result.get()    # Waits for wpi to finish.

        return result

    # -------------------------------------------------------------------------
    # processPair
    # -------------------------------------------------------------------------
    @staticmethod
    @app.task(serializer='pickle')
    def _processPair(pairName, dgScenes, outDir, testMode, unusedLogger):

        if logger:
            logger.info('In DemCreatorCelery._processPair')

        DemCreator._processPair(pairName, dgScenes, outDir, testMode, logger)
