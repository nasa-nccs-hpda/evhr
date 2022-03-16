
from celery import group

from core.model.CeleryConfiguration import app
from evhr.model.DemCreator import DemCreator


# -----------------------------------------------------------------------------
# class DemCreatorCelery
# -----------------------------------------------------------------------------
class DemCreatorCelery(DemCreator):

    # -------------------------------------------------------------------------
    # __init__
    # -------------------------------------------------------------------------
    def __init__(self, outDir, logger=None):

        if logger:
            logger.info('In DemCreatorCelery.__init__')
            
        # Initialize the base class.
        super(DemCreatorCelery, self).__init__(outDir, logger)

    # -------------------------------------------------------------------------
    # processPairs
    #
    # run{env|scenes} -> getPairs -> processPairs
    # -------------------------------------------------------------------------
    def processPairs(self, pairs):
        
        if self._logger:
            self._logger.info('In DemCreatorCelery.processPairs')
            
        wpi = group(DemCreatorCelery._processPair.s(
                        key,
                        pairs[key],
                        self._outDir,
                        self._logger) for key in pairs)

        result = wpi.apply_async()
        result.get()    # Waits for wpi to finish.

        return result

    # -------------------------------------------------------------------------
    # processPair
    # -------------------------------------------------------------------------
    @staticmethod
    @app.task(serializer='pickle')
    def _processPair(pairName, dgScenes, outDir, logger):
        
        if logger:
            logger.info('In DemCreatorCelery._processPair')
            
        DemCreator._processPair(pairName, dgScenes, outDir, logger)
