import logging
import pathlib
from typing import Union

from osgeo import gdal


class EvhrUtils(object):

    @staticmethod
    def createCloudOptimizedGeotiff(destPath: Union[str, pathlib.Path],
                                    srcPath: Union[str, pathlib.Path],
                                    logger: logging.Logger = None) -> None:
        translate_options = {
            'format': 'COG',
            'creationOptions': ['BIGTIFF=YES',
                                'COMPRESS=LZW'],
        }

        if logger:
            logger.info(f'Translating {srcPath} to COG {destPath}')

        options = gdal.TranslateOptions(**translate_options)

        ds = gdal.Translate(str(destPath), str(srcPath), options=options)

        ds = None
