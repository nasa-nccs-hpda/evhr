import os
import pathlib
import sys

import unittest
from unittest.mock import Mock

sys.modules['core.model.Envelope'] = Mock()
sys.modules['core.model.SystemCommand'] = Mock()

from evhr.model.InputDem import InputDem


# ----------------------------------------------------------------------------
# InputDemTestCase
# ----------------------------------------------------------------------------
class InputDemTestCase(unittest.TestCase):

    CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))

    DEM_EXISTS = os.path.join(CURRENT_DIR,
                              '../ASTERGDEM/astergdem.shp')
    DEM_EXISTS = pathlib.Path(DEM_EXISTS)

    DEM_NOT_EXISTS = pathlib.Path('DEM_NOT_EXIST.shp')

    DEM_DIR_EXISTS = pathlib.Path('.')

    DEM_DIR_NOT_EXISTS = pathlib.Path('.DEM_DIR_NOT_EXIST')

    # ------------------------------------------------------------------------
    # test_init
    # ------------------------------------------------------------------------
    def test_init(self):

        with self.assertRaises(FileNotFoundError):

            testDem = InputDem(self.DEM_NOT_EXISTS, self.DEM_DIR_EXISTS, None)
            testDem._validateDem(self.DEM_NOT_EXISTS)

        with self.assertRaisesRegex(FileNotFoundError, 'is expected to be'):

            InputDem(self.DEM_EXISTS, self.DEM_DIR_NOT_EXISTS, None)

        try:

            testDem = InputDem(self.DEM_EXISTS, self.DEM_DIR_EXISTS, None)
            testDem._validateDem(self.DEM_EXISTS)

        except Exception as e:

            print('Unittest not being run on ADAPT')
            print(e)
