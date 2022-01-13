#!/usr/bin/python

import argparse
import logging
import pathlib
import sys

from osgeo.osr import SpatialReference

from core.model.Envelope import Envelope
from evhr.model.DemCreator import DemCreator
from evhr.model.ILProcessController import ILProcessController

# -----------------------------------------------------------------------------
# main
#
# evhr/view/demCreatorCLV.py -o /att/nobackup/rlgill/SystemTesting/testToA/ --scenes /css/nga/WV03/1B/2015/219/WV03_104001000F2D9E00_X1BS_500495393030_01/WV03_20150807213524_104001000F2D9E00_15AUG07213524-M1BS-500495393030_01_P001.ntf
# -----------------------------------------------------------------------------
def main():

    # Process command-line args.
    desc = 'Use this application to produce EVHR DEMs.'
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('--celery',
                        action='store_true',
                        help='Use Celery for distributed processing.')

    parser.add_argument('-e',
                        nargs='*',
                        help='ulx uly lrx lry')

    parser.add_argument('--epsg',
                        type=int,
                        help='EPSG code')

    parser.add_argument('-o',
                        default='.',
                        help='Path to output directory')

    parser.add_argument('--scenes',
                        type=pathlib.Path,
                        nargs='*',
                        help='Fully-qualified path to scene files')

    args = parser.parse_args()

    # Envelope
    if args.epsg:

        srs4326 = SpatialReference()
        srs4326.ImportFromEPSG(4326)
        srs = SpatialReference()
        srs.ImportFromEPSG(args.epsg)

        if srs.IsSame(srs4326):
            srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    env = None

    if args.e:

        env = Envelope()
        env.addPoint(float(args.e[0]), float(args.e[1]), 0, srs)
        env.addPoint(float(args.e[2]), float(args.e[3]), 0, srs)

    # Logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)

    if args.celery:

        with ILProcessController() as processController:
            raise NotImplementedError()
            
    else:

        dc = DemCreator(args.o, logger)
        dc.run(args.scenes)


# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())



