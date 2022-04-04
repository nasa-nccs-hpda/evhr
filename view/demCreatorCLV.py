#!/usr/bin/python

import argparse
import logging
import pathlib
import sys

from osgeo import osr
from osgeo.osr import SpatialReference

from core.model.Envelope import Envelope
from evhr.model.DemCreator import DemCreator
from evhr.model.DemCreatorCelery import DemCreatorCelery
from evhr.model.ILProcessController import ILProcessController

# -----------------------------------------------------------------------------
# main
#
# evhr/view/demCreatorCLV.py -o /att/nobackup/rlgill/SystemTesting/testDEM/ --scenes  '/css/nga/WV02/1B/2018/278/WV02_10300100889D0300_X1BS_502602073050_01/WV02_20181005212447_10300100889D0300_18OCT05212447-P1BS-502602073050_01_P002.ntf'
#
# evhr/view/demCreatorCLV.py -o /att/nobackup/rlgill/SystemTesting/testDEM/ -e -148 65 -147.5 64.5 4326
# -----------------------------------------------------------------------------
def main():

    # Process command-line args.
    desc = 'Use this application to produce EVHR DEMs.'
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('--celery',
                        action='store_true',
                        help='Use Celery for distributed processing.')

    parser.add_argument('--log-to-file',
                        default=False,
                        help='Print messages to a file, instead of the screen.')

    parser.add_argument('-o',
                        default='.',
                        help='Path to output directory')

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('-e',
                       nargs=5,
                       help='ulx uly lrx lry epsg-code')

    group.add_argument('--scenes',
                       type=pathlib.Path,
                       nargs='*',
                       help='Fully-qualified path to scene files')

    args = parser.parse_args()

    # ---
    # Logging
    # ---
    logger = None
    
    if args.log-to-file:
        
        logFile = os.path.join(args.o, 'dem.out')
        logger = logger.basicConfig(logFile)
        
    else:
        logger = logging.getLogger()
    
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)

    # ---
    # Ensure an envelope or scenes were provided.
    # ---
    if args.e and args.scenes:

        logger.info('Both a scene list and an envelope were provided.  \
                    Using the scene list.')

        args.e = None

    # ---
    # Envelope
    # ---
    env = None

    if args.e:

        epsgCode = int(args.e[4])
        srs = SpatialReference()
        srs.ImportFromEPSG(epsgCode)

        srs4326 = SpatialReference()
        srs4326.ImportFromEPSG(4326)

        if srs.IsSame(srs4326):
            srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

        env = Envelope()
        env.addPoint(float(args.e[0]), float(args.e[1]), 0, srs)
        env.addPoint(float(args.e[2]), float(args.e[3]), 0, srs)

    # ---
    # Instantiate the process object and run.
    # ---
    if args.celery:

        with ILProcessController() as processController:

            dc = DemCreatorCelery(args.o, logger)
            
            if env:
                dc.runEnv(env)
            
            elif args.scenes:
                dc.runScenes(args.scenes)

            else:
                raise RuntimeError('Scenes or an envelope must be provided.')

    else:

        dc = DemCreator(args.o, logger)
        
        if env:
            dc.runEnv(env)
            
        elif args.scenes:
            dc.runScenes(args.scenes)

        else:
            raise RuntimeError('Scenes or an envelope must be provided.')

# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())



