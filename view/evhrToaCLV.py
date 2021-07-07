#!/usr/bin/python
import argparse
import logging
import sys

from osgeo import gdal
from osgeo import osr
from osgeo.osr import SpatialReference

from core.model.Envelope import Envelope
from evhr.model.EvhrToA import EvhrToA


# -----------------------------------------------------------------------------
# main
#
# export PYTHONPATH=`pwd`:`pwd`/core:`pwd`/evhr
# evhr/view/evhrToaCLV.py -e -148 65 -147.5 64.5 --epsg 4326 -o /att/nobackup/rlgill/SystemTesting/testToA/
# -----------------------------------------------------------------------------
def main():

    # Process command-line args.
    desc = 'Use this application to produce top of atmosphere images.'
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-e',
                        required=True,
                        nargs='+',
                        help='ulx uly lrx lry')

    parser.add_argument('--epsg',
                        required=True,
                        type=int,
                        help='EPSG code')

    parser.add_argument('-o',
                        default='.',
                        help='Path to output directory')

    args = parser.parse_args()

    print ('GDAL version:', gdal.__version__)

    # Envelope
    srs4326 = SpatialReference()
    srs4326.ImportFromEPSG(4326)
    srs = SpatialReference()
    srs.ImportFromEPSG(args.epsg)
    
    if srs.IsSame(srs4326):
        srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    
    env = Envelope()
    env.addPoint(float(args.e[0]), float(args.e[1]), 0, srs)
    env.addPoint(float(args.e[2]), float(args.e[3]), 0, srs)

    # Logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)
    
    # ToA
    toa = EvhrToA(env, args.o, logger)
    toa.run()

# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())
