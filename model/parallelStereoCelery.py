#!/usr/bin/python

import argparse
import sys


# -------------------------------------------------------------------------
# main
# -------------------------------------------------------------------------
def main():

    # parOpts, sgmOpts, stereoOpts, leftTifGlob, rightTifGlob, leftXmlGlob, rightXmlGlob, rpcDem 
    desc = 'Use this application to produce EVHR DEMs.'
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument('-e',
                        help='a string set of arguments from dg_stereo.sh.')

    parser.add_argument('--parOpts',
                        help='a string set of arguments from dg_stereo.sh.')

    parser.add_argument('--sgmOpts',
                        help='a string set of arguments from dg_stereo.sh.')

    parser.add_argument('--stereoOpts',
                        help='a string set of arguments from dg_stereo.sh.')

    parser.add_argument('--leftTifs',
                        help='a string set of arguments from dg_stereo.sh.')

    parser.add_argument('--rightTifs',
                        help='a string set of arguments from dg_stereo.sh.')

    parser.add_argument('--leftXmls',
                        help='a string set of arguments from dg_stereo.sh.')

    parser.add_argument('--rightXmls',
                        help='a string set of arguments from dg_stereo.sh.')

    parser.add_argument('--rpcDem',
                        help='a string set of arguments from dg_stereo.sh.')

    args = parser.parse_args()

    # cmd = 'parallel_stereo' + \
    #       ' -e ' + args.e + \
    #       ' ' + args.parOpts + \
    #       ' ' + args.stereoOpts + \
    #       ' ' + args.leftTifs + \
    #       ' ' + args.leftXmls + \
    #       ' ' + args.rightTifs + \
    #       ' ' + args.rightXmls + \
    #       ' ' + args.rpcDem
    #
    # import os
    # os.system(cmd)
    
    print('lefts:', args.leftTifs)
    

# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())

