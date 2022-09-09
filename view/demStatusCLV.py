#!/usr/bin/python

import argparse
import glob
import logging
import os
import sys

from evhr.model.DemCreator import DemCreator


# -----------------------------------------------------------------------------
# main
#
# evhr/view/demStatusCLV.py -i /adapt/nobackup/people/rlgill/SystemTesting/testDEM
# -----------------------------------------------------------------------------
def main():

    # Process command-line args.
    desc = 'This application reports the status of a DEM working directory.'
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-i',
                        default='.',
                        required=True,
                        help='Path to DEM directory')

    args = parser.parse_args()

    # Initialize the logger.
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)

    # Find the pair directories.
    pairDirs = glob.glob(os.path.join(args.i, 'WV0*'))

    # Check each pair directory.
    complete = []
    incomplete = []

    for pairDir in pairDirs:

        isComplete = DemCreator.demComplete(pairDir, logger)

        if isComplete:
            complete.append(pairDir)

        else:
            incomplete.append(pairDir)

    # Print results.
    print('\n', len(complete), 'complete DEMs:', complete,
          '\n\n', len(incomplete), 'incomplete DEMs:', incomplete)


# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())
