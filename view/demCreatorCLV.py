#!/usr/bin/python

import csv
import argparse
import logging
import os
import pathlib
import sys

from core.model.ILProcessController import ILProcessController
from evhr.model.DemCreator import DemCreator
from evhr.model.DemCreatorCelery import DemCreatorCelery


# -----------------------------------------------------------------------------
# main
#
# evhr/view/demCreatorCLV.py -t \
#    -o /explore/nobackup/projects/ilab/scratch/SystemTesting/evhr/testDEM \
#    --pairs_in_file /explore/nobackup/projects/ilab/scratch/SystemTesting/evhr/testDEM/demcreator_test_1_pair.csv \
#    --logToFile
# -----------------------------------------------------------------------------
def main():

    # Process command-line args.
    desc = 'Use this application to produce EVHR DEMs.'
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('--pairs_in_file',
                        type=pathlib.Path,
                        help='Fully-qualified path to CSV file containing a '
                        'list of scene files')

    parser.add_argument('--celery',
                        action='store_true',
                        help='Use Celery for distributed processing.')

    parser.add_argument('--logToFile',
                        action='store_true',
                        help='Print messages to a file, \
                              instead of the screen.')

    parser.add_argument('-o',
                        default='.',
                        help='Path to output directory')

    parser.add_argument('-t',
                        action='store_true',
                        help='Run in test mode for speed.')

    args = parser.parse_args()

    # ---
    # Logging
    # ---
    if args.logToFile:

        logFile = os.path.join(args.o, 'dem.out')
        logging.basicConfig(filename=logFile)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)

    # ---
    # Scene pair list file
    #
    # Previous version of DemCreator took various identifiers which
    # were then used to identify pairs through the footprints database.
    #
    # Now, users should directly supply the pairs which they want the
    # dem created for in the form of a CSV where each row represents a
    # pair. Col 1 is reserved for the pair name and the subsequent columns
    # should be paths to scenes which make up the pair to be processed.
    # ---
    pairs = {}

    if args.pairs_in_file:

        with open(args.pairs_in_file, newline='') as csvFile:
            reader = csv.reader(csvFile)
            for row in reader:
                key = row[0]
                pairs[key] = row[1:]

    # ---
    # Instantiate the process object and run.
    # ---
    if args.celery:

        with ILProcessController('evhr.model.CeleryConfiguration') as \
             processController:

            dc = DemCreatorCelery(args.o, logger, args.t)

            dc.runPairs(pairs)

    else:

        dc = DemCreator(args.o, logger, args.t)

        dc.runPairs(pairs)


# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())
