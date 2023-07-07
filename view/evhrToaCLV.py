#!/usr/bin/python
import argparse
import csv
import logging
import pathlib
import sys

from osgeo import gdal
from osgeo import osr
from osgeo.osr import SpatialReference

from core.model.DgFile import DgFile
from core.model.Envelope import Envelope
from core.model.ILProcessController import ILProcessController

from evhr.model.EvhrToA import EvhrToA
from evhr.model.EvhrToaCelery import EvhrToaCelery


# -----------------------------------------------------------------------------
# main
#
# evhr/view/evhrToaCLV.py -e -148 65 -147.5 64.5 4326 -o /explore/nobackup/people/rlgill/SystemTesting/testEVHR/ --pan_res 0.5
#
# -----
# evhr/view/evhrToaCLV.py -o /adapt/nobackup/people/rlgill/SystemTesting/testEVHR/ --scenes_in_file /adapt/nobackup/people/rlgill/SystemTesting/scene-list-file.csv
# 
# -----
# The time it takes to complete the scenes below:
#
# real    280m0.978s
# user    521m44.397s
# sys     36m7.956s
#
# evhr/view/evhrToaCLV.py -o /explore/nobackup/people/rlgill/SystemTesting/testEVHR --scenes /css/nga/WV03/1B/2015/219/WV03_104001000F2D9E00_X1BS_500495393030_01/WV03_20150807213524_104001000F2D9E00_15AUG07213524-M1BS-500495393030_01_P001.ntf /css/nga/WV03/1B/2015/219/WV03_104001000F2D9E00_X1BS_500495393030_01/WV03_20150807213525_104001000F2D9E00_15AUG07213525-M1BS-500495393030_01_P002.ntf /css/nga/WV02/1B/2014/220/WV02_1030010035034E00_X1BS_500138791160_01/WV02_20140808212142_1030010035034E00_14AUG08212142-M1BS-500138791160_01_P001.ntf /css/nga/WV02/1B/2012/226/WV02_103001001B775200_X1BS_052903554020_01/WV02_20120813213619_103001001B775200_12AUG13213619-M1BS-052903554020_01_P002.ntf /css/nga/WV02/1B/2012/226/WV02_103001001B775200_X1BS_052903554020_01/WV02_20120813213618_103001001B775200_12AUG13213618-M1BS-052903554020_01_P001.ntf /css/nga/WV02/1B/2010/215/WV02_103001000621E500_X1BS_052807177090_01/WV02_20100803220600_103001000621E500_10AUG03220600-M1BS-052807177090_01_P002.ntf /css/nga/WV02/1B/2010/215/WV02_103001000621E500_X1BS_052807177090_01/WV02_20100803220558_103001000621E500_10AUG03220558-M1BS-052807177090_01_P001.ntf /css/nga/WV02/1B/2010/226/WV02_1030010006673100_X1BS_052807135060_01/WV02_20100814220608_1030010006673100_10AUG14220608-M1BS-052807135060_01_P002.ntf /css/nga/WV02/1B/2010/226/WV02_1030010006673100_X1BS_052807135060_01/WV02_20100814220607_1030010006673100_10AUG14220607-M1BS-052807135060_01_P001.ntf /css/nga/WV02/1B/2018/223/WV02_1030010080D1AB00_X1BS_502522743050_01/WV02_20180811214527_1030010080D1AB00_18AUG11214527-M1BS-502522743050_01_P002.ntf /css/nga/WV02/1B/2018/223/WV02_1030010080D1AB00_X1BS_502522743050_01/WV02_20180811214526_1030010080D1AB00_18AUG11214526-M1BS-502522743050_01_P001.ntf /css/nga/WV02/1B/2013/214/WV02_10300100252FDE00_X1BS_500122657100_01/WV02_20130802222758_10300100252FDE00_13AUG02222758-M1BS-500122657100_01_P001.ntf /css/nga/WV02/1B/2016/226/WV02_103001005B3B4A00_X1BS_501569343080_01/WV02_20160813220348_103001005B3B4A00_16AUG13220348-M1BS-501569343080_01_P002.ntf /css/nga/WV02/1B/2016/226/WV02_103001005B3B4A00_X1BS_501569343080_01/WV02_20160813220347_103001005B3B4A00_16AUG13220347-M1BS-501569343080_01_P001.ntf /css/nga/WV03/1B/2016/244/WV03_1040010021722200_X1BS_501553711060_01/WV03_20160831213307_1040010021722200_16AUG31213307-M1BS-501553711060_01_P002.ntf /css/nga/WV03/1B/2016/244/WV03_1040010021722200_X1BS_501553711060_01/WV03_20160831213305_1040010021722200_16AUG31213305-M1BS-501553711060_01_P001.ntf /css/nga/WV03/1B/2015/233/WV03_10400100100CD800_X1BS_500496271060_01/WV03_20150821215730_10400100100CD800_15AUG21215730-M1BS-500496271060_01_P002.ntf /css/nga/WV03/1B/2016/231/WV03_1040010021946400_X1BS_501553705070_01/WV03_20160818212236_1040010021946400_16AUG18212236-M1BS-501553705070_01_P001.ntf /css/nga/WV03/1B/2018/223/WV03_1040010040D77A00_X1BS_502531383100_01/WV03_20180811212909_1040010040D77A00_18AUG11212909-M1BS-502531383100_01_P003.ntf /css/nga/WV03/1B/2016/231/WV03_1040010021946400_X1BS_501553705070_01/WV03_20160818212238_1040010021946400_16AUG18212238-M1BS-501553705070_01_P002.ntf /css/nga/WV02/1B/2014/220/WV02_10300100353E2600_X1BS_500138791160_01/WV02_20140808212038_10300100353E2600_14AUG08212038-M1BS-500138791160_01_P001.ntf /css/nga/WV02/1B/2010/213/WV02_10300100061E1100_X1BS_052807177050_01/WV02_20100801213937_10300100061E1100_10AUG01213937-M1BS-052807177050_01_P002.ntf /css/nga/WV02/1B/2010/243/WV02_1030010006BB7000_X1BS_052807142030_01/WV02_20100831214833_1030010006BB7000_10AUG31214833-M1BS-052807142030_01_P002.ntf /css/nga/WV02/1B/2011/242/WV02_103001000D356800_X1BS_052807801080_01/WV02_20110830221615_103001000D356800_11AUG30221615-M1BS-052807801080_01_P002.ntf /css/nga/WV02/1B/2015/233/WV02_1030010047CC2B00_X1BS_500534852100_01/WV02_20150821220448_1030010047CC2B00_15AUG21220448-M1BS-500534852100_01_P002.ntf /css/nga/WV02/1B/2013/213/WV02_1030010026BA7A00_X1BS_500122086120_01/WV02_20130801212448_1030010026BA7A00_13AUG01212448-M1BS-500122086120_01_P001.ntf /css/nga/WV02/1B/2013/213/WV02_1030010026BA7A00_X1BS_500122086120_01/WV02_20130801212449_1030010026BA7A00_13AUG01212449-M1BS-500122086120_01_P002.ntf /css/nga/WV02/1B/2010/227/WV02_1030010006D70C00_X1BS_052807095100_01/WV02_20100815213045_1030010006D70C00_10AUG15213045-M1BS-052807095100_01_P001.ntf /css/nga/WV02/1B/2010/227/WV02_1030010006D70C00_X1BS_052807095100_01/WV02_20100815213046_1030010006D70C00_10AUG15213046-M1BS-052807095100_01_P002.ntf /css/nga/WV02/1B/2013/217/WV02_1030010024148600_X1BS_500122111110_01/WV02_20130805221743_1030010024148600_13AUG05221743-M1BS-500122111110_01_P002.ntf /css/nga/WV02/1B/2013/217/WV02_1030010024148600_X1BS_500122111110_01/WV02_20130805221742_1030010024148600_13AUG05221742-M1BS-500122111110_01_P001.ntf /css/nga/WV02/1B/2011/237/WV02_103001000C493400_X1BS_052807790060_01/WV02_20110825215943_103001000C493400_11AUG25215943-M1BS-052807790060_01_P002.ntf /css/nga/WV02/1B/2013/237/WV02_103001002562AA00_X1BS_500122494010_01/WV02_20130825213959_103001002562AA00_13AUG25213959-M1BS-500122494010_01_P002.ntf /css/nga/WV02/1B/2013/237/WV02_103001002562AA00_X1BS_500122494010_01/WV02_20130825213959_103001002562AA00_13AUG25213959-M1BS-500122494010_01_P001.ntf /css/nga/WV02/1B/2014/222/WV02_1030010034909B00_X1BS_500204185020_01/WV02_20140810214704_1030010034909B00_14AUG10214704-M1BS-500204185020_01_P002.ntf /css/nga/WV02/1B/2020/219/WV02_10300100AB1CDD00_M1BS_504712909030_01/WV02_20200806212340_10300100AB1CDD00_20AUG06212340-M1BS-504712909030_01_P002.ntf /css/nga/WV02/1B/2020/219/WV02_10300100AB1CDD00_M1BS_504712909030_01/WV02_20200806212340_10300100AB1CDD00_20AUG06212340-M1BS-504712909030_01_P001.ntf /css/nga/WV02/1B/2020/218/WV02_10300100ABC88800_M1BS_504712904080_01/WV02_20200805215944_10300100ABC88800_20AUG05215944-M1BS-504712904080_01_P001.ntf /css/nga/WV02/1B/2020/230/WV02_10300100AB5D9200_M1BS_504712939060_01/WV02_20200817211756_10300100AB5D9200_20AUG17211756-M1BS-504712939060_01_P002.ntf /css/nga/WV02/1B/2020/218/WV02_10300100ABC88800_M1BS_504712904080_01/WV02_20200805215944_10300100ABC88800_20AUG05215944-M1BS-504712904080_01_P002.ntf /css/nga/WV02/1B/2020/230/WV02_10300100AB5D9200_M1BS_504712939060_01/WV02_20200817211756_10300100AB5D9200_20AUG17211756-M1BS-504712939060_01_P001.ntf
# -----------------------------------------------------------------------------
def main():

    # Process command-line args.
    desc = 'Use this application to produce top of atmosphere images.'
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('--celery',
                        action='store_true',
                        help='Use Celery for distributed processing.')

    parser.add_argument('-o',
                        default='.',
                        help='Path to output directory')

    parser.add_argument('--pan_res',
                        type=float,
                        default=1,
                        choices=[0.5, 1],
                        help='The resolution, in meters, of panchromatic '
                             'output images')

    parser.add_argument('--pan_sharpen',
                        action='store_true',
                        help='Apply panchromatic sharpening to the output '
                             'ToA images.')

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('-e',
                       nargs=5,
                       help='ulx uly lrx lry epsg-code')

    group.add_argument('--scenes',
                       type=pathlib.Path,
                       nargs='*',
                       help='Fully-qualified path to scene files')

    group.add_argument('--scenes_in_file',
                       type=pathlib.Path,
                       help='Fully-qualified path to CSV file containing a '
                            'list of scene files')

    args = parser.parse_args()

    print('GDAL version:', gdal.__version__)

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
    # Scene list file
    # ---
    scenes = args.scenes

    if args.scenes_in_file:

        with open(args.scenes_in_file, newline='') as csvFile:
            reader = csv.reader(csvFile)
            scenes = [scene[0] for scene in reader]

    # ---
    # Logging
    # ---
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)

    # ---
    # Make DgFiles.
    # ---
    dgScenes = [DgFile(str(s), logger) for s in scenes] if scenes else None

    if args.celery:

        with ILProcessController('evhr.model.CeleryConfiguration') as \
                processController:

            toa = EvhrToaCelery(args.o, args.pan_res, args.pan_sharpen, logger)
            toa.run(env, dgScenes)

    else:

        toa = EvhrToA(args.o, args.pan_res, args.pan_sharpen, logger)
        toa.run(env, dgScenes)


# -----------------------------------------------------------------------------
# Invoke the main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())
