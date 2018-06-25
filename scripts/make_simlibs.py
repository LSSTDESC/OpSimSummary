"""
Script to make Simlib files :
    To get usage : python $0 -h
"""
import os
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import healpy as hp
import opsimsummary as oss
from opsimsummary import Simlibs, OpSimOutput



def genericSimlib(simlibFilename, summary, minVisits, maxVisits, numFields,
                  mapFile='ddf_minion_map.csv', rng=np.random.RandomState(0),
                  mwebv=0., raCol='ditheredRA', decCol='ditheredDec',
                  angleUnit='degrees', opsimversion='lsstv3',
                  indexCol='obsHistID', nside=256, fieldType='DDF',
                  opsimoutput='minion_1016_sqlite.db',
                  opsimsummary_version=oss.__version__):
    """
    Write out simlibs from a summary dataFrame
    """
    simlibs = Simlibs(summary, usePointingTree=True, raCol=raCol,
                      decCol=decCol, angleUnit=angleUnit,
                      indexCol=indexCol, opsimversion=opsimversion)
    surveydf = simlibs.observedVisitsinRegion(minVisits=minVisits,
                                              maxVisits=maxVisits,
                                              writeFile=False,
                                              nside=nside)
    totalfields = len(surveydf)
    surveyPix = simlibs.get_surveyPix(surveydf, numFields=numFields, rng=rng)
    fields = simlibs.simlibs_for_fields(surveyPix, mwebv=mwebv)
    
    area = hp.nside2pixarea(nside, degrees=True) * np.float(totalfields)
    comment = '# Total area corresponding to this simlib is {0:.1f} sq degrees\n'.format(area) 
    comment += '# This is a simlib corresponding to {0} in the OpSim Output {1} using the script OpSimSummary/script/makesimlibs.py in version {2}\n'.format(fieldType, opsimOutput, OpSimSummary_version) 
    simlibs.writeSimlib(simlibFilename, fields, mwebv=mwebv, comments=comment)
    surveyPix = surveyPix.reset_index().query('simlibId > -1').set_index('simlibId')
    surveyPix = surveyPix.reset_index().sort_values(by='simlibId').set_index('simlibId')
    surveyPix.to_csv(mapFile)
    return surveyPix

parser = ArgumentParser(description='write out simlibs from an OpSim Database')
parser.add_argument('--dbname', help='absolute path to sqlite database output from OpSim',
                    default='/Users/rbiswas/data/LSST/OpSimData/minion_1016_sqlite.db')
parser.add_argument('--write_ddf_simlib', help='Whether to write out DDF simlib',
                    type=bool, default=True)
parser.add_argument('--write_wfd_simlib', help='Whether to write out WFD simlib',
                    type=bool, default=True)
parser.add_argument('--opsimversion', help='version of opsim used lsstv3|lsstv4',
                    default='lsstv3')
parser.add_argument('--summaryTableName', help='name of Summary Table Summary|SummaryAllProps',
                    default='Summary')
parser.add_argument('--ddf_simlibfilename', help='absolute path to DDF simlib file to write out',
                    default='ddf_minion_1016_sqlite.simlib')
parser.add_argument('--wfd_simlibfilename', help='absolute path to WFD simlib file to write out',
                    default='wfd_minion_1016_sqlite.simlib')

args = parser.parse_args()

dbname = args.dbname
summaryTableName = args.summaryTableName
opsimoutput = os.path.basename(dbname)
opsimversion = args.opsimversion
write_wfd_simlib = args.write_wfd_simlib
write_ddf_simlib = args.write_ddf_simlib
wfd_simlibfilename = args.wfd_simlibfilename
ddf_simlibfilename = args.ddf_simlibfilename

print(args)
# read the database into a `pd.DataFrame`
opsout = OpSimOutput.fromOpSimDB(dbname,
                                 opsimversion=opsimversion,
                                 tableNames=(summaryTableName, 'Proposal'),
                                 subset='combined')
summary = opsout.summary
if write_ddf_simlib :
    x = genericSimlib(simlibFilename=ddf_simlibfilename,
                      summary=summary, minVisits=10000, maxVisits=None,
                      numFields=50000, mapFile='ddf_minion_1016_sqlite.csv',
                      fieldType='DDF', opsimoutput=opsimoutput)
if write_wfd_simlib :
    x = genericSimlib(simlibFilename=wfd_simlibfilename,
                      summary=summary, minVisits=500, maxVisits=10000,
                      numFields=50000, mapFile='wfd_minion_1016_sqlite.csv',
                      fieldType='WFD', opsimoutput=opsimoutput)
