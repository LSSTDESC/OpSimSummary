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
                  vetoed_hids=None,
                  opsimsummary_version=oss.__version__):
    """
    Write out simlibs from a summary dataFrame

    Parameters
    ----------
    simlibFilename : string
        absolute path to simlib file to be written out
    summary : `pd.dataFrame`
        summary of observations
    minVisits : int
        minimum number of visits, such that fields with lower numbers of visits
        over the survey are not considered.
    maxVisits : int
        maximum number of visits, such that fields with higher numbers of visits
        over the survey are not considered.
    numFields : int
        number of sample fields requested in the simlib file. If this number
        exceeds the number of fields satisfying the number of visits, all of the
        fields satisfying the constraints are included in the simlib
    vetoed_hids : set of integers, defaults to `None`
       if not `None`, set of integers representing healpix Ids in the same
       (`nest`|`ring`) which should not be used. If `None`, ignore.
    """
    simlibs = Simlibs(summary, usePointingTree=True, raCol=raCol,
                      decCol=decCol, angleUnit=angleUnit,
                      indexCol=indexCol, opsimversion=opsimversion)
    surveydf = simlibs.observedVisitsinRegion(minVisits=minVisits,
                                              maxVisits=maxVisits,
                                              writeFile=False,
                                              nside=nside)
    if vetoed_hids is not None:
        hids = set(surveydf.index.values)
        selected = hids - vetoed_hids
        surveydf = surveydf.loc[selected]
    totalfields = len(surveydf)
    surveyPix = simlibs.get_surveyPix(surveydf, numFields=numFields, rng=rng)
    fields = simlibs.simlibs_for_fields(surveyPix, mwebv=mwebv)
    
    area = hp.nside2pixarea(nside, degrees=True) * np.float(totalfields)
    solidangle = hp.nside2pixarea(nside, degrees=False) * np.float(totalfields)
    comment = 'COMMENT: Total area corresponding to this simlib is {0:.1f} sq degrees or a solid angle of {1:4f} \n'.format(area, solidangle) 
    comment += 'COMMENT: This is a simlib corresponding to {0} in the OpSim Output {1} using the script OpSimSummary/script/makesimlibs.py in version {2}\n'.format(fieldType, opsimoutput, opsimsummary_version) 
    simlibs.writeSimlib(simlibFilename, fields, mwebv=mwebv, comments=comment,
                        numLibId=numFields)
    surveyPix = surveyPix.reset_index().query('simlibId > -1').set_index('simlibId')
    surveyPix = surveyPix.reset_index().sort_values(by='simlibId').set_index('simlibId')
    surveyPix.to_csv(mapFile)
    return surveyPix

parser = ArgumentParser(description='write out simlibs from an OpSim Database')
parser.add_argument('--dbname', help='absolute path to sqlite database output from OpSim',
                    default='/Users/rbiswas/data/LSST/OpSimData/minion_1016_sqlite.db')
# parser.add_argument('--write_ddf_simlib', help='Whether to write out DDF simlib',
#                    dest='write_ddf_simlib', action='store_true')
parser.add_argument('--no_write_ddf_simlib', help='Whether to write out DDF simlib',
                    dest='write_ddf_simlib', action='store_false')
# parser.add_argument('--write_wfd_simlib', help='Whether to write out WFD simlib',
#                    dest='write_wfd_simlib', action='store_true')
parser.add_argument('--no_write_wfd_simlib', help='Whether to write out WFD simlib',
                    dest='write_wfd_simlib', action='store_false')
parser.add_argument('--opsimversion', help='version of opsim used lsstv3|lsstv4',
                    default='lsstv3')
parser.add_argument('--summaryTableName', help='name of Summary Table Summary|SummaryAllProps',
                    default='Summary')
parser.add_argument('--ddf_simlibfilename', help='absolute path to DDF simlib file to write out',
                    default='ddf_minion_1016_sqlite.simlib')
parser.add_argument('--wfd_simlibfilename', help='absolute path to WFD simlib file to write out',
                    default='wfd_minion_1016_sqlite.simlib')
parser.add_argument('--numFields_DDF', help='number of locations in DDF where simlib fields are located',
                    default=133, type=int)
parser.add_argument('--numFields_WFD', help='number of locations in DDF where simlib fields are located',
                    default=50000, type=int)

args = parser.parse_args()

dbname = args.dbname
summaryTableName = args.summaryTableName
opsimoutput = os.path.basename(dbname)
opsimversion = args.opsimversion
write_wfd_simlib = args.write_wfd_simlib
write_ddf_simlib = args.write_ddf_simlib
wfd_simlibfilename = args.wfd_simlibfilename
ddf_simlibfilename = args.ddf_simlibfilename
numFields_DDF = args.numFields_DDF
numFields_WFD = args.numFields_WFD

print(args)
# find ddf healpixels
opsout_ddf = OpSimOutput.fromOpSimDB(dbname, opsimversion=opsimversion,
                                     subset='ddf')
simlib_ddf = Simlibs(opsout_ddf.summary, opsimversion=opsimversion,
                     usePointingTree=True)
ddf_hid = set(simlib_ddf.observedVisitsinRegion().index.values) 
print('There are {} pixels in the ddf fields'.format(len(ddf_hid)))
# read the database into a `pd.DataFrame`
opsout = OpSimOutput.fromOpSimDB(dbname,
                                 opsimversion=opsimversion,
                                 tableNames=(summaryTableName, 'Proposal'),
                                 subset='combined')
summary = opsout.summary
if write_ddf_simlib :
    print('writing out simlib for DDF')
    # 133 random locations is similar density of locations in WFD.
    x = genericSimlib(simlibFilename=ddf_simlibfilename,
                      summary=opsout_ddf.summary, minVisits=500, maxVisits=None,
                      numFields=numFields_DDF, mapFile='ddf_minion_1016_sqlite.csv',
                      fieldType='DDF', opsimoutput=opsimoutput)
if write_wfd_simlib :
    print('writing out simlib for WFD')
    x = genericSimlib(simlibFilename=wfd_simlibfilename,
                      summary=summary, minVisits=500, maxVisits=10000,
                      numFields=numFields_WFD, mapFile='wfd_minion_1016_sqlite.csv',
                      fieldType='WFD', opsimoutput=opsimoutput, 
                      vetoed_hids=ddf_hid)
