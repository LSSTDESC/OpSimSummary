"""
Script to make Simlib files from OpSim DataBases (supported versions include
OpSim v3, OpSim v4, AltSched outputs.
    To get usage : python make_simlibs.py -h
"""
import os
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import healpy as hp
import opsimsummary as oss
from opsimsummary import Simlibs, OpSimOutput
import sys
import time
import datetime


def write_genericSimlib(simlibFilename, summary, minVisits, maxVisits, numFields,
                  mapFile='ddf_minion_map.csv', rng=np.random.RandomState(0),
                  mwebv=0., raCol='ditheredRA', decCol='ditheredDec',
                  angleUnit='degrees', opsimversion='lsstv3',
                  indexCol='obsHistID', nside=256, fieldType='DDF',
                  opsimoutput='minion_1016_sqlite.db',
                  vetoed_hids=None,
                  opsimsummary_version=oss.__version__,
                  script_name=None):
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
    opsimsummary_version :
    script_name : 
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
    
    print('Going to write simlib file {0} for opsim output\n')
    if script_name is None:
        script_name = 'OpSimSummary/scripts/make_simlibs.py'
    ts = datetime.datetime.now().isoformat()

    comment = 'COMMENT: Total area corresponding to this simlib is {0:.1f} sq degrees or a solid angle of {1:4f} \n'.format(area, solidangle) 
    comment += 'COMMENT: This is a simlib corresponding to {0} in the OpSim Output {1} using the script {3} in version {2} at time {4}\n'.format(fieldType, opsimoutput, opsimsummary_version, script_name, ts) 
    simlibs.writeSimlib(simlibFilename, fields, mwebv=mwebv, comments=comment,
                        numLibId=numFields)
    surveyPix = surveyPix.reset_index().query('simlibId > -1').set_index('simlibId')
    surveyPix = surveyPix.reset_index().sort_values(by='simlibId').set_index('simlibId')
    surveyPix.to_csv(mapFile)
    return surveyPix

if __name__ == '__main__':
    parser = ArgumentParser(description='write out simlibs from an OpSim Database')
    parser.add_argument('--data_root', help='absolute path to data directory containing dither files and opsim database, defaults to "/"',
                        default='/')
    parser.add_argument('--dbname', help='path to sqlite database output from OpSim relative to data_root, defaults to "minion_1016_sqlite.db"',
                        default='minion_1016_sqlite.db')
    parser.add_argument('--ditherfiles', help='path to ditherfile relative to data root, defaults to `None`', default=None) 
    parser.add_argument('--no_construct_ditherfiles', help='do not try to construct ditherfiles if ditherfiles is None',
                        dest='No_construct_ditherfiles', action='store_true')
    # parser.add_argument('--write_ddf_simlib', help='Whether to write out DDF simlib',
    #                    dest='write_ddf_simlib', action='store_true')
    parser.add_argument('--no_write_ddf_simlib', help='Whether to write out DDF simlib',
                        dest='write_ddf_simlib', action='store_false')
    # parser.add_argument('--write_wfd_simlib', help='Whether to write out WFD simlib',
    #                    dest='write_wfd_simlib', action='store_true')
    parser.add_argument('--no_write_wfd_simlib', help='Whether to write out WFD simlib, defaults to writing it out',
                        dest='write_wfd_simlib', action='store_false')
    parser.add_argument('--opsimversion', help='version of opsim used lsstv3|lsstv4, defaults to lsstv3',
                        default='lsstv3')
    parser.add_argument('--summaryTableName', help='name of Summary Table Summary|SummaryAllProps, defaults to "Summary"',
                        default='Summary')
    parser.add_argument('--ddf_simlibfilename', help='absolute path to DDF simlib file to write out, defaults to `None`',
                        default=None)
    parser.add_argument('--wfd_simlibfilename', help='absolute path to WFD simlib file to write out, defaults to `None`',
                        default=None)
    parser.add_argument('--numFields_DDF', help='number of locations in DDF where simlib fields are located, defaults to 133',
                        default=133, type=int)
    parser.add_argument('--numFields_WFD', help='number of locations in DDF where simlib fields are located, defaults to 50000',
                        default=50000, type=int)
    
    args = parser.parse_args()
    
    data_root = args.data_root
    dbname = os.path.join(data_root, args.dbname)
    basename = dbname.split('/')[-1].split('.db')[0]

    print("read in command line options and figuring out what to do\n")

    print("Obtaining pointing location \n")
    dithercolumns = None
    if args.No_construct_ditherfiles:
        print('Not constructing ditherfile from dbname')
    else:
        if args.ditherfiles is None:
            print('Constructing ditherfiles')
            ditherfiles = os.path.join(args.data_root, 'descDithers_{}.csv'.format(basename))
            print('abs path to  ditherfiles is {}'.format(ditherfiles))
        else:
            print('obtaining ditherfiles from input')
            ditherfiles = os.path.join(args.data_root, args.ditherfiles)
            print(ditherfiles)
            assert os.path.exists(ditherfiles)
            assert os.path.getsize(ditherfiles) > 0
    
            dithercolumns = pd.read_csv(ditherfiles)
            dithercolumns = dithercolumns.rename(columns=dict(observationId='obsHistID',
                                                               descDitheredRA='ditheredRA',
                                                               descDitheredDec='ditheredDec')
                                                              ).set_index('obsHistID')
    
    sys.stdout.flush()
    summaryTableName = args.summaryTableName
    opsimoutput = os.path.basename(dbname)
    opsimversion = args.opsimversion
    write_wfd_simlib = args.write_wfd_simlib
    write_ddf_simlib = args.write_ddf_simlib
    wfd_simlibfilename = args.wfd_simlibfilename
    ddf_simlibfilename = args.ddf_simlibfilename
    
    if wfd_simlibfilename is None:
        wfd_simlibfilename = basename +'_wfd.simlib'
    if ddf_simlibfilename is None:
        ddf_simlibfilename = basename +'_ddf.simlib'
    
    numFields_DDF = args.numFields_DDF
    numFields_WFD = args.numFields_WFD
    print(args)
    
    sys.stdout.flush()
    # find ddf healpixels
    opsout_ddf = OpSimOutput.fromOpSimDB(dbname, opsimversion=opsimversion,
                                         subset='ddf', dithercolumns=dithercolumns)
    simlib_ddf = Simlibs(opsout_ddf.summary, opsimversion=opsimversion,
                         usePointingTree=True)
    ddf_hid = set(simlib_ddf.observedVisitsinRegion().index.values) 
    print('There are {} pixels in the ddf fields'.format(len(ddf_hid)))
    # read the database into a `pd.DataFrame`
    tstart = time.time()
    print("reading database {0} at time {1}. This can take a while ... ".format(dbname, tstart))
    sys.stdout.flush()
    opsout = OpSimOutput.fromOpSimDB(dbname,
                                     opsimversion=opsimversion,
                                     tableNames=(summaryTableName, 'Proposal'),
                                     subset='combined', dithercolumns=dithercolumns)
    tend = time.time()
    print("finished reading database {0} at time {1}".format(dbname, tend))
    print("reading the db took {} minutes".format((tend-tstart)/60.0))
    summary = opsout.summary
    script_name = os.path.abspath(__file__)
    if write_ddf_simlib:
        print('writing out simlib for DDF')
        # 133 random locations is similar density of locations in WFD.
        x = write_genericSimlib(simlibFilename=ddf_simlibfilename,
                          summary=opsout_ddf.summary, minVisits=500, maxVisits=None,
                          numFields=numFields_DDF, mapFile='ddf_minion_1016_sqlite.csv',
                          fieldType='DDF', opsimoutput=dbname,
                          script_name=script_name)
    if write_wfd_simlib :
        print('writing out simlib for WFD')
        x = write_genericSimlib(simlibFilename=wfd_simlibfilename,
                          summary=summary, minVisits=500, maxVisits=10000,
                          numFields=numFields_WFD, mapFile='wfd_minion_1016_sqlite.csv',
                          fieldType='WFD', opsimoutput=dbname, 
                          vetoed_hids=ddf_hid, script_name=script_name)
    print('finished job')
    sys.stdout.flush()
