#!/usr/bin/env python
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
import subprocess


def write_genericSimlib(simlibFilename, summary, minVisits, maxVisits, numFields,
                  author_name=None,
                  mapFile='ddf_minion_map.csv', rng=np.random.RandomState(0),
                  mwebv=0., raCol='ditheredRA', decCol='ditheredDec',
                  angleUnit='degrees', opsimversion='lsstv3',
                  indexCol='obsHistID', nside=256, fieldType='DDF',
                  opsimoutput='minion_1016_sqlite.db',
                  vetoed_hids=None,
                  opsim_output='opsim_output',
                  script_name=None,
                  surveypix_file=None):
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
    author_name : string, defaults to None
        Author name. If not provided, assumed to be None. The author name will be replaced by 
        user login.
    vetoed_hids : set of integers, defaults to `None`
       if not `None`, set of integers representing healpix Ids in the same
       (`nest`|`ring`) which should not be used. If `None`, ignore.
    opsim_output: str
        name of the opsim_output 'opsim_output'
    script_name : 
    surveypix_file : abspath to surveypix_file, defaults to None
        If not None, uses this file to find the libids to simulate. Should
        have an ordered dataframe of ra, dec, simlibIds
    """
                  
    opsimsummary_version=oss.__version__,
    minMJD = summary.expMJD.min()
    maxMJD = summary.expMJD.max()
    simlibs = Simlibs(summary, usePointingTree=True, raCol=raCol,
                      decCol=decCol, angleUnit=angleUnit,
                      indexCol=indexCol, opsimversion=opsimversion)
    surveydf = simlibs.observedVisitsinRegion(minVisits=minVisits,
                                              maxVisits=maxVisits,
                                              writeFile=False,
                                              nside=nside)

    if surveypix_file is None:
        if vetoed_hids is not None:
            hids = set(surveydf.index.values)
            selected = hids - vetoed_hids
            surveydf = surveydf.loc[selected]
        totalfields = len(surveydf)
        surveyPix = simlibs.get_surveyPix(surveydf, numFields=numFields, rng=rng)
    else:
        print('Reading in surveypix file\n')
        print('You should only be using this if you are studying ToO proposals\n')
        surveyPix = pd.read_csv(surveypix_file)
        totalfields = len(surveyPix)

    fields = simlibs.simlibs_for_fields(surveyPix, mwebv=mwebv)

    area = hp.nside2pixarea(nside, degrees=True) * np.float(totalfields)
    solidangle = hp.nside2pixarea(nside, degrees=False) * np.float(totalfields)
 
    print('Going to write simlib file {0} for opsim output\n')
    if script_name is None:
        script_name = 'OpSimSummary/scripts/make_simlibs.py'
    dt = datetime.datetime.now()
    ts = dt.isoformat()


    # Author name
    if author_name is None:
        author_name = os.getlogin()

    ## git information
    try :
        label = subprocess.check_output(['git', 'describe', '--dirty']).strip().decode('utf-8')
        # Assume that there will be git versions before and after 2.22
        branch = subprocess.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD']).strip().decode('utf-8')
        if label.endswith('dirty'):
            git_sha = label.split('-dirty')[0]
            dirty = 'dirty'
        else:
            git_sha = label
            dirty = ''

        git_info = f'The script was part of a {dirty} git repository in the branch {branch} and commit SHA {git_sha} and was labelled with a version : {oss.__version__}'
    except:
        git_info = 'The script was not detected to be a part of a git repository'
    
        
    # OpSim version information
    version = oss.__version__

    comment = 'DOCUMENTATION:\n'
    comment += f'    PURPOSE: simulate LSST based on mock opsim version {opsim_output}\n'
    comment += '    INTENT:   Nominal\n'
    comment += '    USAGE_KEY: SIMLIB_FILE\n'
    comment += '    USAGE_CODE: snlc_sim.exe\n'
    comment += '    VALIDATION_SCIENCE: \n'
    comment += '    NOTES: \n'
    comment += '        PARAMS MINMJD: {}\n'.format(minMJD)
    comment += '        PARAMS MAXMJD: {}\n'.format(maxMJD)
    comment += '        PARAMS TOTAL_AREA: {}\n'.format(area)
    comment += '        PARAMS SOLID_ANGLE: {}\n'.format(solidangle)
    comment += '    VERSIONS:\n'
    comment += f'    - DATE : {dt.strftime(format="%y-%m-%d")}\n'
    comment += f'    AUTHORS : {author_name}, OpSimSummary version {oss.__version__}\n'
    comment += 'DOCUMENTATION_END:\n'

    doc = comment

    comment = 'COMMENT: Total area corresponding to this simlib is {0:.1f} sq degrees or a solid angle of {1:4f} \n'.format(area, solidangle) 
    comment += 'COMMENT: This is a simlib corresponding to {0} in the OpSim Output {1} using the script {3} in version {2} at time {4}. {5}\n'.format(fieldType, opsimoutput,
            opsimsummary_version, script_name, ts, git_info) 
    comment += '\nCOMMENT: PARAMS MINMJD: {}\n'.format(minMJD)
    comment += 'COMMENT: PARAMS MAXMJD: {}\n'.format(maxMJD)
    comment += 'COMMENT: PARAMS TOTAL_AREA: {}\n'.format(area)
    comment += 'COMMENT: PARAMS SOLID_ANGLE: {}\n'.format(solidangle)

    simlibs.writeSimlib(simlibFilename, fields, mwebv=mwebv, doc=doc, comments=comment,
                        numLibId=numFields)
    surveyPix = surveyPix.reset_index().query('simlibId > -1').set_index('simlibId')
    surveyPix = surveyPix.reset_index().sort_values(by='simlibId').set_index('simlibId')
    surveyPix.to_csv(mapFile)
    return surveyPix, surveydf

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
    parser.add_argument('--no_get_ddf_pixels', help='Whether to obtain ddf pixels from opsim or use null set',
                        dest='No_get_ddf_pixels', action='store_false')
    # parser.add_argument('--write_wfd_simlib', help='Whether to write out WFD simlib',
    #                    dest='write_wfd_simlib', action='store_true')
    parser.add_argument('--no_write_wfd_simlib', help='Whether to write out WFD simlib, defaults to writing it out',
                        dest='write_wfd_simlib', action='store_false')
    parser.add_argument('--opsimversion', help='version of opsim used lsstv3|lsstv4, defaults to lsstv3',
                        default='lsstv3')
    parser.add_argument('--summaryTableName', help='name of Summary Table Summary|SummaryAllProps, defaults to "Summary"',
                        default='Summary')
    parser.add_argument('--ddf_surveypix_file', help='absolute path to DDF surveypix_file, defaults to `None`',
                        default=None)
    parser.add_argument('--wfd_surveypix_file', help='absolute path to WFD surveypix_file, defaults to `None`',
                        default=None)
    parser.add_argument('--ddf_simlibfilename', help='absolute path to DDF simlib file to write out, defaults to `None`',
                        default=None)
    parser.add_argument('--wfd_simlibfilename', help='absolute path to WFD simlib file to write out, defaults to `None`',
                        default=None)
    parser.add_argument('--numFields_DDF', help='number of locations in DDF where simlib fields are located, defaults to 133',
                        default=133, type=int)
    parser.add_argument('--numFields_WFD', help='number of locations in DDF where simlib fields are located, defaults to 50000',
                        default=50000, type=int)
    parser.add_argument('--filterNull', help='if added, then the summary table of the OpSim file will be filtered of rows that appear to have null values',
                        dest='filt_Null', action='store_true')
    parser.add_argument('--authorName', help='Name of the author to be recorded in documentation. If skipped, the login will be used instead', type=str, default=None)
    parser.add_argument('--opsim_output', help='name of the opsim output being used, will be obtained from filename if skipped', default=None, type=str)
    print("read in command line options and figuring out what to do\n")
    print("we are using opsimsummary version {0} and the library is located at {1}".format(oss.__version__, oss.__file__))
    print("we are using the path {}".format(sys.path))
    print("Our python version is {}".format(sys.version))
    sys.stdout.flush()
    args = parser.parse_args()
    
    data_root = args.data_root
    author_name = args.authorName
    dbname = os.path.join(data_root, args.dbname)
    basename = dbname.split('/')[-1].split('.db')[0]

    opsim_output = args.opsim_output
    if opsim_output is None:
        opsim_output = basename

    filternulls = False
    if args.filt_Null :
        print('filtering the raw summary table, with this option, null values in `fiveSigmaDepth` will be skipped')
        print("We do not recommend using this option unless you are sure you want this\n")
        filternulls = True

    get_ddf_pixels = True
    if not args.No_get_ddf_pixels:
        get_ddf_pixels = False

    print("get_ddf_pixels", get_ddf_pixels)
    print("\n\n Task: Obtaining pointing location \n")
    dithercolumns = None
    if args.No_construct_ditherfiles:
        print('Not constructing ditherfile from dbname to search for a csv file')
    else:
        if args.ditherfiles is None:
            print('\n\n Task: Constructing ditherfiles as csv files based on opsim output name')
            ditherfiles = os.path.join(args.data_root, 'descDithers_{}.csv'.format(basename))
            print('abs path to  ditherfiles is {}'.format(ditherfiles))
        else:
            print('\n\n Task: obtaining ditherfiles from input')
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
        wfd_simlibfilename = basename +'_WFD.simlib'
    availwfdFileName = basename + "_WFD_avail.csv"
    selectedwfdFileName = basename + "_WFD_sel.csv"
    print('output file names for wfd are {0}, {1}, {2}'.format(wfd_simlibfilename, availwfdFileName, selectedwfdFileName))

    if ddf_simlibfilename is None:
        ddf_simlibfilename = basename +'_DDF.simlib'
    availddfFileName = basename + "_DDF_avail.csv"
    selectedddfFileName = basename + "_DDF_sel.csv"
    print('output file names for DDF are {0}, {1}, {2}'.format(ddf_simlibfilename, availddfFileName, selectedddfFileName))
    
    numFields_DDF = args.numFields_DDF
    numFields_WFD = args.numFields_WFD
    print(args)
    
    sys.stdout.flush()
    # find ddf healpixels
    if get_ddf_pixels:
        print('Finding the DDF healpixels \n')
        opsout_ddf = OpSimOutput.fromOpSimDB(dbname, opsimversion=opsimversion,
                                             subset='ddf', dithercolumns=dithercolumns,
                                             filterNull=filternulls)
        if len(opsout_ddf.summary) > 0:
            print("writing out ddf pixels\n")
            simlib_ddf = Simlibs(opsout_ddf.summary, opsimversion=opsimversion,
                                 usePointingTree=True)
            ddf_hid = set(simlib_ddf.observedVisitsinRegion().index.values) 
    else:

        print("writing out null set of ddf pixels\n")
        ddf_hid = set([])
        print("written out null set of ddf pixels\n")
    print('There are {} pixels in the ddf fields'.format(len(ddf_hid)))
    # read the database into a `pd.DataFrame`
    tstart = time.time()
    print("\n\n Task: reading database {0} at time {1}. This can take a while ... ".format(dbname, tstart))
    sys.stdout.flush()
    opsout = OpSimOutput.fromOpSimDB(dbname,
                                     opsimversion=opsimversion,
                                     tableNames=(summaryTableName, 'Proposal'),
                                     subset='combined', dithercolumns=dithercolumns,
                                     filterNull=filternulls)
    tend = time.time()
    print("finished reading database {0} at time {1}".format(dbname, tend))
    print("reading the db took {} minutes".format((tend-tstart)/60.0))
    sys.stdout.flush()
    summary = opsout.summary
    script_name = os.path.abspath(__file__)
    if write_ddf_simlib:
        print('\n\n Task: writing out simlib for DDF')
        # 133 random locations is similar density of locations in WFD.
        x, y = write_genericSimlib(simlibFilename=ddf_simlibfilename,
                                   summary=opsout_ddf.summary, minVisits=500, maxVisits=None,
                                   numFields=numFields_DDF, mapFile='ddf_minion_1016_sqlite.csv',
                                   fieldType='DDF', opsimoutput=dbname,
                                   script_name=script_name,
                                   author_name=author_name,
                                   opsim_output=opsim_output,
                                   surveypix_file=args.ddf_surveypix_file)
        print('Finished writing out simlib for DDF')
        print('\n\n Task: write mapping to csv')
        x.to_csv(selectedddfFileName)
        y.to_csv(availddfFileName)
        print('Finished writing mapping to csv')
        sys.stdout.flush()
    sys.stdout.flush()
    if write_wfd_simlib :
        print('\n\n Task: writing out simlib for WFD')
        sys.stdout.flush()
        x, y  = write_genericSimlib(simlibFilename=wfd_simlibfilename,
                                    summary=summary, minVisits=500, maxVisits=10000,
                                    numFields=numFields_WFD, mapFile='wfd_minion_1016_sqlite.csv',
                                    fieldType='WFD', opsimoutput=dbname, 
                                    vetoed_hids=ddf_hid,
                                    script_name=script_name,
                                    surveypix_file=args.wfd_surveypix_file)
        print('Finished writing out simlib for WFD')
        print('\n\n Task: write mapping to csv')
        x.to_csv(selectedwfdFileName)
        y.to_csv(availwfdFileName)
        print('Finished writing mapping to csv')
        sys.stdout.flush()
    print('finished job')
    sys.stdout.flush()
