"""
Test file for calculating simlibs
This file reads a couple of lines of opsim outputs and produces a simlib file
which is compared for values of skysig and zptavg against a simlib file
produced by a different code.
"""
from __future__ import absolute_import, print_function
import os
import opsimsummary as oss
import opsimsummary.summarize_opsim as so
import pandas as pd


def test_simlibValues():
    """
    read a couple of lines from a OpSim file (3168), and calculate the SNANA
    SIMLIB quantities and keep this in a dataFrame. Also, read in the simlib
    file calculated before and check that the quantities skysig and zptavg
    match.
    """
    # test_opsimObsData.dat is a snippet from the OpSim output.
    example_data = os.path.join(os.path.split(oss.__file__)[0], 'example_data')
    fname = os.path.join(example_data, 'test_opsimObsData.dat')
    # read it in
    opsim_test_df = pd.read_csv(fname, delimiter='\s+')
    # Change names so that OpSimSummary recognizes the columns
    aliases = dict(MJD='expMJD', SurveyNight='night')
    opsim_test_df.rename(columns=aliases, inplace=True)
    print(opsim_test_df.columns)
    # Run OpSimSummary to get the summary object
    testopsim = so.SummaryOpsim(opsim_test_df, calculateSNANASimlibs=True)
    # restrict to the data regarding fieldID 519
    sl = testopsim.simlib(519).sort_values(by='expMJD')
    # Read the old simlib file into a dataframe
    oldsimlib = os.path.join(example_data, 'oldSimlib.simlib')
    s_old = oss.simlib.Simlib.fromSimlibFile(oldsimlib)
    # Look at the data corresponding to fieldID 519
    old_data = s_old.simlibData(519).sort_values(by='MJD')
    from numpy.testing import assert_allclose
    print('Running assertions for:')
    print('skysig values, to tolerance of 0.01')
    assert_allclose(old_data.SKYSIG.values, sl.simLibSkySig.values, atol=0.01)
    print('zptAvg values, to tolerance of 0.01')
    assert_allclose(old_data.ZPTAVG.values, sl.simLibZPTAVG.values, atol=0.01)
    print('PSF1 values, to tolerance of 0.01')
    assert_allclose(old_data.PSF1.values, sl.simLibPsf.values, atol=0.01)

def test_simlibWrite():
    """
    read a couple of lines from a OpSim file (3168), and calculate the SNANA
    SIMLIB quantities and keep this in a dataFrame. Also, read in the simlib
    file calculated before and check that the quantities skysig and zptavg
    match.
    """
    # test_opsimObsData.dat is a snippet from the OpSim output.
    example_data = os.path.join(os.path.split(oss.__file__)[0], 'example_data')
    fname = os.path.join(example_data, 'test_opsimObsData.dat')
    # read it in
    opsim_test_df = pd.read_csv(fname, delimiter='\s+')
    # Change names so that OpSimSummary recognizes the columns
    aliases = dict(MJD='expMJD',
                   RA='fieldRA',
                   Dec='fieldDec',
                   SurveyNight='night')
    opsim_test_df.rename(columns=aliases, inplace=True)
    print("HELLO", opsim_test_df.columns)
    # Run OpSimSummary to get the summary object
    testopsim = so.SummaryOpsim(opsim_test_df, calculateSNANASimlibs=True)
    # restrict to the data regarding fieldID 519
    sl = testopsim.simlib(519).sort_values(by='expMJD')
    outfile = os.path.join(example_data, 'test_simlib')
    testopsim.writeSimlib(outfile)
    # Read the old simlib file into a dataframe
    oldsimlib = os.path.join(example_data, 'oldSimlib.simlib')
    s_old = oss.simlib.Simlib.fromSimlibFile(oldsimlib)
    # Look at the data corresponding to fieldID 519
    old_data = s_old.simlibData(519).sort_values(by='MJD')
    from numpy.testing import assert_allclose
    print('Running assertions for:')
    print('skysig values, to tolerance of 0.01')
    assert_allclose(old_data.SKYSIG.values, sl.simLibSkySig.values, atol=0.01)
    print('zptAvg values, to tolerance of 0.01')
    assert_allclose(old_data.ZPTAVG.values, sl.simLibZPTAVG.values, atol=0.01)
    print('PSF1 values, to tolerance of 0.01')
    assert_allclose(old_data.PSF1.values, sl.simLibPsf.values, atol=0.01)
if __name__ == '__main__':

    test_simlibValues()
    test_simlibWrite()
