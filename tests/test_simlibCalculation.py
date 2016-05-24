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
from sqlalchemy import create_engine


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

def test_writeSimlib():
    """
    It is important to verify that the format in which the simlib file is
    output can be read in by SNANA. This is hard to do directly within travis.
    Instead, we create a simlib file and test that SNANA simulations can be run
    with it. This file is kept as a template in the example_data directory. In
    this test, we write out a simlib file to disk, and compare the values of
    strings obtained by reading this and the template simlib file read from disk
    and directly compared.
    """
    pkgDir = os.path.split(oss.__file__)[0]
    dbname = os.path.join(pkgDir, 'example_data', 'enigma_1189_micro.db')
    template_simlib = os.path.join(pkgDir, 'example_data',
                                   'Enigma_1189_micro_main.simlib')

    engineFile = 'sqlite:///' + dbname
    engine = create_engine(engineFile)

    # read the database into a `pd.DataFrame`
    Summary = pd.read_sql_table('Summary', engine)

    EnigmaMain = Summary.query('propID == [364]')
    EnigmaMainSummary = so.SummaryOpsim(EnigmaMain, calculateSNANASimlibs=True,
                                        user='rbiswas', host='time')
    simlibfilename = './Enigma_1189_micro_main.simlib'
    EnigmaMainSummary.writeSimlib(simlibfilename)

    with open(template_simlib) as f:
        template_data = f.read()
    with open(simlibfilename) as f:
        new_data = f.read()
    assert new_data == template_data
    if new_data == template_data:
        os.remove(simlibfilename)

if __name__ == '__main__':
    test_simlibValues()
    test_writeSimlib()
