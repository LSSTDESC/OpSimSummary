""" Tests for the code in `opsimsummary/opsim_out.py`
"""
from __future__ import print_function, division, absolute_import
import os
import pytest
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
import opsimsummary as oss
from opsimsummary import OpSimOutput


@pytest.fixture()
def opsfile(fname):
    yield os.path.join(oss.example_data, fname)
    print('teardown')

testdata_propIDDict = [('enigma_1189_micro.db', 'lsstv3', (364, 366)),
                       ('opsimv4_feat_micro.db','sstf', (0, 1))]
@pytest.mark.parametrize("fname,opsimversion,expected", testdata_propIDDict)
def test_get_propIDDict(fname, opsimversion, expected):
    fname = os.path.join(oss.example_data, fname)
    engine = create_engine('sqlite:///' + fname)
    df = pd.read_sql_table(con=engine, table_name='Proposal')
    mydict = OpSimOutput.get_propIDDict(df, opsimversion=opsimversion)
    assert mydict['ddf'] == expected[1]
    assert mydict['wfd'] == expected[0]

test_fromOpSimDB = [('enigma_1189_micro.db', 'lsstv3', 'Summary',
                     (76665, 5, 'radians')),
                    ('opsimv4_feat_micro.db', 'sstf', 'SummaryAllProps',
                     (90846, 5, 'degrees'))]
@pytest.mark.parametrize("fname,opsimversion,tableName,expected", test_fromOpSimDB)
def test_OpSimOutput_fromOmSimDB(fname, opsimversion, tableName, expected):
    fname = os.path.join(oss.example_data, fname)
    opsout = OpSimOutput.fromOpSimDB(fname, opsimversion=opsimversion,
                                     tableNames=(tableName, 'Proposal'))

    # Check the length of the summary table
    assert len(opsout.summary) == expected[0]

    # Check the length of the summary table
    assert len(opsout.proposalTable) == expected[1]

    # Check that 'wfd', and 'ddf' are in the proposal table
    propNameCol = opsout.opsimVars['propName'] 
    proposalNames = opsout.proposalTable[propNameCol].values
    assert 'wfd' in proposalNames
    assert 'ddf' in proposalNames

    # Check that ra are in the right units
    assert opsout.opsimVars['angleUnits'] == expected[2]
    if opsout.opsimVars['angleUnits'] == 'degrees':
        assert opsout.summary.ditheredRA.max() > 10.
    elif opsout.opsimVars['angleUnits'] == 'radians':
        assert opsout.summary.ditheredRA.max() < 2.* np.pi + 1.0e-5
