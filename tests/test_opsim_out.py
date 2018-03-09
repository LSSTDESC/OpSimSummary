import os
import pytest
import pandas as pd
from sqlalchemy import create_engine
import opsimsummary as oss
from opsimsummary import OpSimOutput

testdata_propIDDict = [('enigma_1189_micro.db', 'lsstv3', (364, 366)),
                       ('opsimv4_feat_micro.db','lsstv4', (0, 1))]

#@pytest.fixture()
#def opsfileopsv4():
#    yield os.path.join(oss.example_data, 'enigma_1189_micro.db')
#    print('teardown')

@pytest.fixture()
def opsfile(fname):
    yield os.path.join(oss.example_data, fname)
    print('teardown')

#def test_propIDs(opsfile):
#    opsout = OpSimOutput.fromOpSimDB(opsfile, subset='ddf')
#    assert opsout.summary.propID.unique()[0] == 366

@pytest.mark.parametrize("fname,opsimversion,expected", testdata_propIDDict)
def test_get_propIDDict(fname, opsimversion, expected):
    fname = os.path.join(oss.example_data, fname)
    print(fname)
    engine = create_engine('sqlite:///' + fname)
    df = pd.read_sql_table(con=engine, table_name='Proposal')
    mydict = OpSimOutput.get_propIDDict(df, opsimversion=opsimversion)
    assert(mydict['ddf'] == expected[1])
    assert(mydict['wfd'] == expected[0])
