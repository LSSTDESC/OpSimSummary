import pandas as pd
import os
import OpSimSummary as oss


example_data =  os.path.join(os.path.split(oss.__file__)[0], 'example_data')
fname = os.path.join(example_data, 'test_opsimObsData.dat')
old_simlibfilename = os.path.join(example_data, 'oldSimlib.simlib')
def test_read():
    """
    check that this simlib file exists
    """
    assert os.path.exists(fname) 
def test_shape_is_expected():
    """
    test that the file is being ingested correctly to the extetnt that the
    number of lines matches the expected number in the known file.
    """
    opsim_test_df = pd.read_csv(fname, delimiter='\s+')
    assert len(opsim_test_df) == 2

def test_readSimlib():
    print('hello')
    oldsimlib = oss.simlib.Simlib.fromSimlibFile(old_simlibfilename)
    print(type(oldsimlib))
    print (oldsimlib.meta)
    print dir(oldsimlib)

