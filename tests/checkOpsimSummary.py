#!/usr/bin/env python

import os
import OpSimSummary as oss
import OpSimSummary.summarize_opsim as so
import pandas as pd
"""
This file reads a couple of lines of opsim outputs and produces a simlib file
which is compared for values of skysig and zptavg against a simlib file produced 
by a different code
"""
example_data =  os.path.join(os.path.split(oss.__file__)[0], 'example_data')
fname = os.path.join(example_data, 'test_opsimObsData.dat')

# fname = 'test_opsimObsData.dat'

opsim_test_df = pd.read_csv(fname, delimiter='\s+')
testopsim = so.SummaryOpsim(opsim_test_df)
sl = testopsim.simlib(519).sort('MJD')
oldsimlib = os.path.join(example_data, 'oldSimlib.simlib')
s_old = oss.simlib.Simlib.fromSimlibFile(oldsimlib)
old_data = s_old.simlibData(519).sort('MJD')

from numpy.testing import assert_allclose
assert_allclose(old_data.SKYSIG.values, sl.simLibSkySig.values, atol=0.01)
assert_allclose(old_data.ZPTAVG.values, sl.simLibZPTAVG.values, atol=0.01)
