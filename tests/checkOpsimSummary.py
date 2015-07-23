#!/usr/bin/env python

import os
import OpSimSummary as oss
import OpSimSummary.summarize_opsim as so
import pandas as pd

example_data =  os.path.join(os.path.split(oss.__file__)[0], 'example_data')
fname = os.path.join(example_data, 'test_opsimObsData.dat')

# fname = 'test_opsimObsData.dat'

opsim_test_df = pd.read_csv(fname, delimiter='\s+')
testopsim = so.SummaryOpsim(opsim_test_df)
sl = testopsim.simlib(519)
print sl.simLibSkySig


