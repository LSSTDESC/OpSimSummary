from __future__ import absolute_import, division, print_function
import opsimsummary as oss
from sqlalchemy import create_engine
import pandas as pd
import time
import os
import numpy as np
from opsimsummary import (OpSimOutput,
                          Simlibs)


logfile = 'feature_sim.simlib.log'
opsim_fname = 'feature_rolling_half_mask_10yrs.db'
simlib_fname = 'feature_rolling_half_mask_10yrs.simlib'
script_start = time.time() 
log_str = 'Running script with opsimsummary version {}\n'.format(oss.__version__)
log_val = 'Starting Calculation at {}\n'.format(script_start)
log_str += log_val

pkgDir = os.path.split(oss.__file__)[0]
dbname = os.path.join('/Users/rbiswas/data/', 'LSST/OpSimData',
                      opsim_fname)
log_val = 'The OpSim DataBase used is {}\n'.format(dbname)
log_str += log_val

# read the database into a `pd.DataFrame`
opsout = OpSimOutput.fromOpSimDB(dbname,
                                 opsimversion='lsstv4',
                                 tableNames=('SummaryAllProps', 'Proposal'))
summary = opsout.summary
log_val = 'dataframe read in from database {}\n'.format(time.time())
log_str += log_val

simlibs = Simlibs(summary, opsimversion='lsstv4', usePointingTree=True)
rng = np.random.RandomState(1)
simlibs.randomSimlibs(numFields=50000, fname=simlib_fname)#, rng=rng)
log_val = 'Done'
log_str += log_val

log_val = 'Writing simlib for input to outfile {0} at time {1}\n'.format(simlib_fname, time.time())
print(log_val)
log_str += log_val

with open(logfile, 'w') as f:
    f.write(log_str)
