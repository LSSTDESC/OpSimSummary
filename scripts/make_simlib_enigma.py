from __future__ import absolute_import, division, print_function
import opsimsummary as oss
import opsimsummary.summarize_opsim as so
from sqlalchemy import create_engine
import pandas as pd
import time
import os

script_start = time.time() 
log_str = 'Running script with opsimsummary version {}\n'.format(oss.__VERSION__)
log_val = 'Starting Calculation at {}\n'.format(script_start)
log_str += log_val

pkgDir = os.path.split(oss.__file__)[0]
dbname = os.path.join('/Users/rbiswas/data/', 'LSST/OpSimData',
                      'enigma_1189_sqlite.db')
log_val = 'The OpSim DataBase used is {}\n'.format(dbname)
log_str += log_val

engineFile = 'sqlite:///' + dbname
engine = create_engine(engineFile)

# read the database into a `pd.DataFrame`
Summary = pd.read_sql_table('Summary', engine)
log_val = 'dataframe read in from database {}\n'.format(time.time())
log_str += log_val

def _writeSimlibFor(propIDList, simlibFileName, description='DDF',
                    log_str= log_str):
    df = Summary.query('propID == @propIDList')
    opSummary = so.SummaryOpsim(df, calculateSNANASimlibs=True,
                                user='rbiswas', host='time')
    log_val = 'The summary has {} entries\n'.format(len(df))
    print(log_val)
    log_str += log_val
    log_val = \
    'The summary has {} unique fields\n'.format(len(df.fieldID.unique()))
    print(log_val)
    log_str += log_val

    log_val = 'Writing simlib for {0} input to outfile {1}\n'.format(description,
            simlibFileName)
    print(log_val)
    log_str += log_val
    opSummary.writeSimlib(simlibFileName)
    log_val = 'Done simlib calculation at {0} and simlib written to {1}\n'.\
            format(time.time(), simlibFileName)
    print(log_val)
    log_str += log_val




WFDSimlib = os.path.join('../opsimsummary/example_data/',
                              'Enigma_1189_WFD.simlib')
DDFSimlib = os.path.join('../opsimsummary/example_data/',
                              'Enigma_1189_DDF.simlib')
CombSimlib = os.path.join('../opsimsummary/example_data/',
                              'Enigma_1189_Combined.simlib')


_writeSimlibFor([366], DDFSimlib, description='DDF')
_writeSimlibFor([364], WFDSimlib, description='WFD')
_writeSimlibFor([364, 366], CombSimlib, description='Combined')
logfile = 'enigma_simlibs.log'

with open(logfile, 'w') as f:
    f.write(log_str)
