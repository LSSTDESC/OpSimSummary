from __future__ import absolute_import, division, print_function
import opsimsummary as oss
import opsimsummary.summarize_opsim as so
from sqlalchemy import create_engine
import pandas as pd

import os

pkgDir = os.path.split(oss.__file__)[0]
dbname = os.path.join(pkgDir, 'example_data', 'enigma_1189_micro.db')
print(dbname)

engineFile = 'sqlite:///' + dbname
engine = create_engine(engineFile)

# read the database into a `pd.DataFrame`
Summary = pd.read_sql_table('Summary', engine)

EnigmaMain = Summary.query('propID == [364]')
EnigmaMainSummary = so.SummaryOpsim(EnigmaMain, calculateSNANASimlibs=True,
                                    user='rbiswas', host='time')
simlibfilename = os.path.join('../opsimsummary/example_data/',
                              'Enigma_1189_micro_main.simlib') 
EnigmaMainSummary.writeSimlib(simlibfilename)
