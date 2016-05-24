from __future__ import absolute_import, division, print_function
import opsimsummary as oss
import opsimsummary.summarize_opsim as so
from sqlalchemy import create_engine
import pandas as pd
import os


def test_writeSimlib():
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
    if new_data == template_data :
        os.remove(simlibfilename)

if __name__ == '__main__':
    test_writeSimlib()
