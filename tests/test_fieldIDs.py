from __future__ import division, print_function, absolute_import
import opsimsummary as oss
from sqlalchemy import create_engine
import os
import numpy as np
import pandas as pd


def test_fieldIDs():
    """
    """
    pkgDir = os.path.split(oss.__file__)[0]
    dbname = os.path.join(pkgDir, 'example_data', 'enigma_1189_micro.db')
    engineFile = 'sqlite:///' + dbname
    engine = create_engine(engineFile)

    # read the database into a `pd.DataFrame`
    Summary = pd.read_sql_table('Summary', engine)
    fieldIDFromRADec = oss.fieldID(Summary, np.radians(53.), np.radians(-28.))
    assert(fieldIDFromRADec == 1427)
    fieldIDFromRADec = oss.fieldID(Summary, np.radians(0.), np.radians(-45.))
    assert(fieldIDFromRADec == 744)
    fieldIDFromRADec = oss.fieldID(Summary,
                                   np.radians(85.7),
                                   np.radians(-14.4))
    assert(fieldIDFromRADec == 2001)
