"""
This script is used to make a standard healpixel database in the
`opsimsummary/example_data`  repository. This database is a coarse grained
NSIDE = 256 healpixelized OpSim created from enigma_1189_micro.db and is used for
testing purposes. 

NOTE: To make this database, it is important to run this from the scripts
directory. The `outfile` variable must be uncommented in the first line of the
code. The `opsimsummary/example_data/healpixels_micro.db` file is meant to be a
standard and not be regenerated. In case this needs to be regenerated (for
example, it is discovered that this file is incorrect, then the setup.py must
be run again to include the new file in the package directories for the tests to
pass.
"""
from __future__ import division
import numpy as np
import time
import sqlite3
import healpy as hp
from healpy import query_disc, query_polygon
import opsimsummary as oss
import pandas as pd
from itertools import repeat
import os
from sqlalchemy import create_engine
import opsimsummary as oss
import gzip
import shutil

outfile = os.path.join('../opsimsummary/example_data', 'healpixels_micro.db')
if os.path.exists(outfile):
    raise ValueError('output file already exists. If you want to overwrite, '
                     'rerun script after doing \n\nrm {}'\
                     .format(os.path.abspath(outfile)))
pkgDir = os.path.split(oss.__file__)[0]
dbname = os.path.join(pkgDir, 'example_data', 'enigma_1189_micro.db')
# engineFile = 'sqlite:///' + dbname
# engine = create_engine(engineFile)

# opsim_hdf = '/Users/rbiswas/data/LSST/OpSimData/minion_1016.hdf'
# OpSim_combined = pd.read_sql_query('SELECT * FROM Summary WHERE PropID is 364',
#                                    con=engine, index_col='obsHistID')

ho = oss.HealPixelizedOpSim.fromOpSimDB(opSimDBpath=dbname, subset='wfd',
                                        propIDs=None, NSIDE=4,
                                        raCol='ditheredRA',
                                        decCol='ditheredDec',
                                        vecColName='vec',
                                        fieldRadius=1.75)

ho.writeToDB(outfile)
with open(outfile, 'rb') as f_in, gzip.open(outfile + '.gz', 'wb') as f_out:
    shutil.copyfileobj(f_in, f_out)
