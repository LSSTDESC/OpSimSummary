"""
This script is used to make a standard healpixel database in the
`opsimsummary/example_data`  repository. This database is a coarse grained
NSIDE = 1 healpixelized OpSim created from enigma_1189_micro.db and is used for
testing purposes. 

NOTE: To make this database, it is important to run this from the scripts
directory. The `outfile` variable must be uncommented in the first line of the
code. The `opsimsummary/example_data/healpixels_micro.db` file is meant to be a
standard and not be regenerated. In case this needs to be regenerated (for
example, it is discovered that this file is incorrect, then the setup.py must
be run again to include the new file in the package directories for the tests to
pass.
"""
# outfile = os.path.join('../opsimsummary/example_data', 'healpixels_micro.db')
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

pkgDir = os.path.split(oss.__file__)[0]
dbname = os.path.join(pkgDir, 'example_data', 'enigma_1189_micro.db')
engineFile = 'sqlite:///' + dbname
engine = create_engine(engineFile)

# opsim_hdf = '/Users/rbiswas/data/LSST/OpSimData/minion_1016.hdf'
OpSim_combined = pd.read_sql_query('SELECT * FROM Summary WHERE PropID is 364',
                                    con=engine, index_col='obsHistID')

def addVec(df, raCol='ditheredRA', decCol='ditheredDec'):
    thetas  = - df[decCol] + np.pi /2.
    phis = df[raCol]
    df['vec'] = list(hp.ang2vec(thetas, phis))

addVec(OpSim_combined)
NSIDE = 1
OpSim_combined['hids'] = [query_disc(NSIDE, vec, np.radians(1.75), inclusive=True, nest=True) for vec in OpSim_combined.vec]
# Note this is the less efficient scheme, but the nest scheme is useful later.
lens = map(len, OpSim_combined.hids.values)

rowdata = []
_ = list(rowdata.extend(repeat(i, lens[i])) for i in xrange(len(OpSim_combined)))
obsHistIDs = OpSim_combined.reset_index().ix[rowdata].obsHistID.values
coldata = np.concatenate(OpSim_combined.hids.values)
if os.path.exists(outfile):
    os.remove(outfile)

conn = sqlite3.Connection(outfile)
cur = conn.cursor()
cur.execute('CREATE TABLE simlib (ipix int, obsHistId int)')
tstart = time.time()
told = tstart
for i in range(len(rowdata)):
    cur.execute('INSERT INTO simlib VALUES ({1}, {0})'.format(obsHistIDs[i], coldata[i]))
    if i % 10000000 == 0:
        conn.commit()
        tat = time.time()
        print('committed at {0} taking time {1}'.format(i, tat-told))
        told = tat

conn.commit()
print('Committed the table to disk\n')
# create index
print('Creteing ipix index\n')
cur.execute('CREATE INDEX {ix} on {tn}({cn})'\
                .format(ix='ipix_ind', tn='simlib', cn='ipix'))
print('Creteing obsHistID index\n')
cur.execute('CREATE INDEX {ix} on {tn}({cn})'\
                .format(ix='obshistid_ind', tn='simlib', cn='obsHistId'))
conn.close()

