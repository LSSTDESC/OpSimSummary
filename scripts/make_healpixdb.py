from __future__ import division
import numpy as np
import time
import sqlite3
import healpy as hp
from healpy import query_disc, query_polygon
import opsimsummary as oss
import pandas as pd
from itertools import repeat

opsim_hdf = '/Users/rbiswas/data/LSST/OpSimData/minion_1016.hdf'
OpSim_combined = pd.read_hdf(opsim_hdf, 'Table')

def addVec(df, raCol='ditheredRA', decCol='ditheredDec'):
    thetas  = - df[decCol] + np.pi /2.
    phis = df[raCol]
    df['vec'] = list(hp.ang2vec(thetas, phis))

addVec(OpSim_combined)
NSIDE = 256
OpSim_combined['hids'] = [query_disc(NSIDE, vec, np.radians(1.75), inclusive=True, nest=True) for vec in OpSim_combined.vec]
# Note this is the less efficient scheme, but the nest scheme is useful later.
lens = map(len, OpSim_combined.hids.values)

rowdata = []
_ = list(rowdata.extend(repeat(i, lens[i])) for i in xrange(len(OpSim_combined)))
coldata = np.concatenate(OpSim_combined.hids.values)

conn = sqlite3.Connection('healpixels.db')
cur = conn.cursor()
cur.execute('CREATE TABLE simlib (ipix int, obsHistId int)')
tstart = time.time()
told = tstart
for i in range(len(rowdata)):
    cur.execute('INSERT INTO simlib VALUES ({0}, {1})'.format(rowdata[i], coldata[i]))
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

