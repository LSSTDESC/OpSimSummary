from __future__ import print_function, absolute_import, division

import numpy as np
import healpy as hp
import pandas as pd
from scipy.sparse import csr_matrix
from itertools import repeat 

__all__  = ['addVec', 'HealPixelizedOpSim']

def addVec(df, raCol='ditheredRA', decCol='ditheredDec'):
    """
    Add a column of vectors to the dataFrame df

    Parameters
    ----------
    df : `pd.DataFrame`
        dataframe with two columns raCol and decCol having ra,
        dec in radians.
    raCol : string, optional, defaults to 'ditheredRA'
        column name for the column which has ra values
    decCol : string, optional, defaults to 'ditheredDec'
        column name for the column which has dec values
    """
    thetas  = - df[decCol] + np.pi /2.
    phis = df[raCol]
    df['vec'] = list(hp.ang2vec(thetas, phis))

class HealPixelizedOpSim(object):

    def __init__(self, opsimDF, raCol='ditheredRA', decCol='ditheredDec',
                 NSIDE=1,
                 vecColName='vec',  fieldRadius=1.75):

        self.raCol = raCol
        self.decCol = decCol
        self.opsimdf = opsimDF.reset_index()
        self.cols = opsimDF.columns
        self.vecColName = vecColName
        self.validateDF()
        self._fieldRadius = np.radians(fieldRadius)

        if vecColName not in self.cols:
            addVec(self.opsimdf)
        self.nside = NSIDE
        self.nested = True
        self._rowdata = None
        self._coldata = None
        self._spmat = None

    def obsHistIdsForTile(self, tileID):
        """
        """
        inds = self.sparseMat.getcol(tileID).nonzero()[0]
        return self.opsimdf.ix[inds, 'obsHistID'].values

    @property
    def sparseMat(self):
        if self._spmat is None:
            shape=(len(self.opsimdf), hp.nside2npix(self.nside))
            # numpy ones like number of intersections
            ones_int = np.ones(len(self.rowdata))
            self._spmat = csr_matrix((ones_int, (self.rowdata, self.coldata)),
                                     shape=shape)
        return self._spmat 

    @property
    def rowdata(self):
        if self._rowdata is None:
            self.doPreCalcs()
        return np.asarray(self._rowdata)


    @property
    def coldata(self):
        if self._coldata is None:
            self.doPreCalcs()
        return self._coldata

    def doPreCalcs(self):
        self.opsimdf['hids'] = [hp.query_disc(self.nside,
                                              vec,
                                              self._fieldRadius,
                                              inclusive=True,
                                              nest=True)
                                for vec in self.opsimdf[self.vecColName]]
        lens = map(len, self.opsimdf.hids.values)
        rowdata = []
        _ = list(rowdata.extend(repeat(i, lens[i]))
                 for i in xrange(len(self.opsimdf)))
        coldata = np.concatenate(self.opsimdf.hids.values)
        self._rowdata = rowdata
        self._coldata = coldata

    def validateDF(self):
        if self.raCol not in self.cols:
            raise ValueError('raCol {} not in OpSim cols'.format(self.raCol))
        if self.decCol not in self.cols:
            raise ValueError('decCol {} not in OpSim cols'.format(self.decCol))
        return

