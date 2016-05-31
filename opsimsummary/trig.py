import numpy as np
import pandas as pd

__all__ = ['fieldID']


def angSep( ra1, dec1, ra2, dec2):
    """
    Angular separation between to points on a sphere with coordinates ra, dec
    in the usual conditions for astronomy.
    """
    th1 = - dec1 + np.pi / 2.
    th2 = - dec2 + np.pi / 2.
    ph1 = ra1
    ph2 = ra2
    # Take the dot product
    cos = np.sin(th1) * np.sin(th2) * np.cos(ph1 - ph2)\
        + np.cos(th1)*np.cos(th2)
    # Calculate arccos differently to avoid loss of precision
    # sinhalftheta = np.sqrt(0.5*(1.0 - cos))
    # return 2.*(np.arcsin(sinhalftheta))
    return np.arccos(cos)


def angDistance(ra, dec, df, raCol='fieldRA', decCol='fieldDec'):
    """
    """
    df['dist'] = angSep(ra, dec, df[raCol], df[decCol])
    idx = df.dist.idxmin()
    rval = df.ix[idx]
    df.drop('dist', axis=1, inplace=True)
    return rval


def fieldID(opsimDF, ra, dec, raCol='fieldRA', decCol='fieldDec'):
    """
    """

    df = opsimDF[['fieldID', raCol, decCol]].copy().drop_duplicates()

    rval = angDistance(ra, dec, df)
    return rval['fieldID']

