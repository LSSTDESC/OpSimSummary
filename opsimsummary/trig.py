import numpy as np

__all__ = ['fieldID', 'angSep', 'overlapSummary',
           'convertToSphericalCoordinates']

def convertToSphericalCoordinates(ra, dec, unit='degrees'):
    """
    Convert ra, dec coordinates to usual spherical coordinates theta, phi.
    The ra, dec coordinates have the ranges (0., 2\pi), (-\pi/2.0, \pi/2.0)

    Parameters
    ----------
    ra : sequence of floats or scalar float, mandatory
        ra coordinates in units of degrees or radians
    dec : sequence of floats or scalar float, mandatory
        dec coordinates in units of degrees or radians
    unit : string, optional, defaults to 'degrees'
        'degrees' or 'radians'

    Returns
    -------
    tuple of `numpy.ndarray` of theta and phi of spherical coordinates in
    radians. If scalars were supplied the return is a tuple of two scalars.

    .. note:: The units of ra, dec must be the same. They should have
    consistent lengths
    """
    # Check that unit type has been implemented
    if unit not in ['degrees', 'radians']:
        raise ValueError('Parameter unit to convertToSphericalCoordinates must'
                         ' be degrees or radians\n')

    ra = np.ravel(ra)
    dec = np.ravel(dec)

    # Check consistency of inputs
    if len(ra) != len(dec):
        raise ValueError('The sequences for ra and dec must have equal lengths'
                         'but have lengths {0} and {1}\n'.format(len(ra),
                                                                 len(dec)))
    if unit == 'degrees':
        ra = np.radians(ra)
        dec = np.radians(dec)

    theta = - dec + np.pi / 2.0
    phi = ra

    return theta, phi


def angSep(ra1, dec1, ra2, dec2):
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

def overlapSummary(ra, dec, df, sep=1.75, 
                    raCol='fieldRA', decCol='fieldDec'):

    sep = np.radians(sep)

    df['overlap'] = angSep(ra, dec, df[raCol], df[decCol]) < sep
    summary = df[df['overlap']].copy()
    summary.drop('dist', axis=1, inplace=True)
    return summary
def angDistance(ra, dec, df, raCol='fieldRA', decCol='fieldDec'):
    """
    """
    df['dist'] = angSep(ra, dec, df[raCol], df[decCol])
    idx = df.dist.idxmin()
    rval = df.ix[idx]
    df.drop('dist', axis=1, inplace=True)
    return rval

def obsIndex(opsimDF, ra, dec, raCol='fieldRA', decCol='fieldDec',
              pointinRadius=1.75):
    """
    """
    df = opsimDF[['fieldID', raCol, decCol]].copy().drop_duplicates()
    df['dist'] = angSep(ra, dec, df[raCol], df[decCol])
    pointingRad = np.radians(pointinRadius)
    df = df.query('dist < {}'.format(pointingRad))
    return df.index



def fieldID(opsimDF, ra, dec, raCol='fieldRA', decCol='fieldDec'):
    """
    """

    df = opsimDF[['fieldID', raCol, decCol]].copy().drop_duplicates()

    rval = angDistance(ra, dec, df)
    return rval['fieldID']

