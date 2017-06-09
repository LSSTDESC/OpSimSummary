"""
Module for trigonometric utilities such as conversion between different
conventions and units, as well as calculating angular separations between
points on a sphere.
"""
__all__ = ['fieldID', 'angSep', 'overlapSummary',
           'convertToSphericalCoordinates',
           'convertToCelestialCoordinates',
           'angToVec', 'pixelsForAng', 'pixelsToAng']

import numpy as np
import healpy as hp


def pixelsForAng(lon, lat, nside, nest=True, convention='celestial', unit='degrees'):
    """
    Return the array of healpixels in which the set of points represented by
    the arrays lon, and lat fall
    """
    vecs = angToVec(lon, lat, convention, unit)
    return hp.vec2pix(nside, vecs[:, 0], vecs[:, 1], vecs[:, 2], nest=nest)


def pixelsToAng(tileID, nside, nest=True, convention='celestial', unit='degrees'):
    """
    Return the array of ra, dec for the healpix tile center for the tile with
    pixel ID tileID, 
    the arrays lon, and lat fall
    """
    theta, phi = hp.pix2ang(nside=nside, ipix=tileID, nest=nest)

    if convention == 'celestial':
        ra, dec = convertToCelestialCoordinates(theta, phi, input_unit='radians',
                                                output_unit=unit)
    else:
        raise NotImplementedError('convention other than celestial not implemented')
    return ra, dec


def angToVec(lon, lat, convention, unit):
    """
    Compute an array of unit vectors corresponding to the array of longitude
    and latitude in either celestial conventions (ra, dec) or usual spherical
    coordinate system conventions.

    Parameters
    ----------
    lon : scalar or `numpy.ndarray` of floats
        co-longitude which may be ra, or phi depending on the convention
    lat : scalar or `numpy.ndarray` of floats
        co-latitude which may be dec, or theta depending on the convention
    convention : {'celestial'|'spherical'}
        whether the input angular positions are in spherical or celestial
        coordinates.
    unit : string
        {'degrees'| 'radians'}
    """
    if unit.lower() not in ('degrees', 'radians'):
        raise ValueError('illegal value for units, '
                         'must be either degrees or radians')
    if convention.lower() not in ('celestial', 'spherical'):
        raise ValueError('illegal value for convention, '
                         'must be either spherical or celestial')
    lon = np.ravel(lon)
    lat = np.ravel(lat)

    if convention.lower() == 'celestial':
        theta, phi = convertToSphericalCoordinates(lon, lat, unit=unit.lower())
    if convention.lower() == 'spherical' and unit.lower() == 'degrees':
        theta = np.radians(lat)
        phi = np.radians(lon)

    return hp.ang2vec(theta, phi)

def convertToCelestialCoordinates(theta, phi, input_unit='radians',
                                  output_unit='degrees'):
    """
    Convert theta, phi in usual spherical coordinates to ra, dec
    values. 

    Parameters
    ----------
    theta : `np.ndarray` or float
        co-latitude
    phi :`np.ndarray` or float
        co longitude 
    input_unit : {'degrees', 'radians'}, defautls to radians
    output_unit : {'degrees', 'radians'}, defaults to degrees

    Return
    ------
    tuple of (ra, dec) where ra, dec are `numpy.ndarray` of type float
    of length equal to theta or phi
    """
    if input_unit != 'degrees':
        theta = np.degrees(np.ravel(theta))
        phi = np.degrees(np.ravel(phi))
    dec = - theta + 90.0
    ra = phi
    if output_unit != 'degrees':
        ra = np.radians(phi)
        dec = np.radians(dec)
    return ra, dec

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
    tuple of `numpy.ndarray` (theta, phi) of spherical coordinates in
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

