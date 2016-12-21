"""
Module dealing with visualizing OpSim results
"""
from __future__ import absolute_import, print_function, division
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon 
import matplotlib.pyplot as plt
from .trig import (pixelsForAng,
                   convertToCelestialCoordinates)
from .healpix import healpix_boundaries                   
from .healpixTiles import HealpixTiles
import healpy as hp

__all__ = ['plot_south_steradian_view', 'HPTileVis']
def plot_south_steradian_view(ra, dec, numPoints=100, radius=1.75, boundary=20.,
                              ax=None, show_frame=True):
    """
    """
    if ax is None:
        fig, ax = plt.subplots()
    m = Basemap(projection='spstere', boundinglat=boundary, lon_0=0., ax=ax)
    if show_frame:
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))
    for ra_val, dec_val in zip(ra, dec):
            m.tissot(ra_val, dec_val, radius, numPoints, ax,
                     **dict(fill=False))

    return fig


class HPTileVis(object):
    """
    Class useful for visualizing a combination of OpSim Pointings and Healpix
    Tiles

    Parameters
    ----------
    hpTile :
    opsout :

    Methods
    -------
    
    """
    def __init__(self, hpTile, opsout):
        """
        """
        self.hpTile = hpTile
        self.opsout = opsout
        self.nside = self.hpTile.nside

    @staticmethod
    def tileIDfromCelestialCoordinates(ra, dec, nside, units='degrees'):
        """
        Parameters
        -----------
        ra : float, mandatory
            ra for the point
        dec : float, mandatory
            dec for the point
        units: {'degrees', 'radians'}
        """
        return pixelsForAng(lon=ra,lat=dec, nside=nside, unit=units)

    def tileCenter(self, tileID):
        theta, phi = hp.pix2ang(self.nside, tileID, nest=True)
        ra, dec = convertToCelestialCoordinates(theta, phi, input_unit='radians',
                                                output_unit='degrees')
        return ra, dec
        
    def pointingSummary(self, tileID=None, ra=None, dec=None,
		       	columns=('ditheredRA', 'ditheredDec'), 
                        allPointings=None):
        obsHistIDs = self.hpTile.pointingSequenceForTile(tileID=tileID,
                                                         allPointings=allPointings)
        return self.opsout.summary.ix[obsHistIDs][list(columns)]
    def pointingCenters(self,
                tileID,
                raCol='ditheredRA',
                decCol='ditheredDec',
                query=None):
        summary = self.pointingSummary(tileID)#, columns=[raCol, decCol])
        
        if query is not None:
            summary = summary.query(query)
        ra = summary[raCol].apply(np.degrees).values
        dec = summary[decCol].apply(np.degrees).values
        return ra, dec
    def plotTilePointings(self, tileID, raCol='ditheredRA', decCol='ditheredDec', radius=1.75,
                          paddingFactors=1, query=None, ax=None, projection='cyl',**kwargs):
        """
        Parameters
        ----------
        
        """
        if ax is None:
            fig, ax = plt.subplots()
        padding = np.degrees(hp.max_pixrad(self.nside)) + radius
        ra_tile, dec_tile = self.tileCenter(tileID)
        
        llcrnrlat = dec_tile - padding * paddingFactors
        urcrnrlat = dec_tile + padding * paddingFactors
        llcrnrlon = ra_tile - padding * paddingFactors
        urcrnrlon = ra_tile + padding * paddingFactors
        
        m = Basemap(llcrnrlat=llcrnrlat, llcrnrlon=llcrnrlon,
                    urcrnrlat=urcrnrlat, urcrnrlon=urcrnrlon,
                    projection=projection, lon_0=ra_tile, lat_0=dec_tile,
                    ax=ax)
        
        parallels = np.linspace(llcrnrlat, urcrnrlat, 3)
        meridians = np.linspace(llcrnrlon, urcrnrlon, 3)
        m.drawparallels(parallels, labels=(1, 0, 0, 0)) #np.ones(len(parallels), dtype=bool))
        m.drawmeridians(meridians, labels=(0, 1, 1, 1)) #np.ones(len(meridians), dtype=bool))
        ra, dec = self.pointingCenters(tileID, raCol=raCol, decCol=decCol, query=query)
        lon, lat = healpix_boundaries(tileID, nside=self.nside, units='degrees',
                                      convention='celestial', step=10,
                                      nest=True)
        x, y = m(lon, lat)
        xy = zip(x, y)
        healpixels = Polygon(xy, facecolor='w',fill=False, alpha=1., edgecolor='k', lw=2)
        for ra, dec in zip(ra, dec):
            m.tissot(ra, dec, radius, 100, **kwargs)
        ax.add_patch(healpixels)
        return fig
