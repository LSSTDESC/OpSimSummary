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
        Obtain the tile id for a Healpix tile in the nested convention from
        the ra and dec in degrees
        Parameters
        -----------
        ra : float, mandatory
            ra for the point
        dec : float, mandatory
            dec for the point
        units: {'degrees', 'radians'}
        """
        return pixelsForAng(lon=ra,lat=dec, nside=nside, unit=units)

    def tileCenter(self, ipix):
        """
        return a tuple of arrays returning the ra, dec in degrees of center of
        the healpix tile with id tileID in the nest convention for NSIDE given
        by `self.nside`

        Parameters
        ----------
        ipix : 
        """
        theta, phi = hp.pix2ang(self.nside, ipix, nest=True)
        ra, dec = convertToCelestialCoordinates(theta, phi, input_unit='radians',
                                                output_unit='degrees')
        return ra, dec

    def pointingSummary(self,
                        tileID=None,
                        ra=None,
                        dec=None,
                        angularUnits='degrees',
                        allPointings=None):
        """
        returns the OpSim Summary table corresponding to the maximal set of
        pointings intersecting with the Healpix tile with tileID.
        Parameters
        -----------
        ra : float, mandatory
            ra for the point
        dec : float, mandatory
            dec for the point
        units: {'degrees', 'radians'}
        """
        if tileID is None:
            if ra is None or dec is None:
                raise ValueError("Both tileID and ra, dec cannot be None")

            # Set tileID from ra, dec
            tileID = self.tileIDfromCelestialCoordinates(ra,
                                                         dec,
                                                         self.nside,
                                                         units=angularUnits)

        obsHistIDs = self.hpTile.pointingSequenceForTile(tileID=tileID,
                                                         allPointings=allPointings)
        s = self.opsout.summary.ix[obsHistIDs]
        return s

    def pointingCenters(self,
                tileID,
                raCol='ditheredRA',
                decCol='ditheredDec',
                query=None,
                groupby=None):
        """
        """
        summary = self.pointingSummary(tileID)#, columns=[raCol, decCol])
        
        if query is not None:
            summary = summary.query(query)
        ra = summary[raCol].apply(np.degrees).values
        dec = summary[decCol].apply(np.degrees).values
        return ra, dec

    def plotTilePointings(self,
                          tileID,
                          raCol='ditheredRA', decCol='ditheredDec',
                          radius=1.75,
                          paddingFactors=1,
                          query=None,
                          ax=None,
                          projection='cyl',
                          tile_centers=None,
                          corners=None,
                          **kwargs):
        """
        Plot the Healpix Tile and the maximal set of pointings overlapping with
        it.
        Parameters
        ----------
        tileID : int, mandatory
            Healpix tileID in the nested scheme
        raCol : string, defaults to 'ditheredRA'
            column name with ra that should be used 
        decCol : string, defaults to 'ditheredDec'
            column name with dec that should be used 
        radius : float, defaults to 1.75, degrees
            radius of the field of view
        paddingFactors: float, defaults to 1.0
            controls the size of the figure wrt angular dimensions of the tile.
            paddingFactors=1 is designed to get all the centers of pointings
            overlapping the tile inside the figure
        query : string
            query for the pandas dataframe to select only some of the observations
        ax : instance of matplotlib.figure.axes 
            axes for figure. New figure created insitue if not provided
        projections: string, defaults to `cyl`
            string to specify the `Basemap.projection` 
        tile_centers : tuples of 2 floats, degrees, defaults to None
            (ra of the center of the figure, dec of the center of the figure)
            in degrees
            if None, the center is set by using the center of the Healpixel of
            tileID and `self.nside`
        corners: tuple of floats, defaults to None 
        kwargs: passed to the plotting of the the field of veiw

        Returns
        -------
        fig, tile_centers, corners

        ..notes: Have not thought about wrapover
        """
        # if axes object not provided setup figure
        if ax is None:
            fig, ax = plt.subplots()

        # size of box
        padding = np.degrees(hp.max_pixrad(self.nside)) + radius

        # center of the figure
        if tile_centers is None:
            ra_tile = self.tileCenter(tileID)[0][0]
            dec_tile = self.tileCenter(tileID)[1][0]
        else:
            ra_tile, dec_tile = tile_centers
        
        # corner of the figure
        if corners is not None:
            llcrnrlat, llcrnrlon, urcrnrlat, urcrnrlon = corners
        else:
            llcrnrlat = dec_tile - padding * paddingFactors
            urcrnrlat = dec_tile + padding * paddingFactors
            llcrnrlon = ra_tile - padding * paddingFactors
            urcrnrlon = ra_tile + padding * paddingFactors
        
        # Instantiate basemap
        m = Basemap(llcrnrlat=llcrnrlat, llcrnrlon=llcrnrlon,
                    urcrnrlat=urcrnrlat, urcrnrlon=urcrnrlon,
                    projection=projection, lon_0=ra_tile, lat_0=dec_tile,
                    ax=ax)

        # Draw some parallels and meridians to get a spatial sense        
        parallels = np.linspace(llcrnrlat, urcrnrlat, 3)
        meridians = np.linspace(llcrnrlon, urcrnrlon, 3)
        m.drawparallels(parallels, labels=(1, 0, 0, 0)) #np.ones(len(parallels), dtype=bool))
        m.drawmeridians(meridians, labels=(0, 1, 1, 1)) #np.ones(len(meridians), dtype=bool))

        # obtain the boundaries of the Healpixels and create the polygon patch
        lon, lat = healpix_boundaries(tileID, nside=self.nside, units='degrees',
                                      convention='celestial', step=10,
                                      nest=True)
        x, y = m(lon, lat)
        xy = zip(x, y)
        healpixels = Polygon(xy, facecolor='w', fill=False, alpha=1.,
                             edgecolor='k', lw=2)

        # Obtain the centers of pointings
        ra, dec = self.pointingCenters(tileID, raCol=raCol, decCol=decCol, query=query)
        for ra, dec in zip(ra, dec):
            m.tissot(ra, dec, radius, 100, **kwargs)

        # Draw patch
        ax.add_patch(healpixels)

        # really important for the case where ax is None
        if ax is None:
            fig = ax.figure

        return fig, (ra_tile, dec_tile), (llcrnrlat, llcrnrlon, urcrnrlat, urcrnrlon)
