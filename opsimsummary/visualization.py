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
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
import pyproj
import healpy as hp

__all__ = ['plot_south_steradian_view', 'HPTileVis', 'split_PolygonSegments',
           'AllSkyMap']


def split_PolygonSegments(lon, lat, lon_split=180., epsilon = 0.000001):
    """split spherical segments of colatitude and colongitude that constitute
    a polygon edge if they go across the edge of a projection map. 
    
    Parameters
    ----------
    lon : `np.array` of floats, degrees
        colongitude of points on the borders of a polygon
    lat : `np.array` of floats, degrees
        colatitude of points on the border of a polygon
    lon_split : float in degrees, defaults to 180
        line along which the polygons are to be split
    epsilon : float, defaults to 1.0e-6
        to control overflows
    
    Returns
    --------
    list of tuples, with each of the tuples containing the segments for a polygon
    after splitting. Each such tuple is tuple(lon, lat) for a single polygon
    
    Examples
    --------
    >>> lon, lat = healpix_boundaries([171], nside=8, step=10, convention='celestial)
    >>> split_segs(lon, lat)
    [(array([ 180.,  180.,  180.,  180.,  180.,  180.,  180.,  180.,  180., 
             180.,  180.]),
      array([54.3409123 ,  53.72610365,  53.11021259,  52.49321549,
             51.87508838,  51.25580695,  50.63534652,  50.01368203,
             49.39078804,  48.76663871,  48.14120779])),
     (array([ 181.26760563,  182.5       ,  183.69863014,  184.86486486,
              186.        ,  187.10526316,  188.18181818,  189.23076923,
              190.25316456,  191.25      ,  191.39240506,  191.53846154,
              191.68831169,  191.84210526,  192.        ,  192.16216216,
              192.32876712,  192.5       ,  192.67605634,  192.85714286,
              191.73913043,  190.58823529,  189.40298507,  188.18181818,
              186.92307692,  185.625     ,  184.28571429,  182.90322581,
              181.47540984]),
       array([ 47.51446861,  46.88639405,  46.25695656,  45.62612809,
          44.99388015,  44.36018373,  43.72500931,  43.08832685,
          42.45010577,  41.8103149 ,  42.45010577,  43.08832685,
          43.72500931,  44.36018373,  44.99388015,  45.62612809,
          46.25695656,  46.88639405,  47.51446861,  48.14120779,
          48.76663871,  49.39078804,  50.01368203,  50.63534652,
          51.25580695,  51.87508838,  52.49321549,  53.11021259,  53.72610365]))
    ]
    """
    lon = np.asarray(lon)
    lat = np.asarray(lat)
    mask = lon <= lon_split + epsilon
    lon0 = lon[mask]
    lat0 = lat[mask]
            
    lon1 = lon[~mask]
    lat1 = lat[~mask]
    return list(x for x in ((lon0, lat0), (lon1, lat1)))


class AllSkyMap(Basemap):
    
    def __init__(self, *args, **kwargs):
        Basemap.__init__(self, *args, **kwargs)
        
    def tissot(self, lon_0, lat_0,radius_deg,npts,ax=None, epsilon=1.0e-6, **kwargs):
        """
        Draw a polygon centered at ``lon_0,lat_0``.  The polygon
        approximates a circle on the surface of the earth with radius
        ``radius_deg`` degrees latitude along longitude ``lon_0``,
        made up of ``npts`` vertices.
        The polygon represents a Tissot's indicatrix
        (http://en.wikipedia.org/wiki/Tissot's_Indicatrix),
        which when drawn on a map shows the distortion
        inherent in the map projection.
        .. note::
         Cannot handle situations in which the polygon intersects
         the edge of the map projection domain, and then re-enters the domain.
        Extra keyword ``ax`` can be used to override the default axis instance.
        Other \**kwargs passed on to matplotlib.patches.Polygon.
        returns a matplotlib.patches.Polygon object."""
        ax = kwargs.pop('ax', None) or self._check_ax()
        g = pyproj.Geod(a=self.rmajor,b=self.rminor)
        az12,az21,dist = g.inv(lon_0,lat_0,lon_0,lat_0+radius_deg)
        seg = [self(lon_0,lat_0+radius_deg)]
        delaz = 360./npts
        az = az12
        for n in range(npts):
            az = az+delaz
            lon, lat, az21 = g.fwd(lon_0, lat_0, az, dist)
            seg.append((lon, lat))
        ra, dec = zip(*seg)
        # go to 0, 360.
        ra = np.asarray(ra)
        dec = np.asarray(dec)
        ra[ra < 0] += 360.
        #if np.any(np.abs(ra - lon_0) < radius_deg - epsilon):
        polyseg_list = split_PolygonSegments(ra, dec, lon_split=self.lonmax)
        #else:
        #    polyseg_list = list((ra, dec))
        split_poly_normed = list()
        for radec in polyseg_list:
            if len(radec) > 0:
                ra, dec = radec
                ra = np.asarray(ra)
                mask = ra < self.lonmax + 1.0e-6
                ra[~mask] -= 360.
                split_poly_normed.append((ra, dec))
        pts = list()
        
        for poly_segs in split_poly_normed:
            lon, lat = poly_segs
            if len(lon) > 0:
                pt = self.polygonize(lon, lat, **kwargs)
                ax.add_patch(pt)
                pts.append(pt)
        return split_poly_normed
    
    def polygonize(self, lon, lat, **kwargs):
        x, y = self(lon, lat)
        x = np.asarray(x)
        y = np.asarray(y)
        mask = np.logical_and(x < 1.0e20,  y < 1.0e20)
        xy = zip(x[mask], y[mask])
        pt = Polygon(xy, **kwargs)
        return pt


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
    hpTile : instance of `HealpixTiles`
    opsout : instance of `OpSimOutput`

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
        nside : int, mandatory
            NSIDE for `self.hpTile`
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
        ipix : int
            integer for healpix Tile whose center is requested
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
        tileID : int
            index for healpix tile
        ra : float, mandatory
            ra for the point
        dec : float, mandatory
            dec for the point
        angularUnits: {'degrees', 'radians'}
        allPointings : sequence of int
            obsHistIDs to consider 
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
                opsimAngularUnits='radians',
                outputAngularUnits='degrees',
                groupby=None):
        """
        tileID : int
            index for healpix tile
        raCol : string, defaults to `ditheredRA`
            column in an OpSim output to be used for center
            of pointing in terms of ra
        decCol : string, defaults to `ditheredRA`
            column in an OpSim output to be used for center
            of pointing in terms of dec
        query : string, defaults to None
            if not None, passed as a pandas query to restrict the search size
            in plotting.
        groupby : Not implemented
        """
        summary = self.pointingSummary(tileID)#, columns=[raCol, decCol])
        
        if query is not None:
            summary = summary.query(query)
        if opsimAngularUnits !='radians':
            raise ValueError('not implemented')

        if outputAngularUnits == 'degrees':
            ra = summary[raCol].apply(np.degrees).values
            dec = summary[decCol].apply(np.degrees).values
        else:
            raise ValueError('Notimplemented')
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
                          drawPointings=True,
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
        drawPointings : Bool, defaults to True
            if False, draws only Tiles

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
        m.drawparallels(parallels, labels=(1, 0, 0, 0))
        m.drawmeridians(meridians, labels=(0, 1, 1, 1))

        # obtain the boundaries of the Healpixels and create the polygon patch
        lon, lat = healpix_boundaries(tileID, nside=self.nside, units='degrees',
                                      convention='celestial', step=10,
                                      nest=True)
        x, y = m(lon, lat)
        xy = zip(x, y)
        healpixels = Polygon(xy, facecolor='w', fill=False, alpha=1.,
                             edgecolor='k', lw=2)

        if drawPointings:
            # Obtain the centers of pointings
            ra, dec = self.pointingCenters(tileID, raCol=raCol, decCol=decCol, query=query)
            for ra, dec in zip(ra, dec):
                m.tissot(ra, dec, radius, 100, **kwargs)

        # Draw patch
        ax.add_patch(healpixels)

        # really important for the case where ax is None
        if ax is not None:
            fig = ax.figure

        return fig, (ra_tile, dec_tile), (llcrnrlat, llcrnrlon, urcrnrlat, urcrnrlon)
