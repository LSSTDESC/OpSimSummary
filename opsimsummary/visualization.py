"""
Module dealing with visualizing OpSim results
"""
from __future__ import absolute_import, print_function, division
import abc
from future.utils import with_metaclass
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
import pyproj
import healpy as hp
from .trig import (pixelsForAng,
                   convertToCelestialCoordinates)
from .healpix import healpix_boundaries

from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches

__all__ = ['plot_south_steradian_view', 'HPTileVis', 'split_PolygonSegments',
           'AllSkyMap', 'ObsVisualization', 'AllSkySNVisualization',
           'MilkyWayExtension']

class HandlerEllipse(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = mpatches.Ellipse(xy=center, width=width + xdescent,
                             height=height + ydescent)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]


class MilkyWayExtension(object):
    """Class to describe the part of the sky obscured by the MW brightness

    Attributes
    ----------
    m :
    ax : defaults to None
    color :
    alpha :
    rmat :
    mw_polygon :

    .. note:: Based on code from @ufeind, and Yannick Copin (ycopin@ipnl.in2p3.fr)
    Original warning from authors (This routine is only roughly accurate,
    probably at the arcsec level, and therefore not to be used for astrometric
    purposes. For most accurate conversion, use dedicated
    `kapteyn.celestial.sky2sky` routine.) We mostly use the rotation matrix as
    default rotation matrix and invert it..
    - http://www.dur.ac.uk/physics.astrolab/py_source/conv.py_source
    - Rotation matrix from
      http://www.astro.rug.nl/software/kapteyn/celestialbackground.html
    """
    def __init__(self, m, ax=None, color='y', alpha=1.0,
                 edgecolor='y', lw=0., fill=True, zorder=15,
                 rmat=None, min_lon=-30.):
        self.color = color
        self.alpha = alpha
        self.edgecolor = edgecolor
        self.fill = fill
        self.lw=lw
        self.ax = ax
        self.m = m
        self.zorder = zorder
        self.min_lon = min_lon
        rmat_def = np.array([[-0.054875539396, -0.873437104728, -0.48383499177],
                             [0.494109453628, -0.444829594298, 0.7469822487],
                             [-0.867666135683, -0.198076389613, 0.455983794521]])
        if rmat is None:
            rmat = rmat_def
        self.rmat_inv = np.linalg.inv(rmat)

    def gc2radec(self, gc_phi, gc_theta, min_lon=None):
        """
	convert coordinates in galactic coordinate system to equatorial
        ra, dec.

	Parameters
	----------
        gc_phi : `np.ndarray`, float
            angular coordinate in Galactic coordinate system, where the MW
            disk is at theta=0.
        gc_theta : `np.ndarray`, float 
       	    azimuthal coordinate in Galactic coordinate system.
        min_lon : degrees, defaults to None
            min longitude beyond which the galaxy is cut off. If None,
            `self.min_lon` is used.
        max_lon : None, Not implemented yet	
        """
        vec = hp.ang2vec(gc_theta, gc_phi)
        vec_radec = np.asarray(list(np.dot(self.rmat_inv, v) for v in vec))
        theta, phi = hp.vec2ang(vec_radec)
        ra, dec = convertToCelestialCoordinates(theta, phi)
        if min_lon is None:
            min_lon = self.min_lon
        mask  = dec > min_lon
        return ra[mask], dec[mask]

    @property
    def mw_boundaries(self):
        """tuple of points in ra, dec in degrees wihich show the +20, -20 degree
           boundariess of the MW disk. The return value is
           (ra_h, dec_h), (ra_l, dec_l) with all of the parameters in degrees.
        """
        phi = np.arange(0., 2.0*np.pi, 0.1)
        theta_l = np.ones_like(phi)* 110 * np.pi / 180.
        theta_h = np.ones_like(phi)* 70 * np.pi / 180.
        ra_l, dec_l = self.gc2radec(phi, theta_l)
        ra_h, dec_h = self.gc2radec(phi, theta_h)
        return (ra_h, dec_h), (ra_l, dec_l)

    @property
    def mw_polygon(self):
        (ra_h, dec_h), (ra_l, dec_l) = self.mw_boundaries
        x_l, y_l = self.m(ra_l, dec_l)
        x_h, y_h = self.m(np.roll(ra_h, 3), np.roll(dec_h, 3))
        x, y = (np.concatenate((x_l, x_h[::-1])),
                np.concatenate((y_l, y_h[::-1])))
        p = Polygon(zip(x, y),
                    color=self.color,
                    fill=self.fill,
                    alpha=self.alpha,
                    lw=self.lw,
                    edgecolor=self.edgecolor,
                    zorder=self.zorder)
        return p

    def add_polygons(self, ax):
        p = self.mw_polygon
        _ = self.ax.add_patch(p)
        return self.ax

class ObsVisualization(with_metaclass(abc.ABCMeta, object)):
    """Abstract base class for describing simulations"""
    def __init__(self):
        pass

    @abc.abstractproperty
    def colorCodeRedshifts(self):
        pass

    @abc.abstractproperty
    def show_var_scatter(self):
        pass

    @abc.abstractproperty
    def show_mw(self):
        pass

    @abc.abstractproperty
    def show_visible_fields(self):
        pass

    @abc.abstractproperty
    def radius_deg(self):
        pass

    @abc.abstractproperty
    def band_color_dict(self):
        pass

    @abc.abstractmethod
    def generate_image_bg(self):
        pass

    @abc.abstractmethod
    def generate_mw_polygons(self):
        pass

    @abc.abstractmethod
    def generate_camera(self):
        pass

    @abc.abstractmethod
    def get_visible_field_polygons(self):
        """get a set of polygons to represent the visible fields. The minimal
           input is a time specified as mjd, and a location (which is fixed for
           a ground based observatory. The return values is a list of `Polygons`
        """
        pass

    @abc.abstractmethod
    def generate_var_scatter(self):
        """should generate ra, dec, mag and point size for SN at the time of
        observation"""
        pass

    @abc.abstractmethod
    def generate_image(self):
        """Use methods above to create an image of the sky and optionally save
        it
        """
        pass

    @abc.abstractmethod
    def label_time_image(self, mjd, surveystart):
        pass

    @abc.abstractmethod
    def generate_images_from(self):
        pass


class AllSkySNVisualization(ObsVisualization):
    """ Class implementing simplest ObsVisualization)"""
    def __init__(self, bandColorDict,  radius_deg=4.,
                 showMW=True,
                 colorCodeRedshifts=True,
                 showVisibleFields=False,
                 showVarScatter=False):
        self._radiusDegree = radius_deg
        self._bandColorDict = bandColorDict
        self._show_visible_fields = showVisibleFields
        self._show_var_scatter = showVarScatter
        self._show_mw = showMW
        self._colorCodez = colorCodeRedshifts

    @property
    def colorCodeRedshifts(self):
        return self._colorCodez

    @property
    def show_mw(self):
        return self._show_mw

    @property
    def show_visible_fields(self):
        return self._show_visible_fields

    @property
    def show_var_scatter(self):
        return self._show_var_scatter

    @property
    def radius_deg(self):
        return self._radiusDegree

    @property
    def band_color_dict(self):
        """dictionary definining bands and color reperesentation for bands
        in the survey"""
        return self._bandColorDict

    def generate_mw_polygons(self, m, fill=True,
                             color='y', alpha=0.1,
                             lw=0., edgecolor='y',
                             zorder=15,
                             min_lon=-90.,
                             ax=None):
        """obtain polygons representing the extent of the MW"""
        mwext = MilkyWayExtension(m=m, fill=fill, color=color, alpha=alpha,
                                  lw=lw, edgecolor=edgecolor, ax=ax,
                                  zorder=zorder, min_lon=min_lon)
        polygons = mwext.mw_polygon
        boundaries = mwext.mw_boundaries
        (ra_h, dec_h), (ra_l, dec_l) = boundaries
        x_h, y_h = m(ra_h, dec_h)
        x_l, y_l = m(ra_l, dec_l)
        bounds = (x_h, y_h), (x_l, y_l)
        return polygons, bounds


    def _hack_legend(self, ax, colors, labels, bbox=(1, 1), loc='best'):
        """ hack legend """
        cs, texts = [], []    
        for (colors, text) in zip(('g', 'r','y'),
                                 ('ZTF g band', 'ZTF r band', 'ZTF i band')):
            c = mpatches.Circle((0.5, 0.5), 0.25, facecolor="w",
                                edgecolor=colors, linewidth=1)
            cs.append(c)
            texts.append(text)

        # ax.plot([], [], label='MW', color='w', ls='--')
        l = ax.legend(cs, texts,
                      handler_map={mpatches.Circle: HandlerEllipse()},
                      loc=loc, bbox_to_anchor=bbox)
        #x = []
        #for (c, l) in zip(colors, labels):
        #    ax.hist(x, color=c, label=l)
        #    l  = ax.legend(loc=loc, bbox_to_anchor=bbox)
        return l

    def _breakpoint(self, x):
        ind = np.diff(x).argmax()
        return ind + 1

    def generate_image_bg(self, projection='moll', drawmapboundary=True,
                          bg_color='k', mwcolor='w', mwfill=False,
                          mw_alpha=1.0, mw_edgecolor='w', mw_lw=1.,
                          figsize=(12, 6),
                          **kwargs):
        """Generate a figure axis, and a Basemap child instance"""

        fig, ax = plt.subplots(figsize=figsize)
        m = AllSkyMap(projection=projection, lon_0=0., lat_0=0.,
                      ax=ax, celestial=True)
        _ = m.drawparallels(np.arange(-91., 91., 60.))
        _ = m.drawmeridians(np.arange(-180., 181., 60.))
        _ = m.drawmapboundary(fill_color=bg_color)

        if self.show_mw:
            if mwfill:
                min_lon=-30
            else:
                min_lon=-90
            polygons, bounds = self.generate_mw_polygons(m, color=mwcolor,
                                                         alpha=mw_alpha,
                                                         fill=mwfill,
                                                         lw=mw_lw,
                                                         min_lon=min_lon,
                                                         edgecolor=mw_edgecolor)
            if mwfill:
                _ = ax.add_patch(polygons)
            else:
                (x_h, y_h), (x_l, y_l) = bounds
                break_h = len(x_h) - self._breakpoint(x_h)
                break_l = len(x_l) - self._breakpoint(x_l)
                _ = ax.plot(np.roll(x_h, break_h) , np.roll(y_h, break_h),
                                    ls='dashed', color=mw_edgecolor, lw=mw_lw,
                                    label='MW')
                _ = ax.plot(np.roll(x_l, break_l), np.roll(y_l, break_l),
                                    ls='dashed', color=mw_edgecolor, lw=mw_lw)
        return fig, ax, m

    def generate_camera(self, lon_0, lat_0, m, ax, band='g', radius_deg=4.,
                        default_color='k'):
        """Generate an image of a circular field of view representing the
        camera on the projection with a color representing the bandpass
        filter
        """
        try:
            c = self.band_color_dict[band]
        except (TypeError, KeyError):
            print('band_color_dict found and band are \n',
                  self.band_color_dict, band, len(band))
            print('setting default color\n')
            c = default_color
        camera_polygons = m.tissot(lon_0=lon_0, lat_0=lat_0, radius_deg=radius_deg,
                                   npts=100, ax=ax, add_patch=True,
                                   **dict(fill=False, edgecolor=c,
                                          lw=2))
        return camera_polygons

    def label_time_image(self, mjd, surveystart=None):
        label = '{:0.5f}'.format(mjd)
        return label

    def get_visible_field_polygons(self, mjd, m, facecolor, alpha, **kwargs):
        pass

    def generate_var_scatter(self):
        pass

    def generate_image(self, ra, dec, radius_deg, mjd=None, npts=100, band='g',
                       projection='moll', drawmapboundary=True, mwColor='w',
                       mwAlpha=1.0, mwEdgeColor='w', mwLw=1., mwFill=False,
                       bg_color='k', alpha=0.5, vfcolor='k',
                       cmap=plt.cm.Reds, sndf=None,
                       zlow=0., zhigh=0.2, surveystart=None,
                       bbox=(1, 1), loc=None,
                       **kwargs):
        """
        Use methods above to create an image of the sky and optionally save
        it.
        """
        fig, ax, m = self.generate_image_bg(projection=projection,
                                            drawmapboundary=drawmapboundary,
                                            mwcolor=mwColor, mw_alpha=mwAlpha,
                                            mw_edgecolor=mwEdgeColor,
                                            mwfill=mwFill,
                                            mw_lw=mwLw, bg_color=bg_color,
                                            **kwargs)
        if self.show_visible_fields:
            visible_polygons = self.get_visible_field_polygons(mjd, m,
                                                               facecolor=vfcolor,
                                                               alpha=alpha,
                                                               **kwargs)
            for poly in visible_polygons:
                _ = ax.add_patch(poly)

        camera_polygons = self.generate_camera(lon_0=ra, lat_0=dec, m=m, ax=ax,
                                               band=band)

        xx = None
        if self.show_var_scatter:
            simdf = self.generate_var_scatter(mjd, band, sndf)
            ra = simdf.ra.values
            dec = simdf.dec.values
            s = simdf.area
            x, y = m(ra, dec)
            if self.colorCodeRedshifts:
                z = simdf.z
                normalize = matplotlib.colors.Normalize(vmin=zlow, vmax=zhigh)
                xx = ax.scatter(x, y,  lw=0., s=s,
                                c=z, zorder=10, norm=normalize, cmap=cmap)
                cb = plt.colorbar(xx, ax=ax, norm=normalize,
                                  **dict(orientation='horizontal'))
                cb.set_label('redshift')
            else:
                ax.scatter(x, y, s=s.values, c='w', edgecolors='w',
                           zorder=10)


        label = self.label_time_image(mjd, surveystart)
        ax.set_title(label)
        # add legend
        cvals = (mwColor, vfcolor)
        names = ('MW Region', 'Desired SN fields')
        legend = self._hack_legend(ax, colors=cvals, labels=names, bbox=bbox,
                                  loc=loc) 
        return fig, ax, m, xx

    def generate_images_from(self):
        pass


def split_PolygonSegments(lon, lat, lon_split=180., epsilon=0.000001):
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
        self.lon_split = 180.
        self.epsilon=1.0e-10

    def tissot(self, lon_0, lat_0, radius_deg, npts, ax=None, epsilon=1.0e-10,
               add_patch=True, **kwargs):
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
        g = pyproj.Geod(a=self.rmajor, b=self.rminor)
        az12, az21, dist = g.inv(lon_0, lat_0, lon_0, lat_0 + radius_deg)
        seg = [self(lon_0, lat_0 + radius_deg)]
        delaz = 360. / npts
        az = az12
        for n in range(npts):
            az = az + delaz
            lon, lat, az21 = g.fwd(lon_0, lat_0, az, dist)
            seg.append((lon, lat))
        ra, dec = zip(*seg)
        # go to 0, 360.
        ra = np.asarray(ra)
        dec = np.asarray(dec)
        ## ra[ra < 0] += 360.
        ra = ra % 360.
        lon_0 = lon_0 % 360.
        polyseg_list = split_PolygonSegments(ra, dec,
                                             lon_split=self.lon_split)
        # if np.any(np.abs(ra - lon_0) < radius_deg - epsilon):
        #    print('split_polygon being called')
        #    polyseg_list = split_PolygonSegments(
        #                        ra, dec, lon_split=self.lonmax)
        # else:
        #    polyseg_list = list((ra, dec))
        split_poly_normed = list()
        for radec in polyseg_list:
            if len(radec) > 0:
                ra, dec = radec
                ra = np.asarray(ra)
                mask = ra > self.lon_split
                ra[mask] -= 360.
                # Leave out points on the limb
                points_on_split = np.logical_or(ra > 180. - epsilon,
                                                ra < -180. + epsilon)
                split_poly_normed.append((ra[~points_on_split],
                                          dec[~points_on_split]))

        pts = list()
        for poly_segs in split_poly_normed:
            lon, lat = poly_segs
            if len(lon) > 1:
                pt = self.polygonize(lon, lat, zorder=20., **kwargs)
                if add_patch:
                    ax.add_patch(pt)
                pts.append(pt)
        return pts

    def polygonize(self, lon, lat, zorder=1, **kwargs):
        x, y = self(lon, lat)
        x = np.asarray(x)
        y = np.asarray(y)
        mask = np.logical_and(x < 1.0e20, y < 1.0e20)
        xy = zip(x[mask], y[mask])
        pt = Polygon(xy,  zorder=zorder, **kwargs)
        return pt


def plot_south_steradian_view(ra, dec, numPoints=100, radius=1.75, boundary=20.,
                              ax=None, show_frame=True):
    """
    """
    if ax is None:
        fig, ax = plt.subplots()
    m = Basemap(projection='spstere', boundinglat=boundary, lon_0=0., ax=ax)
    if show_frame:
        m.drawparallels(np.arange(-80., 81., 20.))
        m.drawmeridians(np.arange(-180., 181., 20.))
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
        return pixelsForAng(lon=ra, lat=dec, nside=nside, unit=units)

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
        summary = self.pointingSummary(tileID)  # , columns=[raCol, decCol])

        if query is not None:
            summary = summary.query(query)
        if opsimAngularUnits != 'radians':
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
            ra, dec = self.pointingCenters(
                tileID, raCol=raCol, decCol=decCol, query=query)
            for ra, dec in zip(ra, dec):
                m.tissot(ra, dec, radius, 100, **kwargs)

        # Draw patch
        ax.add_patch(healpixels)

        # really important for the case where ax is None
        if ax is not None:
            fig = ax.figure

        return fig, (ra_tile, dec_tile), (llcrnrlat, llcrnrlon, urcrnrlat, urcrnrlon)
