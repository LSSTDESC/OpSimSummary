"""
Class for describing tesselations of the sphere

Following the usual usage, a tesselation or tiling of the sphere is a
covering of the sphere by a set of gemetrical shapes called tiles, such that
each point of the sphere belongs to one and only one tile (we will be negligent
about edges of tiles). Thus, no two tiles overlap, and there is no point on the
sphere which does not belong to a tile.

Such tilings provide a natural spatial grouping of supernovae on the sky. Such
spatial groupings lead to similar (if not same) observational properties in
terms of the set of observational pointings, as well as the sky properties (eg.
seeing, psf, airmass, MW extinction etc.). The similarities increase as the size of the tiles
decrease.  Therefore, a sensible way to distribute SN simulations is by grouping
the simulations on small tiles together.

Allowed tilings must satisfy certain properties as encoded in the abstract base
class below as well as listed here. Additionally, there are desired properties
that are hepful (for example in speeding up the process), but not essential to
the computations. Such methods may involve properties like equal area of tiles,
hierarchical structure of times, efficient queries, or precomputed results for
OpSim outputs.

Methods and Attributes Required and Intended Usage:
---------------------------------------------------
    - tileIDSequence : Obtain the entire set of tiles which cover the unit
        sphere
    - pointingSequenceForTile :  A method to obtain the maximal set of
        pointings that may overlap a point on a tile. Obviously, there can be
        points in a tile which do not overlap with subsets of such pointings.
        The intended usage is to simulate all SN in a tile together by using
        the set of maximal pointings. This can be done in modes
"""

from __future__ import absolute_import, print_function
from future.utils import with_metaclass
import abc
import numpy as np

__all__ = ["Tiling"]

class Tiling(with_metaclass(abc.ABCMeta, object)):
    """
    Abstract Base Class for tilings specifying the methods that are mandatory

    Attributes
    ----------
    tileIDSequence : sequence
        sequence of IDs indexing the tiles. The IDs may be integers or strings.

    Methods
    -------
    pointingSequenceForTile :
    """
    @abc.abstractmethod
    def __init__(self):
        pass

    @abc.abstractproperty
    def tileIDSequence(self):
        pass


    @abc.abstractmethod
    def pointingSequenceForTile(self, tileID, allPointings=None, columns=None,
                                **kwargs):
        """
        Return a sequence of IDs identifying the maximal set of OpSim pointings
        (obsHistID) that intersect with the tile having ID tileID.

        Parameters
        ----------
        tileID : int, or string, mandatory
            Index for desired tile (should not be for a sequence)
        allPointings : instance of {string, DBConnection, `pandas.DataFrame`}
            Information about a set of pointings we will worry about. The set
            of pointings may be in a database or a dataFrame, and different ways
            of connecting to them is ideally supported.
        columns : tuple of strings, defaults to None
            if None returns only obsHistIDs. Otherwise returns the columns
            listed
        kwargs : extra parameters
            specify method to use optional precomputations to speed up this
            function
        Returns
        -------
        `numpy.ndarray` of obsHistIDs

        .. notes: This can be a crude method returning all of the pointings in
            allPointings if one runs in a pedantic mode later to do a more
            careful filtering. Even in those cases, this may be helpful in
            reducing the nimbers
        """
        pass

    @abc.abstractmethod
    def area(self, tileID):
        """
        return the area of the tile with ID tileID

        Parameters
        ----------
        tileID : int or string, or `numpy.ndarray` thereof, mandatory
            Index for desired tile

        Returns
        -------
        area : `numpy.ndarray` of dtype float
        """
        pass

    @abc.abstractmethod
    def tileIDsForSN(self, ra, dec):
        """
        return a numpy array of tileIDs for point sources located at ra, dec
        where ra, dec are `numpy.ndarrays` each of size equal to number of
        point sources.

        Parameters
        ----------
        ra : `numpy.ndarray` of dtype float, mandatory
            ra values 
        dec : `numpy.ndarray` of dtype float, mandatory
            dec values

        Returns
        -------
        `numpy.ndarray` of tileIDs. So the dtype is probably integer or string
        """
        pass

    @abc.abstractmethod
    def positions(self, tileID, numSamples):
        """
        return a tuple of numpy arrays theta and phi, each of size numSamples
        """

    @staticmethod
    def samplePatchOnSphere(phi, theta, delta, size, rng, degrees=True):
        """
        Uniformly distributes samples on a patch on a sphere between
        phi \pm delta, and theta \pm delta on a sphere. Uniform distribution
        implies that the number of points in a patch of sphere is proportional
        to the area of the patch. Here, the coordinate system is the usual
        spherical coordinate system with the azimuthal angle theta going from
        0 degrees at the North Pole, to 90 degrees at the South Pole, through
        0. at the equator.

        This function is not equipped to handle wrap-around the ranges of theta
        phi and therefore does not work at the poles.
   
        Parameters
        ----------
        phi: float, mandatory, degrees
            center of the spherical patch in ra with range 
        theta: float, mandatory, degrees
        delta: float, mandatory, degrees
        size: int, mandatory
            number of samples
        seed : int, optional, defaults to 1
            random Seed used for generating values
        degrees : bool, optional, defaults to True
            if True, returns angles in degrees, else in
            radians

        Returns
        -------
        tuple of (phivals, thetavals) where phivals and thetavals are arrays of 
            size size in degrees.
        """
        u = rng.uniform(size=size)
        v = rng.uniform(size=size)
        phi = np.radians(phi)
        theta = np.radians(theta)
        delta = np.radians(delta)

        phivals = 2. * delta * u + (phi - delta)
        phivals = np.where(phivals >= 0., phivals, phivals + 2. * np.pi)

        # use conventions in spherical coordinates
        # theta = np.pi/2.0 - theta
        thetamax = theta + delta
        thetamin = theta - delta

        # if thetamax > np.pi or thetamin < 0. :
        #    raise ValueError('Function not implemented to cover wrap around poles')

        # Cumulative Density Function is cos(thetamin) - cos(theta) / cos(thetamin) - cos(thetamax)
        a = np.cos(thetamin) - np.cos(thetamax)
        thetavals = np.arccos(-v * a + np.cos(thetamin))

        if degrees:
            return np.degrees(phivals), np.degrees(thetavals)
        else:
            return phivals, thetavals
