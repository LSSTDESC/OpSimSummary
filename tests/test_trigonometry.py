from __future__ import print_function, absolute_import, division
from opsimsummary import convertToSphericalCoordinates
import numpy as np


def test_convertToSphericalCoordinates():
    """
    Test the workings of convertToSphericalCoordinates by
    - 1. checking a simple scalar input, and comparing the output to expected
        values in degrees
    - 2. Check that the same outputs come out when unit is not specified. ie
        defaults
    - 3. Check that exceptions are raised when the unit is not in degrees or
        radians
    - 4. Check that 
    """
    ra = 24.0
    dec = 90.0

    theta, phi = convertToSphericalCoordinates(ra=ra, dec=dec, unit='degrees')
    assert theta == np.array([0.0])
    assert phi == np.array([0.41887902047863912])
