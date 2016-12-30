from __future__ import print_function, absolute_import, division
from opsimsummary import convertToSphericalCoordinates, convertToCelestialCoordinates
import numpy as np
from numpy.testing import assert_allclose

def test_convertToCelestialCoordinates():
    """
    Test that some limiting conditions work 
    - Check that this works correctly for float inputs in default units
    - Check that this works correctly for numpy inputs in default units

    """
    phi = np.radians(30.)
    theta = np.pi/2.0 
    theta_arr = np.repeat(np.pi/2., 5)
    phi_arr = np.repeat(np.radians(30.), 5)

    ra, dec = convertToCelestialCoordinates(theta, phi) 
    ra_arr, dec_arr = convertToCelestialCoordinates(theta_arr, phi_arr) 
    assert_allclose(dec_arr, np.array([0.0, 0.0, 0.0, 0.0, 0.0]),
                       13)
    assert_allclose(ra_arr, np.array([30.0, 30.0, 30.0, 30.0, 30.0]),
                        13)



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
