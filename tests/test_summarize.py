""" Tests for the code in `opsimsummary/opsim_out.py`
"""
from __future__ import print_function, division, absolute_import
import os
import pytest
import numpy as np
import pandas as pd
import opsimsummary as oss
from opsimsummary import (OpSimOutput, SynOpSim, PointingTree)
import healpy as hp
from numpy.testing import assert_allclose, assert_array_equal


@pytest.fixture()
def opsfile(fname):
    yield os.path.join(oss.example_data, fname)
    print('teardown')


testdata_synopsiminit = [('enigma_1189_micro.db', 'lsstv3',
                         ('Summary', 'Proposal')),
                         ('opsimv4_feat_micro.db', 'lsstv4',
                         ('SummaryAllProps', 'Proposal'))]


@pytest.mark.parametrize("fname,opsimversion,tableNames",
                         testdata_synopsiminit)
def test_SynOpSim_init(fname, opsimversion, tableNames):
    fname = os.path.join(oss.example_data, fname)
    ops_out = OpSimOutput.fromOpSimDB(fname, opsimversion=opsimversion,
                                      tableNames=tableNames)
    num_visits = len(ops_out.summary)
    synopsim = SynOpSim(ops_out.summary)
    assert synopsim.pointings.night.size == num_visits

@pytest.mark.parametrize("fname,opsimversion,tableNames",
                         testdata_synopsiminit)
def test_fromOpSimDB(fname, opsimversion, tableNames):
    fname = os.path.join(oss.example_data, fname)
    ops_out = OpSimOutput.fromOpSimDB(fname, opsimversion=opsimversion,
                                      tableNames=tableNames)
    synopsim = SynOpSim(ops_out.summary)
    synopsimd = SynOpSim.fromOpSimDB(fname, opsimversion=opsimversion,
                                          tableNames=tableNames)
    assert synopsim.pointings.equals(synopsimd.pointings)


@pytest.mark.parametrize("fname,opsimversion,tableNames",
                         testdata_synopsiminit)
def test_pointingTree(fname, opsimversion, tableNames):
    """test that `PointingTree` is set up correctly to handle
    pointings provided by `OpSimOutput.summary` by checking that
    distance calculations for the 10 nearest ones match
    """ 
    fname = os.path.join(oss.example_data, fname)
    ops_out = OpSimOutput.fromOpSimDB(fname, opsimversion=opsimversion,
                                      tableNames=tableNames)
    pointings = ops_out.summary.iloc[:100]

    # Find distances to a specific location ra=54.0 deg, dec=-27.5 deg
    radeg = 54.0
    decdeg = -27.5
    ra = np.radians(radeg)
    dec = np.radians(decdeg)

    ptree = PointingTree(pointings, raCol='_ra', decCol='_dec')

    # Find the 10 nearest pointings 
    dists, idxs = ptree.tree.query([(dec, ra)], k=10)

    # Now calculate the distances more directly
    # order the 10 pointings as above
    obs = pointings.iloc[idxs[0]]
    vec = hp.ang2vec(radeg, decdeg, lonlat=True)
    vecs = hp.ang2vec(np.degrees(obs._ra),
                      np.degrees(obs._dec),
                      lonlat=True)
    assert_allclose(np.dot(vecs, vec), np.cos(dists[0]))


@pytest.mark.parametrize("fname,opsimversion,tableNames",
                         testdata_synopsiminit)
def test_synopsimPTree(fname, opsimversion, tableNames):
    """
    test that the pointings for a particular point at ra = 54.0 deg and
    dec = -27.5 deg are found through the pointing tree and through a
    calculation through the vectors and dot products match
    """
    fname = os.path.join(oss.example_data, fname)
    synopsim = SynOpSim.fromOpSimDB(fname,
                                    opsimversion=opsimversion,
                                    tableNames=tableNames,
                                    usePointingTree=True)


    # Find distances to a specific location ra=54.0 deg, dec=-27.5 deg
    radeg = 54.0
    decdeg = -27.5

    pts = synopsim.pointingTree.pointingsEnclosing(ra=radeg,
                                                   dec=decdeg,
                                                   circRadius=0.,
                                                   pointingRadius=1.75)
    res = np.sort(pts[0])

    # calculate differently
    X = synopsim.pointings[['_dec', '_ra']].apply(np.degrees)

    dec = X.values[:, 0]
    ra = X.values[:, 1]
    vecs = hp.ang2vec(ra, dec, lonlat=True)
    vec = hp.ang2vec(radeg, decdeg, lonlat=True)

    X['dist'] = np.dot(vecs, vec)

    costhresh = np.cos(np.radians(1.75))
    assert_array_equal(np.sort(X.query('dist > @costhresh').index),
                       res)


@pytest.mark.parametrize("fname,opsimversion,tableNames",
                         testdata_synopsiminit)
def test_synopsimPointings(fname, opsimversion, tableNames):
    """
    test that the lengths of pointings for a random set of 100 positions
    in the approximate LSST region match when calculated directly or with
    the pointing tree
    """
    fname = os.path.join(oss.example_data, fname)
    synopsim = SynOpSim.fromOpSimDB(fname,
                                    opsimversion=opsimversion,
                                    tableNames=tableNames,
                                    usePointingTree=True)


    # Find distances to a specific location ra=54.0 deg, dec=-27.5 deg
    num = 100
    radeg = np.random.uniform(0., 360., size=num)
    decdeg = np.random.uniform(-66., 0., size=num)

    # Use the PointingTree
    pts = synopsim.pointingsEnclosing(ra=radeg, dec=decdeg,
                                      circRadius=0., pointingRadius=1.75, 
                                      usePointingTree=True, subset=[])
    ptsh = synopsim.pointingsEnclosing(ra=radeg, dec=decdeg,
                                       circRadius=0., pointingRadius=1.75, 
                                       usePointingTree=False, subset=[])

    ptslens = np.array(list(len(pt) for pt in pts))
    ptshlens = np.array(list(len(pt) for pt in ptsh))
    np.testing.assert_array_equal(ptslens, ptshlens)
    assert max(ptslens) > 0
    
