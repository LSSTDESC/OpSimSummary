""" Tests for the code in `opsimsummary/opsim_out.py`
"""
from __future__ import print_function, division, absolute_import
import os
import pytest
import numpy as np
import pandas as pd
import opsimsummary as oss
from opsimsummary import OpSimOutput, SynOpSim


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
