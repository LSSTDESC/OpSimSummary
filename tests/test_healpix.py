from __future__ import division, print_function, absolute_import
import opsimsummary as oss
import os
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
import unittest
import healpy as hp

class Test_obsHistIDsFortileID(unittest.TestCase):
    """
    """
    def setUp(self):
        # Use example opsim output in example files
        # The test db is from enigma 1189 (current sqlite format),
        pkgDir = os.path.split(oss.__file__)[0]
        dbname = os.path.join(pkgDir, 'example_data', 'enigma_1189_micro.db')
        engineFile = 'sqlite:///' + dbname
        engine = create_engine(engineFile)
    
        # read the WFD fields into a `pd.DataFrame` for test opsim db
        # 364 is hard coded rather than reading the table on OpSim
        opsimdf = pd.read_sql_query('SELECT * FROM Summary WHERE PropID is 364',
                                    con=engine,
                                    index_col='obsHistID') 
        
        # Extremely coarse grained pixels (12 covering the sphere) for tests
        self.nside = 1
        h = oss.HealPixelizedOpSim(opsimDF=opsimdf, NSIDE=self.nside)
        h.doPreCalcs()
        self.hpOps = h

    def test_obsHistIDsForfield8(self):
        """
        Check that a healpix ID p is in the list of hids associated with 
        each of the pointings intersecting with it.
        """
        h = self.hpOps
        hids = np.arange(hp.nside2npix(self.nside))
        l = list(all(h.opsimdf.query('obsHistID in @h.obsHistIdsForTile(@p)')\
                     .hids.apply(lambda x: p in x).values) for p in hids)
        # For p in any of the healix ids:
        #    get set of obsHistIds associated with p using h.obsHistIdsForTile(p)
        #       obtain the healpix id list associated with that pointing through
        #       query disc (precalculated as hids)
        #    Check that p is in hids
        assert all(l)

