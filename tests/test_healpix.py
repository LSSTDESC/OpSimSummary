"""
Tests associated with the `opsimsummary/healpix.py` module
"""
from __future__ import division, print_function, absolute_import
import opsimsummary as oss
import os
import pytest
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
import unittest
import healpy as hp
import sqlite3

class Test_obsHistIDsFortileID(unittest.TestCase):
    """
    Tests associated with obtaining obsHistID values for tileIDs using healpix
    in  the nest scheme.
    """
    @classmethod
    def setUpClass(cls):
        # Use example opsim output in example files
        # The test db is from enigma 1189 (current sqlite format),
        pkgDir = os.path.abspath(os.path.split(oss.__file__)[0])
        dbname = os.path.join(pkgDir, 'example_data', 'enigma_1189_micro.db')
        engineFile = 'sqlite:///' + dbname
        engine = create_engine(engineFile)
    
        # read the WFD fields into a `pd.DataFrame` for test opsim db
        # 364 is hard coded rather than reading the table on OpSim
        sql_query = 'SELECT * FROM Summary WHERE PropID == 364'
        print(sql_query)
        opsimdf = pd.read_sql_query(sql_query,
                                    con=engine,
                                    index_col='obsHistID') 
        
        # Extremely coarse grained pixels (12 covering the sphere) for tests
        cls.nside = 4
        # Standard DB for integration tests
        stdDBname = os.path.join(pkgDir, 'example_data', 'healpixels_micro.db')
        cls.stdDBConn = sqlite3.Connection(stdDBname)
        cls.newDB = 'healpix_micro_new.db'
        if os.path.exists(cls.newDB):
            os.remove(cls.newDB)
        cls.newconn = sqlite3.Connection(cls.newDB)
        h = oss.HealPixelizedOpSim(opsimDF=opsimdf, NSIDE=cls.nside, source=dbname)
        h.doPreCalcs()
        try:
            version = oss.__version__
            h.writeToDB(cls.newDB, version=version)
        except:
            cls.tearDownClass()
            raise Warning('Had to erase teardown the class to set it up')
        cls.hpOps = h

    def test_obsHistIDsForfields(self):
        """
        Check that a healpix ID p is in the list of hids associated with 
        each of the pointings intersecting with it.
        """
        h = self.hpOps
        hids = np.arange(hp.nside2npix(self.nside))
        l = list(all(h.opsimdf.query('obsHistID in @h.obsHistIdsForTile(@p)')\
                     .hids.apply(lambda x: p in x).values) for p in hids)
        # Here is what is going on above
        # For p in any of the healix ids:
        #    get set of obsHistIds associated with p using h.obsHistIdsForTile(p)
        #       obtain the healpix id list associated with that pointing through
        #       query disc (precalculated as hids)
        #    Check that p is in hids
        assert all(l)

    def test_writemethod(self):
        """
        Sanity checks on the output database from writing out ipix, obsHistID database
        """
        # import sqlite3
        # newconn = sqlite3.Connection(self.newDB)
        newcursor = self.newconn.cursor()
        newcursor.execute('SELECT COUNT(*) FROM simlib')
        x = newcursor.fetchone()[0]
        # Hard coded value for enigma_micro database and NSIDE =1
        self.assertEqual(x, 159608)
        newcursor.execute('SELECT MIN(ipix) FROM simlib')
        y = newcursor.fetchone()
        self.assertEqual(y[0], 0)
        x = newcursor.execute('SELECT MAX(ipix) FROM simlib')
        y = x.fetchone()
        self.assertEqual(y[0], 191) 
        x = newcursor.execute('SELECT MIN(ipix) FROM simlib')
        y = x.fetchone()
        self.assertEqual(y[0], 0) 
        newcursor = self.newconn.cursor()
        x = newcursor.execute('SELECT * FROM metadata')
        y = x.fetchall()
        self.assertEqual(len(y), 1)
        self.assertEqual(len(y[0]), 10)
        # Check the version
        version = oss.__version__
        self.assertEqual(y[0][4], version)
        self.assertEqual(np.int(y[0][5]), self.nside)
        self.assertEqual(np.int(y[0][6]), 4)

    @pytest.mark.skip(reason='skipped before')
    def test_compareWithOldDB(self):
        """
        Compare any new database with a stored version of the old database
        """
        newcursor = self.newconn.cursor()

        npix = hp.nside2npix(1)
        ipixvalues = np.arange(npix) 

        stdCursor = self.stdDBConn.cursor()
        for hid in ipixvalues[:11]:
            newcursor.execute('SELECT obsHistID FROM simlib WHERE ipix == {}'\
                              .format(hid))
            _new = newcursor.fetchall()
            new = np.asarray(list(xx[0] for xx in _new))
            stdCursor.execute('SELECT obsHistID FROM simlib WHERE ipix == {}'\
                              .format(hid))
            _std = stdCursor.fetchall()
            std = np.asarray(list(xx[0] for xx in _std))
            np.testing.assert_equal(std, new, verbose=True, err_msg='std = {0} and new ={1}'\
                                  .format(std, new))
        

    def test_compareFunctionWithDB(self):
        """
        Test comparing values in written database to values calculated in-situ
        for 12 values of ipix
        """
        # import sqlite3
        # newconn = sqlite3.Connection(self.newDB)
        newcursor = self.newconn.cursor()

        npix = hp.nside2npix(self.nside)
        ipixvalues = np.arange(npix) 

        for hid in ipixvalues[:11]:
            newcursor.execute('SELECT obsHistID FROM simlib WHERE ipix == {}'\
                              .format(hid))
            x = newcursor.fetchall()
            y = np.asarray(list(xx[0] for xx in x))
    
            h = self.hpOps
            z = h.obsHistIdsForTile(hid)
            np.testing.assert_equal(y, z, verbose=True, err_msg='x = {0} and y ={1}'.format(y, z))
        



    @classmethod
    def tearDownClass(cls):
        import os
        if os.path.exists(cls.newDB):
            os.remove(cls.newDB)




