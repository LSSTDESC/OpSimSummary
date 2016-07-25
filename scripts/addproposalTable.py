""" 
Script to generate truncated output of the OpSim output enigma_1189_sqlite.db.
An example of such output is stored in example_data/enigma_1189_micro.db which
is a sqlite database that has the same structure
structure as the SUMMARY table in enigma_1189_sqlite.db


USAGE:
    form of a sqlite database.  Replace 'path' and __OPSIMDB_
    to the opsim output on your machine, and the path to it.
"""
from __future__ import absolute_import, division, print_funciton
import os.path
import opsimsummary as oss
from .make_smallOpSim import (query_for_schema,
                              db_interact,
                              insert_allvals_statement,
                              insertfromdata)


if __name__ == '__main__':

    import sqlite3



    path = '/Users/rbiswas/data/LSST/OpSimData'
    __OPSIMDB__ = 'enigma_1189_sqlite.db'
    opsimdb = os.path.join(path, __OPSIMDB__)
    test_opsimDB = 'enigma_1189_micro.db'
    test_opdb = os.path.join('../opsimsummary', 'example_data',
                             test_opsimDB)

    proposal_selection = """SELECT * FROM Proposal"""
    schema_query = query_for_schema('Proposal')
    #print(schema_query)

    # get the schema, this is equivalen to .schema on the sqlite command line
    # Since this is part of fetchall, this is a list of tuples (in this case a
    # single tuple. We pull this out and cast it to string from unicode
    schema = str(db_interact(db=opsimdb, query=schema_query, fetchall=True)[0][0])
    print(schema)


    # Clone the database structure
    db_interact(query=schema, db=test_opdb) 


    selected_records = db_interact(db=opsimdb, query=proposal_selection,
                                   fetchall=True)
    
    for record in selected_records:
        s = insert_allvals_statement('Proposal', record)
        db_interact(query=s, db=test_opdb, commit=True)

    db = sqlite3.connect(test_opdb)
    db.commit()
    db.close()



