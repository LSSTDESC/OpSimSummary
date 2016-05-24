""" Script to generate truncated output of the OpSim output enigma_1189_sqlite.db.
An example of such output is stored in example_data/enigma_1189_micro.db which
is a sqlite database that has the same structure as the SUMMARY table in enigma_1189_sqlite.db


USAGE:
    to use it point to  replace opsimdb by an opsim with a valid output in the
    form of a sqlite database.  Replace 'path' and __OPSIMDB__ for your the path
    to the opsim output on your machine, and the path to it.
"""
import os.path
import opsimsummary as oss

def query_for_schema(tableName):
    schema_query = """SELECT sql FROM sqlite_master WHERE tbl_name = '{}'""".format(tableName)
    return schema_query

def db_interact(db, query, cursor=None, fetchall=False, commit=False):

    if cursor is None:
        conn = sqlite3.connect(db)
        cursor = conn.cursor()
    cursor.execute(query)
    if fetchall:
        return cursor.fetchall()
    elif commit:
        conn.commit()

def insert_allvals_statement(tableName, record):
    s = 'INSERT INTO {0} VALUES {1}'.\
        format(tableName, tuple(str(elem) for elem in record))
    return s
        
def insertfromdata(tablename, records, multiple=True):
          """
          construct string to insert multiple records into sqlite3 database
          args:
              tablename: str, mandatory
                  Name of table in the database.
              records: set of records
              multiple:
          returns:
          """
          if multiple:
              lst = records[0]
          else:
              lst = records
          s = 'INSERT INTO ' + str(tablename) + ' VALUES '
          s += "( " + ", ".join(["?"]*len(lst)) + ")"
          return s
if __name__ == '__main__':

    import sqlite3



    path = '/Users/rbiswas/data/LSST/OpSimData'
    __OPSIMDB__ = 'enigma_1189_sqlite.db'
    opsimdb = os.path.join(path, __OPSIMDB__)
    test_opsimDB = 'enigma_1189_micro.db'
    test_opdb = os.path.join('../opsimsummary', 'example_data',
                             test_opsimDB)

    wfd_selection = """ SELECT * FROM Summary WHERE night > 216 and night < 366 and PROPID is 364"""
    ddf_selection = """ SELECT * FROM Summary WHERE night > 300 and night < 366 and PROPID is 366"""
    schema_query = query_for_schema('Summary')
    #print(schema_query)

    # get the schema, this is equivalen to .schema on the sqlite command line
    # Since this is part of fetchall, this is a list of tuples (in this case a
    # single tuple. We pull this out and cast it to string from unicode
    schema = str(db_interact(db=opsimdb, query=schema_query, fetchall=True)[0][0])
    print(schema)


    # Clone the database structure
    db_interact(query=schema, db=test_opdb) 


    selected_records = db_interact(db=opsimdb, query=wfd_selection,
                                   fetchall=True)
    
    selected_ddfrecs = db_interact(db=opsimdb, query=ddf_selection,
                                   fetchall=True)
    for record in selected_records:
        s = insert_allvals_statement('SUMMARY', record)
        db_interact(query=s, db=test_opdb, commit=True)
    for rec in selected_ddfrecs:
        s = insert_allvals_statement('SUMMARY', rec)
        db_interact(query=s, db=test_opdb, commit=True)


    db = sqlite3.connect(test_opdb)
    db.commit()
    db.close()



