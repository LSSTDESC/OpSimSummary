"""
Script to generate truncated output of the OpSim output enigma_1189_sqlite.db.
An example of such output is stored in example_data/enigma_1189_micro.db which
is a sqlite database that has the same structure
structure as the SUMMARY table in enigma_1189_sqlite.db


USAGE:
    form of a sqlite database.  Replace 'path' and __OPSIMDB_
    to the opsim output on your machine, and the path to it.
"""
import os.path
import sqlite3

def query_for_schema(tableName):
    """
    Generate a string that can be used as a query for a sqlite database to get
    the equivalent of the command-line .schema command for a particular table.
    This details how the table was created, and the same structure may be clone
    by executing this query on a separate database.
    Parameters
    ----------
    tableName : string, mandatory
        Name of the table on which this information is desired

    Returns
    -------
        string, query  that when executed on a cursor of a sqlite database
        without the table, will create a table called tableName with the same
        structure (but not data).

    Examples
    --------
    >>> schema_query = query_for_schema('Summary')
    >>> pkgDir = os.path.split(oss.__file__)[0]
    >>> dbname = os.path.join(pkgDir, 'example_data', 'enigma_1189_micro.db')
    >>> x = db_interact(dbname, schema_query, fetchall=True) 
    >>> assert(x[0][0] ==  oldstring) # SKIP 
    True
    """
    schema_query = \
    """SELECT sql FROM sqlite_master WHERE tbl_name = '{}'""".format(tableName)
    return schema_query

def db_interact(db, query, cursor=None, fetchall=False, commit=False):
    """
    interact with a sqlite connection through a cursor

    Parameters
    ----------
    """

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
    __OPSIMDB__ = 'feature_rolling_half_mask_10yrs.db'
    summaryTableName = 'SummaryAllProps'
    wfdPropID = 0
    ddfPropID = 1
    propIDinSummaryTable = 'proposalId'

    opsimdb = os.path.join(path, __OPSIMDB__)
    test_opsimDB = 'opsimv4_feat_micro.db'
    test_opdb = os.path.join('../opsimsummary', 'example_data',
                             test_opsimDB)

    wfd_selection = """ SELECT * FROM {0} WHERE night > 216 and night < 366 and {2} is {1}""".format(summaryTableName, wfdPropID, propIDinSummaryTable)
    ddf_selection = """ SELECT * FROM {0} WHERE night > 300 and night < 366 and {2} is {1}""".format(summaryTableName, ddfPropID, propIDinSummaryTable)
    schema_query = query_for_schema(summaryTableName)
    #print(schema_query)

    # get the schema, this is equivalen to .schema on the sqlite command line
    # Since this is part of fetchall, this is a list of tuples (in this case a
    # single tuple. We pull this out and cast it to string from unicode
    print(db_interact(db=opsimdb, query=schema_query, fetchall=True))
    schema = str(db_interact(db=opsimdb, query=schema_query, fetchall=True)[0][0])
    print(schema)


    # Clone the database structure
    db_interact(query=schema, db=test_opdb) 


    selected_records = db_interact(db=opsimdb, query=wfd_selection,
                                   fetchall=True)
    
    selected_ddfrecs = db_interact(db=opsimdb, query=ddf_selection,
                                   fetchall=True)
    for record in selected_records:
        s = insert_allvals_statement(summaryTableName, record)
        db_interact(query=s, db=test_opdb, commit=True)
    for rec in selected_ddfrecs:
        s = insert_allvals_statement(summaryTableName, rec)
        db_interact(query=s, db=test_opdb, commit=True)


    db = sqlite3.connect(test_opdb)
    db.commit()
    db.close()



