#!/usr/bin/env python 
from __future__ import absolute_import, division, print_function
import os
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
import matplotlib.pyplot as plt

__all__ = ['OpSimOutput']

class OpSimOutput(object):
    def __init__(self, summary=None, propIDDict=None, proposalTable=None,
                 subset=None, propIDs=None):
	self.summary = summary
	self.propIDDict = propIDDict
        self.proposalTable = proposalTable
	self.allowed_subsets = self.get_allowed_subsets()
        self.subset = subset
        self._propID = propIDs
    @property
    def propIds(self):
        if self._propID is not None:
            return self._propID
        elif self.subset is not None and self.propIDDict is not None:
            return self.propIDvals(self.propIDDict, self.subset)

    @classmethod
    def fromOpSimHDF(cls, hdfName, subset='combined',
                     tableNames=('Summary', 'Proposal'), propIDs=None):
        """
        """
	allowed_subsets = cls.get_allowed_subsets()
	subset = subset.lower()
	if subset not in allowed_subsets:
	    raise NotImplementedError('subset {} not implemented'.\
				      format(subset))
        # The hdf representation is assumed to be a faithful representation of
        # the OpSim output
        summarydf = pd.read_hdf(hdfName, key='Summary')

        if 'obsHistID' not in summarydf.colums:
            summarydf.reset_index(inplace=True)
            if 'obsHistID' not in summarydf.colums:
                raise NotImplementedError('obsHistID is not in columns')

        try:
            proposals = pd.read_hdf(hdfName, key='Proposals')
            propDict = cls.get_propIDDict(proposals)
            _propIDs = cls.propIDvals(subset, propDict)
        except:
            pass
        if propIDs is None:
            propIDs = _propIDs
        else:
            if np.asarray(propIDs).sort() != np.asarray(_propIDs).sort():
                raise ValueError('argument propIDs and subset donot match')

        if subset != '_all':
            if propIDs is None:
                if subset != 'unique_all':
                    raise ValueError('propID {0}, subset {1} combination'
                                     'does not seem to make sense.'\
                                     .format(propIDs, subset))
                summary  = summarydf
            else:
                # propIDStrings = ', '.\
                #    join('{}'.format(cls.propIDvals(subset, propDict)))
                summary = summarydf.query('propID == @propIDs')

            summary.drop_duplicates(subset='obsHistID', inplace=True)
        elif subset == '_all':
            summary = summarydf
        else :
            raise NotImplementedError('subset {} not implemented'\
                                      .format(subset))

            summary = summarydf.query('propID == @propIDs')
        summary.set_index('obsHistID', inplace=True)
        return cls(propIDDict=propDict, summary=summary,
                   proposalTable=proposals, subset=subset)

    def _validatePropIDs(self, propIDs, _propIDs):
        if propIDs is None:
            propIDs = _propIDs
        else:
            if np.asarray(propIDs).sort() != np.asarray(_propIDs).sort():
                raise ValueError('argument propIDs and _propIDs do not match')
        return propIDs

    @classmethod
    def fromOpSimDB(cls, dbname, subset='combined', propIDs=None):
	"""
	Class Method to instantiate this from an OpSim sqlite
	database output

	Parameters
	----------
	dbname :
	subset :
	"""
	allowed_subsets = cls.get_allowed_subsets()
	subset = subset.lower()
	if subset not in allowed_subsets:
	    raise NotImplementedError('subset {} not implemented'.\
				      format(subset))
        if not dbname.startswith('sqlite'):
            dbname =  'sqlite:///' + dbname
        print(' reading from database {}'.format(dbname))
        engine = create_engine(dbname, echo=False)

	# Read the proposal table to find out which propID corresponds to
        proposals = pd.read_sql_table('Proposal', con=engine)
        propDict = cls.get_propIDDict(proposals)
        _propIDs = cls.propIDvals(subset, propDict)

        propIDs = cls._validatePropIDs(propIDs, _propIDs)

        # Do the actual sql queries or table reads
        if subset in ('_all', 'unique_all'):
            # In this case read everything (ie. table read)
	    summary = pd.read_sql_table('Summary', con=engine)
        elif subset in ('ddf', 'wfd', 'combined')
	    sql_query = 'SELECT * FROM Summary WHERE PROPID'
	    sql_query += ' in {0}'.format(propDict['ddf'])
	    if subset == 'wfd':
            # _all will be used only to write out other serialized versions
            # of OpSim. Do not drop duplicates, so that different subsets can
            # be constructed from the same hdf file
	    if subset == 'unique_all':
	       summary.drop_duplicates(subset='obsHistID', inplace=True)	
            summary.set_index('obsHistID', inplace=True)
	    return cls(propIDDict=propDict, summary=summary,
                       proposalTable=proposals, subset=subset)
	else:
	    sql_query = 'SELECT * FROM Summary WHERE PROPID'
	    if subset == 'ddf':
		sql_query += ' == {0}'.format(propDict['ddf'])
	    if subset == 'wfd':
		sql_query += ' == {0}'.format(propDict['wfd'])
	    if subset == 'combined':
		sql_query += ' in [{0}, {1}]'.format(propDict['wfd'],
                                                     propDict['ddf'])
        # Read the summary table 
        summary = pd.read_sql_query(sql_query, con=engine)
	summary.drop_duplicates(subset='obsHistID', inplace=True)
	summary.set_index('obsHistID', inplace=True)
        return cls(propIDDict=propDict, summary=summary,
                   proposalTable=proposals, subset=subset)

    @staticmethod
    def get_allowed_subsets():
	return ('_all', 'ddf', 'wfd', 'combined', 'unique_all')
    @staticmethod
    def get_propIDDict(proposalDF):
	"""
        Return a dictionary with keys 'ddf', ad 'wfd' with the proposal IDs
        corresponding to deep drilling fields (ddf) and universal cadence (wfd) 

        Parameters
        ----------
        proposalDF : `pd.DataFrame`, mandatory
            a dataframe with the Proposal Table of the OpSim Run.
        Returns
        -------
        dictionary with keys 'wfd' and 'ddf' with values given by integers
            corresponding to propIDs for these proposals
	"""
	df = proposalDF
	mydict = dict()
	for i, vals in enumerate(df.propConf.values):
	    if 'universal' in vals.lower():
		if 'wfd' in mydict:
		    raise ValueError('Multiple propIDs for WFD found')
		mydict['wfd']  = df.propID.iloc[i]
	    elif 'ddcosmology' in vals.lower():
		if 'ddf' in mydict:
		    raise ValueError('Multiple propIDs for DDF found')
		mydict['ddf']  = df.propID.iloc[i] 
            else:
                mydict[vals.lower()] = df.propID.iloc[i]
	if len(mydict.items()) != len(df):
	    raise ValueError('Unexpected length of dictionary')
	return mydict

    @staticmethod
    def propIDVals(subset, propIDDict, proposalTable):
        """
        Parameters: 
        ----------
        subset : string
            must be member of OpSimOutput.allowed_subsets()
        propIDDict: dictionary, mandatory
            must have subset as a key, and an integer or seq of ints
            as values

        Returns:
        -------
        list of propID values (integers) associated with the subset
        """
        if subset.lower() in ('ddf', 'wfd'):
            return [propIDDict[subset.lower()]]
        elif subset.lower() == 'combined':
            return [propIDDict['ddf'], propIDDict['wfd']] 
        elif subset.lower() in ('_all', 'unique_all'):
            return None
        else:
            raise NotImplementedError('value of subset Not recognized')
        
def OpSimDfFromFile(fname, ftype='hdf', subset='Combined'):
    """
    read a serialized form of the OpSim output into `pd.DataFrame`
    and return a subset of interest

    Parameters
    ----------
    fname : string, mandatory
        absolute path to serialized form of the OpSim database
    ftype : {'sqliteDB', 'ASCII', 'hdf'}
        The kind of serialized version being read from.
            'sqliteDB' : `LSST` project supplied OpSim output format for
                baseline cadences (eg. enigma_1189, minion_1016, etc.) 
            'ASCII' : `LSST` project supplied OpSim output format used in
                older OpSim outputs eg. OpSim v 2.168 output
            'hdf' : `hdf` files written out by `OpSimSummary`
    subset : {'Combined', 'DDF', 'WFD' , 'All'}, defaults to 'Combined' 
        Type of OpSim output desired in the dataframe
        'Combined' : unique pointings in WFD + DDF 
        'WFD' : Unique pointings in WFD
        'DDF' : Unique pointings in DDF Cosmology
        'All' : Entire Summary Table From OpSim
    """
    if ftype == 'sqlite':
        dbname = 'sqlite:///' + fname
        engine = create_engine(dbname, echo=False)
        proposalTable =  pd.read_sql_table('Proposal', con=engine)

        if subset == 'DDF':
            sql


    elif ftype == 'hdf' :
        pass
    elif ftype == 'ASCII':
        pass
    else:
        raise NotImplementedError('ftype {} not implemented'.format(ftype))


