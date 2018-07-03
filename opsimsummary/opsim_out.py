"""
This module deals with representing the data in an OpSim output (to the extent
we will care about it). A description of the OpSim output can be found at
(opsim description)[https://www.lsst.org/scientists/simulations/opsim/summary-table-column-descriptions-v335]

In brief, we will use two tables from the OpSim output:
    - A Summary Table which has the desired information
    - A Proposal Table which contains a dictionary to interpreting the `propID` column
        of Summary.
"""
from __future__ import division, print_function, unicode_literals
__all__ = ['OpSimOutput']
import sys
import traceback
import numpy as np
import pandas as pd
import pdb
from sqlalchemy import create_engine


class OpSimOutput(object):
    """
    Class representing a subset of the output of the OpSim including
    information from the Summary and Proposal Tables with the subset taken over
    the proposals


    Attributes
    ----------
    opsimversion: {'lsstv3'|'sstf'|'lsstv4'}
        version of OpSim corresponding to the output format.
    summary : `pd.DataFrame`
        selected records from the Summary Table of pointings
    propIDDict : dict
        dictionary with strings as keys and integers used in the Summary
        Table to denote these proposals
    proposalTable : `pd.DataFrame`
        the propsal table in the output
    subset: string
        subset of proposals included in this class
    propIDs : list of integers
        integers corresponding to the subset selected through proposals
    zeroDDFDithers : bool, defaults to True
        if True, set dithers in DDF to 0, by setting ditheredRA,
        ditheredDec to fieldRA, fieldDec. This should only be used for
        opsimversion='lsstv3'. For opsimversion='sstf' or 'lsstv4', this
        will be set to False despite inputs, since this is already done, and
        cannot be done with the inputs.
    """
    def __init__(self, summary, propIDDict=None, proposalTable=None,
                 subset=None, propIDs=None, zeroDDFDithers=True,
                 opsimversion='lsstv3'):
        """
        Constructor for the `OpSimOutput` class 

        Parameters
        ----------
        summary: `pd.DataFrame
            a table describing a summary of observations which should have the
            following minimal columns (`ditheredRA`, `ditheredDec`, `expMJD`)
            as columns holding information on the pointing RA and dec of the
            telescope and the mjd of the observation and `obsHistID` a unique 
            index for each entry the name of the index. Other columns
            which are essential in the context of use in simulations are 
            (`FWHMeff`, `filtSkyBrightness`, `fiveSigmaDepth`, `propID`)
            describing the FWHM seeing, sky brightness in the filter of
            observation, and the five sigma depth of the observation and
            an integer index for different programs that an observation.
        propIDDict: dict, defaults to None
            dictionary with keys giving integers that identify proposals and
            values in the form of strings describing the program
        subset: {'wfd'|'ddf'|'combined'| '_all' | 'unique_all'}, defaults to
            `combined` a string that defines a subset of observations to be
            chosen based on the choice of proposals. `wfd` is the LSST `WFD`,
            `ddf` is the LSST `ddf`, `combined` is the combination of 
            `WFD` and `DDF` while, `unique_all` keeps all of the unique
            observations. `_all` is the entire table of of observations from
            OpSim outputs.
        proposalTable: `pd.DataFrame`
            modified table of proposals from the OpSim output
        propIDs: 
        zeroDDFDithers:
        opsimversion: {'lsstv3'|'sstf'|'lsstv4'} , defaults to 'lsstv3'
            a string to denote the version of OpSim outputs. `lsstv3`
            refers to outputs from OpSim version 3 (eg. enigma_1189, minion_1016).
            `sstf` refers to outputs from created at the start of the Observing
            Strategy Task Force available at the following [website](http://altsched.rothchild.me:8080),
            and `lsstv4` refers to outputs available from 
        """
        self.opsimversion = opsimversion
        self.allowed_subsets = self.get_allowed_subsets()
        self.subset = subset
        self.propIDDict = propIDDict
        self.proposalTable = proposalTable

        if opsimversion in ('sstf', 'lsstv4'):
            zeroDDFDithers = False
            ss = 'Warning: Input is zeroDDFDithers = True. But opsimversion is'
            ss += '{} for which this must be False. Setting to False and proceeding\n'.format(opsimversion)
            print(ss) 


        # Check `summary` does not have `nan`s
        if not self.validate_pointings(summary, opsimVars=None):
            print('summary table has nans, exiting\n')
            sys.exit(1)

        if zeroDDFDithers:
            ddfPropID = self.propIDDict['ddf']
            ddfidx = summary.query('propID == @ddfPropID').index
            summary.loc[ddfidx, 'ditheredRA'] = summary.loc[ddfidx, 'fieldRA']
            summary.loc[ddfidx, 'ditheredDec'] = summary.loc[ddfidx, 'fieldDec']
        self._opsimvars = None

        # Have a clear unambiguous ra, dec in radians following LSST convention
        if self.opsimVars['angleUnits'] == 'degrees':
            summary['_ra'] = np.radians(summary['ditheredRA'])
            summary['_dec'] = np.radians(summary['ditheredDec'])
            print('Changing units for {0} from {1}'.format(opsimversion, 'degrees'))
        elif self.opsimVars['angleUnits'] == 'radians':
            print('Keeping units for {0} from {1}'.format(opsimversion, 'radians'))
            summary['_ra'] = summary['ditheredRA']
            summary['_dec'] = summary['ditheredDec']
        else:
            raise ValueError('angle unit of ra and dec Columns not recognized\n')


        print('check that the format matches the expectation')
        if self.validate_pointings(summary, self.opsimVars):
            self.summary = summary
        else:
            raise AssertionError('Pointings are not in required format')
        self._propID = propIDs

    @staticmethod
    def validate_pointings(summary, opsimVars=None):
        """
        Validate a dataframe of pointings for further use. If
        `opsimVars` is `None` then only check that there are no `no.nan`s,
        else check that the table of pointings has the necessary format and
        units by checking that required columns indicated by `opsimVars` exist
        and have sensible values.

        Parameters
        ----------
        summary: `pd.DataFrame` of pointings
        opsimVars: dictionary, defaults to None 
            should be dictionary for each supported OpSim version availble
            from `OpSimOutput.get_opsimVariablesForVersion(opsimversion)`

        Returns
        -------
        Bool (True|False) But exits on False.
        """
        try:
            if opsimVars is not None:
                assert '_ra' in summary.columns
                assert '_dec' in summary.columns

                assert np.fabs(summary['_ra'].max()) <= 2.0 * np.pi
                assert np.fabs(summary['_dec'].min()) >= -1.0 * np.pi

            assert summary.isnull().values.any() == False
        except AssertionError:
            _, _, tb = sys.exc_info()
            traceback.print_tb(tb) # Fixed format
            tb_info = traceback.extract_tb(tb)
            filename, line, func, text = tb_info[-1]

            print ('pointings are not in required format')
            print('An error occurred on line {} in statement {}'.format(line, text))
            sys.exit(1)
        return True

    @property
    def opsimVars(self):
        """
        a set of opsim version dependent variables
        """
        if self._opsimvars is None:
            self._opsimvars = self.get_opsimVariablesForVersion(self.opsimversion)
        return self._opsimvars

    @staticmethod
    def get_dithercolumns(summary,
                          opsimversion='lsstv3',
                          method='default'):
        """
        get dithers
        Parameters
        ----------
        summary : `pd.DataFrame`
           indexed by `obsHistID` and having the columns `fieldRA`, `fieldDec`
        opsimversion : string, defaults to `lsstv3`
           version of the OpSim producing the database.
        method : string
           {'default'} only implemented
        """
        if method == 'default':

           df = summary[['fieldRA', 'fieldDec']]
           # ditheredRA = summary.fieldRA.values
           # ditheredDec = summary.fieldDec.values
           # df = pd.DataFrame(dict(obsHistID=summary.index.values,
           #                       ditheredRA=ditheredRA,
           #                       ditheredDec=ditheredDec))
           df.rename(columns=dict(fieldRA='ditheredRA',
                                  fieldDec='ditheredDec'), inplace=True)
        else:
            raise NotImplementedError('Only the default identity method has been implemented yet\n')
        return df

    @classmethod
    def fromOpSimDB(cls, dbname, subset='combined',
                    tableNames=('Summary', 'Proposal'),
                    propIDs=None, zeroDDFDithers=True,
                    opsimversion='lsstv3', dithercolumns=None,
                    add_dithers=False):
        """
        Class Method to instantiate this from an OpSim sqlite
        database output

        Parameters
        ----------
        dbname : string
            absolute path to database
        subset : string, optional, defaults to 'combined'
            one of {'_all', 'unique_all', 'wfd', 'ddf', 'combined'}
            determines a sequence of propIDs for selecting observations
            appropriate for the OpSim database in use
        propIDs : sequence of integers, defaults to None
            proposal ID values. If present, overrides the use of subset
        tableNames : tuple of strings, defaults to ('Summary', 'Proposal')
            names of tables read from the OpSim database
        zeroDDFDithers : bool, defaults to True
            if True, set dithers in DDF to 0, by setting ditheredRA,
            ditheredDec to fieldRA, fieldDec
        opsimversion: {'lsstv3'|'sstf'|'lsstv4'}
            version of OpSim corresponding to the output format.
        dithercolumns: `pd.DataFrame`, defaults to None
            a pandas dataframe with the columns `ditheredRA`, `ditheredDec` and
            index `obsHistID`
        add_dithers : Bool, defaults to False
            if add_dithers is True, the `ditheredRA`, `ditheredDec`
            columns will be removed and additional dithered columns
            either by `dithercolumns` or `get_dithercolumns` will
            be used.

        """
        print(tableNames)# = ('summary', 'Proposal')
        if opsimversion in ('sstf', 'lsstv4'):
            tableNames=('summaryAllProps', 'Proposal')

        # Because this is in the class method, I am using the staticmethod
        # rather than a property
        opsimVars = cls.get_opsimVariablesForVersion(opsimversion)

        # Check that subset parameter is legal
        allowed_subsets = cls.get_allowed_subsets()
        subset = subset.lower()
        if subset not in allowed_subsets:
            raise NotImplementedError('subset {} not implemented'.\
                                      format(subset))

        # Prepend the abs path with sqlite for use with sqlalchemy
        if not dbname.startswith('sqlite'):
            dbname = 'sqlite:///' + dbname
        print(' reading from database {}'.format(dbname))
        engine = create_engine(dbname, echo=False)

        # Read the proposal table to find out which propID corresponds to
        # the subsets requested
        proposals = pd.read_sql_table(tableNames[1], con=engine)
        propDict = cls.get_propIDDict(proposals, opsimversion=opsimversion)

        # Seq of propIDs consistent with subset
        _propIDs = cls.propIDVals(subset, propDict, proposals)
        # If propIDs and subset were both provided, override subset propIDs
        propIDs = cls._overrideSubsetPropID(propIDs, _propIDs)

        # Do the actual sql queries or table reads
        summaryTableName = opsimVars['summaryTableName']


        # Note OpSim version 4 has different names for the same variable
        # in the Proposal Table and Summary Table.
        propIDNameInSummary = opsimVars['propIDNameInSummary']
        if subset in ('_all', 'unique_all'):
            # In this case read everything (ie. table read)
            summary = pd.read_sql_table(summaryTableName, con=engine)

        elif subset in ('ddf', 'wfd', 'combined'):
            print('Doing combined subset here')
            # In this case use sql queries rather than reading the whole table
            # obtain propIDs in strings for sql queries
            pidString = ', '.join(list(str(pid) for pid in propIDs))
            print(pidString, subset)
            sql_query = 'SELECT * FROM {0} WHERE {1}'.format(summaryTableName,
                                                             propIDNameInSummary
                                                            )
            sql_query += ' in ({})'.format(pidString)

            # If propIDs were passed to the method, this would be used
            print(sql_query)
            summary = pd.read_sql_query(sql_query, con=engine)
        else:
            raise NotImplementedError()


        if cls.validate_pointings(summary, opsimVars=None):
            print('Raw readin works')


        # Standardize names
        replacedict = dict()
        replacedict[opsimVars['obsHistID']] = 'obsHistID'
        replacedict[opsimVars['propIDNameInSummary']] = 'propID'
        replacedict[opsimVars['expMJD']] = 'expMJD'
        replacedict[opsimVars['FWHMeff']] = 'FWHMeff'
        replacedict[opsimVars['filtSkyBrightness']] = 'filtSkyBrightness'
        summary = summary.rename(columns=replacedict)

        if cls.validate_pointings(summary, opsimVars=None):
            print('replacing names works')

        # Drop Duplicates
        if subset != '_all':
            # Drop duplicates unless this is to write out the entire OpSim
            summary = cls.dropDuplicates(summary, propDict, opsimversion)

        # Set Standard Index
        summary.set_index('obsHistID', inplace=True)

        if cls.validate_pointings(summary, opsimVars=None):
            print('dropping duplicates works')

        # At this stage the summary table is read in,
        # and the standard index is set

        # In `lsstv3` minion like baselines, the pointingRA are `ditheredRA` etc.
        # In `sstf` versions, the pointing coordinates are `fieldRA` etc.
        # in `lsstv4`, the pointing coordinates are unsupplied but `ditheredRA` etc.

        if 'ditheredra' not in list(x.lower() for x in summary.columns):
            # eg. has to be done in `lsstv4` and `sstf` unless supplied
            add_dithers = True
        if add_dithers:
            if dithercolumns is not None:
                # If provided with dithers in a dataFrame, use them
                # Check that dithercolumns are available in input
                assert 'ditheredRA' in dithercolumns.columns
                assert 'ditheredDec' in dithercolumns.columns
                assert 'obsHistID' == dithercolumns.index.name
                
                # if the column names already exist in the table remove them
                if 'ditheredra' in list(x.lower() for x in summary.columns):
                    del summary['ditheredRA'] 
                    del summary['ditheredDec'] 

                # Assumption : I have the dither columns in a `pd.DataFrame`
                # with minimal columns `ditheredRA` and `ditheredDec` and
                # index name `obsHistID` which indexes the visits in the 
                # Summary Table

                summary = summary.join(dithercolumns)

            elif add_dithers:
                # No dither column provided
                dithercolumns = cls.get_dithercolumns(summary[['fieldRA',
                                                           'fieldDec']],
                                                      opsimversion=opsimversion,
                                                      method='default')
                if cls.validate_pointings(dithercolumns, opsimVars=None):
                    print('dithercolumns good!')
                try:
                    print('before join number of rows is {}'.format(len(summary)))
                    print('before join, if nulls existed {}'.format(any(summary.isnull())))
                    summary = summary.join(dithercolumns)
                    print(len(summary), len(dithercolumns))
                    if cls.validate_pointings(summary, opsimVars=None):
                        print('join good!')
                except:
                    pass
            else:
                raise NotImplementedError('What did you do ?????')
        else:
            # let pass without further action
            pass

        
        if cls.validate_pointings(summary, opsimVars=None):
            print('joining dithers works')

        return cls(propIDDict=propDict,
                   summary=summary,
                   zeroDDFDithers=zeroDDFDithers,
                   proposalTable=proposals, subset=subset,
                   opsimversion=opsimversion)

    @staticmethod
    def dropDuplicates(df, propIDDict, opsimversion):
        """
        drop duplicates ensuring keeping identity of ddf visits

        Parameters
        ----------
        df : `pd.DataFrame`
        propIDDict : dict

        Returns
        -------
        `pd.DataFrame` with the correct propID and duplicates dropped
        """
        if opsimversion == 'sstf':
            return df

        # As duplicates are dropped in order, reorder IDs so that
        # DDF is lowest, WFD next lowest, everything else as is
        minPropID = df.propID.min()
        ddfID = propIDDict['ddf']
        wfdID = propIDDict['wfd']
        ddfPropID = minPropID - 1
        wfdPropID = minPropID - 2

        ddfmask = df.propID == ddfID
        wfdmask = df.propID == wfdID
        df.loc[ddfmask, 'propID'] = ddfPropID
        df.loc[wfdmask, 'propID'] = wfdPropID

        # drop duplicates keeping the lowest transformed propIDs so that all
        # WFD visits remain, DDF visits which were duplicates of WFD visits are
        # dropped, etc.

        df = df.drop_duplicates(subset='obsHistID', keep='first', inplace=False)
        #df = df.drop_duplicates(subset='obsHistID',
        #                                      keep='first')#.set_index('obsHistID')

        # reset the propIDs to values in the OpSim output
        ddfmask = df.propID == ddfPropID
        wfdmask = df.propID == wfdPropID
        df.loc[ddfmask, 'propID'] = ddfID
        df.loc[wfdmask, 'propID'] = wfdID
        df.sort_values(by='expMJD', inplace=True)
        return df


    @classmethod
    def fromOpSimHDF(cls, hdfName, subset='combined',
                     tableNames=('Summary', 'Proposal'),
                     propIDs=None):
        """
        Construct an instance of a subset of the OpSim
        Output from a serialization in the format of hdf

        Parameters
        ----------
        hdfName :
        subset :
        tableNames :
        propIDs :
        """
        raise NotImplementedError('Not quite working at this moment')
        allowed_subsets = cls.get_allowed_subsets()
        subset = subset.lower()
        if subset not in allowed_subsets:
            raise NotImplementedError('subset {} not implemented'.\
                      format(subset))
        # The hdf representation is assumed to be a faithful representation of
        # the OpSim output
        summarydf = pd.read_hdf(hdfName, key='Summary')

        if 'obsHistID' not in summarydf.columns:
            summarydf.reset_index(inplace=True)
            if 'obsHistID' not in summarydf.columns:
                raise NotImplementedError('obsHistID is not in columns')

        try:
            proposals = pd.read_hdf(hdfName, key='Proposal')
            print('read in proposal')
            propDict = cls.get_propIDDict(proposal)
            print('read in proposal')
            print(subset, propDict)
            _propIDs = cls.propIDVals(subset, propDict, proposals)
        except:
            print('Proposal not read')
            pass

        propIDs = cls._overrideSubsetPropID(propIDs, _propIDs)

        if propIDs is not None:
            if not isinstance(propIDs, list):
                propIDs = propIDs.tolist()
            print('propIDs', propIDs, type(propIDs), type(propIDs[0]))
            print('summarydf cols', summarydf.columns)
            query_str = 'propID == @propIDs'
            print('query_str', query_str)
            print(' Num entries ', len(summarydf))
            summary = summarydf.query(query_str)
        else:
            summary = summarydf
        if propIDs is None and subset not in ('_all', 'unique_all'):
            raise ValueError('No sensible propID and subset combination found')

        if subset != '_all':
            # Usually drop the OpSim duplicates
            summary.drop_duplicates(subset='obsHistID', inplace=True)

        summary.set_index('obsHistID', inplace=True)
        return cls(propIDDict=propDict, summary=summary,
                   proposalTable=proposals, subset=subset)

    @property
    def propIds(self):
        """
        list of values in propID Column of the Summary Table of OpSim
        to be considered for this class, either because they were directly
        provided or through the subset argument.
        """
        if self._propID is not None:
            return self._propID
        elif self.subset is not None and self.propIDDict is not None:
            return self.propIDVals(self.subset, self.propIDDict, self.proposalTable)

    def writeOpSimHDF(self, hdfName):
        """
        Serialize the OpSim output to hdf format in a welldefined way
        The output hdf file has two keys: 'Summary' and 'Proposal'
        """
        if self.subset != '_all':
            raise ValueError('Should be Done only for self.subset == _all')
        self.summary.to_hdf(hdfName, key='Summary', append=False)
        self.proposalTable.to_hdf(hdfName, key='Proposal', append=False)

    @staticmethod
    def _overrideSubsetPropID(propIDs, _propIDs):
        if propIDs is None:
            propIDs = _propIDs
        else:
            if np.asarray(propIDs).sort() != np.asarray(_propIDs).sort():
                raise Warning('argument propIDs and _propIDs do not match')
        return propIDs

    @staticmethod
    def get_allowed_subsets():
        return ('_all', 'ddf', 'wfd', 'combined', 'unique_all')

    @staticmethod
    def get_propIDDict(proposalDF, opsimversion='lsstv3'):
        """
        Return a dictionary with keys 'ddf', ad 'wfd' with the proposal IDs
        corresponding to deep drilling fields (ddf) and universal cadence (wfd) 

        Parameters
        ----------
        proposalDF : `pd.DataFrame`, mandatory
            a dataframe with the Proposal Table of the OpSim Run.
        opsimversion: {'lsstv3'|'sstf'|'lsstv4'}, defaults to 'lsstv3'
            version of opsim from which output is drawn
        Returns
        -------
        dictionary with keys 'wfd' and 'ddf' with values given by integers
            corresponding to propIDs for these proposals
        """
        oss_wfdName = 'wfd'
        oss_ddfName = 'ddf'

        df = proposalDF
        mydict = dict()
        if opsimversion == 'lsstv3':
            propName = 'propConf'
            propIDName = 'propID'
            ops_wfdname = 'conf/survey/Universal-18-0824B.conf'
            ops_ddfname = 'conf/survey/DDcosmology1.conf'
        elif opsimversion == 'sstf':
            propName = 'propName'
            propIDName = 'propId'
            ops_wfdname = 'WideFastDeep'
            ops_ddfname = 'Deep Drilling'
        elif opsimversion == 'lsstv4':
            propName = 'propName'
            propIDName = 'propId'
            ops_wfdname = 'WideFastDeep'
            ops_ddfname = 'DeepDrillingCosmology1'
        else:
            raise NotImplementedError('`get_propIDDict` is not implemented for this `opsimversion`')

        # Rename version based proposal names to internal values
        for idx, row in df.iterrows():
            # remember in enigma outputs, these came with `..` in the beginning
            if ops_wfdname in row[propName]:
                df.loc[idx, propName] = oss_wfdName
            elif ops_ddfname in row[propName]:
                df.loc[idx, propName] = oss_ddfName
            else:
                pass
        return dict(df.set_index(propName)[propIDName])

    @staticmethod
    def get_opsimVariablesForVersion(opsimversion='lsstv3'):
        if opsimversion == 'lsstv3':
            x = dict(summaryTableName='Summary',
                     obsHistID='obsHistID',
                     propName='propConf',
                     propIDName='propID',
                     propIDNameInSummary='propID',
                     ops_wfdname='conf/survey/Universal-18-0824B.conf',
                     ops_ddfname='conf/survey/DDcosmology1.conf',
                     expMJD='expMJD',
                     FWHMeff='FWHMeff',
                     pointingRA='ditheredRA',
                     pointingDec='ditheredDec',
                     filtSkyBrightness='filtSkyBrightness',
                     angleUnits='radians')
        elif opsimversion == 'sstf':
            x = dict(summaryTableName='SummaryAllProps',
                     obsHistID='observationId',
                     propName='propName',
                     propIDName='propId',
                     propIDNameInSummary='proposalId',
                     ops_wfdname='WideFastDeep',
                     ops_ddfname='Deep Drilling',
                     expMJD='observationStartMJD',
                     FWHMeff='seeingFwhmEff',
                     pointingRA='fieldRA',
                     pointingDec='fieldDec',
                     filtSkyBrightness='skyBrightness',
                     angleUnits='degrees')
        elif opsimversion == 'lsstv4':
            x = dict(summaryTableName='SummaryAllProps',
                     obsHistID='observationId',
                     propName='propName',
                     propIDName='propId',
                     propIDNameInSummary='proposalId',
                     ops_wfdname='WideFastDeep',
                     ops_ddfname='DeepDrillingCosmology1',
                     expMJD='observationStartMJD',
                     FWHMeff='seeingFwhmEff',
                     pointingRA='ditheredRA',
                     pointingDec='ditheredDec',
                     filtSkyBrightness='skyBrightness',
                     angleUnits='degrees')
        else:
            raise NotImplementedError('`get_propIDDict` is not implemented for this `opsimversion`')
        return x

    @staticmethod
    def propIDVals(subset, propIDDict, proposalTable):
        """
        Parameters: 
        ----------
        subset : string
            must be member of OpSimOutput.allowed_subsets()
        propIDDict : dictionary, mandatory
            must have subset as a key, and an integer or seq of ints
            as values
        proposalTable : `pd.DataFrame`
            Dataframe representing the proposal table in the OpSim datbase
            output

        Returns:
        -------
        list of propID values (integers) associated with the subset
        """
        if subset is None:
            raise ValueError('subset arg in propIDVals cannot be None')

        if subset.lower() in ('ddf', 'wfd'):
            return [propIDDict[subset.lower()]]
        elif subset.lower() == 'combined':
            return [propIDDict['ddf'], propIDDict['wfd']] 
        elif subset.lower() in ('_all', 'unique_all'):
            if proposalTable is not None:
                return proposalTable.propID.values
            else:
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
    print('This seems to have changed since first written, fixing not a priority')
    raise NotImplementedError('This seems to have changed since first written')
    if ftype == 'sqlite':
        dbname = 'sqlite:///' + fname
        engine = create_engine(dbname, echo=False)
        proposalTable =  pd.read_sql_table('Proposal', con=engine)

        # if subset == 'DDF':
        #    sql

    elif ftype == 'hdf' :
        pass
    elif ftype == 'ASCII':
        pass
    else:
        raise NotImplementedError('ftype {} not implemented'.format(ftype))
