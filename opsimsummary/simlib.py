#!/usr/bin/env python

"""
Module with functionality to represent SNANA simlib data. 
"""
from __future__ import division, print_function, unicode_literals
__all__ = ['SimlibMixin']
import os
import numpy as np
import subprocess
from io import StringIO, BytesIO
from collections import OrderedDict
import pandas as pd

class SimlibMixin(object):
    """
    Mixin for `SummaryOpsim` to provide the following additional functionality
    geared towards creating simlibs for SNANA.
    - Calculate additional columns for simlib either on a complete
        `OpSimOutput.summary` dataframe, or for a dataframe for a particular
        `patch` of sky (LIBID in the SNANA language).
    - Calculate variables required for SNANA simlib outside the Opsim data

    The parent class must have the following attributes:
        - user (default can be None)
        - host (default can be None)
        - telescope
        - survey 
        - pixelSize 
    In order to be able to write out the simlibs to disk, it should also have
    a method to provide a sequence (may be a generator) of fields, and
    `opsimtables`. The fields are instances of a class which has the following
    information

    `fieldID`
    `ra`
    `dec`
    `mwebv` (which may be set to a default value). The responsibility of
    selection of such fields and sorting out the correct requirements is
    of the parent class.
    """
    # Bad patch
    pixelSize = 0.2

    @property
    def simlibVars(self):
        simlibVars = self.get_simlibVars(user=self.user,
                                         host=self.host,
                                         pixelSize=self.pixelSize,
                                         telescope=self.telescope,
                                         survey=self.survey)
        self.user = simlibVars['user']
        self.host = simlibVars['host']
        self.pixelSize = simlibVars['pixelSize']
        self.telescope = simlibVars['telescope']
        self.survey = simlibVars['survey']

        return simlibVars

    @staticmethod
    def get_simlibVars(user=None, host=None, pixelSize=0.2, telescope='LSST',
                       survey='LSST'):
        """ Computes quantities only necessary for SNANA Simlib Calculation
        Parameters
        ----------
        user: string, optional, defaults to None
            user running the program, used in writing out SNANA simlibs only
            if None, the login name of the user is used.
        host: string, optional, defaults to None
            name of host machine, used only in writing out SNANA simlibs
            default of None assigns the output of `hostname` to this variable.
        pixelSize: float, optional, defaults to 0.2
            size of the pixel in arcseconds, defaults to 0.2 as appropriate
            for LSST
        survey: string, optional, defaults to 'LSST'
            name of survey, required only for writing out SNANA simlibs
        telescope: string, optional, defaults to 'LSST'
            name of survey, required only for writing out SNANA simlibs
        """
        # report a user name, either from a constructor parameter, or login name
        if user is None:
            user = os.getlogin()

        # report a host on which the calculations are done. either from
        # constructor parameters or from the system hostname utility 
        if host is None:
            proc = subprocess.Popen('hostname', stdout=subprocess.PIPE)
            host, err = proc.communicate()
        x = (('user', user),
             ('host', host),
             ('pixelSize', pixelSize),
             ('survey', survey),
             ('telescope', telescope))
        return OrderedDict(x)

    def _capitalizeY(self, x):
        """private helper method
        """
        # SNANA has y filter deonoted as Y. Can change in input files to SNANA
        # but more bothersome.
        if 'y' in x:
            return u'Y'
        else:
            return x

    def preprocess_lib(self, opsimtable):
        """
        preprocess the dataframe with data for a single SNANA simlib
        field (ie. data corresponding to a single libid) if necessary.

        Parameters
        ----------
        opsimtable : `pd.DataFrame` with required data from OpSim corresponding
            to a single field.
        """
        # reasonable guess that columns have not been added
        if 'simLibSkySig' not in opsimtable.columns:
            df  = self.add_simlibCols(opsimtable, pixelSize=self.pixelSize)
            df['filter'] = list(map(self._capitalizeY, opsimtable['filter']))
        else:
            df = opsimtable
        return df

    @staticmethod
    def add_simlibCols(opsimtable, pixelSize=0.2):
        """
        Parameters
        ----------
        opsimtable: `~pandas.DataFrame` object, mandatory
            containing an opsim Output of version X. The main requirements here
            are that the columns 'finSeeing', 'fiveSigmaDepth', and
            'filtSkyBrightness' are defined. If the opsim output has differently
            named variables or transformed variables, these should be changed to
            meet the criteria.
        pixelSize: float, units of arc sec, defaults to LSST value of 0.2
            pixel Size
    
    
        Returns
        -------
        DataFrame with additional columns of 'simLibPsf', 'simLibZPTAVG', and
        'simLibSkySig' 
    
        .. note :: This was written from a piece of f77 code by David
            Cinabro sent by email on May 26, 2015. 
        """
        if 'finSeeing' in opsimtable.columns:
            psfwidth = 'finSeeing'
        else:
            psfwidth = 'FWHMeff'
    
        opsim_seeing = opsimtable[psfwidth] # unit of arc sec sq
        # magsky is in units of mag/arcsec^2
        # opsim_maglim is in units of mag
        opsim_maglim = opsimtable['fiveSigmaDepth']
        opsim_magsky = opsimtable['filtSkyBrightness']
    
        # Calculate two variables that come up in consistent units:
        # term1  = 2.0 * opsim_maglim - opsim_magsky
    
        # Area of pixel in arcsec squared
        pixArea = pixelSize * pixelSize
    
        term1 = 2.0 * opsim_maglim - opsim_magsky # * pixArea
        # term2 = opsim_maglim - opsim_magsky
        term2 = - (opsim_maglim - opsim_magsky) # * pixArea
    
    
        # Calculate SIMLIB PSF VALUE
        opsimtable['simLibPsf'] = opsim_seeing /2.35 /pixelSize
       
        # 4 \pi (\sigma_PSF / 2.35 )^2
        area = (1.51 * opsim_seeing)**2.
        
        opsim_snr = 5.
        arg = area * opsim_snr * opsim_snr
    
        # Background dominated limit assuming counts with system transmission only
        # is approximately equal to counts with total transmission
        zpt_approx = term1 + 2.5 * np.log10(arg)
        # zpt_approx = 2.0 * opsim_maglim - opsim_magsky + 2.5 * np.log10(arg)
        # ARG again in David Cinabro's code
    
        val = -0.4 * term2
        # val = -0.4 * (opsim_magsky - opsim_maglim)
    
        tmp = 10.0 ** val
        # Additional term to account for photons from the source, again assuming
        # that counts with system transmission approximately equal counts with total
        # transmission.
    
        zpt_cor = 2.5 * np.log10(1.0 + 1.0 / (area * tmp))
        simlib_zptavg = zpt_approx + zpt_cor
        # ZERO PT CALCULATION 
        opsimtable['simLibZPTAVG'] = simlib_zptavg

    
        # SKYSIG Calculation
        npix_asec = 1. / pixelSize**2.
        opsimtable['simLibSkySig'] = np.sqrt((1.0 / npix_asec) \
        * 10.0 **(-0.4 * (opsim_magsky - simlib_zptavg)))
        return opsimtable

    def fieldheader(self, fieldID, ra, dec, opsimtable, mwebv=0.01):
        """
        Parameters
        ----------
        ra : degrees
        dec : degrees
        """
        # ra = np.degrees(self.ra(fieldID))
        # dec = np.degrees(self.dec(fieldID))
        nobs = len(opsimtable)
        s = '# --------------------------------------------' +'\n' 
        s += 'LIBID: {0:10d}'.format(fieldID) +'\n'
        tmp = 'RA: {0:+10.6f} DECL: {1:+10.6f}   NOBS: {2:10d} MWEBV: {3:5.2f}'
        tmp += ' PIXSIZE: {4:5.3f}'
        s += tmp.format(ra, dec, nobs, mwebv, self.pixelSize) + '\n'
        # s += 'LIBID: {0:10d}'.format(fieldID) + '\n'
        s += '#                           CCD  CCD         PSF1 PSF2 PSF2/1' +'\n'
        s += '#     MJD      IDEXPT  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG' + '\n'
        return s

    @staticmethod
    def fieldfooter(fieldID):
        
        s = 'END_LIBID: {0:10d}'.format(fieldID)
        s += '\n'
        return s
      
    def formatSimLibField(self, fieldID, opsimtable, sep=' '):

        opsimtable = self.preprocess_lib(opsimtable)
        y = ''
        for row in opsimtable.iterrows():
            data = row[1] # skip the index
            # print(data['filter'], type(data['filter']), list(data['filter']))
            #       MJD EXPID FILTER
            lst = ['S:',
                   "{0:5.4f}".format(data.expMJD),
                   "{0:10d}".format(data.obsHistID),
                   data['filter'],
                   "{0:5.2f}".format(1.),                  # CCD Gain
                   "{0:5.2f}".format(0.25),                # CCD Noise
                   "{0:6.2f}".format(data.simLibSkySig),   # SKYSIG
                   "{0:4.2f}".format(data.simLibPsf),      # PSF1
                   "{0:4.2f}".format(0.),                  # PSF2
                   "{0:4.3f}".format(0.),                  # PSFRatio
                   "{0:6.2f}".format(data.simLibZPTAVG),   # ZPTAVG
                   "{0:6.3f}".format(0.005),               # ZPTNoise 
                   "{0:+7.3f}".format(-99.)]               # MAG
            s = sep.join(lst)
            y += s + '\n'
        return y
    
    def simlibFieldasString(self, fh, fieldID, ra, dec, opsimtable, mwebv=0.1):


        #raise NotImplementedError("Has not been checked")
        # Write out the header for each field
        s = self.fieldheader(fieldID, ra, dec, opsimtable, mwebv=mwebv)
        # Write out the actual field
        s += self.formatSimLibField(fieldID, opsimtable, sep=' ')
        # Write out the footer for each field
        s += self.fieldfooter(fieldID)
        return s

    def simLibheader(self):
        sv = self.simlibVars
        user = sv['user']
        host = sv['host']
        telescope = sv['telescope']
        survey = sv['survey']
        # comment: I would like to generalize ugrizY to a sort but am not sure
        # of the logic for other filter names. so ducking for now
        s = 'SURVEY: {0:}    FILTERS: ugrizY  TELESCOPE: {1:}\n'.format(survey, telescope)
        s += 'USER: {0:}     HOST: {1:}\n'.format(user, host) 
        s += 'BEGIN LIBGEN\n'
        return s
    
    def simLibFooter(self, numFields):
        """
        """
        s = 'END_OF_SIMLIB:    {0:10d} ENTRIES'.format(numFields)
        return s


    def writeSimlib(self, filename, fields, opsimtables, comments='\n'):
            
        num_fields = 0
        with open(filename, 'w') as fh:
            # Write out the header to the simlib file
            simlib_header = self.simLibheader()
            fh.write(simlib_header)
            fh.write(comments)

            # Now write the actual simlib data to file
            for field, opsimtable in zip(fields, opsimtables):

                # obtain the set of field dependent parameters from `SynOpSim`
                fieldID = field.fieldID
                ra = field.ra
                dec = field.dec
                mwebv = field.mwebv

                fh.write(self.simlibFieldasString(self, fieldID, ra, dec,
                                                  opsimtable, mwebv=0.1))

                # Write out the header for each field
                # fh.write(self.fieldheader(fieldID, ra, dec, opsimtable,
                #                          mwebv=mwebv))
                # fh.write(self.formatSimLibField(fieldID, opsimtable))

                # Write out the footer for each field
                # fh.write(self.fieldfooter(fieldID))
                num_fields += 1

            # Now write out the footer to the entire simlib file 
            simlib_footer = self.simLibFooter(num_fields)
            fh.write(simlib_footer)

class FieldSimlib(object):
    """
    Class to hold data corresponding to a particular fieldID (LIBID) of a
    SNANA SIMLIB file and methods. The fieldSimlib class for a particular field
    has the following attributes, and may be instantiated by supplying these, 
    or can be conveniently constructed from the string corresponding to the
    description of this data in an SNANA simlib file using the constructor
    fromSimlibString:


    Parameters
    ----------
    simlibdata : a pandas dataFrame
    simlib_meta : a dictionary


    Attributes
    ----------
    fieldID : int
        a unique integer identifying the field of view. Different pointings
        of the same field at different times (but possibly with dithers) are
        associated with the same fieldID
    meta : dict
        metadata associated with the field, which has at least the following
        keys:
        LIBID, RA, DECL, MWEBV, NOBS, PIXSIZE
    data : `~pd.DataFrame` object with the observations and having at least the
        following columns: 'MJD', 'IDEXPT', 'FLT', 'GAIN', 'NOISE', 'SKYSIG',
        'PSF1', 'PSF2', 'PSFRatio', 'ZPTAVG', 'ZPTERR', 'MAG']. The meanings of
        these columns are discussed in the SNANA manual in the sub-section
        'The 'SIMLIB' Observing file (4.7)
    """
    def __init__(self, simlibdata, simlib_meta):
        """
        Instantiate the class from the basic data
        """ 

        self.data = simlibdata 
        self.meta = simlib_meta
        self.fieldID = self.meta['LIBID']

    @classmethod
    def fromSimlibString(cls, simlibstring):
        '''
        Basic constructor method to take a string corresponding to a
        field simlib data corresponding to a single LIBID and parse it to
        metadata containing the properties of the field, a
        `~pandas.DataFrame` containing the data, and the string after the
        data

        Parameters
        ----------
        simlibstring : string, mandatory
        '''

        # split into three parts
        header, data, footer = cls.split_simlibString(simlibstring)

        # create the DataFrame
        clsdata = cls.simlibdata(data)

        # parse header to get header metadata and header fields
        header_metadata, header_fields = cls.split_header(header)
        clsmeta = cls.libid_metadata(header_metadata)

        # Instantiate the class and make sure it works
        myclass = cls(simlibdata=clsdata, simlib_meta=clsmeta) 
        myclass.validate(footer)
        return myclass 

    def validate(self, validate_string):
        """
        Validate the interpretation of the field simlib data from a field simlib
        string by checking 1. the LIBID at the end of the string matches the
        one at the beginnin (ie. somehow multiple fields have not been read in)
        2. the number of rows of the data for this field simlib matches the
        number of observations recorded in the metadata as NOBS

        Parameters
        ----------
        validate_string : string, mandatory
            footer obtained by splitting the simlib corresponding to the field
            usually of the form 
        """

        val = eval(validate_string.split()[-1])
        if int(self.meta['LIBID']) != val:
            print('LIBID value at beginning: ', self.meta['LIBID'])
            print('LIBID value at the end', val)
            raise ValueError('the LIBID values do not match')

        if len(self.data) != self.meta['NOBS']:
            print('NOBS :', self.meta['NOBS'])
            print('len(data) :', len(self.data))
            raise ValueError('the number of observations recorded does not'
                             'match size of data')

    @staticmethod
    def split_simlibString(simlibString):
        '''
        split the string corresponding to a simlib file into header, footer,
        and data pieces

        Parameters
        ----------
        simlibString : string
        '''
        lst = simlibString.split('MAG')
        header = lst[0]
        data, val = lst[1].split('END_LIBID')
        index = data.index('\n')
        return header, data[index+1:], val

    @staticmethod
    def simlibdata(data):
        '''
        manipulate string in the simlibstring to form pandas DataFrame object


        Parameters
        ----------

        data : data
        '''
        fhandle = StringIO(data)
        df = pd.read_csv(fhandle, delimiter="\s+",
                         names=['trash', 'MJD', 'IDEXPT', 'FLT', 'GAIN',
                                'NOISE', 'SKYSIG', 'PSF1', 'PSF2',
                                'PSFRatio', 'ZPTAVG', 'ZPTERR', 'MAG'])
        del df['trash']
        return df

    @staticmethod
    def split_header(header):
        '''
        split header string into metadata and field names

        Parameters
        ----------

        header : header
        '''

        lines = header.split('\n')
        header_metadata = []
        header_fields = []
        for line in lines:
            if line.startswith('#'):
                header_fields.append(line[1:])
            else:
                header_metadata += line.split()
        return header_metadata, header_fields

    @staticmethod
    def libid_metadata(header_metadata):
        '''
        parse header metadata string into a dictionary

        Parameters
        ----------
        header_metadata : header
        '''

        # Even index values 0, 2, 4 are keys
        # remove ':' char at end
        keys = list(map(lambda x: x[:-1], header_metadata[0::2]))

        # odd index values are floats or ints
        vals = list(map(eval, header_metadata[1::2]))

        return dict(zip(keys, vals))

class Simlib(object):

    def __init__(self, simlibDict, simlibMetaData=None):


        self.simlibDict = simlibDict
        self.fieldIDs = self.simlibDict.keys()
        if simlibMetaData is not None:
            self.meta = simlibMetaData

    @classmethod
    def fromSimlibFile(cls, simlibFileName):
        '''
        Constructor for class using an ASCII 

        Parameters
        ----------
        simlibFileName: string, mandatory
            absolute path to SNANA simlib file

        Returns
        -------
        dictionary with LIBIDs (fieldIDs) of type int as keys, and the
        corresponding FieldSimlib objects as values

        Examples
        --------
        >>> sl = Simlib.fromSimlibFile(simlibFileName)
        '''

        file_header, file_data, file_footer = cls.read_simlibFile(simlibFileName)
        mydict = cls.getSimlibs(file_data)
        meta = cls.simlibMetaData(file_header) 
        cls = cls(simlibDict=mydict, simlibMetaData=meta)
        cls.validate(file_footer)

        return cls

    def validate(self, file_footer):
        '''
        '''
        numberlist = list(filter(lambda x: x.isdigit(), file_footer.split()))
        if len(numberlist) !=1:
            raise ValueError('There should only be one integer in the footer')
        numLibId = int(numberlist[0])

        if numLibId != len(self.fieldIDs):
            raise ValueError('The number of fieldIDs is in the simlib does not match the number stated in footer')

        return

    @staticmethod
    def simlibMetaData(simlibHeader):
        '''
        parse the string corresponding to the header of a SNANA simlib file to
        get the simlib MetaData, stored in the form of a string valued
        dictionary with the following keys:
        'SURVEY', 'FILTERS', 'HOST', 'USER', 'COMMENT'

        Parameters
        ----------
        simlibHeader: string corresponding the header of an SNANA simlib file as
            parsed by cls.read_simlibFile.
        Returns
        -------
        dictionary of keys above and values.
        '''

        comments = []
        fields = []

        lines = simlibHeader.split('\n')
        for line in lines:
            if line.startswith('COMMENT') or line.startswith('#'):
                comments.append(line)
            else:
                fields.append(line)
        ss = ' '.join(fields)
        words = ss.split()
        keys = list(map(lambda x: x[:-1], words[0::2]))
        vals = words[1::2]
        if len(keys) != len(vals):
            raise ValueError('the numberof fields in dict should match vals')

        meta = dict(zip(keys, vals))
        meta['COMMENTS'] = '\n'.join(comments)

        return meta


    def simlibData(self, fieldID):
       return self.simlibDict[fieldID].data

    @classmethod
    def getSimlibs(cls, file_data):
    # def getSimlibs(cls, simlibFile):
        # file_header, file_data, file_footer = cls.read_simlibFile(simlibFile)
        simlibStrings = cls.split_simlibStrings(file_data)
        mydict = dict()
        for strings in simlibStrings:
            s = FieldSimlib.fromSimlibString(strings)
            mydict[s.fieldID] = s

        return mydict 
            

    @staticmethod
    def read_simlibFile(simlibfile):
    
        # slurp into a string
        with open(simlibfile) as f:
            ss = f.read()
        
        # split into header, footer and data
        fullfile = ss.split('BEGIN LIBGEN')
        file_header = fullfile[0]
        if 'END_OF_SIMLIB' in ss:
            data, footer = fullfile[1].split('END_OF_SIMLIB')
        else:
            data = fullfile[1]
            footer = ''

        return file_header, data, footer

    @staticmethod
    def split_simlibStrings(simlibStrings):
        simlibs = simlibStrings.split('\nLIBID')[1:]
        simlibs = map(lambda x: 'LIBID' + x.split('# -')[0], simlibs)
        return simlibs
