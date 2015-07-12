#!/usr/bin/env python
import pandas as pd
from cStringIO import StringIO


class fieldSimLib(object):

    from cStringIO import StringIO
    '''
    Class to hold data corresponding to a particular fieldID (LIBID) of a
    SNANA SIMLIB file and methods. The fieldSimlib class for a particular field
    has the following attributes, and may be instantiated by supplying these, 
    or can be conveniently constructed from the string corresponding to the
    description of this data in an SNANA simlib file using the constructor
    fromSimlibString:

    Attributes
    ----------
    fieldID: int
        a unique integer identifying the field of view. Different pointings
        of the same field at different times (but possibly with dithers) are
        associated with the same fieldID
    meta: dict
        metadata associated with the field, which has at least the following
        keys:
        LIBID, RA, DECL, MWEBV, NOBS, PIXSIZE
    data: `~pd.DataFrame` object with the observations and having at least the
        following columns: 'MJD', 'IDEXPT', 'FLT', 'GAIN', 'NOISE', 'SKYSIG',
        'PSF1', 'PSF2', 'PSFRatio', 'ZPTAVG', 'ZPTERR', 'MAG']. The meanings of
        these columns are discussed in the SNANA manual in the sub-section
        'The 'SIMLIB' Observing file (4.7)
    '''
    def __init__(self, simlibdata, simlib_meta):

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
        '''

        # split into three parts
        header, data, footer = cls.split_simlibString(simlibstring)

        # create the DataFrame
        clsdata = cls.simlibdata(data)

        # parse header to get header metadata and header fields
        header_metadata, header_fields = cls.split_header(header)
        clsmeta = cls.libid_metadata(header_metadata)
        # cls.validate_string = footer
        # clsfieldID = cls.meta['LIBID']
        myclass = cls(simlibdata=clsdata, simlib_meta=clsmeta) 
        myclass.validate(footer)
        return myclass 

    def validate(self, validate_string):
        '''
        Validate the interpretation of the field simlib data from a field simlib
        string by checking 
        1. the LIBID at the end of the string matches the one at the beginning
            (ie. somehow multiple fields have not been read in)
        2. the number of rows of the data for this field simlib matches the
        number of observations recorded in the metadata as NOBS

        Parameters
        ----------
        validate_string: string, mandatory
            footer obtained by splitting the simlib corresponding to the field
            usually of the form 
        '''
        val = eval(validate_string.split()[-1])
        if int(self.meta['LIBID']) != val:
            print 'LIBID value at beginning: ', self.meta['LIBID']
            print 'LIBID value at the end', val
            raise ValueError('the LIBID values do not match')

        if len(self.data) != self.meta['NOBS']:
            print 'NOBS :', self.meta['NOBS']
            print 'len(data) :', len(self.data)
            raise ValueError('the number of observations recorded does not match size of data')

    @staticmethod
    def split_simlibString(simlibString):
        '''
        split the string corresponding to a simlib file into header, footer,
        and data pieces
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
        '''

        # Even index values 0, 2, 4 are keys
        # remove ':' char at end
        keys = map(lambda x: x[:-1], header_metadata[0::2])

        # odd index values are floats or ints
        vals = map(eval, header_metadata[1::2])

        return dict(zip(keys, vals))

class simlib(object):

    def __init__(self, simlibFile):
        self.simlibdicts = self.getSimlibs(simlibFile)
        self.fieldIDs = self.simlibdicts.keys()
    
    def simlibData(self, fieldID):
       return self.simlibdicts[fieldID] 

    def getSimlibs(self, simlibFile):
        file_header, file_data, file_footer = self.read_simlibFile(simlibFile)
        simlibStrings = self.split_simlibStrings(file_data)
        mydict = dict()
        for strings in simlibStrings:
            s = fieldSimLib.fromSimlibString(strings)
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
        data, footer = fullfile[1].split('END_OF_SIMLIB')

        return file_header, data, footer

    @staticmethod
    def split_simlibStrings(simlibStrings):
        simlibs = simlibStrings.split('\nLIBID')[1:]
        simlibs = map(lambda x: 'LIBID' + x.split('# -')[0], simlibs)
        return simlibs
