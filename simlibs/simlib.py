#!/usr/bin/env python
import pandas as pd
from cStringIO import StringIO


class fieldSimLib(object):

    from cStringIO import StringIO
    '''
    Class to hold data corresponding to a particular fieldID (LIBID) of a
    SNANA SIMLIB file and methods
    '''

    def __init__(self, simlibstring):
        self.getFieldSimlib(simlibstring)
        # self.meta = _tup[0]
        # self.data = _tup[1]
        # self._val = _tup[2]
        # self.fieldID = self.meta['LIBID']
        self.validate()

    @classmethod
    def getFieldSimlib(cls, simlibstring):
        '''
        Basic constructor method to take a string corresponding to a
        simlib data corresponding to a single LIBID and parse it to
        metadata containing the properties of the field, a
        `~pandas.DataFrame` containing the data, and the string after the
        data
        '''

        # split into three parts
        header, data, footer = cls.split_simlibString(simlibstring)

        # create the DataFrame
        cls.data = cls.simlibdata(data)

        # parse header to get header metadata and header fields
        header_metadata, header_fields = cls.split_header(header)
        cls.meta = cls.libid_metadata(header_metadata)
        cls.validate_string = footer
        cls.fieldID = cls.meta['LIBID']

        # return meta, simlibdat, footer

    @classmethod
    def validate(cls):
        '''
        '''
        val = eval(cls.validate_string.split()[-1])
        if int(cls.meta['LIBID']) != val:
            print 'LIBID value at beginning: ', cls.meta['LIBID']
            print 'LIBID value at the end', val
            raise ValueError('the LIBID values do not match')

        if len(cls.data) != cls.meta['NOBS']:
            print 'NOBS :', cls.meta['NOBS']
            print 'len(data) :', len(cls.data)
            raise ValueError('the number of observations recorded does not match size of data')

    @staticmethod
    def split_simlibString(simlibString):
        '''
        split the string into header, footer, and data pieces
        '''
        lst = simlibString.split('MAG')
        header = lst[0]
        data, val = lst[1].split('END_LIBID')
        index = data.index('\n')
        return header, data[index+1:], val

    @staticmethod
    def simlibdata(data):
        '''
        manipulate data piece to form pandas DataFrame object
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
    
    # @property
    def simlibData(self, fieldID):
       return self.simlibdicts[fieldID] 

    def getSimlibs(self, simlibFile):
        file_header, file_data, file_footer = self.read_simlibFile(simlibFile)
        simlibStrings = self.split_simlibStrings(file_data)
        mydict = dict()
        for strings in simlibStrings:
            s = fieldSimLib(strings)
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
