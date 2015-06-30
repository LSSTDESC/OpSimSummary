#!/usr/bin/env python
import pandas as pd


class fieldSimLib(object):

    '''
    Class to hold data corresponding to a particular fieldID (LIBID) of a
    SNANA SIMLIB file and methods
    '''

    def __init__(self, simlibstring):
        _tup = self.getFieldSimlib(simlibstring)
        self.meta = _tup[0]
        self.data = _tup[1]
        self._val = _tup[2]
        self.fieldID = self.meta['LIBID']
        self.validate()

    def getFieldSimlib(self, simlibstring):
        '''
        Basic constructor method to take a string corresponding to a
        simlib data corresponding to a single LIBID and parse it to
        metadata containing the properties of the field, a
        `~pandas.DataFrame` containing the data, and the string after the
        data
        '''

        # split into three parts
        header, data, footer = self.split_simlibString(simlibstring)

        # create the DataFrame
        simlibdat = self.simlibdata(data)

        # parse header to get header metadata and header fields
        header_metadata, header_fields = self.split_header(header)
        meta = self.libid_metadata(header_metadata)

        return meta, simlibdat, footer

    def validate(self):
        '''
        '''
        val = eval(self._val.split()[-1])
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
        split the string into header, footer, and data pieces
        '''
        lst = simlibString.split('MAG\n')
        header = lst[0]
        data, val = lst[1].split('END_LIBID')
        return header, data, val

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
