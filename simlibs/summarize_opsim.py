#!/usr/bin/env python 
import os
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
import matplotlib.pyplot as plt


def add_simlibCols(opsimtable, pixSize=0.2):
    '''
    Parameters
    ----------
    opsimtable: `~pandas.DataFrame` object, mandatory
        containing an opsim Output of version X. The main requirements here
        are that the columns 'finSeeing', 'fiveSigmaDepth', and
        'filtSkyBrightness' are defined. If the opsim output has differently
        named variables or transformed variables, these should be changed to
        meet the criteria.
    pixSize: float, optional, defaults to LSST value of 0.2
        pixel Size in units of arc seconds.


    Returns
    -------
    DataFrame with additional columns of 'simLibPsf', 'simLibZPTAVG', and
    'simLibSkySig' 
    '''
    
    opsim_seeing = opsimtable['finSeeing']
    opsim_maglim = opsimtable['fiveSigmaDepth']
    opsim_magsky = opsimtable['filtSkyBrightness']

    # Calculate SIMLIB PSF VALUE
    opsimtable['simLibPsf'] = opsim_seeing /2.35 /pixSize
   
    area = (1.51 * opsim_seeing)**2.
    
    opsim_snr = 5.
    arg = area * opsim_snr * opsim_snr
    zpt_approx = 2.0 * opsim_maglim - opsim_magsky + 2.5 * np.log10(arg)
    # ARG again in David Cinabro's code
    val = -0.4 * (opsim_magsky - opsim_maglim)
    tmp = 10.0 ** val
    zpt_cor = 2.5 * np.log10(1.0 + 1.0 / (area * tmp))
    simlib_zptavg = zpt_approx + zpt_cor
    # ZERO PT CALCULATION 
    opsimtable['simLibZPTAVG'] = simlib_zptavg
    
    #SKYSIG Calculation
    npix_asec = 1./ pixSize**2.
    opsimtable['simLibSkySig'] = ((1.0/ npix_asec)*10.0 **(-0.4 * (opsim_maglim - simlib_zptavg)))**0.5
    return opsimtable

class SummaryOpsim(object):
    
    
    def __init__(self, summarydf, user=None, host=None, survey='LSST',
                 telescope='LSST', pixSize=0.2):

        import os 
        import subprocess

        self.df = summarydf.copy(deep=True)
        if 'simLibSkySig' not in self.df.columns:
            self.df  = add_simlibCols(self.df)

        # SNANA has y filter deonoted as Y. Can change in input files to SNANA
        # but more bothersome.
        def capitalizeY(x):
            if 'y' in x:
                return u'Y'
            else:
                return x

        self.df['filter'] = map(capitalizeY, self.df['filter']) 
        self._fieldsimlibs = self.df.groupby(by='fieldID')
        self.fieldIds = self._fieldsimlibs.groups.keys()

        
        # report a user name, either from a constructor parameter, or login name
        if user is None:
            user = os.getlogin()
        self.user = user

        # report a host on which the calculations are done. either from
        # constructor parameters or from the system hostname utility 
        if host is None:
            proc = subprocess.Popen('hostname', stdout=subprocess.PIPE)
            host, err = proc.communicate()
        self.host = host

        self.telescope = telescope
        self.pixelSize = pixSize
        self.survey = survey

    @classmethod
    fromOpSimDB(cls, opSimDB, sql_query):

        from sqlalchemy import create_engine
        import pandas as pd

        engine = create_engine(opSimDB)
        summary = pd.read_sql_table('Summary', engine)
        selected  = summary.query(sql_query)

        return cls(selected)
        

        
    # @property
    def coords(self):

        ra = map(lambda x: self.ra(x), self.fieldIds)
        dec = map(lambda x: self.dec(x), self.fieldIds)

        return ra, dec
    def cadence_Matrix(self, fieldID, sql_query='night < 366',
                     Filters=[u'u', u'g', u'r', u'i', u'z', u'Y'],
                     nightMin=0, nightMax=365, observedOnly=False):
    
        # group on filter and timeindex (night)
        grouping_keys = ['filter', 'night']
        grouped = self.simlib(fieldID).query(sql_query).groupby(grouping_keys)

        # tuples of keys
        filts, nights = zip( *grouped.groups.keys())

        # number of Observations in each group
        numObs = grouped.apply(len).values

        # Create a new dataFrame with nights, Filters, numObs as cols
        cadence_dict = dict()
        cadence_dict['Filters'] = list(filts)
        cadence_dict['night'] = list(nights)

        # If observedOnly: set values above 1 to 1 
        if observedOnly:
            numObs = np.where(np.array(list(numObs)), 1, 0)
        
        cadence_dict['numObs'] = list(numObs)

        # pivot dataFrame to occupation numbers
        Matrix = pd.DataFrame(cadence_dict).pivot('night', 'Filters', 'numObs')


        # First make sure all filters are represented
        for filt in Filters:
            if filt not in Matrix.columns:
                Matrix[filt] = np.nan

        # reorder filters to u,g,r,i,z,y
        M = Matrix[Filters]
        # Extend to all values in plot
        ss = pd.Series(np.arange(nightMin, nightMax))
        Matrix = M.reindex(ss, fill_value=np.nan)

        return Matrix

    def cadence_plot(self, fieldID, sql_query='night < 366',
                     Filters=[u'u', u'g', u'r', u'i', u'z', u'Y'],
                     nightMin=0, nightMax=365, deltaT=5., observedOnly=False,
                     title=True, title_text=None, colorbar=True,
                     colorbarMin=0.):
        '''
        produce a cadence plot that shows the filters and nights observed in
        some subset of the opsim output time span for a field.


        Parameters
        ----------
        fieldID: integer, mandatory
            ID corresponding to the LSST Field as recorded as fieldID in OpSIM
            output
        sql_query: string, optional, defaults to obtaining the first season
            a sql query to select observations within the opsim summary object 
            associated with the field 
        Filters: list of strings, optional, defaults to LSST ugrizY
            a list of strings corresponding to filter names.

        '''


        Matrix = self.cadence_Matrix(fieldID, sql_query=sql_query,
                                Filters=Filters, nightMin=nightMin,
                                nightMax=nightMax, observedOnly=observedOnly)

        if observedOnly:
            axesImage = plt.matshow(Matrix.transpose(), aspect='auto',
                                    cmap=plt.cm.gray_r, vmin=colorbarMin,
                                    vmax=1.)
        else:
            axesImage = plt.matshow(Matrix.transpose(), aspect='auto',
                                    cmap=plt.cm.gray_r, vmin=colorbarMin)


        # setup matplotlib figure and axis objects to manipulate
        ax = axesImage.axes 

        # yticks: annotate with filter names
        # Note that it is also possible to get this from the DataFrame
        # by Matrix.columns, but the columns are sorted according to the order
        # in Filters
        
        ax.set_yticklabels(['0'] +Filters, minor=False)

        # Positiion x ticks at the bottom rather than top
        ax.xaxis.tick_bottom()

        # Add a grid 
        minorxticks = ax.set_xticks(np.arange(0, nightMax - nightMin,
                                              deltaT), minor=True)
        # Hard coding this
        minoryticks = ax.set_yticks(np.arange(-0.5,5.6,1), minor=True)
        ax.set_adjustable('box-forced')
        ax.grid(which='minor')


        # If nightMin is different from 0, translate xticks
        # xticks = ax.get_xticks()
        # xticklabels = ax.get_xticklabels()
        # tmp = [item.get_text() for item in xticklabels]
        # xticklabel = np.array(map(eval, tmp[1:-1])) + nightMin
        # ax.set_xticks(xticks[1:-1])
        # ax.set_xticklabels(xticklabel)

        # Set a title with information about the field
        
        if title:
            # Format field Info from attributes
            t_txt = 'fieldID: {:0>2d} (ra: {:+3f} dec: {:+3f})'
            t_txt = t_txt.format(fieldID, self.ra(fieldID), self.dec(fieldID))

            # if title_text is supplied use that instead
            if title_text is not None:
                t_txt = title_text
            ax.set_title(t_txt)

        # Set a colorbar
        if colorbar:
            plt.colorbar(orientation='horizontal')

        # Get the figure object
        fig = ax.figure

        return fig


    def showFields(self, ax=None, marker=None, **kwargs):
        ra = np.degrees(self.coords()[0])
        dec = self.coords()[1]

        x = np.remainder(ra + 360., 360.)
        ind  = x > 180.
        x[ind] -=360.

        ra = np.radians(x)
        # print ra
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='mollweide')
        if marker is None:
            marker = 'o'
        ax.scatter(ra, dec, marker=marker, **kwargs)
        ax.grid(True)
        fig  = ax.figure
        return fig


    def simlib(self, fieldID):

        return self._fieldsimlibs.get_group(fieldID)
    def ra(self, fieldID):
        ravals = np.unique(self.simlib(fieldID).fieldRA.values)
        if len(ravals)==1:
            return ravals[0]
        else:
            raise ValueError('The fieldDec of this group seems to not be unique\n')
        
    def dec(self, fieldID):
        decvals = np.unique(self.simlib(fieldID).fieldDec.values)
        if len(decvals)==1:
            return decvals[0]
        else:
            raise ValueError('The fieldDec of this group seems to not be unique\n')
    def meta(self, fieldID):
        meta = {}
        meta['LIBID'] = fieldID
        meta['dec'] = self.dec(fieldID)
        return meta
    
    def fieldheader(self, fieldID):
        
        ra = self.ra(fieldID)
        dec = self.dec(fieldID)
        mwebv = 0.01
        pixSize = self.pixelSize 
        nobs = len(self.simlib(fieldID))
        s = '# --------------------------------------------' +'\n' 
        s += 'LIBID: {0:10d}'.format(fieldID) +'\n'
        tmp = 'RA: {0:+10.6f} DECL: {1:+10.6f}   NOBS: {2:10d} MWEBV: {3:5.2f}'
        tmp += ' PIXSIZE: {4:5.3f}'
        s += tmp.format(ra, dec, nobs, mwebv, pixSize) + '\n'
        # s += 'LIBID: {0:10d}'.format(fieldID) + '\n'
        s += '#                           CCD  CCD         PSF1 PSF2 PSF2/1' +'\n'
        s += '#     MJD      IDEXPT  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG' + '\n'
        return s
    
    def fieldfooter(self, fieldID):
        
        s = 'END_LIBID: {0:10d}'.format(fieldID)
        s += '\n'
        return s
        
    def formatSimLibField(self, fieldID, sep=' '):
    
        opSimSummary = self.simlib(fieldID)
        y =''
        for row in opSimSummary.iterrows():
            data = row[1] # skip the index
            #       MJD EXPID FILTER 
            lst = ['S:',
                   "{0:5.3f}".format(data.expMJD),
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
    
    def writeSimLibField(self, fieldID):
        s = self.fieldheader(fieldID)
        s += self.formatSimLibField(fieldID, sep=' ')
        s += self.footer(fieldID)
        return s

    def simLibheader(self): #, user=None, host=None, survey='LSST', telescope='LSST'):
        # comment: I would like to generalize ugrizY to a sort but am not sure
        # of the logic for other filter names. so ducking for now
        s = 'SURVEY: {0:}    FILTERS: ugrizY  TELESCOPE: {1:}\n'.format(self.survey, self.telescope)
        s += 'USER: {0:}     HOST: {1:}\n'.format(self.user, self.host) 
        s += 'BEGIN LIBGEN\n'
        return s
    
    def simLibFooter(self):
        """
        """
        s = 'END_OF_SIMLIB:    {0:10d} ENTRIES'.format(len(self.fieldIds))
        return s


    def writeSimlib(self, filename, comments=''):
        with open(filename, 'w') as fh:
            simlib_header = self.simLibheader()
            simlib_footer = self.simLibFooter()
            fh.write(simlib_header)
            fh.write(comments)
            # fh.write('BEGIN LIBGEN\n')
            # fh.write('\n')
            for fieldID in self.fieldIds:
                fh.write(self.fieldheader(fieldID))
                fh.write(self.formatSimLibField(fieldID))
                fh.write(self.fieldfooter(fieldID))
            fh.write(simlib_footer)


### class simlibInfo(object):
### 
### 
###     def __init__(self, simlibFile):
###         self._file = simlibFile
###         _tup = self.read_simlibFile()
###         self.fileHeader = _tup[0]
###         self.fileData = _tup[1]
###         self.fileFooter = _tup[2]
###         self._simlibs = self.split_simlibFields()
###         self.simlibs = self.get_data()
###     
### 
###     def split_simlibFields(self):
###         simlibs = self._tup[1].split('\n# --------------------------------------------\n')[1:]
###         return simlibs
### 
###     def read_simlibFile(self):
###     
###         # slurp into a string
###         with open(self._file) as f:
###             ss = f.read()
###             
###         # split into header, footer and data
###         fullfile = ss.split('BEGIN LIBGEN')
###         file_header = fullfile[0]
###         data, footer = fullfile[1].split('END_OF_SIMLIB')
###         return file_header, data, footer
### 
### class SummaryOpsimG(object):
###     
###     
###     def __init__(self, summarydf=None, simlibFile=None, user=None, host=None,
###             survey='LSST', telescope='LSST', pixSize=0.2):
### 
###         import os 
###         import subprocess
### 
###         def guessContext(summarydf, simlibFile ):
###             if summarydf is None and simlibFile is None:
###                 raise ValueError('summarydf or simlibFile have to be supplied\n')
###         
###             if not (summarydf is None or simlibFile is None):
###                 raise ValueError('summarydf and simlibFile supplied together\n')
### 
###             if summarydf is None:
###                 self.context = 'SNANAsimlib'
###             if simlibFile is None:
###                 self.context = 'opSimOut'
### 
###             return 
### 
###         if summarydf:
###             self._df = summarydf.copy(deep=True)
###             if 'simLibSkySig' not in self.df.columns:
###                 self._df  = add_simlibCols(self.df)
### 
###         # SNANA has y filter deonoted as Y. Can change in input files to SNANA
###         # but more bothersome.
###         def capitalizeY(x):
###             if 'y' in x:
###                 return u'Y'
###             else:
###                 return x
### 
###         self.df['filter'] = map(capitalizeY, self.df['filter']) 
###         self._fieldsimlibs = self.df.groupby(by='fieldID')
###         self.fieldIds = self._fieldsimlibs.groups.keys()
### 
###         
###         # report a user name, either from a constructor parameter, or login name
###         if user is None:
###             user = os.getlogin()
###         self.user = user
### 
###         # report a host on which the calculations are done. either from
###         # constructor parameters or from the system hostname utility 
###         if host is None:
###             proc = subprocess.Popen('hostname', stdout=subprocess.PIPE)
###             host, err = proc.communicate()
###         self.host = host
### 
###         self.telescope = telescope
###         self.pixelSize = pixSize
###         self.survey = survey
###     def simlib(self, fieldID):
### 
###         return self._fieldsimlibs.get_group(fieldID)
###     def ra(self, fieldID):
###         ravals = np.unique(self.simlib(fieldID).fieldRA.values)
###         if len(ravals)==1:
###             return ravals[0]
###         else:
###             raise ValueError('The fieldDec of this group seems to not be unique\n')
###         
###     def dec(self, fieldID):
###         decvals = np.unique(self.simlib(fieldID).fieldDec.values)
###         if len(decvals)==1:
###             return decvals[0]
###         else:
###             raise ValueError('The fieldDec of this group seems to not be unique\n')
###     def meta(self, fieldID):
###         meta = {}
###         meta['LIBID'] = fieldID
###         meta['dec'] = self.dec(fieldID)
###         return meta
###     
###     def fieldheader(self, fieldID):
###         
###         ra = self.ra(fieldID)
###         dec = self.dec(fieldID)
###         mwebv = 0.01
###         pixSize = self.pixelSize 
###         nobs = len(self.simlib(fieldID))
###         s = '# --------------------------------------------' +'\n' 
###         s += 'LIBID: {0:10d}'.format(fieldID) +'\n'
###         tmp = 'RA: {0:+10.6f} DECL: {1:+10.6f}   NOBS: {2:10d} MWEBV: {3:5.2f}'
###         tmp += ' PIXSIZE: {4:5.3f}'
###         s += tmp.format(ra, dec, nobs, mwebv, pixSize) + '\n'
###         # s += 'LIBID: {0:10d}'.format(fieldID) + '\n'
###         s += '#                           CCD  CCD         PSF1 PSF2 PSF2/1' +'\n'
###         s += '#     MJD      IDEXPT  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG' + '\n'
###         return s
###     
###     def fieldfooter(self, fieldID):
###         
###         s = 'END_LIBID: {0:10d}'.format(fieldID)
###         s += '\n'
###         return s
###         
###     def formatSimLibField(self, fieldID, sep=' '):
###     
###         opSimSummary = self.simlib(fieldID)
###         y =''
###         for row in opSimSummary.iterrows():
###             data = row[1] # skip the index
###             #       MJD EXPID FILTER 
###             lst = ['S:',
###                    "{0:5.3f}".format(data.expMJD),
###                    "{0:10d}".format(data.obsHistID),
###                    data['filter'], 
###                    "{0:5.2f}".format(1.),                  # CCD Gain
###                    "{0:5.2f}".format(0.25),                # CCD Noise 
###                    "{0:6.2f}".format(data.simLibSkySig),   # SKYSIG
###                    "{0:4.2f}".format(data.simLibPsf),      # PSF1 
###                    "{0:4.2f}".format(0.),                  # PSF2 
###                    "{0:4.3f}".format(0.),                  # PSFRatio 
###                    "{0:6.2f}".format(data.simLibZPTAVG),   # ZPTAVG
###                    "{0:6.3f}".format(0.005),               # ZPTNoise 
###                    "{0:+7.3f}".format(-99.)]               # MAG
###             s = sep.join(lst)
###             y += s + '\n'
###         return y
###     
###     def writeSimLibField(self, fieldID):
###         s = self.fieldheader(fieldID)
###         s += self.formatSimLibField(fieldID, sep=' ')
###         s += self.footer(fieldID)
###         return s
### 
###     def simLibheader(self): #, user=None, host=None, survey='LSST', telescope='LSST'):
###         # comment: I would like to generalize ugrizY to a sort but am not sure
###         # of the logic for other filter names. so ducking for now
###         s = 'SURVEY: {0:}    FILTERS: ugrizY  TELESCOPE: {1:}\n'.format(self.survey, self.telescope)
###         s += 'USER: {0:}     HOST: {1:}\n'.format(self.user, self.host) 
###         s += 'BEGIN LIBGEN\n'
###         return s
###     
###     def simLibFooter(self):
###         """
###         """
###         s = 'END_OF_SIMLIB:    {0:10d} ENTRIES'.format(len(self.fieldIds))
###         return s
### 
### 
###     def writeSimlib(self, filename, comments=''):
###         with open(filename, 'w') as fh:
###             simlib_header = self.simLibheader()
###             simlib_footer = self.simLibFooter()
###             fh.write(simlib_header)
###             fh.write(comments)
###             # fh.write('BEGIN LIBGEN\n')
###             # fh.write('\n')
###             for fieldID in self.fieldIds:
###                 fh.write(self.fieldheader(fieldID))
###                 fh.write(self.formatSimLibField(fieldID))
###                 fh.write(self.fieldfooter(fieldID))
###             fh.write(simlib_footer)
### 
### 
