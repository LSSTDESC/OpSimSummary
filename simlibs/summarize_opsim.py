#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
from sqlalchemy import create_engine


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

class SummaryOpsimG(object):
    
    
    def __init__(self, summarydf=None, simlibFile=None, user=None, host=None, survey='LSST',
                 telescope='LSST', pixSize=0.2):

        import os 
        import subprocess

        if summarydf is None and simlibFile is None:
            raise ValueError('summarydf or simlibFile have to be supplied\n')
        
        if not (summarydf is None or simlibFile is None):
            raise ValueError('Both summarydf and simlibFile cannot be supplied\n')

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


