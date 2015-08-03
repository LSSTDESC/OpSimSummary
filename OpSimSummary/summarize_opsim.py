#!/usr/bin/env python 
import os
import numpy as np
import pandas as pd
from sqlalchemy import create_engine
import matplotlib.pyplot as plt
# __all__ = ['SummaryOpsim']
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

    .. note :: This was written from a piece of f77 code by David
        Cinabro sent by email on May 26, 2015. 
    '''
    
    opsim_seeing = opsimtable['finSeeing'] # unit of arc sec sq
    # magsky is in units of mag/arcsec^2
    # opsim_maglim is in units of mag
    opsim_maglim = opsimtable['fiveSigmaDepth']
    opsim_magsky = opsimtable['filtSkyBrightness']

    # Calculate two variables that come up in consistent units:
    # term1  = 2.0 * opsim_maglim - opsim_magsky

    # Area of pixel in arcsec squared
    pixArea = pixSize * pixSize

    term1 = 2.0 * opsim_maglim - opsim_magsky # * pixArea
    # term2 = opsim_maglim - opsim_magsky
    term2 = - (opsim_maglim - opsim_magsky) # * pixArea


    # Calculate SIMLIB PSF VALUE
    opsimtable['simLibPsf'] = opsim_seeing /2.35 /pixSize
   
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
    
    #SKYSIG Calculation
    npix_asec = 1./ pixSize**2.
    opsimtable['simLibSkySig'] = np.sqrt((1.0/ npix_asec)*10.0 **(-0.4 * (opsim_magsky - simlib_zptavg)))
    return opsimtable

class SummaryOpsim(object):
    
    
    def __init__(self, summarydf, user=None, host=None, survey='LSST',
                 telescope='LSST', pixSize=0.2):
        '''
        Create a summary of the OpSim output. 

        Parameters
        ----------
        summarydf: mandatory, `~pandas.DataFrame` (alternative constructors too)
            DataFrame corresponding to the opsim output (with appropriate
            selections) to be summarized.
        user: string, optional, defaults to None
            user running the program, used in writing out SNANA simlibs only 
            if None, the login name of the user is used.
        host: string, optional, defaults to None
            name of host machine, used only in writing out SNANA simlibs
            default of None assigns the output of `hostname` to this variable.
        survey: string, optional, defaults to 'LSST'
            name of survey, required only for writing out SNANA simlibs
        telescope: string, optional, defaults to 'LSST'
            name of survey, required only for writing out SNANA simlibs
        pixSize: float, optional, defaults to 0.2
            size of the pixel in arcseconds, defaults to 0.2 as appropriate
            for LSST
        '''
        import os 
        import subprocess

        self.df = summarydf.copy(deep=True)
        self.df['MJDay'] = np.floor(self.df.expMJD.values)
        self.df['MJDay'] = self.df['MJDay'].astype(int)
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
    def fromOpSimASCII(cls, opSimFlatFile, **kwargs):
        """
        instantiates the sumamry table from opSim ASCII outputs like
        the outputs for version 2.168


        Parameters
        ---------
        opSimFlatFile : string, mandatory
            absolute path to ASCII file

        """
        import pandas as pd


        summary = pd.read_csv(opSimFlatFile, delimiter=r"\s+")
        return cls(summary, **kwargs)


    @classmethod
    def fromOpSimDB(cls, opSimDBfilename, tablename=None, sql_query=None,
                    dialect='sqlite', **kwargs):
        '''
        used to instantiate the summary from the opsim database


        Parameters
        ----------
        opSimDBfilename : str, mandatory
            absolute path to opsim sqlite database
        tablename : string, optional, defaults to None
            name of table in case it is not in sql_query
        sql_query : strng, mandatory
            sql_query to get the required OpSim visits
        dialect : string, defaults to 'sqlite'
            Only dialect supported at present

        Returns
        -------

        '''

        from sqlalchemy import create_engine
        import pandas as pd

        if dialect == 'sqlite':
            opSimDB = 'sqlite:///' + opSimDBfilename
        else:
            raise ValueError('other dialects not supported yet\n')
        engine = create_engine(opSimDB)
        if sql_query is None:
            summary = pd.read_sql_table(tablename, engine, **kwargs)
        else:
            summary = pd.read_sql_query(sql_query, engine, **kwargs)

        return cls(summary, **kwargs)
        
    def coords(self):

        ra = map(lambda x: self.ra(x), self.fieldIds)
        dec = map(lambda x: self.dec(x), self.fieldIds)

        return ra, dec


    def cadence_Matrix(self, fieldID, sql_query='night < 366',
                       mjd_center=None, mjd_range=[-50., 50.],
                       Filters=[u'u', u'g', u'r', u'i', u'z', u'Y'],
                       nightMin=0, nightMax=365, observedOnly=False):
    
        timeMin = nightMin
        timeMax = nightMax
        timeIndex = 'night'
        if mjd_center is not None:
            timeIndex = 'MJDay'
            if 'mjd' not in sql_query.lower():
                timeMin = mjd_center + mjd_range[0]
                timeMax = mjd_center + mjd_range[1]
                sql_query = 'MJDay > ' + str(timeMin) 
                sql_query += ' and MJDay < ' + str(timeMax)

        ss = pd.Series(np.arange(timeMin, timeMax))

        # group on filter and timeIndex (night)
        grouping_keys = ['filter', timeIndex]
        print timeIndex, sql_query
        print grouping_keys

        queriedOpsim = self.simlib(fieldID).query(sql_query)

        if queriedOpsim.size == 0 :
            Matrix = pd.DataFrame(index=ss, columns=Filters)
            Matrix.fillna(0., inplace=True)
            return Matrix

        
        grouped = queriedOpsim.groupby(grouping_keys)
 
        print grouped.groups.keys()
        # tuples of keys
        filts, times = zip( *grouped.groups.keys())
        print times

        # number of Observations in each group
        numObs = grouped.apply(len).values

        # Create a new dataFrame with nights, Filters, numObs as cols
        cadence_dict = dict()
        cadence_dict['Filters'] = list(filts)
        cadence_dict[timeIndex] = list(times)

        # If observedOnly: set values above 1 to 1 
        if observedOnly:
            numObs = np.where(np.array(list(numObs)), 1, 0)
        
        cadence_dict['numObs'] = list(numObs)

        # pivot dataFrame to occupation numbers
        X = pd.DataFrame(cadence_dict)
        Matrix = X.pivot(timeIndex, 'Filters', 'numObs')

        # First make sure all filters are represented
        for filt in Filters:
            if filt not in Matrix.columns:
                Matrix[filt] = np.nan

        # reorder filters to u,g,r,i,z,y
        M = Matrix[Filters]
        
        # Extend to all values in plot
        print ss.size, timeMin, timeMax
        Matrix = M.reindex(ss, fill_value=np.nan)

        return Matrix #, X


    @staticmethod
    def mjdvalfornight(night):
        return night + (49561 - 208)

    @staticmethod
    def nightformjd(mjd) :
        return mjd - (49561 - 208)

    def cadence_plot(self, fieldID, sql_query='night < 366', mjd_center=None,
                     mjd_range = [-50, 50],
                     Filters=[u'u', u'g', u'r', u'i', u'z', u'Y'],
                     nightMin=0, nightMax=365, deltaT=5., observedOnly=False,
                     title=True, title_text=None, colorbar=True,
                     colorbarMin=0., showmjd=True):
        """
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

        """


        Matrix = self.cadence_Matrix(fieldID, sql_query=sql_query,
                                mjd_center=mjd_center, mjd_range=mjd_range,
                                Filters=Filters, nightMin=nightMin,
                                nightMax=nightMax, observedOnly=observedOnly)

        if mjd_center is not None:
            timeMin = mjd_center + mjd_range[0]
            timeMax = mjd_center + mjd_range[1]
        else:
            timeMin = nightMin
            timeMax = nightMax

        if observedOnly:
            axesImage = plt.matshow(Matrix.transpose(), aspect='auto',
                                    cmap=plt.cm.gray_r, vmin=colorbarMin,
                                    vmax=1., extent=(timeMin - 0.5,
                                                     timeMax + 0.5,
                                                     -0.5, 5.5))
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
        # if mjd_center is not None:
        #    nightMin = mjd_center + mjd_range[0]
        #    nightMax = mjd_center + mjd_range[1]
        minorxticks = ax.set_xticks(np.arange(timeMin, timeMax,
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

        return fig, Matrix


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
        
        ra = np.degrees(self.ra(fieldID))
        dec = np.degrees(self.dec(fieldID))
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


