import pandas as pd
from dataclasses import dataclass
import glob
import numpy as np
from abc import ABC, abstractmethod
import healpy as hp
from sn_tools.sn_io import loopStack
from sn_tools.sn_utils import multiproc


@dataclass
class Simu:
    type: str
    num: str
    dir: str
    list: str
    nside: int


class Infos:
    """
    class to build a dataframe
    with requested infos to make plots

    Parameters
    ---------------
    simu: dataclass of type Simu

    """

    def __init__(self, simu, ip):

        self.simu = simu
        self.ip = ip
        self.families = []
        self.colors = ['b', 'k', 'r', 'g', 'm', 'c']
        self.markers = [".", "o", "v", "^", "<", ">", "1", "2", "3", "4", "8",
                        "s", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d", "|", "_"]
        self.markers += ["x", "X", "D", "d", "|", "_"]

        dbList = pd.read_csv(simu.list, comment='#')
        print(dbList)

        # self.rlist = self.cleandbName(dbList)
        self.rlist = dbList['dbName'].to_list()
        self.resdf = self.dfInfos()

    def cleandbName(self, dbList):
        """
        Method to clean dbNames by removing all char after version number (included): v_..

        Parameters
        ---------------
        dbList: pandas df
          containing the list of dbNames

        Returns
        ----------
        list of 'cleaned' dbNames

        """
        r = []
        # get 'cleaned' dbName
        for i, val in dbList.iterrows():
            spl = val['dbName'].split('v{}'.format(self.simu.num))[0]
            if spl[-1] == '_':
                spl = spl[:-1]
            r.append(spl)

        return r

    def clean(self, fam):
        """
        Method to clean a string

        Parameters
        ----------------
        fam: str
          the string to clean

        Returns
        -----------
        the final string

        """
        if fam[-1] == '_':
            return fam[:-1]

        if fam[-1] == '.':
            return fam[:-2]

        return fam

    def family(self, dbName):
        """
        Method to get a family from a dbName

        Parameters
        --------------
        dbName: str
           the dbName to process

        Returns
        ----------
        str: the name of the 'family'


        """

        ro = []
        fam = dbName
        for i in range(len(dbName)):
            stre = dbName[:i+1]
            num = 0
            for kk in self.rlist:
                if stre == kk[:i+1]:
                    num += 1
            # print(stre, num)
            ro.append(num)
            if i > 5 and ro[-1]-ro[-2] < 0:
                fam = dbName[:i]
                break

        return self.clean(fam)

    def dfInfos(self):
        """
        Method to build a pandas df
        with requested infos for the plotter
        """

        resdf = pd.DataFrame()
        # get families and fill infos
        for va in self.rlist:
            resdf = pd.concat((resdf, self.getInfos(va)))

        return resdf

    def getInfos(self, dbName):
        """
        Method to build a df with infos for plotter
        for a single dbName

        Parameters
        ---------------
        dbName: str
          dbName to process

        Returns
        -----------
        pandas df with infos as cols.


        """
        fam = self.family(dbName)

        if fam not in self.families:
            self.families.append(fam)
        #print('infos', fam, len(self.markers), self.families.index(fam))
        imark = self.families.index(fam)
        print(self.simu.type, self.simu.dir, dbName, fam,
              self.colors[self.ip], self.markers[self.families.index(fam)])

        return pd.DataFrame({'simuType': [self.simu.type],
                             'simuNum': [self.simu.num],
                             'dirFile': [self.simu.dir],
                             'dbName': [dbName],
                             'family': [fam],
                             'color': [self.colors[self.ip]],
                             'marker': [self.markers[self.families.index(fam)]]})


class ProcessData:

    def __init__(self, nside, metricName, fieldType):

        self.nside = nside
        self.metricName = metricName
        self.fieldType = fieldType

    def processMulti(self, toproc, outFile, process_class, nproc=1):
        """
        Method to analyze metric output using multiprocesses
        The results are stored in outFile (npy file)

        Parameters
        --------------
        toproc: pandas df
         data to process
        outFile: str
         output file name
        nproc: int, opt
         number of cores to use for the processing

        """
        params = {}
        params['fieldType'] = self.fieldType
        params['metricName'] = self.metricName
        params['nside'] = self.nside
        params['process_class'] = process_class

        print('multiprocessing', nproc)
        resdf = multiproc(toproc, params, self.processLoop, nproc)
        np.save(outFile, resdf.to_records(index=False))

    def processLoop(self, toproc, params, j=0, output_q=None):
        """
        Function to analyze a set of metric result files

        Parameters
        --------------
        toproc: pandas df
          data to process
        j: int, opt
          internal int for the multiprocessing
        output_q: multiprocessing.queue
         queue for multiprocessing

        Returns
        -----------
        pandas df with the following cols:
        zlim, nsn, sig_nsn, nsn_extra, dbName, plotName, color,marker
        """

        fieldType = params['fieldType']
        metricName = params['metricName']
        nside = params['nside']
        process_class = params['process_class']
        npixels = -1

        # this is to get summary values here
        resdf = pd.DataFrame()
        for index, val in toproc.iterrows():
            metricdata = process_class(
                val, self.metricName, self.fieldType, self.nside, npixels)

            # metricdata.plot()
            # plt.show()
            if metricdata.data_summary is not None:
                resdf = pd.concat((resdf, metricdata.data_summary))

        print('end of proc', j)
        if output_q is not None:
            output_q.put({j: resdf})
        else:
            return resdf


class ProcessFile(ABC):

    def __init__(self, info, metricName, fieldType, nside, npixels):
        """
        class to analyze results from NSN metric

        Parameters
        ---------------
        info: array
          various infos (dirfile, dbname, ...)
        metricName: str
          metric name
        fieldType: str
          type of field to process
        nside: int
           healpix nside parameter
        npixels: int
          total number of pixels processed

        """
        self.info = info
        self.metricName = metricName
        self.fieldType = fieldType
        self.nside = nside
        self.npixels = npixels

        self.processFiles()

    def processFiles(self):
        """
        Method to process a set of files


        """

        search_path = '{}/{}/{}/*{}Metric_{}*_nside_{}_*.hdf5'.format(
            self.info['dirFile'], self.info['dbName'], self.metricName, self.metricName, self.fieldType, self.nside)
        print('looking for', search_path)

        fileNames = glob.glob(search_path)

        print(fileNames)
        if len(fileNames) > 0:
            self.data_summary = self.process(fileNames)
        else:
            print('Missing files for', self.info['dbName'])
            self.data_summary = None

    @abstractmethod
    def process(self, fileNames):
        """
        Abstract method to process metric values from files

        Parameters
        ---------------
        fileNames: list(str)
          list of files to process


        """
        pass


class MetricValues:
    def __init__(self, dirFile, dbNames, metricName='NSN',
                 fieldType='DD', fieldNames=['COSMOS'], nside=128, nproc=8):
        """
        Class to transform input data and match to DD fieldsgrab metric values

        Parameters
        ---------------
        dirFile: str
          directory of the files to process
        dbNames: list(str) 
           list of OS to process
        metricName: str, opt
          name of the metric to consider (default: NSN)
        fieldType: str,opt
          field type to consider (default: DD)
        fieldNames: list(str), opt
          fieldNames to process (default: ['COSMOS'])
        nside: int, opt
          nside healpix parameter (default: 128)
        nproc: int, opt
          number of procs to use (default: 8)

        Returns
        ----------


        """

        # get pixelArea
        self.pixArea = hp.nside2pixarea(nside, degrees=True)
        x1_colors = [(-2.0, 0.2), (0.0, 0.0)]
        self.corr = dict(zip(x1_colors, ['faint', 'medium']))
        # self.data = self.process_loop(dirFile, metricName, fieldType, fieldNames,
        #                              nside, forPlot).to_records()

        params = {}
        params['dirFile'] = dirFile
        params['metricName'] = metricName
        params['fieldType'] = fieldType
        params['fieldNames'] = fieldNames
        params['nside'] = nside

        self.data = multiproc(dbNames, params, self.process, nproc)

    def process(self, dbNames, params, j=0, output_q=None):

        dirFile = params['dirFile']
        metricName = params['metricName']
        fieldType = params['fieldType']
        fieldNames = params['fieldNames']
        nside = params['nside']

        restot = pd.DataFrame()
        for dbName in dbNames:
            for fieldName in fieldNames:
                res = self.process_field(
                    dirFile, dbName, metricName, fieldType, fieldName, nside)
                restot = pd.concat((restot, res))

        if output_q is not None:
            return output_q.put({j: restot})
        else:
            return restot

    def process_field(self, dirFile, dbName, metricName, fieldType, fieldName, nside):
        """
        Single file processing
        This method load the files corresponding to dbName and transform it
        so as to have all infos on one line.


        Parameters
        ----------------
        dirFile: str
         directory where the files are located
        dbName: str
          name of the cadence to processe
        metricName: str
          name of the metric of interest
        fieldType: str
          field type: DD or WFD
        nside: int
          nside for healpix tessallation


        Returns
        -----------
         pandas df with metric values and additional info (dbName, fieldName, ...)

        """
        if fieldType == 'DD':
            search_path = '{}/{}/{}_{}/*{}Metric_{}*_nside_{}_*.hdf5'.format(
                dirFile, dbName, metricName,  fieldName, metricName, fieldType, nside)
        if fieldType == 'WFD':
            search_path = '{}/{}/{}/*{}Metric_{}*_nside_{}_*.hdf5'.format(
                dirFile, dbName,metricName, metricName, fieldType, nside)
        print('looking for', search_path)
        vars = ['pixRA', 'pixDec', 'healpixID', 'season', 'status']
        # vars = ['healpixID', 'season']
        fileNames = glob.glob(search_path)
        print(fileNames)
        finaldf = pd.DataFrame()
        if fileNames:
            # plt.plot(metricValues['pixRA'],metricValues['pixDec'],'ko')
            # plt.show()
            metricValues = loopStack(fileNames, 'astropyTable').to_pandas()
            metricValues = metricValues.round({'pixRA': 3, 'pixDec': 3})
            metricValues['dbName'] = dbName
            metricValues['fieldname'] = fieldName
            metricValues['pixArea'] = self.pixArea
            metricValues['filter'] = 'grizy'
            dbName_split = dbName.split('_')
            n = len(dbName_split)
            metricValues['dbName_plot'] = '_'.join(dbName_split[0:n-2])

            print(metricValues.columns)
            return metricValues


def get_dist(data, pixRA_mean=-1, pixDec_mean=-1):
    """
    Function to estimate the distance dist = sqrt((deltaRA*cos(Dec))**2+deltaDec**2)

    Parameters
    ---------------
    data: pandas df
      data to process

    Returns
    ----------
    pandas df with dist col

    """
    if pixRA_mean == -1:
        pixRA_mean = np.mean(data['pixRA'])
        pixDec_mean = np.mean(data['pixDec'])
    data['dist'] = np.sqrt(((data['pixRA']-pixRA_mean)*np.cos(np.deg2rad(data['pixDec'])))**2
                           + (data['pixDec']-pixDec_mean)**2)
    data['pixRA_mean'] = pixRA_mean
    data['pixDec_mean'] = pixDec_mean

    return data


def dumpcsv_medcad(metricTot,prefix='metric_summary_DD'):
    """
    Function to dump metric results in csv file

    Parameters
    --------------
    metricTot: pandas df
      data to process

    """

    data = pd.DataFrame(metricTot)

    summary = data.groupby(['dbName','family']).agg({'nsn': 'sum',
                                            'zcomp': 'median',
                                            }).reset_index()

    summary_fields = data.groupby(['dbName','family','fieldname']).agg({'nsn': 'sum',
                                                               'zcomp': 'median',
                                                                }).reset_index()
    summary_fields_season = data.groupby(['dbName', 'family','fieldname', 'season']).agg({'nsn': 'sum',
                                                                                'zcomp': 'median',
                                                                                 }).reset_index()
    if 'healpixID' in data.columns and 'gap_max' in data.columns:
        summary_fields_pixels = data.groupby(['dbName', 'family','fieldname', 'healpixID','pixRA','pixDec','season']).agg({'nsn': 'sum',
                                                                                        'zcomp': 'median',
                                                                                        'cadence': 'median',
                                                                                        'gap_max': 'median',
                                                                                        'season_length': 'median'
                                                                                        }).reset_index()
        
    if 'healpixID' in data.columns :
        summary_fields_pixels = data.groupby(['dbName', 'family','fieldname', 'healpixID','pixRA','pixDec','season']).agg({'nsn': 'sum',
                                                                                        'zcomp': 'median'
                                                                                        }).reset_index()    
    print(summary)
    print(summary_fields)
    summary.to_csv('{}.csv'.format(prefix), index=False)
    summary_fields.to_csv('{}_fields.csv'.format(prefix), index=False)
    summary_fields_season.to_csv('{}_fields_season.csv'.format(prefix), index=False)
    if 'healpixID' in data.columns:
        summary_fields_pixels.to_csv(
            '{}_fields_pixels.csv'.format(prefix), index=False)
