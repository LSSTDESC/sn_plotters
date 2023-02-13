from sn_tools.sn_io import Read_LightCurve
from sn_fitter.fit_sn_cosmo import Fit_LC


class VisuLC:
    def __init__(self, metaFileInput, metaDirInput):
        """
        Class to visualize (and fit) LCs

        Parameters
        ----------
        metaFileInput : str
            metadata file.
        metaDirInput : str
            location dir of meta data file.

        Returns
        -------
        None.

        """

        meta = Read_LightCurve(file_name=metaFileInput, inputDir=metaDirInput)
        metaTable = meta.get_table(path='meta')

        metadata = metaTable.meta

        # get lc
        lcDir = metadata['lc_dir']
        lcName = metadata['lc_fileName']

        self.lcs = Read_LightCurve(file_name=lcName, inputDir=lcDir)

        # print SNIDS

        print(metaTable['SNID'])

        # fit instance
        self.fit = Fit_LC(model='salt3', version='2.0', display=True)

    def plot(self, lcpath):

        lc = self.lcs.get_table(lcpath)

        # trying to fit here
        outfit = self.fit(lc, plot=True)
        #result, fitted_model = fit('sncosmo', plot=True)
