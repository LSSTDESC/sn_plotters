from sn_tools.sn_io import Read_LightCurve
from sn_fitter.fit_sn_cosmo import Fit_LC
import sncosmo


class VisuLC:
    def __init__(self, metaFileInput, metaDirInput):

        meta = Read_LightCurve(file_name=metaFileInput, inputDir=metaDirInput)
        metaTable = meta.get_table(path='meta')

        metadata = metaTable.meta
        # get lc
        lcDir = metadata['directory']
        lcName = metadata['file_name']

        self.lcs = Read_LightCurve(file_name=lcName, inputDir=lcDir)

        # print SNIDS

        print(metaTable['path'])

    def plot(self, lcpath):

        lc = self.lcs.get_table(lcpath)

        # trying to fit here
        fit = Fit_LC(lc)
        result, fitted_model = fit('sncosmo', plot=True)
