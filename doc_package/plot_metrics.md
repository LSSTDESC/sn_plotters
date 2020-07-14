# sn_plotter_metrics

## cadencePlot.py ##
|name | type | task|
|----|----|----|
|Lims | class | to handle light curve of SN|
|plotCadence| function|display two plots: cadence plot and histogram of redshift limits|
|plotMollview|function|Mollweid or Cart view(multiple)|
|plotViewIndiv|function|Mollweid or Cart view (single)|
|plotDDLoop|function|multiple plots for DDF|
|plotDDCadence_barh|function|multiple barh-plots for DDF|
|plotDDCorrel|function|plot display for DD correlation|
|plotDDFit|function|Function to perform some linear fit for two variables|
|plotDDLoopCorrel|function|multiple plots for DDF correlation study|
|plotDDCadence|function|DDF cadence plot|


## snrPlot.py ##
|name | type | task|
|----|----|----|
|SNRPlot|function|Signal-to-Ratio vs MJD plot for one field and one season |
| | |SNR of  a SN with T0=MJD-10 days|
|detecFracPlot|function|Plot Mollweid view of detection rates|
|detecFracHist|function|Plot histogram of detection rates|
|detecFracHist_bandseason|function|Plot histogram of detection rates per band and per season|
|GetHealpix|function|Get Healpix map of data|

## globalPlot.py ##
|name | type | task|
|----|----|----|
|PlotHist|class|plot histograms for a set of OS|
|PlotTime|class|plot set of variables for a given OS|
|PlotStat|class|display median values of the Global Metric|

## nsnPlot.py ##
|name | type | task|
|----|----|----|
|NSNAnalysis |class|analyze results from NSN metric|
|plot_DDSummary|function|to display NSN results for DD fields|
|plotNSN|function| Plot NSN vs redshift limit|