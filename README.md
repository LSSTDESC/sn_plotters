# sn_plotters
A set of python scripts to display results of the Survey Strategy Support pipeline (sn_pipe)

```
This software was developed within the LSST DESC using LSST DESC resources, and so meets the criteria 
given in, and is bound by, the LSST DESC Publication Policy for being a "DESC product".
We welcome requests to access code for non-DESC use; if you wish to use the code outside DESC please contact the developers.

```
## Release Status

|Release|Date|packages|
|---|---|---|
|v1.0.0|2020/07/15|sn_plotter_fitlc_v1.0.0,sn_plotter_metrics_v1.0.0,sn_plotter_simu_v1.0.0|

## Feedback, License etc

If you have comments, suggestions or questions, please [write us an issue](https://github.com/LSSTDESC/sn_plotters/issues).

This is open source software, available for re-use under the modified BSD license.

```
Copyright (c) 2020, the sn_plotters contributors on GitHub, https://github.com/LSSTDESC/sn_plotters/graphs/contributors.
All rights reserved.
```

## Content of sn_plotters ##

 * **docs**
 * **\_\_init\_\_.py**
 * **version.py**: package version
 * **LICENCE**
 * **setup.py**
 * [**sn_plotter_fitlc**](doc_package/plot_fitlc.md) : set of scripts to plot results from fit
 * [**sn_plotter_simu**](doc_package/plot_simu.md) : set of scripts to plot results from simulation
 * [**sn_plotter_metrics**](doc_package/plot_metrics.md)  : set of scripts to plot results from metrics
 * **README.md** : this readme

## Complete tree ##
```bash
|-- docs
|   |-- api
|   |   |-- sn_plotter_fitlc.fitlcPlot.rst
|   |   |-- sn_plotter_fitlc.rst
|   |   |-- sn_plotter_metrics.cadencePlot.rst
|   |   |-- sn_plotter_metrics.globalPlot.rst
|   |   |-- sn_plotter_metrics.nsnPlot.rst
|   |   |-- sn_plotter_metrics.rst
|   |   |-- sn_plotter_metrics.snrPlot.rst
|   |   |-- sn_plotter_simu.rst
|   |   |-- sn_plotter_simu.simuPlot.rst
|   |-- conf.py
|   |-- index.rst
|   |-- make.bat
|   |-- Makefile
|-- __init__.py
|-- version.py
|-- LICENCE
|-- README.md
|-- setup.py
|-- sn_plotter_fitlc
|   |-- fitlcPlot.py
|   |-- __init__.py
|   |-- version.py
|-- sn_plotter_metrics
|   |-- cadencePlot.py
|   |-- globalPlot.py
|   |-- __init__.py
|   |-- nsnPlot.py
|   |-- snrPlot.py
|   |-- version.py
|-- sn_plotter_simu
    |-- __init__.py
    |-- simuPlot.py
    |-- version.py
```