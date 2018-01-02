# OpSimSummary

Tools for summarizing the LSST Operations Simulations outputs for use with catalog simulation tools or other analyses. 

[![Build Status](https://travis-ci.org/rbiswas4/OpSimSummary.svg?branch=master)](https://travis-ci.org/rbiswas4/OpSimSummary)
- documentation is hosted at [github_pages](https://rbiswas4.github.io/OpSimSummary)
# Requirements
- python 2.7 or 3.5
- pandas 
- numpy
- sqlalchemy
- ipython notebook (only to run example files)
- healpy (also installed by lsst.sims.maf)
- lsst.sims.maf (The example can be run, if you already know the propID corresponding to the proposal you want) 
- basemap: can be installed from conda https://anaconda.org/conda-forge/basemap

Everything other than lsst.sims.maf is part of the anaconda python distribution and available on pip. For lsst.sims.maf installation instructions, see the confluence page. 
# Running
- Assuming you have miniconda installed, 
you can install the code and its dependencies for the current user by 
```
./install/install_all.sh
```
If you have the dependecies installed, you can install the code only by running 
```
.install/install_opsimsummary.sh
```


(The --user argument only installs this for the user only, leaving that out will install for all users on the computer if you have appropriate permissions) 

- open the examples/ExampleSimlib.ipynb notebook and execute cells with shift+enter
