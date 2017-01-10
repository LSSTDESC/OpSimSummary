# Generating simlib Files

Tools for studying the LSST OpSim outputs largely based on `pandas`. This can be used for example to generates
 simlib files for SNANA simulations or similar products that are useful for generating SN simulations.

[![Build Status](https://travis-ci.org/rbiswas4/OpSimSummary.svg?branch=master)](https://travis-ci.org/rbiswas4/OpSimSummary)
# Requirements

- pandas 
- numpy
- sqlalchemy
- ipython notebook (only to run example files)
- healpy (also installed by lsst.sims.maf)
- lsst.sims.maf (The example can be run, if you already know the propID corresponding to the proposal you want) 
- basemap: can be installed from conda https://anaconda.org/conda-forge/basemap

Everything other than lsst.sims.maf is part of the anaconda python distribution and available on pip. For lsst.sims.maf installation instructions, see the confluence page. 
# Running
- Install the package by 

```
python setup.py install --user
```


(The --user argument only installs this for the user only, leaving that out will install for all users on the computer if you have appropriate permissions) 

- open the examples/ExampleSimlib.ipynb notebook and execute cells with shift+enter
