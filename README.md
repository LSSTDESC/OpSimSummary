# OpSimSummary

Tools for summarizing the LSST Operations Simulations outputs for use with catalog simulation tools or other analyses. 

[![Build Status](https://travis-ci.org/rbiswas4/OpSimSummary.svg?branch=master)](https://travis-ci.org/rbiswas4/OpSimSummary)
- documentation is hosted at [github_pages](https://rbiswas4.github.io/OpSimSummary)

# Installation 
`opsimsummary` runs on either python 2.7 or 3.5+ . The list of required software to run `opsimsummary` is listed [here](./install/requirements.md). Assuming miniconda is already installed, `opsimsummary` can be installedalong with its dependencies for the current user by 
```
./install/install_all.sh
```
and it will install dependencies from conda channels or use pip.
If the dependecies [listed](./install/requirements.md) are already installed, `opsimsummary` can be installed by running 
```
./install/install_opsimsummary.sh
```

Obviously, running the commands in the scripts will also allow the user to achieve the same results.

# Running
- open the `examples/Demo_OpSimOutput.ipynb` and execute the cells.
