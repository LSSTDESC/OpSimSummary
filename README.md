# OpSimSummary

[![Build Status](https://travis-ci.org/LSSTDESC/OpSimSummary.svg?branch=master)](https://travis-ci.org/LSSTDESC/OpSimSummary)

`OpSimSummary` is a codebase developed to interact with the LSST Operations Simulator outputs. Currently they are used for catalog Time Domain Simulations. 
This includes a library that can be called by simulation codes to obtain the set of LSST pointings observing a particular point, as well as a script which uses
the library and precomputes such pointings and store them in an observation library. This storage is in a format specific to [`SNANA`](http://snana.uchicago.edu/)

- Documentation is hosted at [github_pages](https://lsstdesc.github.io/OpSimSummary/build/html/index.html)

## Installation  and Software Requirements

`opsimsummary` runs on either python 2.7 or 3.6+ . The list of required software to run `opsimsummary` is listed [here](./install/requirements.md). For installation methods, please see the [documentation](https://lsstdesc.github.io/OpSimSummary/build/html/Installation.html)

