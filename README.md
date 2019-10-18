# OpSimSummary

[![Build Status](https://travis-ci.org/LSSTDESC/OpSimSummary.svg?branch=master)](https://travis-ci.org/LSSTDESC/OpSimSummary)
[![DOI](https://zenodo.org/badge/37937479.svg)](https://zenodo.org/badge/latestdoi/37937479)

`OpSimSummary` is a codebase developed to interact with the LSST Operations Simulator outputs. Currently they are used for catalog Time Domain Simulations. 
This includes a library that can be called by simulation codes to obtain the set of LSST pointings observing a particular point, as well as a script which uses
the library and precomputes such pointings and store them in an observation library. This storage is in a format specific to [`SNANA`](http://snana.uchicago.edu/)

# Using `OpSimSummary`
`OpSimSumamry` is open source and licensed under [BSD 3-clause](./LICENSE). While a release is available at [zenodo](https://zenodo.org/record/2671955#.XNPhvi2ZM1g) the software will continue to be developed as needed in this github repository. If you plan to use a released version on zenodo, you can cite the specific release (for example the current release using the doi (for example see the `Export` section of the [zenodo](https://zenodo.org/record/2671955#.XNPhvi2ZM1g) page for the release.) Alternatively, you can cite the a specific non-release version on the github repository through the git SHA.

Additionally, in both cases, please cite the [code paper](https://arxiv.org/abs/1905.02887). If you are using bibtex, [NASA ADS](http://adsabs.harvard.edu) provides a [bibtex](https://ui.adsabs.harvard.edu/abs/2019arXiv190502887B/exportcitation)  for this reference (as also other formatted forms).

- Documentation is hosted at [github_pages](https://lsstdesc.github.io/OpSimSummary/build/html/index.html)

## Installation  and Software Requirements

`opsimsummary` runs on either python 2.7 or 3.6+ . The list of required software to run `opsimsummary` is listed [here](./install/requirements.md). For installation methods, please see the [documentation](https://lsstdesc.github.io/OpSimSummary/build/html/installation.html)

