Overview
========
`OpSimSummary` was developed to interact with outputs of the [LSST Operations Simulator](https://www.lsst.org/scientists/simulations/opsim) for the purposes of transient and variable (mostly
 supernova cosmology) studies. There are other tools that can interact with such outputs (downloadable from [here](https://www.lsst.org/scientists/simulations/opsim/summary-table-column-descriptions-v335) and currently formatted as sqlite databases with schema described [here](https://www.lsst.org/scientists/simulations/opsim/summary-table-column-descriptions-v335).

The main functionality of `OpSimSummary` is to find out the set of LSST visits in an `OpSim` database that will observe a non-moving object at a given location of the sky. This can be done in
two ways:

1. In memory : where the visit properties are passed as a `pandas.dataFrame` using a generator.  This is demonstrated in this [notebook](https://github.com/rbiswas4/OpSimSummary/blob/master/example/Demo_SynOpSim.ipynb)
2. Written out to a text file: This is to provide a `simlib` file to [SNANA](http://snana.uchicago.edu) to enable `SNANA` simulations of LSST. To use this, one could (after installing `OpSimSummary`) run to get the options for this [script](https://github.com/rbiswas4/OpSimSummary/blob/master/scripts/make_simlibs.py)
.. code-block:: console
    $python scripts/make_simlibs.py -h 
