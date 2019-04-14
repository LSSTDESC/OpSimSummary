Overview
========


``OpSimSummary`` was developed to interact with outputs of the `LSST Operations Simulator`_  for the purposes of transient and variable (mostly supernova cosmology) studies. Examples of such outputs are downloadable from `the opsim webpage`_ formatted as sqlite databases with fixed schema_. 

.. _`the opsim webpage`: https://www.lsst.org/scientists/simulations/opsim/opsim-v335-benchmark-surveys
.. _schema: https://www.lsst.org/scientists/simulations/opsim/summary-table-column-descriptions-v335
.. _`LSST Operations Simulator`: https://www.lsst.org/scientists/simulations/opsim


The key functionality of ``OpSimSummary`` is to find the set of LSST visits in an ``OpSim`` database that will observe a given location of the sky. This set of visits and the properties of these visits (also obtained from ``OpSim`` databases) is a useful piece of information for the simulations of observations of Time Domain Astronomy Sources (TDAS). The functionality is available in two ways:

1. As a library: ``OpSimSummary`` can used as a library by a simulation program. In this case, the simulator can use ``OpSimSummary`` to obtain the set of visits and their properties in the form of a `pandas.DataFrame` using a generator.  The basic code snippet to be used is demonstrated in the code block below, and a usable example is in the this demo notebook_

.. _notebook: https://github.com/rbiswas4/OpSimSummary/blob/master/example/Demo_SynOpSim.ipynb

.. code-block:: console

    from opsimsummary import SynOpSim
    
    synopsim = SynOpSim.fromOpSimDB(myopsimv3,
                                    opsimversion='lsstv3',
                                    usePointingTree=True)
    
    # If we want to obtain visits which observe two transients : ta at  (54., -27.5)
    # and tb at (52., -23.0)
    gen = synopsim.pointingsEnclosing(ra=(54.0, 52.),
                                      dec=(-27.5, -23.0),
                                      circRadius=0.,
                                      pointingRadius=1.75,
                                      usePointingTree=True)
    
    df_ta = next(gen)
    df_tb = next(gen)
    # df_ta is a `pandas.dataFrame` with all of the opsim columns for the visits
    # observing transient t_a (assumed to be fixed at (54., -27.5) with the
    # observationId or obsHistId of the visit as the index of the data frame

2. Pre-computed files: Alternatively, such sets and visit properties (or derived products) may be obtained written out to a text file for a set of discrete (but possibly large number of locations). Currently, we do this to provide an observation library (formerly called ``simlib``) for  SNANA_  to enable SNANA_ simulations of LSST. To use this, one could (after installing ``OpSimSummary``) use a `script to make the simlib files`_. The options for the script can be found by using the following code block.


.. _SNANA: http://snana.uchicago.edu 
.. _`script to make the simlib files`: https://github.com/rbiswas4/OpSimSummary/blob/master/scripts/make_simlibs.py

.. code-block:: console

    $python scripts/make_simlibs.py -h 



.. _codestructure: ./codestructure.rst

A further description of the main building blocks of this code are found in the `codestructure`_.


