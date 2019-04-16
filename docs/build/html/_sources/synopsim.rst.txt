SynOpSim
========

This is the key object in `OpSimSummary`. This has a memory footprint of the size of the opsim database,
and contains all the information about the database. It can be constructed in the following way:


.. autoclass:: SynOpSim
    :members


.. code-block:: console

    from opsimsummary import SynOpSim
    
    synopsim = SynOpSim.fromOpSimDB(myopsimv3,
                                    opsimversion='lsstv3',
                                    usePointingTree=True)
    

.. code-block:: console
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


