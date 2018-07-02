Operation Simulator Outputs
===========================

`OpSimSummary` can be used with Operation Simulator Outputs for OpSim version 3, and 4 as well as a set of OpSim outputs released for the work of the DESC Survey Strategy Task Force listed here_ . The outputs of these versions are different, and as of now, `OpSimSummary` requires an input to know which of these versions are being used.

.. _here: http://altsched.rothchild.me:8080/


The main class of `OpSimSummary` are required to read in `OpSim` outputs is `OpSimOutputs`. Note that the class `SynOpSim` also uses `OpSimOutputs` to
 read in `OpSim` databases. `OpSimOutput` can be instantiated using a `summary` table and `proposal` table read in from any source. The `summary` table has a set of unique pointings of the telescope with observational characteristics of each of the pointings (such as found in the OpSim `summary` or `summaryAllProps` tables). The `Proposal` Table has a list of different `proposals` with an identifying integer index and descriptive names.


.. code-block:: python

    opsimout = OpSimOutput(summary=summary, proposalTable=proposalTable,
                           zeroDDFDithers=False)

It is, however, far easier and strongly recommended to use this directly using a class method and an absolute path to the `OpSim` database `dbname`. The following is the example code for `lsstv3` (for example `minion_1016_sqlite.db`. The choices for `opsimversion` and `tableNames` for the names of the `summary` and `proposal` table in OpSim v3 are important: 

.. code-block:: python

    $opsout = OpSimOutput.fromOpSimDB(dbname, subset='combined',
                                      tableNames=('Summary', 'Proposal'),
                                      propIDs=None, zeroDDFDithers=True,
                                      opsimversion='lsstv3')

For the set of outputs for `OpSim` versions used by the DESC Survey Strategy Task Force: these variables are :

.. code-block:: python

    $opsout = OpSimOutput.fromOpSimDB(dbname, subset='combined',
                                      tableNames=('summaryAllProps', 'Proposal'),
                                      propIDs=None, zeroDDFDithers=True,
                                      opsimversion='sstf')
