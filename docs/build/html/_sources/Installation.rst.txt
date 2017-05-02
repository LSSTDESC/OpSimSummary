Installation
============

The engine of `OpSimSummary` depends on commonly used python packages has both essential and and optional requirements. 

Currently, all of the requirements can be installed with `anaconda`

Install Dependencies
====================
There are scripts in the `./install` directory which can be used to install the dependencies:

From the root level of the package:

.. code-block:: console

   $./install/install_pip_requirements.sh
   $./install_conda_requirements.sh

installs the dependencies aside from `lsst.sims.maf`. Currently, this is not an essential requirement, but his might become the case. 

To install `lsst.sims.maf` , please look at the `instructions here <https://confluence.lsstcorp.org/display/SIM/Catalogs+and+MAF>`_ .

Finally, to install `OpSimSummary` please clone the the `repository<https://github.com/rbiswas4/OpSimSummary>`_ and run the setup:

.. code-block:: console

    $git clone git@github.com:rbiswas4/OpSimSummary.git
    $cd OpSimSummary
    $python setup.py install --user


