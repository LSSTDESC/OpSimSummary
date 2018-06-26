Installation
============

The engine of `OpSimSummary` depends on commonly used python packages. Currently, all of the
requirements can be installed with `anaconda` and `pip`. To use `OpSimSummary` as a library or
the scripts in `OpSimSummary`, you will need the requirements listed [here](https://github.com/rbiswas4/OpSimSummary/blob/fix_release_help/install/requirements.md)

Installing Dependencies
========================
There are scripts in the `./install` directory which can be used to install the dependencies. To use these, you must have a conda python
installation. If you don't have a such a python installation, you can use [a script](https://github.com/rbiswas4/install_utils/tree/master/scripts) to install a
miniconda distribution. 

To install the dependenencies and `OpSimSummary`, from the root level of the package:

.. code-block:: console

   $bash ./install/install_all.sh

installs the dependencies and `OpSimSummary`.

To just install `OpSimSummary` if you already have taken care of the requirements, use :

.. code-block:: console 

    $bash install/install_opsimsummary

To install other dependencies that are useful (for example for notebooks, and utlities for plotting) use:

.. code-block:: console

    $conda install jupyter seaborn


