Installation
============

The engine of ``OpSimSummary`` depends on commonly used python packages. Currently, all of the
requirements can be installed with `anaconda` and `pip`. To use ``OpSimSummary`` as a library or
the scripts in ``OpSimSummary``, you will need the requirements listed here_ .

.. _here: https://github.com/LSSTDESC/OpSimSummary/blob/master/install/requirements.md

Installing Dependencies
========================
There are scripts in the ``./install`` directory which can be used to install the dependencies. To use these, you must have a conda python
installation. If you don't have a such a python installation, you can use a script_ to install a miniconda distribution. 

.. _script: https://github.com/rbiswas4/install_utils/blob/master/scripts/install_python.sh
using the code
```
bash install_python.sh path_to_desired_python_location
```
where `path_to_desired_python_location` must be a directory that exists.

To install the dependenencies and ``OpSimSummary``, from the root level of the package:

.. code-block:: console

   $bash ./install/install_all.sh

installs the dependencies and ``OpSimSummary``.

To just install ``OpSimSummary`` if you already have taken care of the requirements, use :

.. code-block:: console 

    $bash install/install_opsimsummary

To install other dependencies that are useful (for example for notebooks, and utlities for plotting) use:

.. code-block:: console

    $conda install jupyter seaborn


