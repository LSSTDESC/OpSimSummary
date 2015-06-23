# Gnerating simlib Files

Generates simlib files for SNANA simulations.

# Requirements

- pandas 
- mumpy
- sqlalchemy
- ipython notebook (only to run example files)
- lsst.sims.maf (The example can be run, if you already know the propID corresponding to the proposal you want) 

Everything other than lsst.sims.maf is part of the anaconda python distribution and available on pip. For lsst.sims.maf installation instructions, see the confluence page. 
# Running
- Install the package by 

```
python setup.py install --user
```


(The --user argument only installs this for the user only, leaving that out will install for all users on the computer if you have appropriate permissions) 

- open the examples/ExampleSimlib.ipynb notebook and execute cells with shift+enter
