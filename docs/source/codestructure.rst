Package Structure
================

The key components are in three modules. 

1. `opsim_out.py`: This module contains a class `OpSimOutput`. Its function is to read in the OpSim database, irrespective of the OpSim version, and create the same object for the rest of the code.


2. `summarize_opsim.py`: The module contains two classes of interest: `SynOpSim` which is central object in `OpSimSummary`, and `PointingTree`. `SynOpSim` can communicate directly with `OpSim` outputs using the `OpSimOutput` class mentioned above. It also uses `PointingTree` to calculate the `OpSim` pointings that will observe a certain transient.

3. `simlib.py`:
- `summarize_opsim`: The major operations of interest happen here. The user facing object is the class `SynOpSim` which can read an OpSim database directly through the convenience method `fromOpSimDB` which uses the code in `OpSimOutput` for this purpose.







The first module opsim_out.py_  deals with inputs, being designed to ingest two pandas dataframes (tables) including the `summaryAllProps` table with unique observations, the proposal table, and the \verb OpSim  version as instantiation parameters of the class \verb OpSimOutput . It is expected that users will want to read in \verb OpSim  databases, and hence this has a method to instantiate directly from an OpSim output database through \verb OpSimOutput.fromOpSimDB  method. The use of this class is also meant to guard against changes in the `OpSim` output in between versions, providing a stable interface for transient/variable simulation software across multiple versions simultaneously. This is done through a dictionary that keeps track of important conceptual variables as keys and their values across versions, remapping the values (names of quantites and also units) to a standard choice to fulfill our objective of version stability. The functionalities required for obtaining the re-organziation and re-summarizing of the database content are in the module \verb summarize_opsim . The principal object here is the class \verb SynOpSim , which can be instantiated using two pandas dataframes representing the `summaryAllProps' table and the `proposal' table, or directly from the outputs of LSST opsim (in both version 3 and 4) through methods that essentially call \verb OpSimOutput.fromOpSimDB . Reading in the OpSim database into pandas dataframes uses a memory of $\sim 1\rm{GB}$ wich is usually acceptable as a memory load on modern computers. At this point, \verb SynOpSim  builds a tree of the centers of the focal plane in each visit using methods in \verb scikit-learn (The \verb BallTree  class to be particular). Given a sequence of ra, and dec values in degrees describing the location of astrophysical sources, the \verb SynOpSim  method \verb pointingsEnclosing  takes a sequence of source location as input, and returns a sequence of quantities, each quantity being for a source. Each quantity (ie. corresponing to a source) includes all the visits observing the location of that source. While the size of this object may be large, a source may be observed about a thousand times in the wide fast deep survey (and $\sim 20, 000 $ times in the DDF) and the length of transients can be extremely long, the use of a generator to work on this one transient at a time helps in managing memory usage.

In order to create simlib files \verb SNANA  using this library, it is recommended that users use the simlib creation script provided in \verb scripts/make_simlibs.py , or write a simliar script. The arguments to the script are best ascertained by running the examples in Listing.~\ref{code:runsimlib} and further details of what this script does, and therefore of the data products created using this script are in the next subsection.
The library module used for the production of simlibs is in the module \verb simlibs  , after using the functionality of `SynOpSim` above is `simlibs.py`. This provides a \verb Simlib  mixin for use with \verb SynOpSim . The quantities required in a simlib file are transformations of the quantities available in the \verb OpSim  database. These conversions are performed using the method \verb SimlibMixin.add_simlibCols  . 
