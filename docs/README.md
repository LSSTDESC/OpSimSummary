# Procedure to generate the documentation.

- clean the source directory of `opsimsummary.*.rst` files. `opsimsummary.rst` and `modules.rst`. 
- generate these files with ```sphinx-apidoc ../opsimsummary/ -o source/ --separate``` from the `docs` folder.
- remove the file documenting version.py
- manually add the line ```.. include:: codestructure.rst``` to `source/module.rst` file before the table.
- run the Makefile

All the steps can be run using the `makdocs.sh`. If you are on a mac, please check that your `sed` is linked to the `gsed`
