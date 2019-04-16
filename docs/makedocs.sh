rm source/{opsimsummary.rst,opsimsummary.*.rst,modules.rst}
sphinx-apidoc ../opsimsummary/ -o source/ --separate
rm source/opsimsummary.version.rst
sed -i '4i   .. include:: codestructure.rst' source/modules.rst
make clean; make html
