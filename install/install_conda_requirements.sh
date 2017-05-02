conda config --add channels anaconda

conda install -q --yes --file ./install/conda_requirements.txt
conda list --explicit > spec-file.txtconda list --explicit > ./install/spec-file.txt;
