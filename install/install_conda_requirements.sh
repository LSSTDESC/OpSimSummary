conda config --add channels anaconda
conda config --add channels openastronomy

conda install --no-update-dependencies --yes -v --file ./install/conda_requirements.txt
conda list --explicit > ./install/spec-file.txt;
