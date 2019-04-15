conda install -c anaconda -v --yes --file ./install/conda_requirements_anaconda.txt
echo "Done install conda requirements from anaconda"
conda install -c conda-forge -v --yes --file ./install/conda_requirements_conda_forge.txt
echo "Done install conda requirements from conda forge"
conda list --explicit > ./install/spec-file.txt;
