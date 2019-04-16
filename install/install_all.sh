#!/usr/bin/env bash
echo "installing dependencies from pip"
./install/install_pip_requirements.sh
echo "installing dependencies from conda"
conda update -n base -c defaults conda
conda clean --all
conda update --all
./install/install_conda_requirements.sh
echo "Done installing repositories"
echo "Install OpSimSummary"
./install/install_opsimsummary.sh

echo "run install test script"
# python ./install/check_install.py


