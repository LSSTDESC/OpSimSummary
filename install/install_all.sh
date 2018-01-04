#!/usr/bin/env bash
echo "installing dependencies from pip"
./install/install_pip_requirements.sh
echo "installing dependencies from conda"
./install/install_conda_requirements.sh
echo "Done installing repositories"
echo "Install OpSimSummary"
./install/install_opsimsummary.sh


