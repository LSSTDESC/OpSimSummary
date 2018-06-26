### This script will install miniconda packages to ${HOME}/astroml_miniconda3 by default or to a user specifed
### location specified through ```$0 path/mydir``` , where path must already exist.

## Directory to install miniconda to
miniconda_dir=${HOME}/miniconda_oss
if [ $# -gt 1 ];
    then echo "illegal number of parameters"
fi
if [ $# -eq 1 ]; then
    miniconda_dir=$@
    echo ${miniconda_dir}
    if [ -d ${miniconda_dir} ]; then
        echo "install directory ${miniconda_dir} already exists, please delete first to clobberr"
        exit
    fi
    dir=$(dirname ${miniconda_dir})
    if [ ! -d "${dir}" ]; then
        echo "path to miniconda install does not exist"
        exit
    fi
fi
echo " installing python to $miniconda_dir"

## Download miniconda 
# miniconda version
MINICONDA_VERSION="latest"
# platform linux/osx
system=`uname -s`
echo ${system}
if [ "$system" = "Linux" ]; then
   platform="Linux"
elif [ "$system" = "Darwin" ]; then
   platform="MacOSX" 
else
    echo "Platform unknown\n"
fi
echo "platform = ${platform}\n"
fname="Miniconda3-${MINICONDA_VERSION}-${platform}-x86_64.sh"
url="https://repo.continuum.io/miniconda/${fname}"
curl -OL ${url}

## Install miniconda
bash ${fname} -b -p ${miniconda_dir}
printf "You can clean up by typing \n rm -f $fname\n"
echo "In order to use this python, please set your path to include ${miniconda_dir}/bin, for example in bash:"
echo "export PATH=${miniconda_dir}"'/bin:$PATH'
echo "export PATH=${miniconda_dir}"'/bin:$PATH'>./install/setup.sh
export PATH=${miniconda_dir}/bin:$PATH
#conda env create -n astroml --file environment.yml
#echo "source activate astroml">>./install/setup.sh
echo "You can either run ./install/setup.sh everytime you want to use this, or copy the lines to your environment to use this python all the time"
