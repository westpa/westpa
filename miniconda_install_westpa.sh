#!/bin/bash

if ! echo "$0" | grep '\.sh$' > /dev/null; then
   echo 'Please execute this script directly rather than sourcing the script.'
   echo "Usage: ./miniconda_install_westpa.sh <Miniconda Installation Prefix>"
   return 1
fi

if [[ "$#" -ne 1 || $1 == "-h" || $1 == "--help" ]]; then
   echo "Please specify the installation prefix for miniconda."
   echo -e "Usage: ./miniconda_install_westpa.sh <Miniconda Installation Prefix>"
   echo ""
   exit
fi

PREFIX="$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"

if [ -d $PREFIX ]; then
   echo "$PREFIX already exists."
   echo "Please specify a different directory or remove the current one and try again."
   exit
else
   mkdir -p $(dirname $PREFIX)
fi

### +------------------------------------------------+ ########################################
### | Step 1.  Install Miniconda                     | ########################################
### +------------------------------------------------+ ########################################

# Determine platform of machine (currently only for Linux or Mac) 
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     platform=Linux;;
    Darwin*)    platform=Mac;;
    *)          platform="UNKNOWN:${unameOut}"
esac
echo ${platform}

# Assign download target based on platform
if [ ${platform} == "Linux" ]; then
   download=Miniconda2-latest-Linux-x86_64.sh
elif [ ${platform} == "Mac" ]; then
   download=Miniconda2-latest-MacOSX-x86_64.sh
else
   echo "Installation script only supports Linux or Mac OS X."
   echo "Exiting installation."
   exit
fi

echo "download = " ${download}

# Download Anaconda Python 2.7 from https://www.continuum.io
if command -v wget >&/dev/null; then
   wget https://repo.continuum.io/miniconda/$download -P $(dirname $PREFIX)
else
    curl -o $(dirname $PREFIX)/$download https://repo.continuum.io/miniconda/$download
fi

# Batch install to specified PREFIX
bash $(dirname $PREFIX)/$download -b -p $PREFIX

# Prepend the Miniconda2 install location to PATH
export PATH="$PREFIX/bin:$PATH"

### +------------------------------------------------+ ########################################
### | Step 2.  Install WESTPA in virtual environment | ########################################
### +------------------------------------------------+ ########################################

# Specify conda environment name
conda_env=westpa-2019.8

# Install WESTPA in virtual environment
conda create --yes -c conda-forge -n $conda_env westpa

. $PREFIX/etc/profile.d/conda.sh

conda activate $conda_env

# Place WESTPA environment variables inside conda env
export ENV_PREFIX="$(dirname $(dirname `which python3`))"
mkdir -p $ENV_PREFIX/etc/conda/activate.d
mkdir -p $ENV_PREFIX/etc/conda/deactivate.d
touch $ENV_PREFIX/etc/conda/activate.d/env_vars.sh
cat << EOC >> $ENV_PREFIX/etc/conda/activate.d/env_vars.sh 
. $(dirname $(dirname `which python3`))/$conda_env/westpa.sh
EOC
touch $ENV_PREFIX/etc/conda/deactivate.d/env_vars.sh

conda deactivate

# Reminder to add Miniconda2 install location to PATH
echo "Be sure to add the following lines to your .bashrc startup file"
echo "export PATH="$PREFIX/bin:'$PATH'""
echo ". $PREFIX/etc/profile.d/conda.sh"

