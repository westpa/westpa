#!/bin/bash

PREFIX=$1

if [ -d $PREFIX ]; then
   echo "$PREFIX already exists."
   echo "Please specify a different directory or remove the current one and try again."
   return
else
   mkdir -p $(dirname $PREFIX)
fi

# Specify conda environment name
conda_env=westpa-2017.10

# Determine whether Linux or Mac
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
   return
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

# Install WESTPA in virtual environment
conda create --yes -n $conda_env westpa

# Place WESTPA environment variables inside conda env
mkdir -p $PREFIX/envs/$conda_env/etc/conda/activate.d
mkdir -p $PREFIX/envs/$conda_env/etc/conda/deactivate.d
touch $PREFIX/envs/$conda_env/etc/conda/activate.d/env_vars.sh
cat <<EOC >> $PREFIX/envs/$conda_env/etc/conda/activate.d/env_vars.sh
. $(dirname $(dirname `which python2.7`))/envs/$conda_env/$conda_env/westpa.sh
EOC

touch $PREFIX/envs/$conda_env/etc/conda/deactivate.d/env_vars.sh

# Reminder to add Miniconda2 install location to PATH
echo "Be sure to add the following to your .bashrc startup file"
echo "export PATH="$PREFIX/bin:'$PATH'""
