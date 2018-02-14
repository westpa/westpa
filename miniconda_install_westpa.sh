#!/bin/bash

PREFIX=$1

if [ -d $PREFIX ]; then
   echo "$PREFIX already exists."
   echo "Please specify a different directory or remove the current one and try again."
   return
else
   mkdir $(dirname $PREFIX)
fi

# Download Anaconda Python 2.7 from https://www.continuum.io
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -P $(dirname $PREFIX)

# Batch install to specified PREFIX
bash $(dirname $PREFIX)/Miniconda2-latest-Linux-x86_64.sh -b -p $PREFIX

# Prepend the Miniconda2 install location to PATH
export PATH="$PREFIX/bin:$PATH"

# Add conda-forge channel
conda config --add channels conda-forge

# Install WESTPA in virtual environment
conda create --yes -n westpa-2017.10 westpa

# Reminder to add Miniconda2 install location to PATH
echo "Be sure to add the following to your .bashrc startup file"
echo "export PATH="$PREFIX/bin:'$PATH'""
