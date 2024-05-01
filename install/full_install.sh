#!/bin/bash

# -----------------------------------------------------------------------------
# Determine OS Platform
OS="`uname`"
case $OS in
  'Linux')
    OS='Linux'
    ;;
  'Darwin')
    OS='MacOS'
    ;;
  *) ;;
esac


# -----------------------------------------------------------------------------
if command -v conda > /dev/null; then
  echo "Conda already installed."
else
  echo "Installing Miniconda"
  # Specify Miniconda version
  MINICONDA_VERSION="Miniconda3-latest"
  if [ "$OS" = "Linux" ]; then
    MINICONDA_SH="${MINICONDA_VERSION}-Linux-x86_64.sh"
  elif [ "$OS" = "MacOS" ]; then
    MINICONDA_SH="${MINICONDA_VERSION}-MacOSX-x86_64.sh"
  fi

  # Download Miniconda
  wget https://repo.anaconda.com/miniconda/$MINICONDA_SH

  # Install Miniconda
  bash $MINICONDA_SH -b -p $HOME/miniconda

  # Initialize Conda
  $HOME/miniconda/bin/conda init
  source $HOME/.bashrc

  echo "Miniconda installed."
  rm Miniconda*.sh
fi

# -----------------------------------------------------------------------------
# Activate Miniconda
# echo "Installing Mamba into base environment"
# source "$HOME/miniconda3/bin/activate"

# # Install Mamba from Conda-Forge
# conda install mamba -n base -c conda-forge
# echo "Mamba installed."

# -----------------------------------------------------------------------------
# Activate conda
echo "Installing the seidart environment and all dependencies."
conda init
source $HOME/.bashrc
conda activate
# source $HOME/miniconda3/bin/activate

# Create conda environment from YAML
conda env create -f seidart-environment.yml

echo 'Done! Type "conda activate seidart" to access the seidart environment.'

