#!/bin/sh

# Default values for the variables
FULL=false
ENVUPGRADE=false

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --full) FULL=true ;;
        --envupgrade) ENVUPGRADE=true ;;
        *) echo "Unknown option: $1" ;;
    esac
    shift
done

# --------------------------------------------------------------------------
# Function to log messages with a timestamp
log_with_timestamp() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" > install_log.txt
}

# --------------------------------------------------------------------------
# Check if Anaconda/Miniconda are installed
# Determine OS Platform
OS="$(uname)"
case $OS in
  'Linux') OS='Linux' ;;
  'Darwin') OS='MacOS' ;;
  *) ;;
esac

# Miniconda installation and environment creation
if [ "$FULL" = "true" ]; then 
    if command -v conda > /dev/null; then
        echo "Conda already installed."
    else
        echo "Installing Miniconda"
        MINICONDA_VERSION="Miniconda3-latest"
        if [ "$OS" = "Linux" ]; then
            MINICONDA_SH="${MINICONDA_VERSION}-Linux-x86_64.sh"
        elif [ "$OS" = "MacOS" ]; then
            MINICONDA_SH="${MINICONDA_VERSION}-MacOSX-x86_64.sh"
        fi

        wget https://repo.anaconda.com/miniconda/$MINICONDA_SH
        bash $MINICONDA_SH -b -p $HOME/miniconda
        $HOME/miniconda/bin/conda init
        source $HOME/.bashrc
        echo "Miniconda installed."
        rm Miniconda*.sh
    fi

    conda init
    if [ "$SHELL" = "/bin/bash" ]; then
        source $HOME/.bashrc
    elif [ "$SHELL" = "/bin/zsh" ]; then
        source $HOME/.zshrc
    elif [ "$SHELL" = "/bin/csh" ] || [ "$SHELL" = "/bin/tcsh" ]; then
        source $HOME/.cshrc
    fi

    # Environment setup logic
else 
    if conda env list | grep -w "^seidart\s" > /dev/null; then
        echo "Conda environment 'seidart' exists."
        conda init
        conda activate seidart
        echo Uninstalling SeidarT if it has been installed and removing any previous wheels.

        # Uninstall the existing seidart package and clean build directories
        pip uninstall seidart 

        rm -rf _skbuild/ dist/ 
    else
        echo "Conda environment 'seidart' does not exist. Creating..."
        conda create -n seidart -f seidart-environment.yml
        conda init
        conda activate seidart
    fi

    if [ "$ENVUPGRADE" = "true" ]; then
        conda env update -n seidart -f seidart-environment.yml --prune
    fi
fi

# --------------------------------------------------------------------------
# Build and install the package
echo "Building the wheels..."
build_result=$(python setup.py bdist_wheel 2>&1)
if [ $? -ne 0 ]; then
    log_with_timestamp "Failed to build wheel package."
    log_with_timestamp "$build_result"
    echo "Failed to build wheel package. See install_log.txt for details."
    exit 1
fi

# --------------------------------------------------------------------------
# Install the wheel package
if [ "$SHELL" = "/bin/csh" ] || [ "$SHELL" = "/bin/tcsh" ]; then
    set wheel=`ls dist/*.whl`
    set install_result=`pip install $wheel 2>&1`
else
    wheel=$(ls dist/*.whl)
    install_result=$(pip install "$wheel" 2>&1)
fi

echo "Installing SeidarT..."

# Check installation success
echo "$install_result" | grep -q "Successfully installed"
if [ $? -eq 0 ]; then
    echo "seidart installed successfully."
    rm -rf _skbuild/ dist/
else
    log_with_timestamp "Failed to install seidart."
    log_with_timestamp "$install_result"
    echo "Failed to install seidart. See install_log.txt for details."
fi
