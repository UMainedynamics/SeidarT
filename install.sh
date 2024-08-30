#!/bin/sh

# Function to log messages with a timestamp
log_with_timestamp() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >> install_log.txt
}

# Uninstall the existing seidart package and clean build directories
pip uninstall seidart 
rm -rf _skbuild/ dist/ 

# Build the wheel package
build_result=$(python setup.py bdist_wheel 2>&1)
if [ $? -ne 0 ]; then
    log_with_timestamp "Failed to build wheel package."
    log_with_timestamp "$build_result"
    echo "Failed to build wheel package. See install_log.txt for details."
    exit 1
fi

# Install the wheel package
if [ "$SHELL" = "/bin/csh" ] || [ "$SHELL" = "/bin/tcsh" ]; then
    set wheel=`ls dist/*.whl`
    set install_result=`pip install $wheel 2>&1`
else
    wheel=`ls dist/*.whl`
    install_result=$(pip install "$wheel" 2>&1)
fi

# Check if installation was successful
echo "$install_result" | grep -q "Successfully installed"
if [ $? -eq 0 ]; then
    echo "seidart installed successfully."
else
    log_with_timestamp "Failed to install seidart."
    log_with_timestamp "$install_result"
    echo "Failed to install seidart. See install_log.txt for details."
fi
