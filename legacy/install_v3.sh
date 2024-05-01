#!/bin/bash

# Install script for SEIDART toolbox

# ----------------------------- Anaconda install ------------------------------
VER="v2"
echo "--------------------------------------------
seidart Anaconda-based installer $VER
Univ. of Maine / Univ. of Washington, 2020
--------------------------------------------
This installer will check for Anaconda/Miniconda
and install a seidart environment prior to compiling from
source.
You will have the option to install Miniconda
if no existing conda is found.

Please follow instructions in prompts.
"
read -rp $'Press Enter to continue...\n'


echo '
Do you wish to delete source files after installation?
This will not affect how the program runs, but you will
not be able to edit source code.
'
# ask if we should enable pure end-user mode and delete the source scripts
read -rp $'Type "yes" and press Enter to delete source files. (default: no)\n' EUMODE

if [[ "$EUMODE" == "yes" || "$EUMODE" == "Yes" || "$EUMODE" == "YES" ]] ; then
        # end-user mode
        EUMODE="yes"
else
        # developer mode
        unset EUMODE
        echo "Developer mode enabled. Source scripts will not be deleted."
fi

bash conda_deps.sh ||
echo "Conda installation failed. Try installing dependencies the run the noconda_install script." ||
exit 1

`grep etc/profile.d/conda.sh ~/.bashrc` &&
conda activate seidart &&
echo "conda activate seidart : Successfully activated seidart environment." ||
echo "Could not find seidart conda environment. Exiting." ||
exit 1

echo ""
echo "Starting compiling/installation process..."

# -----------------------------------------------------------------------------
# Make sure we have a folder to put everything in
if [ ! -d "bin" ]; then
	mkdir bin
fi

# --------------------------------- Clean Up ----------------------------------
# Typically, install will overwrite everything but when there is a compile
# error this can make debugging difficult if it goes unnoticed. Not everyone
# can read the f2py or gfortran outputs

# Just clean up if specified
CLEAN=${1:-"none"}
if [[ $CLEAN == "clean" ]] ; then
        # Clear the bin folder
        rm -rf bin/*
        # Remove .mod and .o (if any) files generated during the fortran compile
        rm source/*.mod
        rm source/*.o
        exit 1
fi


if [ ! -z ${EUMODE+x} ]; then
        echo "Deleting source files..."
        rm -rfv exe materials vis io survey_wrappers
        echo '
[end user mode]:
Source files deleted. You will need to download the
software again if you would like to edit source files.
Developers should not commit changes on this branch as
this will cause unwanted consequences in source control.
Developers can also use the git command "git restore ."
to restore all deleted files to source control.'
        unset EUMODE
fi

echo '
Done!
Type "conda activate seidart" to access the seidart environment.'
