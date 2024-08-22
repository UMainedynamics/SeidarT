SeidarT
=======

.. <!-- ### Table of Contents -->
.. <!-- [Introduction](#introduction)  
.. [Installation](#install)  
.. [Auto-Install](#auto-installation) 
.. [Manual Install](#manual-installation)  
.. [Hardware Requirements](#hardware-requirements)  
.. [Operating System Requirements](#operating-system-requirements)   -->

Full documentation appears in the docs folder. We have also created a separate GitHub repo for examples called `SeidarT-Recipes <https://github.com/UMainedynamics/SeidarT-Recipes>`_ which can be easily installed via pip. 

..  ======================================================================

Introduction
------------
The Seismic and Radar Toolbox (SeidarT) is a collaboration between researchers at the Universities of Maine and Washington to provide an open source platform for forward modeling mechanical and electromagnetic wave propagation. The major objective of the project is to easily and quickly implement isotropic and anisotropic complex geometries and/or velocity structures to develop prior constraints for - not limited too - investigating, estimating, and imaging englacial ice structure, sub-glacial boundary conditions on the sub-regional scale. Larger problems would require the curvature of the Earth to be taken into consideration, but many glacier seismic and radar experiments do not expand into regional parameter estimation and velocity modeling.

Much of the Staggered Grid FDTD code has been adopted from the *SEISMIC_CPML* software provided by  `Computational Infrastucture for Geophysics (CIG) <https://geodynamics.org/cig/software/>`_. Further details to the backend numerical code can be found in the [References](#references) section below.


Installation
------------

The dynamical programming language of **Python3** is used as a command line interface to run the more computationally extensive modeling schemes in **Fortran**. There are a number of ways to install this software depending on your desired level of control over the process. Most users should be fine with the "automatic installation" in the section below.

SeidarT package binaries are publicly availble on the `PyPi repository <https://pypi.org/project/seidart/>`_ and `source code <https://github.com/UMainedynamics/SeidarT>`_ can be found on GitHub. 


.. ============================================================================

"Auto" installation  
^^^^^^^^^^^^^^^^^^^

For Windows users, see the VM setup. Unix/Linux users can download the tarball `with this link <https://github.com/UMainedynamics/SeidarT/blob/main/install.tar.xz?raw=1>`_. 

Extract the *install* directory from the *install.tar.gz* which includes an install script, *full_install.sh*, and the *seidart-environment.yml*. It's not necessary to know or do much more than execute a few command line entries via a bash terminal or powershell terminal. The install script checks for and installs if necessary the Anaconda/Miniconda package manager. A virtual environment is created to avoid causing system incompatibilities and complicated software dependencies. If Anaconda/Miniconda is not defined in your 'PATH' variable than it will be installed using the default install location. After installing the Conda package, the *seidart* environment is built using pre-defined dependencies in a YAML (Yet Another Markup Language) file. Both *Bash* executables and Python modules are built during install. In order to use either, the environment must be active. This can be easily done from the *Bash* command line interface (CLI) using the command::

    conda activate seidart

Documentation for managing conda environments with Miniconda or Anaconda can be found `here <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_. The full Anaconda release has a GUI called Navigator for managing environments. This can be found on the `webpage <https://docs.anaconda.com/free/navigator/tutorials/manage-environments/>`_.  

.. -----------------------------------------------------------------------------

VM Install 
^^^^^^^^^^

There are a few different options for virtual machine software, but *VirtualBox* is a robust and free VM management software that can be installed on Windows, MacOS, and Linux. Follow the directions to download and install the software from the `website <https://www.virtualbox.org/>`_. 

There are a couple options for creating a virtual machine for *SeidarT*. The simplest is to download the `files <https://drive.google.com/drive/folders/1zVzlKLug95wfy6NCwYGtsbD_cJK8CW1S?usp=drive_link>`_ (~15GB) for the VM clone from Google Drive. The VM clone requires 4 GB of RAM. Currently, the clone is setup as 15 GB of hard disk space which, even though it is a large file, is still limited in capacity. It is recommended that an external drive be used for creating and building models. It is common for a few GB of hard disk space to be used up when running the models which can easily be recovered by deleting the .dat outputs. 

After installing VirtualBox and downloading the VM clone files, launch the VirtualBox software. In the VM manager, click on the *Add* button and you will be prompted to choose a .vbox file. Navigate to the directory with the VM clone and select it. This will launch a Debian Linux clone with a GNOME desktop environment. The username and password are *seidart* which can be changed. `Here <https://reintech.io/blog/managing-users-groups-debian-12>`_ is an example of how to do so. To get started, open up a *Bash* terminal and activate the conda environment (see above) and start a Python session by entering into the command line:;

    python 

For users that would like to build a VM with more control, different preferences (i.e. hard disk space), for a different flavor of Linux, or using a different VM manager you will need to download the .iso file for the desired Linux then create a new VM. This will prompt you with the setup parameters. Following setup, you can launch your VM and open up a terminal. From here you can follow the *Auto Installation* (above) or  *Manual Installation* (below) directions.

.. note::
    Windows offers the Windows Subsystem for Linux (WSL) which allows users to run a Linux environment without the need for a separate virtual machine or dual boot. `Setup and Install <https://learn.microsoft.com/en-us/windows/wsl/install>`_ of WSL is relatively easy, and file systems, desktop environments, and clipboards are shared. If using WSL follow the auto-install or the manual install directions. 

.. -----------------------------------------------------------------------------

Manual installation
^^^^^^^^^^^^^^^^^^^

The full repo can be found on GitHub and is hosted on PyPi. *SeidarT* has been tested on Python 3.11 and is not yet supported with Python 3.12. For users that prefer building virtual environments with Anaconda, the install folder contains the *seidart-environment.yml* or it can be found in the root directory of the GitHub repo. To clone the directory, you will need to setup a GitHub account then open a bash terminal and enter the command:

    git clone git@github.com:UMainedynamics/SeidarT.git

Change directories into the install folder 

    cd SeidarT/install

Create and build the *seidart* environment::

    conda env create -f seidart-environment.yml

This will install all dependencies and the latest *seidart* version found on PyPi. For users who prefer more control in their installation, their is a small list of dependencies that must be met. These are:  *gcc*>10, *gfortran*, *ghostscript*, *imagemagick*, *numpy*, *pandas*, *matplotlib*, *scipy*, *glob2*, *pyevtk*, *mplstereonet*. For Linux users, all of these packages can be installed from Anaconda, however for MacOS users, suitable GCC and GFortran versions aren't supported or up to date in Anaconda. If using MacOS you will need to install a few packages via Homebrew::

    xcode-select --install
    brew install gcc gfortran

Additionally for MacOS users, make sure that the conda environment doesn't default to the older GCC and GFortran versions and instead uses the system wide compiler. To do so remove them from the environment::

    conda remove gfortran gcc 

Following install of all dependencies,:: 

    pip install seidart

will pull the latest version and install the package from PyPi. 

Alternatively, if you run into incompatibility with the operating system and the package wheels, you can install locally. To do this, activate your conda environment then navigate to the SeidarT folder. This is the root directory and from here the command::
    
    pip install . 

will install dependencies, build wheels based on your OS, and install SeidarT into the environment. For development, or for testing specific repo branches, this method works great. If you want to keep up to date on the version and this is the install method that works for you, pull updates from the main branch first::
    
    git pull origin main 

.. =============================================================================

Hardware Requirements
---------------------

*SeidarT* was tested and developed on a quad core 5th gen i7 processor with 16 Gb of RAM without any burden on the system so a typical modern laptop is sufficient for many application. When running models with large domains or a high number of time steps, the computational load is obviously increased, however the storage requirements become more significant. It can be easy to fill up 10's of Gb of storage, but an external drive can resolve that problem. The Apple M-chips may have compatability issues with particular types of software and Python packages, but we have maintained a relatively simple design along with leveraging some of the most commonly used Python packages which should help to mitigate any issues with computing on an M-chip. 

.. =============================================================================

Operating System requirements
-----------------------------

All of the development was carried out on a Linux operating system and limited to Debian, Ubuntu, Solus 2, and Fedora. No compatibility issues between Linux flavors arose. The binaries are built on Github Actions for Windows 10 and 11 (latest), MacOS 13 and 14 (latest), and most flavors of Linux. Cross-platform usability is one of the core tenets in the development of the software and needs to be maintained in future development. 

.. =============================================================================

Upgrading Versions
