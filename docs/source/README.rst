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


Dependencies
^^^^^^^^^^^^

SeidarT was built with the goal of having minimal system dependencies to avoid system updates causing deprecations and incompatability. First, ensure that the GCC compiler is up to date with version 14 or greater. Second, install Miniconda (recommended) or Anaconda by following the directions posted on their `install page <https://docs.anaconda.com/miniconda/install/>`_.

.. note::

    Documentation for managing conda environments with Miniconda or Anaconda can be found `here <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_. The full Anaconda release has a GUI called Navigator for managing environments. This can be found on the `webpage <https://docs.anaconda.com/free/navigator/tutorials/manage-environments/>`_.


Installation
^^^^^^^^^^^^

*SeidarT* was built to be installed on 64bit Arm and Intel architectures for Linux and Mac operating systems. For Windows users , `Windows Subsystem for Linux <https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux>`_ is recommended. Windows Subsystem for Linux (WSL) which allows users to run a Linux environment without the need for a separate virtual machine or dual boot. `Setup and Install <https://learn.microsoft.com/en-us/windows/wsl/install>`_ of WSL is relatively easy, and file systems, desktop environments, and clipboards are shared.  *SeidarT* has been tested on Fedora, Solus 2, and Debian which covers a variety of flavors of Linux. The dynamical programming language of **Python3** is used as a command line interface to run the more computationally extensive modeling schemes in **Fortran**. There are a number of ways to install this software depending on your desired level of control over the process.

SeidarT package binaries are publicly availble on the `PyPi repository <https://pypi.org/project/seidart/>`_ and `source code <https://github.com/UMainedynamics/SeidarT>`_ can be found on GitHub. To get started, clone or download and extract the Github repo. The dropdown in the green "<> Code" will provide the option to download the ZIP file. If cloning via command line you will need to use the HTTPS link. From the command line::

    git clone https://github.com/UMainedynamics/SeidarT.git 

will create the SeidarT folder in the current working folder. Change directories into the root folder (SeidarT) and you will find the environment YAML file for creating the Miniconda/Anaconda environment and install all dependencies. Create the environment from the Bash terminal with the following commands::

    conda create -f seidart-environment.yml 
    conda activate seidart

To finish installing run::

    pip install . 

This will install the SeidarT package via source code in the *src* folder. Many of the standard routines are setup as entry points and can be called via the command line. To easily test if SeidarT was installed, you can call the help menu of one of the routines. For instance,::

    prjbuild -h 

If a **command not found** output is returned, make sure that the *seidart* environment is activated and was activated during the *pip install .* command. If the problem persists then it is likely an issue with the precompiled binaries. The precompiled binaries are in the *src/seidart/binaries* folder which will be put in the correct path during the *pip* install, however it is possible that your system architecture is incompatible even if you are using Mac or Linux. Please submit a ticket on the repository `issue tracker <https://github.com/UMainedynamics/SeidarT/issues>`_ This will help accommodate more system architectures and operating systems. 

Building the binaries locally should solve any incompatibilities. To do so, you will need to clone the *SeidarT-Backend* repo::

    git clone https://github.com/UMainedynamics/SeidarT-Backend.git

Change directories into the SeidarT-Backend folder and activate the *seidart* environment if it is inactive. You can run the build script via Bash terminal::

    bash build.sh 

This will create the *seidart* executable and put them in the binaries folder associated to the environment. You can check to make sure there is only one *seidartfdtd* executable::
 
    ls $CONDA_PREFIX/bin/seidartfdtd*

If you see multiple *seidartfdtd* files, delete them and recompile again with the *build.sh* script seen above. 


.. =============================================================================

Hardware Requirements
^^^^^^^^^^^^^^^^^^^^^

*SeidarT* was tested and developed on a quad core 5th gen i7 processor with 16 Gb of RAM without any burden on the system so a typical modern laptop is sufficient for many application. When running models with large domains or a high number of time steps, the computational load is increased, however the storage requirements become more significant. It can be easy to fill up 10's of Gb of storage, but an external drive can resolve that problem. We have tested *SeidarT* on a limited number of operating systems due to practical considerations. The Apple M chips and previous Intel based Apple computer chips have not been an issue. Firmware and software variations between operating systems are most likely the cause of incompatibility issues. For animations of large models, there may be some limits from graphics cards that can't sufficiently support them.

.. =============================================================================

Operating System requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All of the development was carried out on a Linux operating system and limited to Debian, Ubuntu, Solus 2, and Fedora. No compatibility issues between Linux flavors arose. The binaries are built on Github Actions for Windows 10 and 11 (latest), MacOS 13 and 14 (latest), and most flavors of Linux. Cross-platform usability is one of the core tenets in the development of the software and needs to be maintained in future development. 

.. =============================================================================

Upgrading Versions
^^^^^^^^^^^^^^^^^^

If you want to keep up to date on the version and this is the install method that works for you, pull updates from the main branch first::
    
    git pull origin main 

