SeidarT
=======

.. <!-- ### Table of Contents -->
.. <!-- [Introduction](#introduction)  
.. [Installation](#install)  
.. [Auto-Install](#auto-installation) 
.. [Manual Install](#manual-installation)  
.. [Hardware Requirements](#hardware-requirements)  
.. [Operating System Requirements](#operating-system-requirements)   -->

SeidarT: Seismic and Radar Toolbox for full-waveform modeling in snow, ice, and other complex materials.

Full documentation appears in the docs folder, but can also be accessed in `Github Pages <https://umainedynamics.github.io/SeidarT/docs/build/html/index.html>`_.  We have also created a separate GitHub repo for examples called `SeidarT-Recipes <https://github.com/UMainedynamics/SeidarT-Recipes>`_
..  ======================================================================

Introduction
------------

Introduction
------------

The Seismic and Radar Toolbox (SeidarT) is an open-source, community-driven software package for full-waveform finite-difference time-domain (FDTD) modeling of elastic and electromagnetic wave propagation in heterogeneous, anisotropic, and attenuating media.  SeidarT targets applications in snow and ice but is general enough for other geological and engineered materials, including porous and partially saturated media.

SeidarT natively supports full-tensor anisotropy (all 21 stiffness coefficients), frequency-independent attenuation via a generalized Q formulation, and unified treatment of both seismic and electromagnetic wave physics on a simple Cartesian grid.  The computational core is implemented in Fortran for performance, while a Python interface provides an accessible command-line and scripting environment for model construction, workflow automation, and post-processing.

Models are defined using an intuitive image-plus-JSON workflow: users draw the model geometry as a PNG image (e.g., in Gimp or Inkscape), where pixel colors map to material IDs, and configure the simulation through a human-readable JSON project file specifying domain size, material properties (stiffness or permittivity tensors, density, conductivity), sources, and boundary conditions.  This approach eliminates complex meshing and makes it easy to implement isotropic and anisotropic complex geometries and velocity structures for building priors, performing parameter studies, and generating synthetic data.

Although the manuscript focuses on cryospheric use cases (e.g., glacier internal structure, temperate layers, and englacial aquifers), SeidarT is designed as a general wave-physics toolbox.  It can be used for problems such as environmental monitoring, subsurface and infrastructure characterization, and other scenarios where fully coupled elastic or electromagnetic wavefields are needed.


Dependencies
^^^^^^^^^^^^


SeidarT was built with the goal of having minimal system dependencies to avoid system updates causing deprecations and incompatibility. First, ensure that the GCC (gfortran) compiler is up to date with version 14 or greater. Second, install Miniconda (recommended) or Anaconda by following the directions posted on their `install page <https://docs.anaconda.com/miniconda/install/>`_.

.. note::

    Documentation for managing conda environments with Miniconda or Anaconda can be found `here <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_. The full Anaconda release has a GUI called Navigator for managing environments. This can be found on the `webpage <https://docs.anaconda.com/free/navigator/tutorials/manage-environments/>`_.

SeidarT uses Python 3 as the high-level interface and Fortran for the time-critical finite-difference solvers.  Runtime dependencies include NumPy and standard scientific Python tools for I/O and visualization; see the environment YAML file in the repository for the exact dependency list. 

Installation
^^^^^^^^^^^^

*SeidarT* was built to be installed on 64bit Arm and Intel architectures for Linux and Mac operating systems. For Windows users , `Windows Subsystem for Linux <https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux>`_ is recommended. Windows Subsystem for Linux (WSL) allows users to run a Linux environment without the need for a separate virtual machine or dual boot. `Setup and Install <https://learn.microsoft.com/en-us/windows/wsl/install>`_ of WSL is relatively easy, and file systems, desktop environments, and clipboards are shared.  *SeidarT* has been tested on Fedora, Solus 2, and Debian which covers a variety of flavors of Linux. The dynamical programming language of **Python3** is used as a command line interface to run the more computationally extensive modeling schemes in **Fortran**. There are a number of ways to install this software depending on your desired level of control over the process.

`Source code <https://github.com/UMainedynamics/SeidarT>`_ can be found on GitHub. To get started, clone or download and extract the Github repo. The dropdown in the green "<> Code" will provide the option to download the ZIP file. If cloning via command line you will need to use the HTTPS link. From the command line::

    git clone https://github.com/UMainedynamics/SeidarT.git 

will create the SeidarT folder in the current working folder. Change directories into the root folder (SeidarT) and you will find the environment YAML file for creating the Miniconda/Anaconda environment and install all dependencies. Create the environment from the Bash terminal with the following commands::

    conda create -f seidart-environment.yml 
    conda activate seidart

To finish installing run::

    pip install . 

This will install the SeidarT package via source code in the *src* folder. To be sure that the installation was successful enter:: 
    
    conda list | grep seidart 

If there are any issues on install please submit a ticket on the repository `issue tracker <https://github.com/UMainedynamics/SeidarT/issues>`_ This will help accommodate more system architectures and operating systems. 

Fortran binaries are built locally to avoid operating system incompatibilities. To do so, you will need to clone the *SeidarT-Backend* repo::

    git clone https://github.com/UMainedynamics/SeidarT-Backend.git

Change directories into the SeidarT-Backend folder and activate the *seidart* environment if it is inactive. You can run the build script via Bash terminal::

    bash build.sh 

This will create the *seidart* executable and put them in the binaries folder associated to the environment. You can check to make sure there is only one *seidartfdtd* executable::
 
    ls $CONDA_PREFIX/bin/seidartfdtd*

If you see multiple *seidartfdtd* files, delete them and recompile again with the *build.sh* script seen above. 


.. =============================================================================

Hardware Requirements
^^^^^^^^^^^^^^^^^^^^^

SeidarT was tested and developed on a quad-core 5th-generation Intel i7 processor with 16 GB of RAM without any noticeable burden on the system, so a typical modern laptop is sufficient for many applications.  For models with large domains, high temporal sampling, or 2.5D surveys, computational load and especially storage requirements increase; full-wavefield output can easily reach tens of gigabytes, so using an external or high-capacity internal drive is recommended. 

The software has been exercised on a limited number of operating systems in practice, but runs without issue on Apple M-series chips and previous Intel-based Apple architectures.  Remaining incompatibilities are most likely related to system-specific firmware, drivers, or library versions rather than to the SeidarT source itself, and issues can usually be resolved by rebuilding the Fortran backend in the target environment. 

For animations or interactive visualization of large models, GPU and display limitations may become the bottleneck rather than CPU or memory; in those cases, subsampling output in space or time and exporting pre-rendered animations can improve usability. 

.. =============================================================================
Operating System requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All core development was carried out on Linux, primarily Debian, Ubuntu, Solus 2, and Fedora, and no compatibility issues were observed between these distributions when using locally compiled Fortran binaries and conda environments.  The SeidarT backend is built and tested via continuous integration on GitHub Actions for Windows 10/11 (latest), macOS 13/14 (latest), and multiple Linux distributions, and binaries are intended to be portable across common 64‑bit Intel and ARM platforms. 

Cross‑platform usability is a core tenet of SeidarT’s design.  For users on native Windows systems, Windows Subsystem for Linux (WSL) is recommended and has been successfully used to run the Linux toolchain and simulations. 

.. =============================================================================

Upgrading Versions
^^^^^^^^^^^^^^^^^^

If you want to keep up to date with the latest version, then in the top folder for the repository pull updates from the main branch first::

    git pull origin main

After pulling changes in the SeidarT-Backend folder, re-run the build of the Fortran backend to ensure that any updates to the numerical kernels are reflected in your local binaries.  When upgrading across major versions, consult the project documentation and changelog for notes on any changes to the JSON project format, material-parameter handling, or default stability and CPML settings. 
