SeidarT
=======

.. <!-- ### Table of Contents -->
.. <!-- [Introduction](#introduction)  
.. [Installation](#install)  
.. [Auto-Install](#auto-installation) 
.. [Manual Install](#manual-installation)  
.. [Hardware Requirements](#hardware-requirements)  
.. [Operating System Requirements](#operating-system-requirements)   -->

SeidarT is a Python and Fortran toolbox for full-waveform seismic and electromagnetic modeling in snow, ice, and other heterogeneous materials.

It uses a PNG-plus-JSON project workflow to define geometry, materials, sources, receivers, and boundary conditions. The repository also includes command-line tools and Python modules for building project files, running simulations, extracting receiver data, and generating plots or animations.

Full documentation is in the `docs <https://umainedynamics.github.io/SeidarT/docs/build/html/index.html>`_ folder and published on `GitHub Pages <https://umainedynamics.github.io/SeidarT/docs/build/html/index.html>`_. Example workflows live in the separate `SeidarT-Recipes <https://github.com/UMainedynamics/SeidarT-Recipes>`_ repository.

-------------------------------------------------------------------------------------

Introduction
------------

SeidarT is an open-source wave-physics toolkit for finite-difference time-domain modeling of elastic and electromagnetic propagation on a Cartesian grid. It is aimed at snow, ice, and similar heterogeneous media, but the workflow is general enough for other geological and engineered materials.

The core solvers are implemented in Fortran, while the Python package provides project-file generation, simulation setup, receiver extraction, plotting, animation, and VTK export. The package supports anisotropic materials, frequency-independent attenuation, and both 2D and 2.5D-style workflows through the project file and companion command-line tools.


Dependencies
^^^^^^^^^^^^

SeidarT is designed to keep the computing and documentation environments separate. For simulation work, use the `seidart-environment.yml` file. For docs generation, use `documentation-environment.yml`.

First, ensure that the GCC `gfortran` compiler is installed. Then install Miniconda or Anaconda using the directions on their `install page <https://docs.anaconda.com/miniconda/install/>`_.

.. note::

    Documentation for managing conda environments with Miniconda or Anaconda can be found `here <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_. The full Anaconda release has a GUI called Navigator for managing environments. This can be found on the `webpage <https://docs.anaconda.com/free/navigator/tutorials/manage-environments/>`_.

SeidarT uses Python 3 as the high-level interface and Fortran for the time-critical finite-difference solvers. Runtime dependencies include NumPy and the scientific Python stack used for project setup, plotting, and post-processing; see the environment YAML files for the exact dependency lists.

Installation
^^^^^^^^^^^^

*SeidarT* is intended for Linux and macOS systems on 64-bit Intel or ARM hardware. On Windows, `Windows Subsystem for Linux <https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux>`_ is the supported path. The Python layer drives the command-line workflows, and the Fortran backend does the heavy numerical work.

`Source code <https://github.com/UMainedynamics/SeidarT>`_ is available on GitHub. To get started, clone the repository and move into the root directory::

    git clone https://github.com/UMainedynamics/SeidarT.git 

will create the SeidarT folder in the current working folder. Change directories into the root folder (SeidarT) and you will find the environment YAML file for creating the Miniconda/Anaconda environment and install all dependencies. Create the environment from the Bash terminal with the following commands::

    conda env create -f seidart-environment.yml
    conda activate seidart

Then install the package from source::

    pip install . 

This installs the SeidarT package from the `src` layout. To confirm the install, run::

    conda list | grep seidart 

If there are any issues on install please submit a ticket on the repository `issue tracker <https://github.com/UMainedynamics/SeidarT/issues>`_.

Fortran binaries are built locally to avoid operating-system incompatibilities. To do so, clone the *SeidarT-Backend* repo::

    git clone https://github.com/UMainedynamics/SeidarT-Backend.git

Change directories into the SeidarT-Backend folder and activate the *seidart* environment if it is inactive. Then run::

    bash build.sh 

This creates the `seidartfdtd` executable in the environment `bin` directory. You can confirm it with::
 
    ls $CONDA_PREFIX/bin/seidartfdtd*

If you see multiple *seidartfdtd* files, delete them and recompile again with the *build.sh* script seen above. 

*SeidarT-Recipes* does not need to be installed, but it contains example workflows and input files. Clone it wherever you keep examples::
    
    git clone https://github.com/UMainedynamics/SeidarT-Recipes.git 

It uses a `src` layout and includes runnable examples under `src/seidart-recipes`.


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
