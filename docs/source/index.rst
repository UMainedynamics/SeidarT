.. SeidarT documentation master file, created by
   sphinx-quickstart on Thu Apr  4 10:10:36 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive. To generate the docs, you need to run 
   from the root directory 'sphinx-apidocs -o docs/source src/seidart -e'. 


Welcome to SeidarT's documentation!
===================================

.. _tutorial:

.. toctree::
   :numbered:
   :maxdepth: 2
   :caption: General Overview
   
   README.md

.. _modules:

.. toctree::
   :maxdepth: 2
   :caption: Module Components:
   
   
   ROUTINES.md
   seidart.routines.definitions
   seidart.routines.materials
   seidart.routines.prjbuild
   seidart.routines.sourcefunction
   seidart.routines.prjrun
   seidart.routines.arraybuild
   seidart.visualization.imgen
   seidart.visualization.im2anim
   seidart.simulations.commond_offset
   seidart.routines.orientsynth
   
.. toctree::
   :maxdepth: 2
   :caption: Command Line Interface
   
   CLI.md
   prjbuild
   sourcefunction
   prjrun 
   arraybuild
   im2anim
   

.. toctree::
    :maxdepth: 2
    :caption: Extra material

    png-creation
    build-docs
    materials_description
    units
    references

Search documentation
======================

Need to look something up?

* :ref:`search`