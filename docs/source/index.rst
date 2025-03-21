.. SeidarT documentation master file, created by
   sphinx-quickstart on Thu Apr  4 10:10:36 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive. To generate the docs, you need to run 
   from the root directory 'sphinx-apidocs -o docs/source src/seidart -e'. 


Welcome to SeidarT's documentation!
===================================

.. _intro: 

.. toctree::
   :maxdepth: 2 
   :caption: Introduction 
   
   README

.. _tutorial:

.. toctree::
   :maxdepth: 2
   :caption: General Overview
   
   getting_started
   initialize_project 
   project_file
   anisotropic_materials

.. _recipes:

.. toctree::
   :maxdepth: 2
   :caption: Recipes
   
   RECIPES

.. _commandline:

.. toctree::
   :maxdepth: 2
   :caption: Command Line Interface
   
   CLI.md
   prjbuild
   arraybuild
   im2anim
   
.. _modules:

.. toctree::
   :maxdepth: 2
   :caption: Module Components:
   
   ROUTINES
   seidart.routines.definitions
   seidart.routines.materials
   seidart.routines.prjbuild
   seidart.routines.arraybuild
   seidart.routines.fabricsynth
   VISUALIZATION.md
   seidart.visualization.imgen
   seidart.visualization.im2anim
   SIMULATIONS.md
   seidart.simulations.common_offset


.. toctree::
    :maxdepth: 2
    :caption: Extra material

    png-creation
    build-docs
    collaboration
    materials_description
    units
    references

Search documentation
======================

Need to look something up?

* :ref:`search`
