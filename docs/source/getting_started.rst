Getting Started
---------------

General workflow
~~~~~~~~~~~~~~~~

SeidarT uses a small number of input files and a repeatable workflow.

#. Prepare a PNG model image. Each unique RGB color maps to one material.
#. Create a receiver file as a CSV or XYZ-style text file with receiver coordinates.
#. If one or more materials are anisotropic, provide the Euler-angle file referenced by the project JSON.
#. Generate a project file with ``prjbuild``.
#. Edit the JSON project file to set the domain, materials, source, and survey settings.
#. Run the survey workflow that matches your experiment, such as single-shot, common-offset, common-midpoint, or multioffset processing.
#. Post-process the results with the plotting and animation utilities such as ``arraybuild``, ``im2anim``, ``implot``, and ``vtkbuild``.

Seismic outputs are written as wavefield and receiver data in units of m/s. Electromagnetic outputs use V/m for the field components.

Files to generate or edit
~~~~~~~~~~~~~~~~~~~~~~~~~

* *PNG image (.png)*

    This defines the geometry of your model system. A good starting size is
    50 to 500 pixels for each direction. Each RGB color represents a different
    material, so **the file must be saved with no antialiasing**. Typically each pixel represents the same distance in x and z (in meters).
    To get started on a new project create a new folder and save the image
    there. From the command line, change directories to the project folder
    then run::

        prjbuild -i /path/to/geometry/image.png -o project_filename.json

    Below, we describe the project-file structure and how to edit it.

* *receiver locations (text file, commonly receivers.xyz)*

    A comma-separated list of X,Y,Z coordinates, one receiver per line. The
    coordinates can be given in pixels or physical units, depending on how the
    project file is configured.

* *Euler-angle file* if using anisotropic materials

    A plain-text file with one Euler-angle triplet per line. The project file
    points to this file for materials marked as anisotropic.


* *project file (.json)*

    This file defines the domain, materials, source parameters, and survey
    settings for both seismic and electromagnetic runs. It is the main file you
    edit after running ``prjbuild``.
