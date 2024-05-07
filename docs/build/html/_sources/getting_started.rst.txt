Getting Started
---------------

General workflow
~~~~~~~~~~~~~~~~

Using SeidarT follows a relatively simple workflow.

#. You need two or three files to start:

  * A 2D image saved in *png* format.
  * A *csv* file listing the X,Y,Z coordinates of receivers for your survey
  * If your material is anisotropic, you need a file in the format delimited file specifying the

  Euler angles for a number of crystals, with one triplet per line. See an example orientation
  file and/or generate one using the ``orientation_tensor`` function.

#. Generate a project file (using ``prjbuild``) and edit that text file to set up your survey.
#. Create files describing the radar or seismic source (``sourcefunction``).
#. Choose the style of survey you want to do [single shot, common offset, common midpoint, or (in development) polarimetric] and run the calculations.
#. For single shot, you can create an animation of the wave propagation (``im2anim`` for 2D or ``vtkbuild`` for 2.5D).
#. Display your results as radar- or seismograms, or wiggle plots. You can also save the timeseries-receiver data in a *csv* file for further processing in different software.

Output from the seismic model is m/s and from the radar model is 

Files to generate or edit
~~~~~~~~~~~~~~~~~~~~~~~~~

* *PNG image (.png)*

    This defines the geometry of your model system. A good starting size is
    50 to 500 pixels for each direction. Each RGB color represents a different
    material, so **the file must be saved with no antialiasing**. Typically each pixel represents the same distance in x and z (in meters).
    To get started on a new project create a new folder and save the image
    to the folder. From the command line, change directories to the
    project folder then enter the following::

        prjbuild -i /path/to/geometry/image.png -p project_filename.prj

    Below, we describe the *prj* file structure and how to edit it.

* *receiver locations (text file, commonly receivers.xyz)*

    A comma separated list of X,Y,Z coordinates (one set per line,
    with X,Y,Z as the first line) for receiver locations. Can use pixels, but
    more typically meters as the units.


* *project file (.prj)*

    This file is the heart of the software. It defines domain values, material properties, and survey conditions for
    electromagnetic and seismic runs. Here, we identify what each line means and which to edit.
    All lines with # are comments. Bold text indicates a line the user should edit.