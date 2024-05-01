arraybuild
###################

*generates a csv file (receiver_array.csv) with the timeseries in columns*
*by receiver*

**Usage**

.. code-block:: bash

    arraybuild [-h] -p PROJECTFILE -r RCXFILE -c CHANNEL -i -z -d -s -P -C -g GAIN -e EXAGGERATION


**Inputs**

* ``-h``, ``--help``

    show this help message and exit

* ``-p PROJECTFILE``, ``--prjfile PROJECTFILE``

    The project file path

* ``-r RCXFILE``, ``--rcxfile RCXFILE``

    The file path for the text file of receiver locations

* ``-c CHANNEL``, ``--channel CHANNEL``

    The channel to query.

    | For radar: Ex, Ey, or Ez
    | For seismic: Vx, Vy, Vz, S1, S2, S3, S4, S5, S6
    | S# represents the stress tensor:

    * S1 - sigma_11
    * S2 - sigma_22
    * S3 - sigma_33
    * S4 - sigma_23
    * S5 - sigma_13
    * S6 - sigma_12

* ``-i``, ``--index``

    Indicate whether the receiver file contains coordinate
    indices or if these are the locations in meters.
    Default is False (meters)

* ``-z``, ``--is_complex``

    Flag if the saved model outputs are complex valued. 

* ``-d``, ``--double_precision``

    Flag if the saved model outputs are double precision. Values can be both complex and double precision.

* ``-s``, ``--save``

    Flag to save the array object. This will be saved as a .pkl (pickle format) file with the project name, channel, and source location in the filename. 

* ``-P``, ``--plot``

    Flag for plot output. This displays a Matplotlib plot that can be saved. 

* ``-C``, ``--plot_complex``

    If the model is output, flag if you want to plot the complex values of the time series.

* ``-g GAIN``, ``--gain GAIN``

    Specify the window length for the auto gain control function. A value of 0 or 1 is equivalent to 1-bit normalization and a value of the number of time steps is no gain applied. Defualt is the number of time steps.

* ``-e EXAGGERATION``, ``--exaggeration EXAGGERATION``

    Specify the float value for the horizontal to vertical ratio. The default value is 0.5. 

**Outputs**

CSV of amplitude values for all pixels surveyed (*receiver_array.csv*)


