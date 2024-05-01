prjrun
##########################

*Performs single-shot wave propagation with receiver locations specified*
*later in postprocessing.*

**Usage**

.. code-block:: bash

    prjrun -p PROJECTFILE -m [n e s] -a [0 1] -d -c -o DATA_OUTPUT_DIRECTORY

**Inputs**

* ``-p PROJECTFILE``, ``--prjfile PROJECTFILE`` .prj file

    The file path for the project file, completely filled in for the model
    type used except for permittivity and stiffness coefficients, and dt

* ``-m [n e s]`` Which model to run

    * ``n`` calculates only timesteps and material tensors, necessary before running sourcefunction
    * ``e`` electromagnetic propagation
    * ``s`` seismic wave propagation

* ``-a`` (optional)

    Flag to append/recompute the coefficients to the permittivity and
    stiffness matrices; Default is False. Do not
    recompute if you have made manual changes to the matrices in the prj file.

* ``-d`` (optional)

    Flag to store output as double precision. Default is single precision. Higher precision produces larger files. 
    
* ``-c`` (optional)

    Flag if the complex equations should be used for the computation of the electromagnetic model. Default is False. 

* ``-o DATA_OUTPUT_DIRECTORY`` (not yet integrated)

    Provide the output data directory if different from the current working directory to store outputs from each time step. 

**Outputs**

.dat files equal in number to the number of time steps specified in the .prj file.

