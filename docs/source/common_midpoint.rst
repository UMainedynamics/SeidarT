Common Midpoint
~~~~~~~~~~~~~~~

Simulating a common midpoint (CMP) survey is relatively easy in SeidarT. A typical CMP is likely to have a fixed linear array and multiple source points along the array, although there are many different techniques for achieving different results. This tutorial creates a single source point, but can easily be expanded to having multiple source points (and/or multiple array layouts) by duplicating the project file for each source location, changing the source location in the project file and rerunning the model, or by looping through an array of source locations in a Python script. The Fortran binary files corresponding to each model run will contain the source location in the filename so that each model doesn't overwrite the previous model. 

The general steps to create a single shot survey are:

1. **Build the project file**

2. **Create the receiver file**. This is a CSV delimited file with header fields X,Y,Z. The values can be the indices of the receiver (i.e. nodal point) or the distance from the origin - :math:`\Delta x_i \cdot \text{Indice}_{\text{rcx}}`. The origin is defined as the top left corner of the model PNG. 

3. **Load the project file**. This creates the model, domain, and material objects required for modeling.

4. **Create source function**. 

5. **Generate tensor coefficients for each material**. The example below corresponds to a project file that already has tensor coefficients. These can be edited and loaded into the model objects. 

6. **Run the model**. 

7. **Create the Array object**. 

8. **Visualize and/or save outputs**.  

Below is a simple template for a single shot survey for a 91 channel linear array at the surface. All of the files including the source code can be found in the `src/seidart/recipes/single_source <https://github.com/UMainedynamics/SeidarT/tree/main/src/seidart/recipes/single_source>`_ folder. 


.. code-block:: python
    :linenos:
    
    import numpy as np 
    from seidart.routines import prjbuild, prjrun, sourcefunction
    from seidart.routines.arraybuild import Array
    from seidart.visualization.im2anim import build_animation

    # Define the necessary files. Adjust the relative paths if necessary. 
    prjfile = 'single_source.prj' 
    rcxfile = 'receivers.xyz'

    # Initiate the model and domain objects
    dom, mat, seis, em = prjrun.domain_initialization(prjfile)

    # Compute the permittivity coefficients and check to make sure the project file has all required values
    prjrun.status_check(
        em, mat, dom, prjfile, seismic = False, append_to_prjfile = True
    )
    
    # Create the source function
    timevec, fx, fy, fz, srcfn = sourcefunction(em, 10, 'gaus1', 'e')
    
    # The non-complex equations aren't necessary but are also a solution to the PDE
    complex_values = False
    prjrun.runelectromag(em, mat, dom, use_complex_equations = complex_values)
    
    # Create the array object
    array_ex = Array('Ex', prjfile, rcxfile, is_complex = complex_values)
    # Add an AGC function for visualization
    array_ex.gain = int(em.time_steps/3)
    # We need to scale the axes
    array_ex.exaggeration = 0.1
    # Create the plot 
    array_ex.sectionplot(
        plot_complex = False
    )
    
    # Create the GIF so that we can view the wavefield
    build_animation(
            prjfile, 
            'Ex', 10, 10, 0.3, 
            is_complex = complex_values, 
            is_single_precision = True
    )

    # --------------------------------------------------------------------------
    # We can do the same for the vertical electric field as above
    array_ez = Array('Ez', prjfile, rcxfile, is_complex = complex_values)
    array_ez.gain = int(em.time_steps/3)
    array_ez.exaggeration = 0.1
    array_ez.sectionplot(
        plot_complex = False
    )
    build_animation(
            prjfile, 
            'Ex', 10, 10, 0.3, 
            is_complex = complex_values, 
            is_single_precision = True,
            plottype = 'energy_density'
    )