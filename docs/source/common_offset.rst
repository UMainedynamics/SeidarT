Common Offset
~~~~~~~~~~~~~

We have wrapped other *SeidarT* functions into a class that will simulate a common offset survey. Along with the project file and the receiver file, a source location file needs to be provided. The source and receiver files are formatted the same with the same length. Each line in the source file corresponds to the same line in the receiver file. This allows the user to build a survey over topography, and even though it is called common offset, this can be very easily adapted for cross borehole surveys. 

The general steps to build a common offset profile are:

1. **Build the project file**. The shot locations set in the .prj file are arbitrary because the common offset simulation requires a separate shot file. 

2. **Create the receiver and source file**. If the source spacing is the same as the source-receiver offset, one way approach is to copy the receiver or source file and delete the first and last set of values for the source and receiver files, respectively. Another approach to creating the source or receiver file is to draw them in the png when building the model. In the *seidart.routines.definitions* module is the *rcxgen* function. The given rgb value will create a CSV file of all of the locations where that color is found. This color will have to be defined in the project file but can be a duplicate material.

3. **Create the object**. Import the *CommonOffset* class and initiate it with the required inputs: project file, source file, receiver file, and the channel. Currently, the survey is for a single channel. If you would like to see the cross polar outputs, rotate the source 90<sup>o</sup>. If you rotate it in the x-y direction (:math:`\phi`)

4. **Run the model**.

.. note::
    Common offset models can take a long time to run. To avoid having to re-run the model because the number of time steps are not sufficient, it is good practice to create a GIF from a single source model to visually inspect if the duration is long enough to capture the desired reflections. 

5. **Display the survey**. There are 

6. **Save outputs**. You can save the object as a `pickle format <https://docs.python.org/3/library/pickle.html>`_ (defualt) or specify a .csv output. There are benefits to either and no information is lost. The pickle format does maintain all of the information in the .prj file and source and receiver files that was loaded when the object was initiated, but it can be easily referenced if using a CSV output. 

.. note::
    For those that prefer to use analysis or visualization programs outside of Python, there are ways to read pickle files into `R <https://stackoverflow.com/questions/35121192/reading-a-pickle-file-pandas-python-data-frame-in-r>`_, `MATLAB <https://www.mathworks.com/matlabcentral/answers/1738695-run-pickle-file-in-matlab>`_, `Julia <https://stackoverflow.com/questions/65720584/how-to-load-python-pickle-from-julia>`_, etc. The hyperlinks will give you common examples of each. 


Below is a template for building and creating a common offset survey. You can find the source code and all files in the GitHub repository folder `src/seidart/recipes/common_offset <https://github.com/UMainedynamics/SeidarT/tree/main/src/seidart/recipes/common_offset>`_.

.. code-block:: python
    :linenos:
    
    # Import modules
    from seidart.simulations.common_offset import CommonOffset
    
    # Define the 
    prjfile = 'common_offset.prj' 
    rcxfile = 'common_offset_receivers.xyz'
    srcfile = 'common_offset_sources.xyz'
    
    # A polarization in the x-direction (phi=0, theta=0) was defined in the project file so we want outputs in the x-direction electric field.
    channel = 'Ex'
    
    # Let's compute the real components only. Change this to True if you would like to use the complex set of equations
    complex = False
    
    # Initiate the CO object
    co = CommonOffset(
        srcfile, 
        channel, 
        prjfile, 
        rcxfile, 
        receiver_indices = False, 
        is_complex = complex,
        single_precision = True
    )
    
    # Run the model. This may take some time
    co.co_run(seismic = False)
    
    # For display purposes lets put an AGC on each of the time series
    co.gain = 800
    
    # Plot axes need to be scaled appropriately. The exaggeration is the vertical/horizontal ratio
    co.exaggeration = 0.05
    
    # Plot all of the receivers. 
    co.sectionplot(plot_complex = complex)
    
    


