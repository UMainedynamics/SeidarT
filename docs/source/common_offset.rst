Common Offset
~~~~~~~~~~~~~

We have wrapped other *SeidarT* functions into a class that will simulate a common offset survey. Along with the project file and the receiver file, a source location file needs to be provided. The source and receiver files are formatted the same with the same length. Each line in the source file corresponds to the same line in the receiver file. This allows the user to build a survey over topography, and even though it is called common offset, this can be very easily adapted for cross borehole surveys. 

The general steps to build a common offset profile are:

1. Create the receiver and source file. If the source spacing is the same as the source-receiver offset, one way approach is to copy the receiver or source file and delete the first and last set of values for the source and receiver files, respectively. Another approach to creating the source or receiver file is to draw them in the png when building the model. In the *seidart.routines.definitions* module is the *rcxgen* function. The given rgb value will create a CSV file of all of the locations where that color is found. This color will have to be defined in the project file but can be a duplicate material.
 
2. Import the *CommonOffset* class from simulations. 
::
    from seidart.simulations.common_offset import CommonOffset
::

