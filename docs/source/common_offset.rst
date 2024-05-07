Common Offset
~~~~~~~~~~~~~

We have wrapped other *SeidarT* functions into a class that will simulate a common offset survey. Along with the project file and the receiver file, a source location file needs to be provided. The source and receiver files are formatted the same with the same length. Each line in the source file corresponds to the same line in the receiver file. This allows the user to build a survey over topography, and even though it is called common offset, this can be very easily adapted for cross borehole surveys. 

::
    from seidart.simulations.common_offset import CommonOffset

