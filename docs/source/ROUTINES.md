## Standard Routines <a name="standard-routines"></a>

<!-- ### Table of Contents
[Standard Routines](#standard-routines)  
[.prj File](#prj-file)  
[Domain Inputs](#domain-inputs)  
[Material Inputs](#material-inputs)  
[Seismic and EM Inputs](#seismic-and-electromagnetic-inputs)  
[Attenuation Inputs](#attenuation-input)  
[Tensor Components](#tensor-components)  
[Class Definitions](#class-definitions)   -->

The *routines* directory contains the primary functionality for building and running models. *Bash* executables and Python modules are built during the install for many of the scripts found in this directory. Those are: 

- prjbuild 
- prjrun 
- sourcefunction 
- arraybuild 

*SeidarT* models are built from a .png file where each unique color corresponds to a different material. There is built in support for a limited number of materials, which can be referenced in the *materials.py* script, however the motivation for building this toolbox is for studying snow and ice so those have much more robust representation. There is room to expand in the material definitions. Building a project starts with *prjbuild* which creates a comma delimited '.prj' file that contains the domain, material, and model configurations for either a seismic or electromagnetic wave propagation model. 

### .prj File <a name="project-file"></a>

In a .prj file, there are 8 categories which can be identified by the line prefix. These are:

* I - PNG image reference
* D - domain inputs 
* M - material inputs
* S - seismic model inputs 
* A - attenuation inputs (seismic only)
* C - seismic tensor components 
* E - electromagnetic model inputs 
* P - permittivity and conductivity tensors 

### Domain Inputs <a name="domain-inputs"></a>

**User defined inputs and defintions**  
*dim* - This is the dimension of the model which can be 2 or 2.5; default = 2.  
*ny* - For 2.5D models, this number must be >=1. Default is 'n/a'.  
*dx* - The spatial step in the x-direction   
*dy* - The spatial step in the y-direction  
*dz* - The spatial step in the z-direction  
*cpml* - the thickness of the Convolutional Perfectly Matched Layer often referred to as the absorbing boundary. Default is 10.  

When building a project from a .png file, the number of spatial steps in the x- and z-directions is autofilled. For a 2.5D model, **dim** needs to be set to 2.5, and **ny** needs a non-zero value. Their are default values that are emplaced but can be changed accordingly. *cpml* defines the absorbing boundary thickness for all 4 (or 6) borders of the model domain. Less power is absorbed when the model is poorly conditioned such that the source function has a corresponding wavelength that is comparable to or smaller than the spatial step - **dx**, **dy**, or **dz**. One solution is to increase the absorbing boundary thickness or decrease the maximum spatial step. 

### Material Inputs <a name="material-inputs"></a>

There are 9 fields for each material. These are:  
1.  *number* - the number ID that relates to the color. This is used to relate info from inputs to the image.   
2.  *id* - this is the string ID of a material found in the material.py dictionary.  
3.  *R/G/B* - This is a unique RGB color to help the user identify which material belongs here.  
4.  *Temperature* - the temperature of the material in Celsius.  
5.  *Density* - the density of the material in g/cm<sup>3. 
6.  *Porosity* - the percentage of pore space in the material.  
7.  *Water_Content* - the percentage of the pore space that is water. All other pore space is assumed to be air.  
8.  *Anisotropic* - A flag whether the material is anisotropic.  
9.  *ANG_File* - If *Anisotropic* is True then a file containing the Euler angles in radians for z-x-z transformation needs to be supplied. This could be a .txt. Columns 1-3 in the ANG_File must contain the values for the z-x-z rotation.  

### Seismic and Electromagnetic Inputs <a name="model-inputs"></a>

*time_steps* - the number of time steps to run the simulation. 
*x, y, z* - the spatial location of the source in the cartesian reference frame of the model.
*f0* - the center frequency of the source function in Hz.
*theta* - source rotation in the x-z plane in degrees. A value of 90<sup>o</sup> is a vertical source (i.e. F<sub>x</sub> = 0)  
*phi* - source rotation in the x-y plance in degrees. A value of 0<sup>o</sup> has a zero component in the y-direction.

### Attenuation Input <a name="attenuation-inputs"></a>

### Tensor Components <a name="tensor-components"></a>



**dt is not passed to the Fortran FDTD code, and is computed from the smallest spatial step and the maximum propagation velocity. Changing this value will alter the plotting axes so it is best to leave it alone. 

#### Class Definitions <a name="classes"></a>


