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


Top-Level Dictionary
====================

This JSON file contains three main sections:

    Domain
    Materials
    Seismic
    Electromagnetic


Domain 
^^^^^^ 


The Domain object describes the computational domain and image parameters.

.. list-table:: :header-rows: 1 :widths: 25 75

        Key
        Description / Value
        dim
        Dimensionality of the domain (e.g., 2 for 2D).
        nx
        Number of cells in the x direction.
        ny
        Number of cells in the y direction.
        nz
        Number of cells in the z direction.
        dx
        Cell size in the x direction.
        dy
        Cell size in the y direction.
        dz
        Cell size in the z direction.
        cpml
        Thickness (in number of cells) of the Convolutional Perfectly Matched Layer (CPML).
        nmats
        Number of distinct material types.
        image_file
        Name of an image file used as a reference or mask (englacialwater.png).

Materials
^^^^^^^^^

The Materials array contains nmats material definitions. Each object has:

.. list-table:: :header-rows: 1 :widths: 15 85

        Key
        Description / Value
        id
        Material identifier (integer).
        name
        Descriptive name of the material (e.g., snow, granite, etc.).
        rgb
        RGB color code in the form R/G/B.
        temperature
        Temperature (in °C or desired unit).
        density
        Density (in kg/m^3 or desired unit).
        porosity
        Porosity (in percent, fraction, or desired unit).
        water_content
        Water content (in percent, fraction, or desired unit).
        is_anisotropic
        Boolean indicating if the material has anisotropic properties.
        euler_angles
        Euler angles for anisotropic materials (null if isotropic).

Example 
^^^^^^^

Below is a partial list of the materials from the JSON:

    Material 0:
        id = 0
        name = snow
        rgb = 0/255/255
        temperature = -5
        density = 910
        porosity = 2
        water_content = 0
        is_anisotropic = false
        euler_angles = null

    Material 1:
        id = 1
        name = granite
        rgb = 26/26/26
        etc...

(.. and so on for all 7 entries.)
Seismic

The Seismic object contains three parts:

    Source
    Attenuation
    Stiffness_Coefficients

Source 
^^^^^^

.. list-table:: :header-rows: 1 :widths: 20 80

        Key
        Description / Value
        dt
        Time step for the seismic simulation.
        time_steps
        Total number of time steps in the seismic simulation.
        x, y, z
        Physical location of the source in the domain (in user-defined units).
        xind, yind, zind
        Discrete grid indices for the source.
        source_frequency
        Dominant frequency of the source wavelet.
        x-z_rotation, x-y_rotation
        Rotation angles of the source in degrees (or radians) around the respective axes.
        amplitude
        Amplitude of the source wavelet.
        source_type
        Wavelet type (e.g., gaus1, Ricker, etc.).

Attenuation ^^^^^^^^^^^

An array of material-specific attenuation parameters. Each object has:

    id: material ID (matching the Materials list).
    gamma_x, gamma_y, gamma_z: attenuation coefficients in the x, y, z directions.
    reference_frequency: the frequency at which attenuation parameters are measured or referenced.

Stiffness_Coefficients ^^^^^^^^^^^^^^^^^^^^^^

An array of stiffness tensors for seismic wave propagation. Each entry corresponds to a material ID and provides all components of the (possibly anisotropic) stiffness matrix (C_ij). For isotropic materials, many of these may be zero. Each object includes:

    id: Material ID.
    c11, c12, c13, c14, c15, c16, c22, c23, ... etc.: Components of the stiffness matrix.
    rho: Density used in the seismic model (may be repeated from the Materials section if needed).

Electromagnetic

The Electromagnetic object also contains:

    Source
    Permittivity_Coefficients
    Conductivity_Coefficients

Source ^^^^^^

.. list-table:: :header-rows: 1 :widths: 20 80

        Key
        Description / Value
        dt
        Time step for the EM simulation.
        time_steps
        Total number of time steps.
        x, y, z
        Physical location of the EM source in the domain.
        xind, yind, zind
        Discrete grid indices for the EM source.
        source_frequency
        Frequency of the EM wave.
        x-z_rotation, x-y_rotation
        Rotation angles of the source in the x–z or x–y planes.
        amplitude
        Amplitude of the EM source.
        source_type
        Type of the wavelet/pulse (e.g., gaus1).

Permittivity_Coefficients ^^^^^^^^^^^^^^^^^^^^^^^^^

An array of tensors for the relative permittivity of each material. Each entry has:

    id: Material ID.
    e11, e12, e13, e22, e23, e33: Tensor components for the permittivity.

Conductivity_Coefficients ^^^^^^^^^^^^^^^^^^^^^^^^^

An array of tensors for the conductivity of each material. Each entry has:

    id: Material ID.
    s11, s12, s13, s22, s23, s33: Tensor components for the electrical conductivity.

.. ======================================================

### Project File <a name="project-file"></a>

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


