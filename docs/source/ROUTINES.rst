

.. _standard-routines: 

Standard Routines 
*****************

.. contents:: 
    :local: 
    :depth: 34
 
.. _project-file: 

Project File
============

The *routines* directory contains the primary functionality for building and running models. *Bash* executables and Python modules are built during the install for many of the scripts found in this directory. Those are: 

- prjbuild 
- prjrun 
- sourcefunction 
- arraybuild 

*SeidarT* models are built from a .png file where each unique color corresponds to a different material. There is built in support for a limited number of materials, which can be referenced in the *materials.py* script, however the motivation for building this toolbox is for studying snow and ice so those have much more robust representation. There is room to expand in the material definitions. 

Building a project starts with *prjbuild* which creates a JSON file that contains the domain, material, and model configurations for either a seismic or electromagnetic wave propagation model. Default values are populated in the fields to prevent any errors during the model run. 

This JSON file contains three main sections:

    - Domain
    - Materials
    - Seismic
    - Electromagnetic

.. _domain-inputs: 

Domain 
------


The Domain object describes the computational domain and image parameters.

.. list-table:: Domain Inputs
    :widths: 25 50 25
    :align: center
    :header-rows: 1 

    * - Key 
      - Description 
      - Value
    * - dim 
      - Dimensionality of the domain 
      - 2 or 2.5.
    * - nx
      - Number of cells in the x direction.
      - integer
    * - ny
      - Number of cells in the y direction.
      - integer (default is 0)
    * - nz
      - Number of cells in the z direction.
      -  integer
    * - dx
      - Cell size in the x direction.
      - float
    * - dy
      - Cell size in the y direction.
      - float (default is 1.0)
    * - dz
      - Cell size in the z direction.
      - float
    * - cpml
      - Thickness (in number of cells) of the CPML.
      - integer
    * - nmats
      - Number of distinct material types.
      - integer (autofilled)
    * - image_file
      - Name of an image file used as a reference or mask.
      - string (autofilled)

.. _material-inputs: 

Materials
---------

The Materials array contains nmats material definitions. Each object has:

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


.. list-table:: Material Inputs 
    :widths: 25 50 25 
    :align: center
    :header-rows: 1 

    * - Key 
      - Description 
      - Value
    * - id
      - Material identifier
      - integer
    * - name
      - Descriptive name of the material (e.g., snow, granite, etc.).
      - string
    * - rgb
      - RGB color code in the form R/G/B.
      - integer/integer/integer
    * - temperature
      - Temperature (in °C or desired unit).
      - float
    * - density
      - Density (in kg/m^3 or desired unit).
      - float
    * - porosity
      - Porosity (in percent, fraction, or desired unit).
      - float
    * - water_content
      - Water content (in percent, fraction, or desired unit).
      - float
    * - is_anisotropic
      - Boolean indicating if the material has anisotropic properties.
      - boolean
    * - euler_angles
      - File path containing the z-x-z angles in radians (null if isotropic).
      - string


Source 
------

Both seismic and electromagnetic sections of the JSON file contain a Source sections. These are identical parameter keys although values/units will differ for some. 


.. list-table:: Source Parameters 
    :widths: 35 50 25 
    :align: center 
    :header-rows: 1 

    * - Key 
      - Description 
      - Value
    * - dt
      - Time interval
      - float (autofilled)
    * - time_steps
      - Total number of time steps of simulation.
      - integer
    * - x, y, z
      - Physical location of the source in the domain (in meters).
      - float
    * - xind, yind, zind
      - Discrete grid indices for the source.
      - integer (autofilled)
    * - source_frequency
      - Dominant frequency of the source wavelet (Hz).
      - float
    * - x-z_rotation, x-y_rotation
      - Rotation angles of the source in degrees around the respective axes.
      - float
    * - amplitude
      - Scalar amplitude of the source wavelet.
      - float
    * - source_type
      - Wavelet type - gaus0, gaus1, gaus2.
      - string


.. note::

    dt is not passed to the Fortran FDTD code, and is computed from the smallest spatial step and the maximum propagation velocity. Changing this value will alter the plotting axes so it is best to leave it alone. 

Seismic
-------

Along with Source, the Seismic object contains two more additional parts:

    Attenuation
    Stiffness_Coefficients

These values can be edited after being calculated internally or input manually. 

Attenuation 
^^^^^^^^^^^

An array of material-specific attenuation parameters. Each object has:

    id: material ID (matching the Materials list).
    gamma_x, gamma_y, gamma_z: attenuation coefficients in the x, y, z directions.
    reference_frequency: the frequency at which attenuation parameters are measured or referenced.

Stiffness_Coefficients 
^^^^^^^^^^^^^^^^^^^^^^

An array of stiffness tensors for seismic wave propagation. Each entry corresponds to a material ID and provides all components of the (possibly anisotropic) stiffness matrix (C_ij). For isotropic materials, many of these may be zero. Each object includes:

    id: Material ID.
    c11, c12, c13, c14, c15, c16, c22, c23, ... etc.: Components of the stiffness matrix.
    rho: Density used in the seismic model (may be repeated from the Materials section if needed).


Electromagnetic
---------------

The Electromagnetic object also contains:

    Permittivity_Coefficients
    Conductivity_Coefficients


Permittivity_Coefficients 
^^^^^^^^^^^^^^^^^^^^^^^^^

An array of tensors for the relative permittivity of each material. Each entry has:

    id: Material ID.
    e11, e12, e13, e22, e23, e33: Tensor components for the permittivity.

Conductivity_Coefficients 
^^^^^^^^^^^^^^^^^^^^^^^^^

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






Class Definitions
=================
