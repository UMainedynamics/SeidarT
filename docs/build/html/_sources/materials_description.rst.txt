Currently available materials
#####################################

*SeidarT* primarily focuses on ice and snow therefore the bulk properties of these materials can be modeled by changing inputs such as pore space, density, liquid water content, and temperature. Permittivity and stiffness tensors for snow are considered isotropic while ice1h can be either isotropic or anisotropic. All other materials are isotropic, and water and air have non-zero stiffness coefficients except for the values :math:`C_{ii}\ \forall i=1,2,3`. 

Isotropic materials beyond ice and snow are also dependent on certain variables. They are computed and bounded by their seismic P and S velocities or their permittivity and conductivity values. We implement an inverse trigonometric function (i.e. :math:`V_{P\text{mod}} = a \tan \rho + V_{P\text{min}}`) to scale the seismic velocities between their upper and lower bounds using overburden pressure as an independent variable prior to computing the Lamé constants. Dynamic stress fields and changes in water pressure are have not been incorporated in the computation of the stiffness coefficients. For permittivity and conductivity coefficients, a porosity correction is applied. 

Any permittivity, conductivity, and stiffness values can be directly edited in the project file prior to modeling. 

The current list of materials are:

* air
* ice1h
* snow
* soil
* water
* oil
* dry_sand
* wet_sand
* granite
* gneiss
* basalt
* limestone
* anhydrite
* coal
* salt

If higher order bulk property modeling is required, users are encouraged to do so using the ThermoElastic and Seismic Analysis (`TESA <https://umaine.edu/mecheng/vel/software/tesa_toolbox/>`_) toolbox. 