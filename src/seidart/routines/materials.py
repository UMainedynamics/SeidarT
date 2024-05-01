import numpy as np
import pandas as pd
from typing import Union, Tuple
from scipy.interpolate import interp1d

__all__ = [
    'pressure_array',
    'anisotropic_boolean',
    'get_seismic',
    'get_perm',
    'rho_water_correction',
    'isotropic_stiffness_tensor',
    'isotropic_permittivity_tensor',
    'porewater_correction',
    'ice_stiffness',
    'ice_permittivity',
    'snow_permittivity',
    'snow_conductivity',
    'read_ang',
    'rotator_zxz',
    'bond',
    'fujita_complex_permittivity'
]

# Global constants
eps0 = 8.85418782e-12 # used for em only
mu0 = 4.0*np.pi*1.0e-7 # used for em only

"""
Seismic values can be found in: 
      Acoustics of Porous Media (1992), Bourbie, Coussy, and Zinszner
      
Permittivity values can be found in:
        Electrical Properties of Rocks and Minerals
    The following values are provided by:
        Seismic Velocities - https://pangea.stanford.edu/courses/gp262/Notes/8.SeismicVelocity.pdf
        Permeabilitues - http://www.geo.umass.edu/faculty/wclement/dielec.html
        Conductivities - Duba et al. (1977), Duba et al. (1978), Watanabe (1970), 
                    Mohammadi and Mohammadi (2016),
                    https://www.nrcs.usda.gov/INTERNET/FSE_DOCUMENTS/NRCS142P2_053280.PDF,
                    https://en.wikipedia.org/wiki/Electrical_resistivity_and_conductivity
    Values are: 
        Vp_min, Vp_max, Vs_min, Vs_max, Rel_Perm_min, Rel_Perm_max, Conductivity_min, Conductivity_max
    
    Permittivity is given as the relative permittivity and for most cases we 
    will assume that relative permeability is unity; however, if we include
    materials that are high in magnetite, hematite, etc. then we will need to
    accomodate for better permeability estimates.
    
    We are given a range of velocities of different materials found empirically. 
    For isotropic materials we can determine the Lame constants from the equations:
        Vp = sqrt( lambda + 2 mu  / rho ),
        Vs = sqrt( mu / rho ),
        c11 = c22 = c33 = lambda + 2 mu,
        c12 = c13 = c23 = lambda,
        c44 = c55 = c66 = mu
        
    Created by Steven Bernsen
"""

# =============================================================================
#                       Define material dictionaries
# =============================================================================

isotropic_materials = {
    "air":np.array([343, 343, 0.0, 0.0, 1.0, 1.0, 1.0e-16, 1.0e-15]),
    "ice1h":np.array([3400, 3800, 1700, 1900, 3.1, 3.22, 1.0e-7, 1.0e-6]),
    "snow":np.array([100, 2000, 50, 500, 1.0, 70, 1.0e-9, 1.0e-4]),
    "soil":np.array([300, 700, 100, 300, 3.9, 29.4, 1.0e-2, 1.0e-1]), # Permittivity estimates are best constructed with the snow_permittivity() function
    "water":np.array([1450, 1500, 0, 0, 80.36, 80.36, 5.5e-6, 5.0e-2]), # This can change drastically depending on the ions in solution
    "oil":np.array([1200, 1250, 0, 0, 2.07, 2.14, 5.7e-8, 2.1e-7]),
    "dry_sand":np.array([400, 1200, 100, 500, 2.9, 4.7, 1.0e-3, 1.0e-3]), # perm porositiy dependence
    "wet_sand":np.array([1500, 2000, 400, 600, 2.9, 105, 2.5e-4, 1.2e-3]), 
    "granite":np.array([4500, 6000, 2500, 3300, 4.8, 18.9, 4.0e-5, 2.5e-4]),
    "gneiss":np.array([4400, 5200, 2700, 3200, 8.5, 8.5, 2.5e-4, 2.5e-3]),
    "basalt":np.array([5000, 6000, 2800, 3400, 12, 12, 1.0e-6, 1.0e-4]),
    "limestone":np.array([3500, 6000, 2000, 3300, 7.8, 8.5, 2.5e-4, 1.0e-3]),
    "anhydrite":np.array([4000, 5500, 2200, 3100, 5, 11.5, 1.0e-6, 1.0e-5]), # permittivity value from Gypsum
    "coal":np.array([2200, 2700, 1000, 1400, 5.6, 6.3, 1.0e-8, 1.0e-3]), # This has a water dependency
    "salt":np.array([4500, 5500, 2500, 3100, 5.6, 5.6, 1.0e-7, 1.0e2]) # This is dependent on ions and water content
}
# =============================================================================
#                                  Functions
# =============================================================================

# -----------------------------------------------------------------------------
def pressure_array(
        im: Union[list, np.ndarray], 
        temp: Union[list, np.ndarray], 
        rho: Union[list, np.ndarray], 
        dz: Union[list, np.ndarray], 
        porosity: Union[list, np.ndarray] = [0], 
        lwc: Union[list, np.ndarray] = [0]
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes the hydrostatic pressure, temperature, and density at each grid 
    point based on material ID, accounting for porosity and water content.

    :param im: An m-by-n array of integer values representing material IDs.
    :param temp: Temperatures at each grid point.
    :param rho: Densities at each grid point.
    :param dz: Vertical spacing between grid points.
    :param porosity: Porosities at each grid point, default is a list of zeros.
    :param lwc: Liquid water contents at each grid point, default is a list 
        of zeros.
    :type im: Union[list, np.ndarray]
    :type temp: Union[list, np.ndarray]
    :type rho: Union[list, np.ndarray]
    :type dz: Union[list, np.ndarray]
    :type porosity: Union[list, np.ndarray], optional
    :type lwc: Union[list, np.ndarray], optional
    :return: A tuple containing arrays for temperature, density, and hydrostatic 
        pressure at each grid point.
    :rtype: Tuple[np.ndarray, np.ndarray, np.ndarray]
    """


    # First match the size of 
    k = np.unique(im)

    m, n = im.shape
    # allocate the pressure, temperature and density
    pressure = np.zeros([m, n])
    density = np.zeros([m, n])

    if not temp.shape == im.shape:
        temperature = np.zeros([m, n])
    for j in range(0,n):
        for i in range(0, m):
            temperature[i,j] = temp[ im[i,j] ]
            density[i,j],_,_ = porewater_correction(temperature[i,j], rho[ im[i,j] ], 
                porosity[ im[i,j]], lwc[ im[i,j]])
            pressure[i,j] = np.mean(density[0:i,j]) * 9.80665 * i * dz
    
    return(temperature, density, pressure)

# -----------------------------------------------------------------------------
def anisotropic_boolean(
        im: Union[np.ndarray, list], 
        matbool: Union[np.ndarray, list], 
        angvect: Union[np.ndarray, list]
    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Determines if materials identified in an image are anisotropic and provides
    corresponding angular file names if applicable.

    :param im: An array representing material IDs in an image.
    :param matbool: Array indicating whether a material is anisotropic.
    :param angvect: Array of angular file names for anisotropic materials.
    :type im: Union[np.ndarray, list]
    :type matbool: Union[np.ndarray, list]
    :type angvect: Union[np.ndarray, list]
    :return: A tuple of two arrays; the first indicates anisotropy (boolean), 
        the second contains angular file names.
    :rtype: Tuple[np.ndarray, np.ndarray]
    """

    m,n = im.shape
    anisotropic = np.zeros([m,n], dtype = bool)
    afile = np.zeros([m,n], dtype = str)

    for i in range(0, m):
        for j in range(0,n):

            # The booolean in anisotropic could be true, True, TRUE     
            anisotropic[i,j] = (
                matbool[ im[i,j] ] == 'true' or \
                matbool[ im[i,j] ] == 'TRUE' or \
                matbool[ im[i,j] ] == 'True'
            )
            
            if anisotropic[i,j]:
                afile[i,j] = angvect[ im[i,j] ]

    return(anisotropic, afile)

# -----------------------------------------------------------------------------
def get_seismic(
        material_name: Union[list, np.ndarray] = [None], 
        temp: Union[list, np.ndarray] = [None], 
        rho: Union[list, np.ndarray] = [None], 
        porosity: Union[list, np.ndarray] = [0], 
        lwc: Union[list, np.ndarray] = [0], 
        anisotropic: Union[list, np.ndarray] = [False], 
        angfile: Union[list, np.ndarray] = [None]
    ) -> np.ndarray:
    """
    Calculates seismic stiffness coefficients based on material properties,
    accounting for anisotropic conditions where applicable.

    :param material_name: Names of materials.
    :param temp: Temperatures associated with each material.
    :param rho: Densities of materials.
    :param porosity: Porosities of materials.
    :param lwc: Liquid water contents of materials.
    :param anisotropic: Indicates if materials are anisotropic.
    :param angfile: Angular files associated with anisotropic materials.
    :type material_name: Union[list, np.ndarray], optional
    :type temp: Union[list, np.ndarray], optional
    :type rho: Union[list, np.ndarray], optional
    :type porosity: Union[list, np.ndarray], optional
    :type lwc: Union[list, np.ndarray], optional
    :type anisotropic: Union[list, np.ndarray], optional
    :type angfile: Union[list, np.ndarray], optional
    :return: A tensor containing seismic stiffness coefficients.
    :rtype: np.ndarray
    """

    m = len(temp)
    tensor = np.zeros([m, 23])

    # Adjust the stiffness tensor according to the pressure and temperature 
    # conditions for ice
    for ind in range(0, m):
        density,_,_ = porewater_correction(
            temp[ind], 
            rho[ind], 
            porosity[ind], 
            lwc[ind] 
        )

        if anisotropic[ind] and material_name[ind] == 'ice1h':
            euler = read_ang(angfile[ind])
            p = len(euler[:,0])

            cvoigt = np.zeros([6,6])
            creuss = np.zeros([6,6])
            C = np.zeros([6,6])
            
            # Assume a constant pressure of 0.1 MPa (Why? because this is 
            # approximately 1 ATM)
            C = ice_stiffness(temp[ind], 0.1)
            S = np.linalg.inv(C)

            for k in range(0, p):
                R = rotator_zxz(euler[k,:] )
                M = bond(R)
                N = np.linalg.inv(M)
                cvoigt = cvoigt + ( np.matmul( M, np.matmul(C, M.T) ) )
                creuss = creuss + ( np.matmul( N, np.matmul(S, N.T) ) )

            cvoigt = cvoigt/p
            creuss = creuss/p 
            creuss = np.linalg.inv(creuss) 

            # Calculate the hill average 
            C = (cvoigt + creuss)/2
        elif not anisotropic[ind] and material_name[ind] == 'ice1h':
            C = ice_stiffness(temp[ind], 0.1)
        else:
            material_limits = isotropic_materials[ material_name[ind] ]
            C = isotropic_stiffness_tensor(0.1, density, material_limits )

        tensor[ind, :] = (
            ind, 
            C[0,0], C[0,1], C[0,2], C[0,3], C[0,4], C[0,5],
            C[1,1], C[1,2], C[1,3], C[1,4], C[1,5],
            C[2,2], C[2,3], C[2,4], C[2,5],
            C[3,3], C[3,4], C[3,5],
            C[4,4], C[4,5],
            C[5,5], 
            density
        )

    return(tensor)

# -----------------------------------------------------------------------------
def get_perm(
        material,
        modelclass
    ) -> np.ndarray:
    """
    Computes the permittivity and conductivity tensors for materials based on attributes
    contained within material and model class instances.

    :param material: An instance of a class containing attributes for materials, including
                     temperature, density, porosity, water content, and anisotropic properties.
    :type material: Material
    :param modelclass: An instance of a class containing modeling parameters, such as the
                       frequency of interest for electromagnetic modeling.
    :type modelclass: Model
    :return: A tensor array containing permittivity and conductivity values for each material.
    :rtype: np.ndarray
    """

    #!!!!! Use the material class as an input since it will have all of the values that you need instead of inputting all of them 
    material.temp
    
    m = len(material.temp)
    # We will always compute the complex tensor. 
    tensor = np.zeros([m, 13], dtype = complex)
    
    # Adjust the stiffness tensor according to the pressure and temperature 
    # conditions for ice
    for ind in range(0, m):
        
        if material.material[ind] == 'ice1h':
            permittivity = ice_permittivity(
                material.temp[ind],
                material.rho[ind],
                center_frequency = modelclass.f0
            )
            conductivity = snow_conductivity(
                permittivity = permittivity, frequency = modelclass.f0
            )
            
        elif material.material[ind] == 'snow':
            permittivity = snow_permittivity(
                temperature = material.temp[ind],
                lwc = material.lwc[ind], 
                porosity = material.porosity[ind]
            )
            conductivity = snow_conductivity(
                permittivity = permittivity, frequency = modelclass.f0
            )
        else:
            permittivity = np.round(
                isotropic_permittivity_tensor(
                    material.temp[ind], 
                    material.porosity[ind], 
                    material.lwc[ind], 
                    material.material[ind])[0], 
                    3
                )
            conductivity = isotropic_permittivity_tensor(
                material.temp[ind], 
                material.porosity[ind], 
                material.lwc[ind], material.material[ind]
            )[1]
        
        if material.is_anisotropic[ind]:
            euler = read_ang(angfile[ind])
            p = len(euler[:,0])
            
            pvoigt = np.zeros([3,3])
            preuss = np.zeros([3,3])
            permittivity = np.zeros([3,3])
            
            # Assume a constant pressure of 0.1 MPa
            
            S = np.linalg.inv(permittivity)
            
            for k in range(0, p):
                R = rotator_zxz(euler[k,:] )
            
                Ri = np.linalg.inv(R)
                #!!!!! We need to do the same for conductivity.  
                pvoigt = pvoigt + ( np.matmul( R, np.matmul(permittivity, R.T) ) )
                preuss = preuss + ( np.matmul( Ri, np.matmul(S, Ri.T) ) )
            
            pvoigt = pvoigt/p
            preuss = preuss/p 
            preuss = np.linalg.inv(preuss) 
            
            # Calculate the hill average 
            permittivity = (pvoigt + preuss)/2
            
        tensor[ind, :] = np.array(
            [
                ind,
                permittivity[0,0], permittivity[0,1], permittivity[0,2],
                permittivity[1,1], permittivity[1,2],
                permittivity[2,2],
                conductivity[0,0], conductivity[0,1], conductivity[0,2],
                conductivity[1,1], conductivity[1,2],
                conductivity[2,2]
            ] 
        )
    
    return(tensor)

# -----------------------------------------------------------------------------
def rho_water_correction(temperature: float = 0.0) -> float:
    """
    Corrects the density of water based on temperature using the empirical 
    formula derived from the Kell equation.

    :param temperature: The temperature at which to compute the water density 
        correction.
    :type temperature: float
    :return: The corrected water density.
    :rtype: float
    """
    rho_water = (
        999.83952 + \
            16.945176 * temperature - \
                7.9870401e-3 * temperature**2 - \
                    46.170461e-6 * temperature**3 + \
                        105.56302e-9 * temperature**4 - \
                            280.54253e-12 * temperature**5
        )/(1 + 16.897850e-3 * temperature)
    return(rho_water)

# -----------------------------------------------------------------------------
def isotropic_stiffness_tensor(
        pressure: float, 
        density: float, 
        material_limits: np.ndarray
        ) -> np.ndarray:
    """
    Computes the isotropic stiffness tensor for a given material based on 
    pressure, density, and predefined material properties.

    :param pressure: The hydrostatic pressure to which the material is 
        subjected.
    :param density: The density of the material.
    :param material_limits: An array containing the material's velocity limits 
        and other relevant properties.
    :type pressure: float
    :type density: float
    :type material_limits: np.ndarray
    :return: The isotropic stiffness tensor for the material.
    :rtype: np.ndarray
    """

    Vp = material_limits[0:2]
    Vs = material_limits[2:4]
    cp = 2*(Vp[1] - Vp[0])/np.pi 
    cs = 2*(Vs[1] - Vs[0])/np.pi
    
    # Correct for pressure
    pvelocity = cp*np.arctan(pressure ) + Vp[0]
    svelocity = cs*np.arctan(pressure) + Vs[0]
    
    # Compute the lame parameters
    mu = density*(svelocity**2)
    lam = density*(pvelocity**2) - 2*mu

    # Assign the matrix
    C = np.zeros([6,6])
    C[0:3,0:3] = lam
    np.fill_diagonal(C, C.diagonal() + mu)
    C[0,0] = lam + 2*mu
    C[1,1]= C[1,1] + mu
    C[2,2] = C[2,2] + mu

    return(C)

# -----------------------------------------------------------------------------
def isotropic_permittivity_tensor(
        temperature: float, 
        porosity: float, 
        water_content: float, 
        material_name: str
    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the isotropic permittivity tensor for a material based on 
    temperature, porosity, water content, and the material's inherent properties.

    :param temperature: The temperature of the material.
    :param porosity: The porosity of the material.
    :param water_content: The water content in the material.
    :param material_name: The name of the material.
    :type temperature: float
    :type porosity: float
    :type water_content: float
    :type material_name: str
    :return: A tuple containing the permittivity and conductivity tensors.
    :rtype: Tuple[np.ndarray, np.ndarray]
    """

    material_limits = isotropic_materials[ material_name ]
    perm0 = material_limits[4]
    perm1 = material_limits[5]

    cond0 = material_limits[6]
    cond1 = material_limits[7]

    # Calculate the slope         
    if material_name == 'ice1h':
        # We'll assume that the maximum porosity of ice (a.k.a. fresh pow pow)
        # is 85%. The porosity is given as percent [0,100]
        perm0 = 3.1884 + 9.1e-4 * temperature
        perm_coef = (perm1 - perm0)/85
        cond_coef = (cond1 - cond0)/85
        permittivity = np.eye(3,3) * (perm_coef*(porosity) + perm0)
        conductivity = np.eye(3,3) * (cond_coef*(porosity) + cond0)
            
    elif material_name == 'soil' or material_name == 'dry sand':
        # The limit of these two materials is around 55%
        perm_coef = (perm1 - perm0)/55
        cond_coef = (cond1 - cond0)/55
        permittivity = np.eye(3,3) * (perm_coef*(porosity) + perm0)
        conductivity = np.eye(3,3) * (cond_coef*(porosity) + cond0)
    
    elif material_name == 'salt':
        # Permittivity will change a little bit but let's neglect it for now
        permittivity = np.eye(3,3) * perm0
        # Let's use a simple linear trend for a maximum of 20%? water content
        cond_coef = (cond1 - cond0)/20
        conductivity = np.eye(3,3) * (cond_coef*(water_content) +  cond0 )
    
    elif material_name == 'water' or material_name == 'oil':
        # Water and oil do not have a porosity.
        permittivity = np.eye(3,3) * material_limits[4]
        conductivity = np.eye(3,3) * material_limits[6]
    else:
        # For other materials we'll assume the limit is around 3%
        perm_coef = (perm1 - perm0)/3
        cond_coef = (cond1 - cond0)/3
        permittivity = np.eye(3,3) * (perm_coef*(porosity) + perm0)
        conductivity = np.eye(3,3) * (cond_coef*(porosity) + cond0)
    
    return(permittivity, conductivity)

# -----------------------------------------------------------------------------
def porewater_correction(
        temperature: float, 
        density: float, 
        porosity: float, 
        liquid_water_content: float
    ) -> Tuple[float, float, float]:
    """
    Applies corrections to the bulk density of a material based on its porosity 
    and water content, considering temperature adjustments to the densities of 
    air and water.

    :param temperature: The temperature of the material.
    :param density: The initial density of the material.
    :param porosity: The porosity percentage of the material.
    :param liquid_water_content: The percentage of the pore space filled with 
        water.
    :type temperature: float
    :type density: float
    :type porosity: float
    :type liquid_water_content: float
    :return: A tuple containing the corrected density, the density contribution 
        from air, and the density contribution from water.
    :rtype: Tuple[float, float, float]
    """

    rho_air = 0.02897/(8.2057338e-5 * (273 + temperature) )
    # Kell Equation; This doesn't incorporate pressure. That would be a nice
    # addition so that we can mimic super cooled water at depth. 
    rho_water = rho_water_correction(temperature)
    
    # rho_water = -4.6074e-7*temperature**4 + \
    #   1.0326e-4*temperature**3 - 1.0833e-2*temperature**2 + \
    #       9.4207e-2*temperature**1 + 999.998

    # There are certain limits such as phase changes so let's put practical 
    # limits on this
    rho_water = np.max( (rho_water, 950) ) # We can't quite accomodate supercooled water density
    rho_water = np.min( (rho_water, rho_water_correction() )) # beyond the freezing and vaporization temperatures, things get wonky
    
    # the water content is the percent of pores that contain water
    grams_air = (1-liquid_water_content/100)*rho_air
    grams_water = (liquid_water_content/100)*rho_water
    rho_wc = grams_air + grams_water
        
    density = (1-porosity/100)*density + (porosity/100)*rho_wc

    return(density, grams_air, grams_water)

# -----------------------------------------------------------------------------
def ice_stiffness(
        temperature: float = None, 
        pressure: float = 0.0
    ) -> np.ndarray:
    """
    Computes the stiffness tensor for ice under specified temperature and 
    pressure conditions based on empirical relationships.

    :param temperature: The temperature at which to compute the stiffness tensor.
    :param pressure: The pressure at which to compute the stiffness tensor.
    :type temperature: float, optional
    :type pressure: float
    :return: The stiffness tensor for ice.
    :rtype: np.ndarray
    """

    # Allocate space for the stiffness tensor
    C = np.zeros([6,6])
    
    C[0,0] = 136.813 - 0.28940*temperature - 0.00178270*(temperature**2) \
      + 4.6648*pressure - 0.13501*(pressure**2) 
    C[0,1] = 69.4200 - 0.14673*temperature - 0.00090362*(temperature**2) \
      + 5.0743*pressure + .085917*(pressure**2)
    C[0,2] = 56.3410 - 0.11916*temperature - 0.00073120*(temperature**2) \
      + 6.4189*pressure - .52490*(pressure**2)
    C[2,2] = 147.607 - 0.31129*temperature - 0.0018948*(temperature**2) \
      + 4.7546*pressure - .11307*(pressure**2)
    C[3,3] = 29.7260 - 0.062874*temperature - 0.00038956*(temperature**2) \
      + 0.5662*pressure + .036917*(pressure**2)
    
    # Fill in the symmetry
    C[1,1] = C[0,0]
    C[1,0] = C[0,1]
    C[2,0] = C[0,2]
    C[1,2] = C[0,2]
    C[2,1] = C[1,2]
    C[4,4] = C[3,3]
    C[5,5] = (C[0,0] - C[0,1] )/2
    
    stiffness = C*1e8

    return(stiffness)

# -----------------------------------------------------------------------------
def ice_permittivity(
        temperature: float, 
        density: float, 
        center_frequency: float = None,
        method: str = "fujita"
    ) -> np.ndarray:
    """
    Computes the complex permittivity of ice given its temperature, density, and
    the frequency of interest. Supports different methods of calculation.

    :param temperature: Temperature of the ice in degrees Celsius.
    :param density: Density of the ice in kg/m^3.
    :param center_frequency: Frequency at which to compute permittivity, in Hz.
    :param method: The method used for calculating permittivity. Supports 
        "kovacs" and "fujita".
    :type temperature: float
    :type density: float
    :type center_frequency: float, optional
    :type method: str
    :return: The complex permittivity tensor for ice.
    :rtype: np.ndarray
    """
    #Allocate
    P = np.zeros([3,3], dtype = complex)

    # The following is for 2-10 GHz. The details can be found in 
    if method == "kovacs":
        perm = (1 + 0.845 * density)**2
    else: # Fujita et al. (2000)
        perm = 3.1884 + 9.1e-4 * temperature
        dP = 0.0256 + 3.57e-5 * (6.0e-6) * temperature
        complex_perm = fujita_complex_permittivity(
            temperature, center_frequency
        )
        perm = complex(perm, complex_perm)
    
    permittivity = np.eye(3,3) * perm 
    if method == 'fujita':
        permittivity[2,2] = perm + dP 

    return(permittivity)

# -----------------------------------------------------------------------------
def snow_permittivity(
        density: float = 917., 
        temperature: float = 0., 
        lwc: float = 0., 
        porosity: float = 50.,
        method: str = "shivola-tiuri"
    ) -> np.ndarray:
    """
    Calculates the complex permittivity of snow based on its density, temperature, liquid water content (LWC), 
    porosity, and the chosen calculation method.

    :param density: Density of the snow in kg/m^3.
    :param temperature: Temperature of the snow in degrees Celsius.
    :param lwc: Liquid water content of the snow in percentage.
    :param porosity: Porosity of the snow in percentage.
    :param method: The method to be used for calculating permittivity. Defaults to "shivola-tiuri".
    :type density: float
    :type temperature: float
    :type lwc: float
    :type porosity: float
    :type method: str
    :return: The complex permittivity tensor for snow.
    :rtype: np.ndarray
    """

    # Temperature equations
    # jones (2005), liebe et al. (1991)
    # Density equations 
    # shivola and tiuri (1986), wise
    
    rho_d,grams_air,grams_water = porewater_correction(
        temperature, density, porosity, lwc
    )
    
    # LWC is in kg/m3 but we need it in g/cm3
    lwc = grams_water / 1000
    rho_d = rho_d / 1000
    # Put temperature in terms of kelvin
    T = temperature + 273.15

    if method == "shivola-tiuri":
        perm = 8.8*lwc + 70.4*(lwc**2) + 1 + 1.17*rho_d + 0.7*(rho_d**2)
    elif method == "wise":
        perm = 1 + 1.202*rho_d + 0.983*(rho_d**2) + 21.3*lwc
    elif method == "jones":
        perm = 78.51 * (
            1 - 4.579 * 1e-3 * (T - 298) + \
                1.19 * 1e-5 * (T - 298)**2 - \
                    2.8*1e-8 * (T - 298)**3
        )
    else: # Liebe et al.
        perm = 77.66 - 103.3 * (1 - (300/(T)) )
    
    # Compute the complex part
    complex_permittivity = 0.8*lwc + 0.72*(lwc**2)
    permittivity = np.eye(3,3) * complex(perm, complex_permittivity)

    return(permittivity)

# -----------------------------------------------------------------------------
def water_permittivity(temperature):
    pass

# -----------------------------------------------------------------------------
def snow_conductivity(
        lwc: float = None, 
        permittivity: np.ndarray = None, 
        frequency: float = None
    ) -> np.ndarray:
    """
    Computes the electrical conductivity of snow given its liquid water content 
    (LWC), permittivity, and frequency of interest.

    :param lwc: Liquid water content of the snow, optionally used if 
        permittivity is not provided.
    :param permittivity: The complex permittivity of snow, if available.
    :param frequency: Frequency at which conductivity is to be calculated, 
        in Hz.
    :type lwc: float, optional
    :type permittivity: np.ndarray, optional
    :type frequency: float, optional
    :return: The conductivity tensor for snow.
    :rtype: np.ndarray
    """

    if np.iscomplexobj(permittivity):
        sigma = permittivity.imag * frequency * eps0
    else:
        # granlund et al. (2010)
        _,_,grams_water = porewater_correction(
            temperature, density, porosity, lwc
        )
        
        # LWC is in kg/m3 but we need it in g/cm3
        lwc = grams_water / 1000
        # granlund provides an equation for micro-S/m so we multiply by 1e-4 to
        # put it in terms of S/m
        sigma = (20 + 3e3 * lwc) * 1e-4 
    
    conductivity = np.eye(3,3) * sigma 
    return(conductivity)

# -----------------------------------------------------------------------------
def read_ang(filepath: str, delimiter = " ") -> np.ndarray:
    """
    Reads Euler angles from a .ang file, typically associated with 
    EBSD (Electron Backscatter Diffraction) data.

    :param filepath: The path to the .ang file.
    :type filepath: str
    :param delimiter: The delimiter for the input file
    :type delimiter: str
    :return: An array of Euler angles extracted from the file.
    :rtype: np.ndarray

    Note:
        The .ang file is expected to contain columns for Euler angles in 
        radians, following Bunge's notation (z-x-z rotation), among other data 
        related to EBSD measurements.
    """

    
    # Load the file in as a data frame
    # if delimiter = ",":
    #     euler = pd.read_csv(filepath).to_numpy()
    # else:
    #     euler = np.genfromtxt(filepath, delimiter = delimiter)
    euler = pd.read_table(filepath, delimiter = delimiter).to_numpy()
    
    # take only the euler angles...for now
    if euler.shape[0] > 3 :
        euler = euler[:,0:3]

    # Unfortunately, the space delimiters are inconsistent :(
    # We know there are 10 columns and the rows shouldn't contain all NA's
    m, n = np.shape(euler)

    # reshape to M-by-1 vector
    euler = euler.reshape(m*n,1)

    # remvoe 'nan'
    euler = euler[~np.isnan(euler)]

    # reshape back to array
    euler = euler.reshape(m, int( len(euler)/m ) )

    # save ferris
    return(euler)

# -----------------------------------------------------------------------------
def rotator_zxz(eul: np.ndarray) -> np.ndarray:
    """
    Generates a rotation matrix from Euler angles using the z-x-z rotation 
    convention.

    :param eul: An array containing the three Euler angles.
    :type eul: np.ndarray
    :return: The 3x3 rotation matrix derived from the Euler angles.
    :rtype: np.ndarray
    """

    # From the 3 euler angles for the zxz rotation, compute the rotation matrix
    R = np.zeros([3,3])
    D = np.zeros([3,3])
    C = np.zeros([3,3])
    B = np.zeros([3,3])

    D[0,:] = [ np.cos( eul[0] ), -np.sin( eul[0] ), 0.0 ]
    D[1,:] = [ np.sin( eul[0] ), np.cos( eul[0] ), 0.0 ]
    D[2,:] = [ 0.0, 0.0, 1.0 ]

    C[0,:] = [ 1.0, 0.0, 0.0 ]
    C[1,:] = [ 0.0, np.cos( eul[1] ), -np.sin( eul[1] ) ]
    C[2,:] = [ 0.0, np.sin( eul[1] ), np.cos( eul[1] ) ]

    B[0,:] = [ np.cos( eul[2] ), -np.sin( eul[2] ), 0.0 ] 
    B[1,:] = [ np.sin( eul[2] ), np.cos( eul[2] ), 0.0 ]
    B[2,:] = [ 0.0, 0.0, 1.0 ]

    R = np.matmul(D, C)
    R = np.matmul(R, B)

    return(R)

# -----------------------------------------------------------------------------
def bond(R: np.ndarray) -> np.ndarray:
    """
    Calculates the 6x6 Bond transformation matrix from a 3x3 rotation matrix, 
    useful for transforming stiffness or compliance matrices in crystallography 
    and materials science.

    :param R: The 3x3 rotation matrix.
    :type R: np.ndarray
    :return: The 6x6 Bond transformation matrix.
    :rtype: np.ndarray
    """

    # From the euler rotation matrix, compute the 6-by-6 bond matrix
    M = np.zeros([6,6])
    M[0,:] = [ 
        R[0,0]**2, R[0,1]**2, R[0,2]**2, 
        2*R[0,1]*R[0,2], 2*R[0,2]*R[0,0], 2*R[0,0]*R[0,1] 
    ]
    M[1,:] = [ 
        R[1,0]**2, R[1,1]**2, R[1,2]**2, 
        2*R[1,1]*R[1,2], 2*R[1,2]*R[1,0], 2*R[1,0] * R[1,1] 
    ]
    M[2,:] = [ 
        R[2,0]**2, R[2,1]**2, R[2,2]**2, 
        2*R[2,1]*R[2,2], 2*R[2,2]*R[2,0], 2*R[2,0] * R[2,1] 
    ]
    M[3,:] = [ 
        R[1,0]* R[2,0], R[1,1] * R[2,1], 
        R[1,2] * R[2,2], R[1,1] * R[2,2] + R[1,2]*R[2,1], 
        R[1,0]*R[2,2] + R[1,2]*R[2,0], R[1,1]*R[2,0] + R[1,0]*R[2,1] 
    ]
    M[4,:] = [ 
        R[2,0]* R[0,0], R[2,1] * R[0,1], 
        R[2,2] * R[0,2], R[0,1] * R[2,2] + R[0,2]*R[2,1], 
        R[0,2]*R[2,0] + R[0,0]*R[2,2], R[0,0]*R[2,1] + R[0,1]*R[2,0] 
    ]
    M[5,:] = [ 
        R[0,0]* R[1,0], R[0,1] * R[1,1], 
        R[0,2] * R[1,2], R[0,1] * R[1,2] + R[0,2]*R[1,1], 
        R[0,2]*R[1,0] + R[0,0]*R[1,2], R[0,0]*R[1,1] + R[0,1]*R[1,0] 
    ]

    return(M)
    
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# The following is for the complex permittivity calculations that were defined
# by Fujita et al. 
T = np.array(
    [190, 200, 220, 240, 248, 253, 258, 263, 265]
)
A = np.array(
    [0.005, 0.010, 0.031, .268, .635, 1.059, 1.728, 2.769, 3.326]
)*10.e-4
B = np.array(
    [1.537, 1.747, 2.469, 3.495, 4.006, 4.380, 4.696, 5.277, 5.646]
)*10.e-5
C = np.array(
    [1.175, 1.168, 1.129, 1.088, 1.073, 1.062, 1.056, 1.038, 1.024]
)
# Interpolation functions for A, B, and C
A_interp = interp1d(T, A, kind='cubic', fill_value='extrapolate')
B_interp = interp1d(T, B, kind='cubic', fill_value='extrapolate')
C_interp = interp1d(T, C, kind='cubic', fill_value='extrapolate')

def fujita_complex_permittivity(temperature: float, frequency: float) -> float:
    """
    Calculates the complex permittivity of ice using Fujita's method, based on 
    the provided temperature and frequency.

    :param temperature: The temperature of ice in degrees Celsius.
    :param frequency: The frequency at which to calculate permittivity, in Hz.
    :type temperature: float
    :type frequency: float
    :return: The complex permittivity value.
    :rtype: float
    """
    # frequency = 1 is equivalent to 1 GHz or 1e9 Hz. The input is in Hz.
    temperature = temperature + 273 # Convert to Kelvins
    frequency = frequency / 1e9
    A_val = A_interp(temperature)
    B_val = B_interp(temperature)
    C_val = C_interp(temperature)
    epsilon_val = A_val/frequency + B_val*(frequency**C_val)
    return epsilon_val
    