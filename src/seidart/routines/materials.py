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
    'astiffness2moduli',
    'bulk_modulus_water',
    'bulk_modulus_air',
    'moduli_basalt',
    'moduli_granite',
    'moduli_gneiss',
    'moduli_sand_silt_clay',
    'snow_stiffness',
    'ice_stiffness',
    'ice_permittivity',
    'ice_density',
    'tensor2velocities',
    'snow_permittivity',
    'snow_conductivity',
    'sand_silt_clay_permittivity_conductivity',
    'read_ang',
    'rotator_zxz',
    'bond',
    'highfreq_stability_conditions',
    'fujita_complex_permittivity'
]

# Global constants
eps0 = 8.85418782e-12 # used for em only
mu0 = 4.0*np.pi*1.0e-7 # used for em only
c_light = 299792458
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
# When adding new materials use:
#   "name":np.array([Vp_min, Vp_max, Vs_min, Vs_max, Perm_min, Perm_max, Cond_min, Cond_max])

isotropic_materials = {
    "air":np.array([343, 343, 0.0, 0.0, 1.0, 1.0, 1.0e-16, 1.0e-15]),
    "ice1h":np.array([3400, 3800, 1700, 1900, 3.1, 3.22, 1.0e-7, 1.0e-6]),
    "snow":np.array([100, 2000, 50, 500, 1.0, 70, 1.0e-9, 1.0e-4]),
    "soil":np.array([300, 700, 100, 300, 3.9, 29.4, 1.0e-2, 1.0e-1]), # Permittivity estimates are best constructed with the snow_permittivity() function
    "water":np.array([1450, 1500, 0, 0, 80.36, 80.36, 5.5e-6, 5.0e-2]), # This can change drastically depending on the ions in solution
    "oil":np.array([1200, 1250, 0, 0, 2.07, 2.14, 5.7e-8, 2.1e-7]),
    "sand":np.array([400, 1200, 100, 500, 2.9, 4.7, 1.0e-3, 1.0e-3]), # perm porositiy dependence
    "silt":np.array([1500, 2000, 400, 600, 2.9, 105, 2.5e-4, 1.2e-3]), 
    "clay":np.array([1500, 2000, 400, 600, 2.9, 105, 2.5e-4, 1.2e-3]), 
    "granite":np.array([4500, 6000, 2500, 3300, 4.8, 18.9, 4.0e-5, 2.5e-4]),
    "gneiss":np.array([4400, 5200, 2700, 3200, 8.5, 8.5, 2.5e-4, 2.5e-3]),
    "basalt":np.array([5000, 6000, 2800, 3400, 12, 12, 1.0e-6, 1.0e-4]),
    "limestone":np.array([3500, 6000, 2000, 3300, 7.8, 8.5, 2.5e-4, 1.0e-3]),
    "anhydrite":np.array([4000, 5500, 2200, 3100, 5, 11.5, 1.0e-6, 1.0e-5]), # permittivity value from Gypsum
    "coal":np.array([2200, 2700, 1000, 1400, 5.6, 6.3, 1.0e-8, 1.0e-3]), # This has a water dependency
    "salt":np.array([4500, 5500, 2500, 3100, 5.6, 5.6, 1.0e-7, 1.0e2]) # This is dependent on ions and water content
}

# ==============================================================================
#                                   General
# ==============================================================================
# ------------------------------------------------------------------------------
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
    point based on material ID, accounting for porosity and water content. The
    output units of pressure are in Pascals. Gammon uses kbar to compute the 
    ice stiffness tensor. 

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
            density[i,j],_,_ = porewater_correction(
                temperature[i,j], 
                rho[ im[i,j] ], 
                porosity[ im[i,j]], 
                lwc[ im[i,j]]
            )
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
def vrh(C, S, eulerangles):
    m, n = eulerangles.shape
    for k in range(m):
        R = rotator_zxz(eulerangles[k,:] )
        M = bond(R)
        N = np.linalg.inv(M)
        cvoigt = cvoigt + ( np.matmul( M, np.matmul(C, M.T) ) )
        creuss = creuss + ( np.matmul( N, np.matmul(S, N.T) ) )
    
    # Normalize the Voigt and Reuss estimates
    cvoigt /= m
    creuss /= m 
    creuss = np.linalg.inv(creuss) 

    # Calculate the hill average 
    hill = (cvoigt + creuss)/2
    return hill

# -----------------------------------------------------------------------------
def read_ang(filepath: str) -> np.ndarray:
    """
    Reads Euler angles from a space delimited .ang file, typically associated 
    with EBSD (Electron Backscatter Diffraction) data. Synthetic data can be 
    built with the fabricsynth module.

    :param filepath: The path to the .ang file.
    :type filepath: str
    :return: An array of Euler angles extracted from the file.
    :rtype: np.ndarray

    Note:
        The .ang file is expected to contain columns for Euler angles in 
        radians, following Bunge's notation (z-x-z rotation), among other data 
        related to EBSD measurements.
    """

    # Load the file in as a data frame
    euler = pd.read_table(filepath, delimiter = " ").to_numpy()
    
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

    R = np.matmul(D, np.matmul(C, B) )


    return(R)

# -----------------------------------------------------------------------------
def bond(R: np.ndarray) -> np.ndarray:
    """
    Calculates the 6x6 Bond transformation matrix from a 3x3 rotation matrix, 
    useful for transforming stiffness or compliance matrices in crystallography 
    and materials science.

    :param R: The 3x3 rotation matrix.
    :type R: np.ndarray
    :return M: The 6x6 Bond transformation matrix.
    :rtype M: np.ndarray
    """

    # Initialize the 6x6 Bond matrix
    M = np.zeros([6,6])

    # Fill the Bond matrix according to the transformation rules
    M[0,:] = [ 
        R[0,0]**2, R[0,1]**2, R[0,2]**2, 
        2*R[0,1]*R[0,2], 2*R[0,2]*R[0,0], 2*R[0,0]*R[0,1] 
    ]
    M[1,:] = [ 
        R[1,0]**2, R[1,1]**2, R[1,2]**2, 
        2*R[1,1]*R[1,2], 2*R[1,2]*R[1,0], 2*R[1,0]*R[1,1] 
    ]
    M[2,:] = [ 
        R[2,0]**2, R[2,1]**2, R[2,2]**2, 
        2*R[2,1]*R[2,2], 2*R[2,2]*R[2,0], 2*R[2,0]*R[2,1] 
    ]
    M[3,:] = [ 
        R[1,0]*R[2,0], R[1,1]*R[2,1], R[1,2]*R[2,2], 
        R[1,1]*R[2,2] + R[1,2]*R[2,1], R[1,2]*R[2,0] + R[1,0]*R[2,2], R[1,0]*R[2,1] + R[1,1]*R[2,0] 
    ]
    M[4,:] = [ 
        R[2,0]*R[0,0], R[2,1]*R[0,1], R[2,2]*R[0,2], 
        R[2,1]*R[0,2] + R[2,2]*R[0,1], R[2,2]*R[0,0] + R[2,0]*R[0,2], R[2,0]*R[0,1] + R[2,1]*R[0,0] 
    ]
    M[5,:] = [ 
        R[0,0]*R[1,0], R[0,1]*R[1,1], R[0,2]*R[1,2], 
        R[0,1]*R[1,2] + R[0,2]*R[1,1], R[0,2]*R[1,0] + R[0,0]*R[1,2], R[0,0]*R[1,1] + R[0,1]*R[1,0] 
    ]

    return M

# ------------------------------------------------------------------------------
def highfreq_stability_conditions(T):
    """
    This is a check to make sure that no numerical instabilities or unreal wave 
    modes are generated. See Bécache (2003) and Komatitsch (2007) for more 
    details.
    
    :param T: This is either a 2nd order (3-by-3) or 4th order (6-by-6) tensor 
    :type T: np.ndarray
    :return: The stability condition values 
    :rtype: tuple
    """
    cond1 = (( T[0,1] + T[2,2] )**2 - T[0,0] * (T[1,1] - T[2,2] )) * \
        (( T[0,1]+T[2,2] )**2 + T[2,2]*(T[1,1] - T[2,2] ))
    cond2 = ( T[0,1] + 2 * T[2,2] )**2 - T[0,0]*T[1,1]
    cond3 = ( T[0,1] + T[2,2] )**2 - T[0,0]*T[2,2] - T[2,2]**2

    if cond1 > 0:
        print('Failed to pass high frequency stability condition 1.')
    if cond2 > 0:
        print('Failed to pass high frequency stability condition 2.')
    if cond3 > 0:
        print('Failed to pass high frequency stability condition 3.')

    return cond1, cond2, cond3 

# -----------------------------------------------------------------------------
def tensor2velocities(T: np.ndarray, rho: float = 910, seismic: bool = True):
    if T.shape == (6,):
        temp = np.zeros([3,3])
        temp[0,:] = T[0:3] 
        temp[1,1:] = T[3:5]
        temp[2,2] = T[5]
        T = np.triu(temp).T + np.triu(temp, 1)
    if T.shape == (21,):
        temp = np.zeros([6,6])
        temp[0,:] = T[0:6] 
        temp[1,1:] = T[6:11]
        temp[2,2:] = T[11:15]
        temp[3,3:] = T[15:18]
        temp[4,4:] = T[18:20]
        temp[5,5] = T[20] 
        T = np.triu(temp).T + np.triu(temp, 1)     
    if seismic:
        vpv = np.sqrt( T[2,2]/rho )
        vph = np.sqrt( T[0,0]/rho )
        vsv = np.sqrt( T[3,3]/rho )
        vsh = np.sqrt( (T[0,0] - T[0,1])/(2*rho) )
        return vpv, vph, vsv, vsh
    else:
        vx = c_light * np.sqrt(1/T[0,0])
        vy = c_light * np.sqrt(1/T[1,1])
        vz = c_light * np.sqrt(1/T[2,2])
        return vx, vy, vz 

# ==============================================================================
#                                   Seismic
# ==============================================================================
# -----------------------------------------------------------------------------
def astiffness2moduli(C):
    """
    Computes the isotropic bulk and shear modulus from an anisotropic stiffness tensor using Hill's average.
    
    Parameters:
        C (numpy array): 6x6 stiffness tensor in Voigt notation (Pa)
    
    Returns:
        (K_iso, mu_iso): Tuple of isotropic bulk and shear moduli (Pa)
    """
    # Voigt Average
    K_V = (1/9) * (C[0,0] + C[1,1] + C[2,2]) + (2/9) * (C[0,1] + C[0,2] + C[1,2])
    mu_V = (1/15) * (C[0,0] + C[1,1] + C[2,2] - C[0,1] - C[0,2] - C[1,2]) + (1/5) * (C[3,3] + C[4,4] + C[5,5])

    # Compliance tensor (inverse of stiffness tensor)
    S = np.linalg.inv(C)

    # Reuss Average
    K_R = 1 / ((1/9) * (S[0,0] + S[1,1] + S[2,2]) + (2/9) * (S[0,1] + S[0,2] + S[1,2]))
    mu_R = 1 / (5 * ((4/15) * (S[0,0] + S[1,1] + S[2,2] - S[0,1] - S[0,2] - S[1,2]) + (3/5) * (S[3,3] + S[4,4] + S[5,5])))

    # Hill Average
    K_iso = (K_V + K_R) / 2
    mu_iso = (mu_V + mu_R) / 2
    
    C_isotropic = np.zeros([6,6])
    lambda_ = K_iso - 2 * mu_iso / 3
    # Fill in the bulk/shear moduli components
    C_isotropic[:3, :3] = lambda_  # λ terms
    np.fill_diagonal(C_isotropic[:3, :3], K_iso + 4 * mu_iso / 3)  # K + 4/3 μ on the diagonal
    np.fill_diagonal(C_isotropic[3:, 3:], mu_iso)  # Shear moduli (μ) for Voigt indices 4,5,6
    return K_iso, mu_iso, C_isotropic
    
# -----------------------------------------------------------------------------
def bulk_modulus_water(T):
    """
    Returns the bulk modulus of water (in GPa) as a function of temperature (°C).
    Includes the supercooled regime down to -30°C.
    
    See the journal article:
    
    Compressibility of water as a function of temperature and pressure
    by
    Fine and Millero 1973
    """
    if T >= 0:
        # Liquid water (0°C to 100°C)
        K = 2.15 - 0.0164 * T + 4.72e-5 * T**2 - 6.5e-7 * T**3
    elif T >= -30:
        # Supercooled water (-30°C to 0°C)
        K = 2.25 - 0.016 * T - 6.9e-5 * T**2
    else:
        raise ValueError("Temperature is too low for known water bulk modulus data.")
    
    K *= 1e9
    C = np.zeros([6,6])
    C[0:3,0:3] = K
    
    return K, C

# -----------------------------------------------------------------------------
def bulk_modulus_air(T, P=101325):
    """
    Computes the bulk modulus of air (Pa) as a function of temperature (°C) and pressure (Pa),
    while accounting for variable air density.

    :param T_C: Temperature in Celsius
    :param P: Pressure in Pascals (default is atmospheric pressure, 101325 Pa)
    :return: Bulk modulus in Pascals (Pa)
    """
    gamma = 1.4  # Adiabatic index for dry air
    R = 8.314  # Universal gas constant (J/(mol·K))
    M = 0.029  # Molar mass of dry air (kg/mol)
    
    T = T + 273.15  # Convert temperature to Kelvin
    rho = (P * M) / (R * T)  # Compute air density (kg/m³)
    
    K = gamma * rho * (R / M) * T  # Compute bulk modulus
    
    C = np.zeros([6,6])
    C[0:3, 0:3] = K
    
    return K, C, rho

# -----------------------------------------------------------------------------
def moduli_basalt(T, P):
    # Pressure needs to be in MPa
    K = (50e9 - 0.02e9 * T + 0.1 * P)  # Pa
    G = (32e9 - 0.015e9 * T + 0.08 * P)  # Pa
    
    return K, G

# -----------------------------------------------------------------------------
def moduli_granite(T, P):
    K = (45e9 - 0.025e9 * T + 0.12 * P)  # Pa
    G = (27e9 - 0.02e9 * T + 0.1 * P)  # Pa

# -----------------------------------------------------------------------------
def moduli_gneiss(T, P):
    K = (50e9 - 0.03e9 * T + 0.15 * P)  # Pa
    G = (25e9 - 0.02e9 * T + 0.1 * P)  # Pa

# -----------------------------------------------------------------------------
def moduli_sand_silt_clay(
        P, T, lwc, porosity = 0.0,
        composition = {"sand": 1.0, "silt": 0.0, "clay": 0.0}
    ):
    """
    Computes the bulk modulus (K) and shear modulus (G) for a sand-silt-clay mixture 
    based on pressure, temperature, and water content using Gassmann's equation.

    :param P: Overburden pressure in Pascals (Pa).
    :param T: Temperature in degrees Celsius.
    :param lwc: Water content fraction (0 to 1).
    :param composition: Dictionary with fractions of sand, silt, and clay (must sum to 1).
    :return: Bulk modulus (Pa), Shear modulus (Pa)
    """
    
    # Porosity and liquid water content are input as percentage
    porosity = porosity/100 
    lwc = lwc/100 
    
    
    # Check that composition sums to 1
    if not np.isclose(sum(composition.values()), 1.0):
        raise ValueError("Composition fractions must sum to 1.")
    
    # Material properties for sand, silt, and clay (from experimental data)
     # Material properties for sand, silt, and clay (from experimental data)
    properties = {
        "sand": {"K_dry": 30e9, "G_dry": 25e9, "K_mineral": 36e9, "G_mineral": 30e9, "rho_solid": 2650},
        "silt": {"K_dry": 20e9, "G_dry": 15e9, "K_mineral": 28e9, "G_mineral": 22e9, "rho_solid": 2700},
        "clay": {"K_dry": 10e9, "G_dry": 5e9, "K_mineral": 20e9, "G_mineral": 15e9, "rho_solid": 2750},
    }
    
    # Bulk modulus of water at given temperature (Fine & Millero, 1973)
    K_water, _ = bulk_modulus_water(T)
    
    # Compute weighted dry bulk modulus and shear modulus
    K_dry_mix = sum(composition[mat] * properties[mat]["K_dry"] for mat in composition)
    G_dry_mix = sum(composition[mat] * properties[mat]["G_dry"] for mat in composition)
    K_mineral_mix = sum(composition[mat] * properties[mat]["K_mineral"] for mat in composition)
    G_mineral_mix = sum(composition[mat] * properties[mat]["G_mineral"] for mat in composition)
    # porosity = sum(composition[mat] * properties[mat]["porosity"] for mat in composition)
    
    # Gassmann’s equation for saturated bulk modulus
    K = K_dry_mix + ((1 - K_dry_mix / K_mineral_mix)**2) / (
        porosity / K_water + (1 - porosity) / K_mineral_mix - K_dry_mix / K_mineral_mix**2)
    
    # Hardin-Drnevich empirical model for shear modulus (adjusted for water content)
    void_ratio = porosity / (1 - porosity)  # Approximate e from porosity
    OCR = 1.0  # Assume normally consolidated for now
    plasticity_index = 0.3  # Representative value for mixed soils
    
    G_max = 14760 * ((2.973 - void_ratio)**2 / (1 + void_ratio)) * (OCR)**plasticity_index * (P)**0.5
    
    # Adjust shear modulus based on water content (more water reduces G)
    G = G_max * (1 - 0.5 * lwc) 
    
    # We need the effective density too 
    # Compute dry density (assuming fully dry state)
    rho_solid_mix = sum(composition[mat] * properties[mat]["rho_solid"] for mat in composition)
    rho_dry = rho_solid_mix * (1 - porosity)
    
    rho_water = rho_water_correction(T) 
    _, _, rho_air = bulk_modulus_air(T) 
    rho_eff = (1 - porosity) * rho_dry + porosity * lwc * rho_water + porosity * (1 - lwc) * rho_air
    
    return K, G, rho_eff

def gassmann_bulk_modulus(K_dry, K_mineral, K_fluid, porosity):
    """
    Computes the bulk modulus of a fluid-saturated porous material using Gassmann's equation.

    :param K_dry: Bulk modulus of the dry matrix (Pa)
    :param K_mineral: Bulk modulus of the solid grains (Pa)
    :param K_fluid: Bulk modulus of the saturating fluid (Pa)
    :param porosity: Porosity (fraction, 0 to 1)
    :return: Saturated bulk modulus (Pa)
    """
    numerator = (1 - K_dry / K_mineral) ** 2
    denominator = porosity / K_fluid + (1 - porosity) / K_mineral - (K_dry / K_mineral**2)

    return K_dry + (numerator / denominator)

def compute_volume_fractions(porosity, lwc):
    """
    Computes volume fractions of ice, air, and water based on porosity and liquid water content.
    
    :param porosity: Porosity fraction (0 to 1)
    :param lwc: Liquid water content fraction (0 to 1)
    :return: (V_ice, V_air, V_water)
    """
    V_air = porosity * (1 - lwc)  # Air fraction
    V_water = porosity * lwc      # Water fraction
    V_ice = 1 - porosity          # Ice fraction
    
    return V_ice, V_air, V_water

# -----------------------------------------------------------------------------
def eshelby_tensor(C_matrix, aspect_ratio=1.0):
    """
    Computes the Eshelby tensor for a spherical pore in an isotropic matrix.
    
    :param C_matrix: 6x6 stiffness tensor of the matrix material (ice).
    :param aspect_ratio: Shape of the inclusion (1.0 for spheres, <1 for oblate).
    :return: 6x6 Eshelby tensor.
    """
    S = np.zeros((6,6))  # Placeholder for Eshelby tensor
    # Simplified Eshelby tensor for spherical pores
    S[:3, :3] = np.eye(3) * (3 * aspect_ratio) / (4 * (1 + aspect_ratio))
    S[3:, 3:] = np.eye(3) * (aspect_ratio) / (2 * (1 + aspect_ratio))
    return S

# -----------------------------------------------------------------------------
def self_consistent_approximation(C, porosity, tol=1e-6, max_iter=100):
    """
    Computes the effective stiffness tensor for highly porous snow using SCA.

    :param C: 6x6 stiffness tensor.
    :param porosity: Volume fraction of air (0 to 1).
    :return: 6x6 effective stiffness tensor of porous media.
    """

    C_eff = C.copy()  # Initial guess
    C_air = np.zeros((6,6))  # Air has negligible stiffness
    S = eshelby_tensor(C)  # Compute Eshelby tensor for spherical pores

    for _ in range(max_iter):
        # Compute the interaction tensor
        P = np.linalg.inv(np.eye(6) + np.matmul(S, (C_air - C_eff)))

        # Update the effective stiffness tensor
        C_new = C_eff + porosity * np.matmul((C_air - C_eff), P)

        # Convergence check
        if np.linalg.norm(C_new - C_eff) < tol:
            break

        C_eff = C_new

    return C_eff

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
        subjected. Units are Pa.
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
def get_seismic(
        self,
        material
    ):
    """
    Calculates seismic stiffness coefficients based on material properties,
    accounting for anisotropic conditions where applicable.

    :param material: The class object that contains material parameters.
    :type material: Material 
    """

    m = len(material.temp)
    stiffness_coefficients = np.zeros([m, 22])

    # Adjust the stiffness tensor according to the pressure and temperature 
    # conditions for ice
    for ind in range(0, m):
        density,_,_ = porewater_correction(
            material.temp[ind], 
            material.rho[ind], 
            material.porosity[ind], 
            material.lwc[ind] 
        )

        if material.material[ind] == 'ice1h':
            print('Computing the stiffness for anisotropic material ice1h.')
            euler = read_ang(material.angfile[ind])
            p = len(euler[:,0])

            cvoigt = np.zeros([6,6])
            creuss = np.zeros([6,6])
            C = np.zeros([6,6])
            
            # Assume a constant pressure of 0.1 MPa (Why? because this is 
            # approximately 1 ATM)
            pressure = 0.1 * 1e-1 # in kbar
            Cice = ice_stiffness(material.temp[ind], pressure)
            Sice = np.linalg.inv(Cice)
            if material.is_anisotropic[ind]:
                C = vrh(Cice, Sice, euler)
            else:
                print('Computing the homogeneous stiffness coefficients for ice1h.')
                _,_, C = astiffness2moduli(Cice)
            # The stiffness matrix needs to be positive definite which means all
            # of the eigenvalues must be positive
            eigenvalues, eigenvectors = np.linalg.eigh(C)
            if np.sum(eigenvalues <= 0) > 0:
                print('Stiffness tensor is not positive definite.')
            else:
                print('Stiffness tensor is positive definite')            
        elif material.material[ind] == 'snow':
            if material.is_anisotropic[ind]:
                C_ice = ice_stiffness(material.temp[ind], pressure)
            else:
                pass
                #k_snow, mu_
        else:
            material_limits = isotropic_materials[ material.material[ind] ]
            C = isotropic_stiffness_tensor(0.1, density, material_limits )
        
        # The high frequency stability conditions will never pass for anisotropic ice
        cond1, cond2, cond3 = highfreq_stability_conditions(C)
        
        # Assign the values        
        C = np.round(C)
        self.stiffness_coefficients.loc[ind] = np.array([
            C[0,0], C[0,1], C[0,2], C[0,3], C[0,4], C[0,5],
                    C[1,1], C[1,2], C[1,3], C[1,4], C[1,5],
                            C[2,2], C[2,3], C[2,4], C[2,5],
                                    C[3,3], C[3,4], C[3,5],
                                            C[4,4], C[4,5],
                                                    C[5,5], 
            density
        ])
    
# -----------------------------------------------------------------------------
def snow_stiffness(
        temperature, lwc, porosity, pressure=0.0,
        use_fabric=False, eulerangles=None, method=None
    ):
    """
    Computes the effective stiffness tensor for snow using Hill, SCA, or Gassmann's equations.

    :param temperature: Temperature in Celsius.
    :param lwc: Liquid water content (% of pore space filled with water).
    :param porosity: Total porosity of the snow (% of volume that is pore space).
    :param pressure: Pressure in Pascals (default: 0 Pa).
    :param use_fabric: If True, fabric-based stiffness is used; otherwise, isotropic stiffness is computed.
    :param eulerangles: (m x 3) array of Euler angles for fabric-based calculations.
    :param method: Override homogenization method; options: "Hill", "SCA", or "Gassmann".
        - Default: "Hill" for porosity ≤ 50%.
        - "SCA" for high porosity (>50%).
        - "Gassmann" if porosity ≤ 50% and lwc ≥ 80%.
    :return: (K_snow, G_snow, C_snow, rho_eff)
        - K_snow: Bulk modulus of snow (Pa).
        - G_snow: Shear modulus of snow (Pa).
        - C_snow: 6x6 effective stiffness tensor (Pa).
        - rho_eff: Effective density (kg/m³).
    """
    # Convert percentage inputs to fractions
    porosity /= 100.0
    lwc /= 100.0

    # Determine default method
    if method is None:
        if porosity <= 0.50:
            method = "Gassmann" if lwc >= 0.80 else "Hill"
        else:
            method = "SCA"

    # Get properties of water, air, and ice
    K_water, C_water = bulk_modulus_water(temperature)
    rho_water = rho_water_correction(temperature)
    K_air, C_air, rho_air = bulk_modulus_air(temperature)
    G_air = 0.0  # Air has no shear strength
    G_water = 0.0  # Water has no shear strength

    # Compute ice stiffness tensor and bulk/shear moduli
    C_ice = ice_stiffness(temperature, pressure)
    K_ice, G_ice, C_ice_iso = astiffness2moduli(C_ice)  # Convert to isotropic equivalent if needed

    # Compute volume fractions of ice, air, and water
    V_ice, V_air, V_water = compute_volume_fractions(porosity, lwc)

    # Compute effective density
    rho_eff = V_ice * ice_density(temperature) + V_air * rho_air + V_water * rho_water

    # Use fabric-based stiffness if specified
    if use_fabric:
        C_ice = vrh(C_ice, np.linalg.inv(C_ice), eulerangles)

    # Compute bulk modulus using the selected method
    if method == "SCA":
        # Self-Consistent Approximation (SCA) for high porosity
        C_snow = self_consistent_approximation(C_ice, porosity)
        K_snow, G_snow, _ = astiffness2moduli(C_snow)

    elif method == "Gassmann":
        # Compute dry bulk modulus
        K_dry = gassmann_bulk_modulus(K_ice, K_ice, K_air, porosity)

        # Compute saturated bulk modulus using Gassmann's equation
        K_snow = gassmann_bulk_modulus(K_dry, K_ice, K_water, porosity)
        G_snow = G_ice  # Gassmann assumes shear modulus does not change with fluid saturation
        C_snow = C_ice_iso.copy()
        C_snow[:3, :3] = K_snow  # Update bulk modulus component

    else:  # Default to Hill estimate
        C_voigt = (1 - porosity) * C_ice + porosity * (1 - lwc) * C_air + porosity * lwc * C_water
        C_reuss = np.linalg.inv((1 - porosity) * np.linalg.inv(C_ice) + porosity * (1 - lwc) * np.linalg.inv(C_air) + porosity * lwc * np.linalg.inv(C_water))
        C_snow = (C_voigt + C_reuss) / 2
        K_snow, G_snow, _ = astiffness2moduli(C_snow)

    return K_snow, G_snow, C_snow, rho_eff

# -----------------------------------------------------------------------------
def ice_stiffness(
        temperature: float = None, 
        pressure: float = 0.0
    ) -> np.ndarray:
    """
    Computes the stiffness tensor for ice under specified temperature and 
    pressure conditions based on empirical relationships.

    :param temperature: The temperature at which to compute the stiffness 
        tensor. Units are in Celsius.
    :param pressure: The pressure at which to compute the stiffness tensor. 
        Units are in kbar.
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
    # C[2,1] = C[1,2]
    C[2,1] = C[0,2]
    C[4,4] = C[3,3]
    C[5,5] = (C[0,0] - C[0,1] )/2
    
    stiffness = C*1e8

    return(stiffness)

# -----------------------------------------------------------------------------
def ice_density(temperature: float, method: str = None):
    '''
    Compute the density of an ice crystal given the temperature.

    :param temperature: The value of the temperature in celsius
    :type temperature: float
    :param method: Specify which method to be used. Default is 
    '''
    rho_o = 917 # Reference density at 0C
    alpha = 51.0e-6
    if method == 'gammon':
        # This is suitable for warmer ice (i.e. > -20 C)
        rho = (1/rho_o) * (
            1+ 1.576e-4*(temperature) - \
                    2.778e-7*(temperature**2) + \
                        8.850e-9*(temperature**3) - \
                            1.778e-10*(temperature**4) 
            )
        rho = 1/rho 
    else:
        # The reference temperature is 273.15K; 0C
        rho = rho_o * ( 1 - alpha * temperature)
    return(rho)

# -----------------------------------------------------------------------------
def rho_water_correction(temperature: float = 0.0) -> float:
    """
    Corrects the density of water based on temperature using the empirical 
    formula derived from the Kell equation for water >= 0C and Wagner et al. 
    (1994) for supercoold water.

    :param temperature: The temperature at which to compute the water density 
        correction.
    :type temperature: float
    :return: The corrected water density.
    :rtype: float
    """
    if temperature >= 0.0: # Use Kell equation for 0-100 C
        rho_water = (
            999.83952 + \
                16.945176 * temperature - \
                    7.9870401e-3 * temperature**2 - \
                        46.170461e-6 * temperature**3 + \
                            105.56302e-9 * temperature**4 - \
                                280.54253e-12 * temperature**5
            )/(1 + 16.897850e-3 * temperature)
    else:
        rho_water = (
            999.85 + 0.07055 * temperature + \
                0.000273 * temperature ** 2 + \
                    2.2e-7 * temperature ** 3 
        )
    
    return rho_water 

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
    # We can't quite accomodate supercooled water density
    rho_water = np.max( (rho_water, 950) ) 
    # beyond the freezing and vaporization temperatures, things get wonky
    rho_water = np.min( (rho_water, rho_water_correction() )) 
    
    # the water content is the percent of pores that contain water
    grams_air = (1-liquid_water_content/100)*rho_air
    grams_water = (liquid_water_content/100)*rho_water
    rho_wc = grams_air + grams_water
        
    density = (1-porosity/100)*density + (porosity/100)*rho_wc

    return(density, grams_air, grams_water)



# ==============================================================================
#                                Permittivity
# ==============================================================================
# -----------------------------------------------------------------------------

def sand_silt_clay_permittivity_conductivity(
        T, lwc, porosity, 
        composition = {"sand": 1.0, "silt": 0.0, "clay": 0.0}
    ):
    """
    Computes the relative permittivity (ε_r) and electrical conductivity (σ) 
    for a sand-silt-clay mixture based on pressure, temperature, and water content.

    :param P: Overburden pressure in Pascals (Pa).
    :param T: Temperature in degrees Celsius.
    :param lwc: Water content fraction (0 to 1).
    :param composition: Dictionary with fractions of sand, silt, and clay (must sum to 1).
    :return: Relative permittivity (dimensionless), Electrical conductivity (S/m)
    """
    
    porosity /= 100.0 
    lwc /= 100.0
    
    # Check that composition sums to 1
    if not np.isclose(sum(composition.values()), 1.0):
        raise ValueError("Composition fractions must sum to 1.")

    # Material properties for sand, silt, and clay (from experimental data)
    properties = {
        "sand": {"epsilon": 4.7, "sigma": 1e-3},  # Relative permittivity, Conductivity (S/m)
        "silt": {"epsilon": 6.0, "sigma": 2e-3},
        "clay": {"epsilon": 12.0, "sigma": 1e-2},
    }

    # Permittivity and conductivity of water at given temperature (Debye model)
    epsilon_water = 80.36 - 0.3 * (T - 25)  # Empirical relationship
    sigma_water = 0.055 * np.exp(0.02 * (T - 25))  # S/m (increases with temp)
    
    # Permittivity and conductivity of air
    epsilon_air = 1.0
    sigma_air = 0.0
    
    # Compute weighted dry permittivity and conductivity
    epsilon_dry_mix = sum(composition[mat] * properties[mat]["epsilon"] for mat in composition)
    sigma_dry_mix = sum(composition[mat] * properties[mat]["sigma"] for mat in composition)
    
    # Mix permittivity using a nonlinear Lichtenecker (logarithmic) model
    epsilon = np.exp(
        (1 - porosity) * np.log(epsilon_dry_mix) +
        (porosity * (1 - lwc)) * np.log(epsilon_air) +
        (porosity * lwc) * np.log(epsilon_water)
    )

    # Mix conductivity using a weighted average model
    sigma = (1 - porosity) * sigma_dry_mix + \
        porosity * (1 - lwc) * sigma_air + \
            porosity * lwc * sigma_water
    
    return epsilon, sigma

# -----------------------------------------------------------------------------
def get_perm(
        self,
        material,
    ):
    """
    Computes the permittivity and conductivity tensors for materials based on attributes
    contained within material and model class instances.

    :param material: An instance of a class containing attributes for materials, including
                     temperature, density, porosity, water content, and anisotropic properties.
    :type material: Material
    """

    #!!!!! Use the material class as an input since it will have all of the values that you need instead of inputting all of them 
    # material.temp
    
    m = len(material.temp)
    # We will always compute the complex tensor. 
    permittivity_coefficients = np.zeros([m, 6], dtype = complex)
    conductivity_coefficients = np.zeros([m, 6])
    
    # Adjust the stiffness tensor according to the pressure and temperature 
    # conditions for ice
    for ind in range(0, m):
        
        if material.material[ind] == 'ice1h':
            permittivity = ice_permittivity(
                material.temp[ind],
                material.rho[ind],
                center_frequency = self.f0
            )
            conductivity = snow_conductivity(
                permittivity = permittivity, frequency = self.f0
            )
            
        elif material.material[ind] == 'snow':
            permittivity, rho_d = snow_permittivity(
                temperature = material.temp[ind],
                lwc = material.lwc[ind], 
                porosity = material.porosity[ind],
                center_frequency = self.f0
            )
            conductivity = snow_conductivity(
                permittivity = permittivity, frequency = self.f0
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
        
        # Convert back to real
        permittivity = permittivity.real 
        conductivity = conductivity.real
        
        if material.is_anisotropic[ind]:
            
            euler = read_ang(material.angfile[ind])
            p = len(euler[:,0])
            
            pvoigt = np.zeros([3,3])
            cvoigt = np.zeros([3,3])
            preuss = np.zeros([3,3])
            creuss = np.zeros([3,3])
            
            Sp = np.linalg.inv(permittivity)
            Sc = np.linalg.inv(conductivity)
            
            for k in range(0, p):
                R = rotator_zxz(euler[k,:] )
            
                Ri = np.linalg.inv(R)
                #!!!!! We need to do the same for conductivity.  
                pvoigt = pvoigt + ( np.matmul( R, np.matmul(permittivity, R.T) ) )
                cvoigt = cvoigt + ( np.matmul( R, np.matmul(conductivity, R.T) ) )
                preuss = preuss + ( np.matmul( Ri, np.matmul(Sp, Ri.T) ) )
                creuss = creuss + ( np.matmul( Ri, np.matmul(Sc, Ri.T) ) )
            
            pvoigt = pvoigt/p
            cvoigt = cvoigt/p
            preuss = preuss/p 
            creuss = creuss/p 
            preuss = np.linalg.inv(preuss) 
            creuss = np.linalg.inv(creuss) 
            
            # Calculate the hill average 
            permittivity = (pvoigt + preuss)/2
            conductivity = (cvoigt + creuss)/2
            
            # Check to make sure the permittivity tensor is positive definite
            eigenvalues, eigenvectors = np.linalg.eigh(permittivity)
            if np.sum(eigenvalues <= 0) > 0:
                print('Permittivity tensor is not positive definite.')
                
        self.permittivity_coefficients.loc[ind] = np.array(
            [
                permittivity[0,0], permittivity[0,1], permittivity[0,2],
                permittivity[1,1], permittivity[1,2],
                permittivity[2,2]
            ]
        )
        self.conductivity_coefficients.loc[ind] = np.array(
            [
                conductivity[0,0], conductivity[0,1], conductivity[0,2],
                conductivity[1,1], conductivity[1,2],
                conductivity[2,2]
            ] 
        )

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


# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
def snow_permittivity(
        density: float = 917., 
        temperature: float = 0., 
        lwc: float = 0., 
        porosity: float = 50.,
        center_frequency: float = 1e7,
        method: str = "shivola-tiuri"
    ) -> np.ndarray:
    """
    Calculates the complex permittivity of snow based on its density, 
    temperature, liquid water content (LWC), porosity, and the chosen 
    calculation method.

    :param density: Density of the snow in kg/m^3.
    :param temperature: Temperature of the snow in degrees Celsius.
    :param lwc: Liquid water content of the pore space of the snow in 
        percentage.
    :param porosity: Porosity of the snow in percentage.
    :param method: The method to be used for calculating permittivity. Defaults 
        to "shivola-tiuri".
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
    if lwc > 0:
        complex_permittivity = 0.8*lwc + 0.72*(lwc**2)
    else:
        # The complex permittivity tensor for ice is anisotropic, but we will
        # assume snow is homogeneous and use the mean of the complex 
        # permittivity
        eps_i = np.diag(ice_permittivity(temperature, 917, center_frequency).imag).mean()
        complex_permittivity = eps_i * (0.52 * rho_d + 0.62 * rho_d**2)
    
    permittivity = np.eye(3,3) * complex(perm, complex_permittivity)

    return(permittivity, rho_d)

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
