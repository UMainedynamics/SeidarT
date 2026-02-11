import numpy as np
import pandas as pd
from typing import Union, Tuple
from scipy.interpolate import interp1d

__all__ = [
    'vrh2',
    'vrh4',
    'pressure_array',
    'anisotropic_boolean',
    'get_seismic',
    'get_perm',
    'rho_water_correction',
    'isotropic_stiffness_tensor',
    'isotropic_permittivity_tensor',
    'porewater_correction',
    'moduli2stiffness',
    'astiffness2moduli',
    'voigt_to_full_tensor',
    'get_christoffel_matrix',
    'get_complex_refractive_index',    
    'anisotropy_parameters',
    'compute_volume_fractions',
    'bulk_modulus_water',
    'bulk_modulus_air',
    'moduli_basalt',
    'moduli_granite',
    'moduli_gneiss',
    'moduli_sand_silt_clay',
    'snow_stiffness',
    'ice_stiffness_petrenko',
    'ice_stiffness_gagnon',
    'ice_permittivity',
    'ice_density',
    'tensor2velocities',
    'snow_permittivity',
    'snow_conductivity',
    'sand_silt_clay_permittivity_conductivity',
    'read_ang',
    'rotator',
    'rotator_euler',
    'trend_plunge_to_rotation_matrix',
    'rotation_matrix_to_trend_plunge_orientation',
    'rotation_matrix_to_euler_angles',
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
def vrh2(permittivity, conductivity, eulerangles):
    """
    
    """
    
    p = len(eulerangles[:,0])
    
    pvoigt = np.zeros([3,3])
    cvoigt = np.zeros([3,3])
    preuss = np.zeros([3,3])
    creuss = np.zeros([3,3])
    
    Sp = np.linalg.inv(permittivity)
    Sc = np.linalg.inv(conductivity)
    
    for k in range(0, p):
        R = rotator_euler(eulerangles[k,:] )
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
    
    return permittivity, conductivity

# -----------------------------------------------------------------------------
def vrh4(C, eulerangles):
    """
    
    """
    S = np.linalg.inv(C)
    m, n = eulerangles.shape
    cvoigt = np.zeros([6,6])
    creuss = np.zeros([6,6])
    for k in range(m):
        R = rotator_euler(eulerangles[k,:] )
        M = bond(R)
        N = np.linalg.inv(M).T
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
def rotator(angle: float, axis: str) -> np.ndarray:
    """
    General 3x3 rotation matrix calculation for any axis. 
    """
    R = np.zeros([3,3])
    if axis == 'x':
        R[0,:] = [ 1.0, 0.0, 0.0 ]
        R[1,:] = [ 0.0, np.cos( angle ), -np.sin( angle ) ]
        R[2,:] = [ 0.0, np.sin( angle ),  np.cos( angle ) ]
    elif axis == 'y':
        R[0,:] = [  np.cos( angle ), 0.0, np.sin( angle )]
        R[1,:] = [ 0.0, 1.0, 0.0 ]
        R[2,:] = [ -np.sin( angle ), 0.0, np.cos( angle ) ]
    elif axis == 'z':
        R[0,:] = [ np.cos( angle ), -np.sin( angle ), 0.0 ] 
        R[1,:] = [ np.sin( angle ),  np.cos( angle ), 0.0 ]
        R[2,:] = [ 0.0, 0.0, 1.0 ]
    else: 
        raise ValueError("axis must be 'x', 'y', or 'z'.")
    
    return R

# This is a generic routine and the name should be changed here and anywhere it is used.
def rotator_euler(eul: np.ndarray, rot: str = 'zxz') -> np.ndarray:
    """
    Generates a rotation matrix from Euler angles for a given rotation sequence.
    Default is a 'zxz' rotation.
    
    :param eul: An array containing the three Euler angles. Angles are in radians
    :type eul: np.ndarray
    :return: The 3x3 rotation matrix derived from the Euler angles.
    :rtype: np.ndarray
    """
    if len(rot) != 3:
        raise ValueError(
             f"Rotation sequence needs three axes specified. For example, 'zxz' or 'xyx', etc. You have specified '{rot}'."
        ) 
     
    # From the 3 euler angles for the zxz rotation, compute the rotation matrix
    D = rotator(eul[0], rot[0])
    C = rotator(eul[1], rot[1]) 
    B = rotator(eul[2], rot[2])
    
    # R = np.matmul(B, np.matmul(C, D) )
    R = B @ C @ D
    return(R)

# --------------------------------------------------------------------------
def trend_plunge_to_rotation_matrix(trend, plunge, orientation): 
    '''
    Compute the rotation matrix from a trend and plunge pair. 
    
    :param trend: The trend value in degrees 
    :type trend: float 
    :param plunge: The plunge value in degrees  
    :type plunge: float 
    :param orientation: The orientation/rotation angle around the axis of 
        the plunge.
    :type orientation: float
    
    :return: rotation_matrix
    :rtype: np.ndarray
    '''
    # Convert trend and plunge to radians
    trend_rad = np.radians(trend)
    plunge_rad = np.radians(plunge)
    orient_rad = np.radians(orientation)
    
    # Rotation matrix for trend (rotation around Z-axis)
    Rt = np.array([
        [np.cos(trend_rad), -np.sin(trend_rad), 0],
        [np.sin(trend_rad), np.cos(trend_rad), 0],
        [0, 0, 1]
    ])
    
    # Rotation matrix for plunge (rotation around X-axis)
    Rp = np.array([
        [1, 0, 0],
        [0, np.cos(plunge_rad), -np.sin(plunge_rad)],
        [0, np.sin(plunge_rad), np.cos(plunge_rad)]
    ])
    
    Ro = np.array([
        [np.cos(orient_rad), -np.sin(orient_rad), 0],
        [np.sin(orient_rad), np.cos(orient_rad), 0],
        [0, 0, 1]
    ])
    
    # Combined rotation matrix
    rotation_matrix = Rt @ Rp @ Ro
    return rotation_matrix

def rotation_matrix_to_trend_plunge_orientation(R):
    """
    Consistent inverse for trendplungetorotationmatrix.
    Returns trend, plunge, orientation in radians.
    """
    # 1) rotation matrix -> Bunge Euler angles (you already have this)
    phi1, Phi, phi2 = rotation_matrix_to_euler_angles(R)

    # 2) Rebuild the same rotation using your trend/plunge/orientation
    #    but we want the c-axis, so just grab it from R directly:
    #    c_local = (0,0,1), so c_global = R @ c_local = column 2
    cx, cy, cz = R @ np.array([0.0, 0.0, 1.0])
    
    # 3) trend and plunge of c-axis
    plunge = np.arcsin(cz)          # radians
    trend  = np.arctan2(cy, cx)     # radians
    
    # 4) orientation: rotation about c-axis.
    #    Use the image of local x as "a-axis" and measure its angle around c.
    ax, ay, az = R @ np.array([1.0, 0.0, 0.0])
    
    c_vec = np.array([cx, cy, cz])
    a_vec = np.array([ax, ay, az])
    
    # Project a_vec into plane perpendicular to c_vec
    a_perp = a_vec - np.dot(a_vec, c_vec) * c_vec
    a_perp /= np.linalg.norm(a_perp)
    
    # Reference direction in that plane: projection of global x
    ref = np.array([1.0, 0.0, 0.0])
    ref_perp = ref - np.dot(ref, c_vec) * c_vec
    ref_perp /= np.linalg.norm(ref_perp)
    
    cos_o = np.clip(np.dot(ref_perp, a_perp), -1.0, 1.0)
    sin_o = np.dot(c_vec, np.cross(ref_perp, a_perp))
    orientation = np.arctan2(sin_o, cos_o)  # radians
    
    return trend, plunge, orientation

# --------------------------------------------------------------------------
def rotation_matrix_to_euler_angles(rotation_matrix):
    """
    Convert rotation matrix to Euler angles
    
    :param rotation_matrix: The 3-by-3 orthogonal rotation matrix
    :type rotation_matrix: np.ndarray 
    
    """
    if rotation_matrix.shape != (3, 3):
        raise ValueError("Input matrix must be 3x3")
    
    # Ensure the matrix is a valid rotation matrix
    if not np.allclose(
            np.dot(rotation_matrix, rotation_matrix.T), np.eye(3)
        ) or not np.isclose(np.linalg.det(rotation_matrix), 1):
        raise ValueError("Input matrix is not a valid rotation matrix")
    
    # Extract the Euler angles
    Phi = np.arccos(rotation_matrix[2, 2])
    phi1 = np.arctan2(rotation_matrix[0, 2], -rotation_matrix[1, 2])
    phi2 = np.arctan2(rotation_matrix[2, 0], rotation_matrix[2, 1])
    
    return phi1, Phi, phi2

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
        if isinstance(T, pd.DataFrame):
            vx = c_light * np.sqrt(1/T['e11'])
            vy = c_light * np.sqrt(1/T['e22'])
            vz = c_light * np.sqrt(1/T['e33'])
            return vx, vy, vz
        else:
            temp = np.zeros([3,3])
            temp[0,:] = T[0:3] 
            temp[1,1:] = T[3:5]
            temp[2,2] = T[5]
            T = np.triu(temp).T + np.triu(temp, 1)
    if T.shape == (21,) or T.shape == (22,): # First case is a numpy array and second is a dataframe
        temp = np.zeros([6,6])
        if isinstance(T, pd.DataFrame):
            rho = T['rho']
            vpv = np.sqrt(T['c33']/rho)
            vph = np.sqrt(T['c11']/rho)
            vsv = np.sqrt(T['c44']/rho)
            vsh = np.sqrt( (T['c11'] - T['c12'] )/(2*rho))
            return vpv, vph, vsv, vsh
        else:
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

# -----------------------------------------------------------------------------
def voigt_to_full_tensor(C_voigt):
    """
    Convert a 6x6 stiffness matrix in Voigt notation to a 3x3x3x3 tensor.
    
    Parameters:
        C_voigt: numpy.ndarray
            6x6 stiffness matrix in Voigt notation.
            
    Returns:
        C_full: numpy.ndarray
            3x3x3x3 stiffness tensor.
    """
    # Initialize full tensor with zeros
    C_full = np.zeros((3, 3, 3, 3))
    
    # Mapping from sorted index pairs to Voigt index
    voigt_map = {
        (0, 0): 0,
        (1, 1): 1,
        (2, 2): 2,
        (1, 2): 3,
        (0, 2): 4,
        (0, 1): 5
    }
    
    # Conversion factor: 1 if i == j, else 1/sqrt(2)
    def factor(i, j):
        return 1.0 if i == j else 1.0/np.sqrt(2)
    
    # Loop over all tensor indices
    for i in range(3):
        for j in range(3):
            # Determine the Voigt index for the (i, j) pair (order doesn't matter)
            I = voigt_map[tuple(sorted((i, j)))]
            for k in range(3):
                for l in range(3):
                    # Determine the Voigt index for the (k, l) pair
                    J = voigt_map[tuple(sorted((k, l)))]
                    # Multiply the appropriate conversion factors and assign
                    C_full[i, j, k, l] = factor(i, j) * factor(k, l) * C_voigt[I, J]
    return C_full

# -----------------------------------------------------------------------------
def voigt_to_mandel(C_voigt):
    """
    Convert a 6x6 fourth-order tensor in Voigt notation to Mandel notation.
    
    Parameters
    ----------
    C_voigt : (6, 6) array_like
        Tensor in Voigt notation.
    
    Returns
    -------
    C_mandel : (6, 6) ndarray
        Tensor in Mandel notation.
    """
    C_voigt = np.asarray(C_voigt, dtype=float)
    
    # Scaling vector: [1, 1, 1, sqrt(2), sqrt(2), sqrt(2)]
    s = np.ones(6)
    s[3:] = np.sqrt(2.0)
    
    # Similarity transform: C_M = S * C_V * S, with S = diag(s)
    # Use outer product to avoid constructing the full diagonal matrix explicitly.
    C_mandel = C_voigt * s[:, None] * s[None, :]
    
    return C_mandel

# -----------------------------------------------------------------------------
def get_christoffel_matrix(
        coefficients: pd.Series, 
        f0: float, 
        direction, 
    ):
    """
    
    """
    direction_map = {
        'x':  np.array([1, 0, 0]),
        'y':  np.array([0, 1, 0]),
        'z':  np.array([0, 0, 1]),
        'xy': np.array([1, 1, 0]) / np.sqrt(2),
        'xz': np.array([1, 0, 1]) / np.sqrt(2),
        'yz': np.array([0, 1, 1]) / np.sqrt(2)
    }
            
    # Handle string input (predefined directions)
    if isinstance(direction, str):
        if direction not in direction_map:
            raise ValueError(f"Invalid direction: {direction}. Choose from {list(direction_map.keys())}.")
        n = direction_map[direction]
    
    # Handle unit vector input
    elif isinstance(direction, np.ndarray) and direction.shape == (3,):
        n = direction
    else:
        raise TypeError("Direction must be either a string ('x', 'y', ...) or a 3-element numpy array.")
    
    # Normalize direction vector
    if np.linalg.norm(n) == 0:
        raise ValueError("Propagation vector cannot be zero.")
    n = n / np.linalg.norm(n)
    
    rho = coefficients['rho']
    C_full = np.zeros((3, 3, 3, 3))
    # Extract stiffness coefficients from the specified row
    C = np.array([
        [coefficients['c11'], coefficients['c12'], coefficients['c13'], coefficients['c14'], coefficients['c15'], coefficients['c16']],
        [coefficients['c12'], coefficients['c22'], coefficients['c23'], coefficients['c24'], coefficients['c25'], coefficients['c26']],
        [coefficients['c13'], coefficients['c23'], coefficients['c33'], coefficients['c34'], coefficients['c35'], coefficients['c36']],
        [coefficients['c14'], coefficients['c24'], coefficients['c34'], coefficients['c44'], coefficients['c45'], coefficients['c46']],
        [coefficients['c15'], coefficients['c25'], coefficients['c35'], coefficients['c45'], coefficients['c55'], coefficients['c56']],
        [coefficients['c16'], coefficients['c26'], coefficients['c36'], coefficients['c46'], coefficients['c56'], coefficients['c66']]
    ])
    
    C_full = voigt_to_full_tensor(C)
    Gamma = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            Gamma[i, j] = np.sum(C_full[i, :, j, :] * np.outer(n, n))/rho
    
    eigenvalues, eigenvectors = np.linalg.eigh(Gamma) 
    velocities = np.sqrt(np.clip(eigenvalues, a_min=0.0, a_max = None) )
    velocities.sort() 
    return Gamma, eigenvalues, eigenvectors, velocities 

# -----------------------------------------------------------------------------
def get_complex_refractive_index(
        coefficients: pd.Series, 
        f0: float, 
        direction, 
        mu = 4e-7 * np.pi 
    ):
    """
    Coefficients are given as relative permittivity but not relative conductivity
    """
    direction_map = {
        'x':  np.array([1, 0, 0]),
        'y':  np.array([0, 1, 0]),
        'z':  np.array([0, 0, 1]),
        'xy': np.array([1, 1, 0]) / np.sqrt(2),
        'xz': np.array([1, 0, 1]) / np.sqrt(2),
        'yz': np.array([0, 1, 1]) / np.sqrt(2)
    }
            
    # Handle string input (predefined directions)
    if isinstance(direction, str):
        if direction not in direction_map:
            raise ValueError(f"Invalid direction: {direction}. Choose from {list(direction_map.keys())}.")
        n = direction_map[direction]
    
    # Handle unit vector input
    elif isinstance(direction, np.ndarray) and direction.shape == (3,):
        n = direction
    else:
        raise TypeError("Direction must be either a string ('x', 'y', ...) or a 3-element numpy array.")
    
    # Normalize direction vector
    if np.linalg.norm(n) == 0:
        raise ValueError("Propagation vector cannot be zero.")
    n = n / np.linalg.norm(n)
    
    # --------------------------------------------------------------------------
        
    # Extract permittivity values from the specified row
    epsilon = np.array([
        [coefficients['e11'], coefficients['e12'], coefficients['e13']],
        [coefficients['e12'], coefficients['e22'], coefficients['e23']],
        [coefficients['e13'], coefficients['e23'], coefficients['e33']]
    ], dtype = float)
    
    # Extract conductivity values from the specified row
    sigma = np.array([
        [coefficients['s11'], coefficients['s12'], coefficients['s13']],
        [coefficients['s12'], coefficients['s22'], coefficients['s23']],
        [coefficients['s13'], coefficients['s23'], coefficients['s33']]
    ], dtype = float)
    
    omega = 2.0 * np.pi * f0  # Angular frequency
    
    # complex permittivity
    epsilon_eff = epsilon - 1j * sigma / (omega*eps0)  # Complex permittivity
    # epsilon_eff *= eps0
    epsilon_eff_inv = np.linalg.inv(epsilon_eff)
    
    # --------------------------------------------------------------------------
    I = np.eye(3, dtype = float) 
    P = I - np.outer(n,n)
    A = P @ epsilon_eff_inv @ P
    
    # --------------------------------------------------------------------------
    e1 = np.cross(n, [0.0, 0.0, 1.0] )
    if np.linalg.norm(e1) < 1e-8:
        e1 = np.cross(n, [0.0, 1.0, 0.0])
    
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    e2 = e2 / np.linalg.norm(e2) 
    
    A_perp = np.array([
        [e1 @ A @ e1, e1 @ A @ e2], 
        [e2 @ A @ e1, e2 @ A @ e2]
    ], dtype = complex) 
    
    lam, vecs = np.linalg.eig(A_perp) 
    
    complex_index = np.empty_like(lam, dtype = complex) 
    for i, val in enumerate(lam):
        if val == 0.0:
            complex_index[i] = np.nan + 1j * np.nan 
        else:
            index_i = 1.0 / np.sqrt(val)
            if index_i.real < 0:
                index_i = -index_i
            
            complex_index[i] = index_i
    
    return lam, vecs, complex_index                 

# -----------------------------------------------------------------------------
def anisotropy_parameters(tensor):
    eigvals, eigvecs = np.linalg.eigh(tensor)
    ratio = eigvals.max() / eigvals.min() 
    birefringence = eigvals.max() - eigvals.min() 
    spread = birefringence / eigvals.mean()
    return ratio, birefringence, spread

# ==============================================================================
#                                   Seismic
# ==============================================================================
def moduli2stiffness(K, G):
    C = np.zeros([6,6])
    C[(0,1,2),(0,1,2)] = K + 4*G / 3 
    C[(0,0,1,1,2,2),(1,2,0,2,0,1)] = K - 2*G/3
    C[(3,4,5),(3,4,5)] = G
    return C

# -----------------------------------------------------------------------------
def astiffness2moduli(C, transversely_isotropic=True):
    """
    Computes the isotropic bulk and shear modulus from an anisotropic stiffness tensor
    using Hill's average for a transversely isotropic stiffness tensor.
    
    Parameters:
        C (numpy.ndarray): 6x6 stiffness tensor in Voigt notation (Pa)
    
    Returns:
        (K_iso, mu_iso, C_isotropic): Tuple of isotropic bulk modulus (Pa), shear modulus (Pa),
                                      and the 6x6 isotropic stiffness tensor (Pa)
    """
    S = np.linalg.inv(C)
    
    if transversely_isotropic:
        # --- Voigt average for bulk modulus ---
        Kv = (C[0,0] + C[0,1] + 2*C[0,2] + 2*C[2,2] + 2*C[3,3]) / 6  # fixed typo in parens

        # --- Reuss average for bulk modulus (approximate form for TI) ---
        # Use compliance-based version to be safe
        Kr = 1 / (S[0,0] + S[1,1] + S[2,2] + 2 * (S[0,1] + S[0,2] + S[1,2]))  # matches isotropic form

        # --- Voigt average for shear modulus ---
        mu_V = (2*C[5,5] + C[3,3] + C[4,4] + C[3,3]) / 5  # general TI, no assumption about c44=c55

        # --- Reuss average for shear modulus ---
        mu_R = 5 / (2/C[5,5] + 1/C[3,3] + 1/C[4,4] + 1/C[3,3])  # same terms in reciprocal
    else:
        # --- Voigt average ---
        K_V = (C[0,0] + C[1,1] + C[2,2] + 2 * (C[0,1] + C[0,2] + C[1,2])) / 9
        mu_V = ((C[0,0] + C[1,1] + C[2,2] - C[0,1] - C[0,2] - C[1,2]) +
                3 * (C[3,3] + C[4,4] + C[5,5])) / 15

        # --- Reuss average ---
        K_R = 1 / (S[0,0] + S[1,1] + S[2,2] + 2 * (S[0,1] + S[0,2] + S[1,2]))
        mu_R = 15 / ((4 * (S[0,0] + S[1,1] + S[2,2] - S[0,1] - S[0,2] - S[1,2])) +
                     3 * (S[3,3] + S[4,4] + S[5,5]))

        Kv = K_V
        Kr = K_R

    # --- Hill averages ---
    K_iso = (Kv + Kr) / 2
    mu_iso = (mu_V + mu_R) / 2

    # --- Construct isotropic stiffness tensor ---
    C_isotropic = np.zeros((6, 6))
    lam = K_iso - 2 * mu_iso / 3

    for i in range(3):
        for j in range(3):
            C_isotropic[i, j] = lam
        C_isotropic[i, i] = lam + 2 * mu_iso

    C_isotropic[3:, 3:] = np.eye(3) * mu_iso

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

def gassmann_bulk_modulus(K_dry, K_mineral, K_fluid, porosity, small_denominator=1e-12):
    """
    Computes the bulk modulus of a fluid-saturated porous material using Gassmann's equation.
    
    :param K_dry: Bulk modulus of the dry matrix (Pa).
    :param K_mineral: Bulk modulus of the solid grains (Pa).
    :param K_fluid: Bulk modulus of the saturating fluid (Pa).
    :param porosity: Porosity (fraction, 0 to 1).
    :param small_denominator: A small threshold to avoid division by zero or negative denominators.
    :return: Saturated bulk modulus (Pa).
    """
    numerator = (1 - K_dry / K_mineral) ** 2
    denominator = (porosity / K_fluid) + ((1 - porosity) / K_mineral) - (K_dry / (K_mineral**2))
    
    # Check if the denominator is too small or negative.
    if denominator <= small_denominator:
        # You can clamp or do some fallback:
        denominator = small_denominator
    
    return K_dry + (numerator / denominator)

def compute_dry_bulk_modulus(K_ice, porosity_percent, exponent=2):
    """
    Estimate the dry frame bulk modulus (K_dry) for porous snow/ice
    using the relationship:
    
        K_dry = K_ice * (1 - porosity_fraction)**exponent
        
    Parameters:
        K_ice (float): Bulk modulus of ice (Pa).
        porosity_percent (float): Porosity as a percentage (0-100).
        exponent (float, optional): Exponent controlling sensitivity. 
                                    Default is 2.
    
    Returns:
        float: Estimated dry frame bulk modulus (Pa).
    """
    porosity_fraction = porosity_percent / 100.0
    return K_ice * (1 - porosity_fraction) ** exponent

def compute_volume_fractions(porosity, lwc):
    """
    Computes volume fractions of ice, air, and water based on porosity and liquid water content.
    
    :param porosity: Porosity fraction (0 to 1)
    :param lwc: Liquid water content fraction (0 to 1)
    :return: (V_ice, V_air, V_water)
    """
    # porosity /= 100 
    # lwc /= 100
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
def self_consistent_approximation(C, C_air, porosity, tol=1e-6, max_iter=100):
    """
    Computes the effective stiffness tensor for highly porous snow using SCA.
    
    :param C: 6x6 stiffness tensor.
    :param porosity: Volume fraction of air (0 to 1).
    :return: 6x6 effective stiffness tensor of porous media.
    """
    
    C_eff = (1 - porosity) * C  # Initial guess
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
        material, snow_method = 'Hill', ice_method = 'Petrenko'
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
            # Assume a constant p ressure of 0.1 MPa (Why? because this is 
            # approximately 1 ATM)
            pressure = 0.1 * 1e-1 # in kbar
            if ice_method == 'Gagnon':
                C_ice = ice_stiffness_gagnon(material.temp[ind], pressure)
            else:
                C_ice = ice_stiffness_petrenko(material.temp[ind], pressure)
            
            if material.is_anisotropic[ind]:
                euler = read_ang(material.angfile[ind])
                C = vrh4(C_ice, euler)
            else:
                print('Computing the homogeneous stiffness coefficients for ice1h.')
                _,_, C = astiffness2moduli(C_ice)
            # The stiffness matrix needs to be positive definite which means all
            # of the eigenvalues must be positive
            eigenvalues, eigenvectors = np.linalg.eigh(C)
            if np.sum(eigenvalues <= 0) > 0:
                print('Stiffness tensor is not positive definite.')
            else:
                print('Stiffness tensor is positive definite')            
        elif material.material[ind] == 'snow':
            if material.is_anisotropic[ind]:
                __, __, C, density = snow_stiffness(
                    material.temp[ind], material.lwc, material.porosity,
                    use_fabric=True, eulerangles=euler, method=snow_method, 
                    exponent=2, epsilon=1e-6, small_denominator=1e-12
                )
            else:
                __, __, C, density = snow_stiffness(
                    material.temp[ind], material.lwc[ind], 
                    material.porosity[ind], method=snow_method, exponent=3
                )
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
        use_fabric=False, eulerangles=None, method=None, exponent=2, 
        epsilon=1e-6, small_denominator=1e-12, ice_method = 'Petrenko'
    ):
    """
    Computes the effective stiffness tensor for snow using Hill, SCA, or Gassmann's equations.
    
    :param temperature: Temperature in Celsius.
    :param lwc: Liquid water content (% of pore space filled with water).
    :param porosity: Total porosity of the snow (% of volume that is pore space).
    :param pressure: Pressure in Pascals (default: 0 Pa).
    :param use_fabric: If True, fabric-based stiffness is used; otherwise, isotropic stiffness is computed.
    :param eulerangles: (m x 3) array of Euler angles for fabric-based calculations.
    :param method: Override homogenization method; options: "Hill", "SCA", "Gassmann", "SCA-Gassmann", etc.
    :param exponent: Exponent for reducing K_ice (and optionally G_ice) by porosity.
    :param epsilon: Small regularization for fluid shear terms in Hill or SCA-Hill mixing.
    :param small_denominator: Threshold to clamp Gassmann's denominator.
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
        if porosity <= 0.50 and lwc >= 0.8:
            method = "Gassmann"
        elif porosity > 0.5 and lwc > 0.0:
            method = "SCA-Gassmann" if lwc >= 0.8 else "SCA-Hill"
        elif porosity > 0.5 and lwc == 0.0:
            method = "SCA"
        else:
            method = "Hill" 
    
    # Get properties of water, air, and ice
    K_water, C_water = bulk_modulus_water(temperature)
    rho_water = rho_water_correction(temperature)
    K_air, C_air, rho_air = bulk_modulus_air(temperature)
    G_air = 0.0  # Air has no shear strength
    G_water = 0.0  # Water has no shear strength

    # Compute ice stiffness tensor and bulk/shear moduli
    if ice_method == 'Gagnon':
        C_ice = ice_stiffness_gagnon(temperature, pressure)
    else: 
        C_ice = ice_stiffness_petrenko(temperature, pressure)
    
    K_ice, G_ice, C_ice_iso = astiffness2moduli(C_ice)  # Convert to isotropic equivalent if needed

    # If using fabric-based stiffness
    if use_fabric:
        # Typically you'd do: S_ice = np.linalg.inv(C_ice)
        # C_ice = vrh4(C_ice, S_ice, eulerangles)
        # but depends on your function signature
        pass

    # Compute volume fractions of ice, air, and water
    V_ice, V_air, V_water = compute_volume_fractions(porosity, lwc)

    # Compute effective density
    rho_eff = V_ice * ice_density(temperature) + V_air * rho_air + V_water * rho_water

    # Empirical "dry" frame moduli for the porous skeleton
    K_dry = K_ice * (1 - porosity)**exponent
    # (Optional) You can do the same for shear:
    G_dry = G_ice * (1 - porosity)**exponent

    # Now compute the final bulk/shear moduli using your chosen method
    if method == "SCA":
        # Self-Consistent Approximation (SCA) for high porosity
        C_snow = self_consistent_approximation(C_ice, C_air, porosity)
        K_snow, G_snow, _ = astiffness2moduli(C_snow)

    elif method == "Gassmann":
        # Single-step fluid substitution with water, skipping air.
        # This is recommended if K_dry is your best estimate of the "dry" skeleton.
        K_snow = gassmann_bulk_modulus(K_dry, K_ice, K_water, porosity, small_denominator)
        # Gassmann doesn't change the shear modulus, so keep the "dry" shear
        G_snow = G_dry
        C_snow = moduli2stiffness(K_snow, G_snow)

    elif method == "SCA-Gassmann":
        # 1) SCA to get a better "dry" framework at high porosity
        C_snow = self_consistent_approximation(C_ice, C_air, porosity)
        K_dry_sca, G_sca, _ = astiffness2moduli(C_snow)
        # 2) Now apply Gassmann with water
        K_snow = gassmann_bulk_modulus(K_dry_sca, K_ice, K_water, porosity, small_denominator)
        # Shear remains the SCA "dry" shear
        G_snow = G_sca
        C_snow = moduli2stiffness(K_snow, G_snow)
    elif method == "SCA-Hill":
        # Example: SCA for the dry skeleton, then Hill average with fluid phases
        C_snow = self_consistent_approximation(C_ice, C_air, porosity)
        # Regularize fluid shear
        C_air_reg = C_air.copy()
        C_water_reg = C_water.copy()
        for i in range(3, 6):
            if C_air_reg[i, i] == 0:
                C_air_reg[i, i] = epsilon
            if C_water_reg[i, i] == 0:
                C_water_reg[i, i] = epsilon
        # Weighted Voigt/Reuss
        C_voigt = (1 - porosity) * C_snow + porosity * (1 - lwc) * C_air + porosity * lwc * C_water
        C_reuss = np.linalg.pinv(
            (1 - porosity) * np.linalg.pinv(C_snow) +
            porosity * (1 - lwc) * np.linalg.pinv(C_air_reg) +
            porosity * lwc * np.linalg.pinv(C_water_reg)
        )
        C_snow = (C_voigt + C_reuss)/2
        K_snow, G_snow, _ = astiffness2moduli(C_snow)
    else:  # Default to Hill estimate
        C_air_reg = C_air.copy()
        C_water_reg = C_water.copy()
        for i in range(3, 6):
            if C_air_reg[i, i] == 0:
                C_air_reg[i, i] = epsilon
            if C_water_reg[i, i] == 0:
                C_water_reg[i, i] = epsilon
        
        C_voigt = (1 - porosity)*C_ice + porosity*(1 - lwc)*C_air + porosity*lwc*C_water
        C_reuss = np.linalg.pinv(
            (1 - porosity)*np.linalg.pinv(C_ice) +
            porosity*(1 - lwc)*np.linalg.pinv(C_air_reg) +
            porosity*lwc*np.linalg.pinv(C_water_reg)
        )
        C_snow = (C_voigt + C_reuss)/2
        K_snow, G_snow, _ = astiffness2moduli(C_snow)
    
    return K_snow, G_snow, C_snow, rho_eff

# -----------------------------------------------------------------------------
def ice_stiffness_petrenko(temperature: float, pressure = 0) -> np.ndarray:
    """
    Calculate the 6×6 stiffness tensor (in Voigt notation) for ice Ih based on
    empirical formulas similar to those given in "Physics of Ice" by Petrenko & Whitworth.
    
    Parameters:
        temperature (float): Temperature in °C.
        pressure (float): Pressure in kbar.
        
    Returns:
        np.ndarray: 6×6 stiffness tensor in Pa.
        
    Note:
        The five independent elastic constants for hexagonal ice are C11, C12, C13, C33, and C44.
        The sixth constant, C66, is given by (C11 - C12)/2.
        The formulas below use linear (or quadratic) pressure and temperature derivatives.
        The coefficients here are only examples—please check with the actual values from
        "Physics of Ice" for your application.
    """
    # Example reference values at T0 = 0 °C and P0 = 0 kbar (in GPa)
    C11_0 = 11.0    # GPa
    C12_0 = 5.0     # GPa
    C13_0 = 4.5     # GPa
    C33_0 = 11.5    # GPa
    C44_0 = 3.2     # GPa
    
    # Temperature derivatives (example values, in GPa/°C)
    dC11_dT = -0.05  
    dC12_dT = -0.03  
    dC13_dT = -0.03  
    dC33_dT = -0.05  
    dC44_dT = -0.02  
    
    # Pressure derivatives (example values, in GPa/kbar)
    dC11_dP = 0.4    
    dC12_dP = 0.3    
    dC13_dP = 0.35   
    dC33_dP = 0.5    
    dC44_dP = 0.1    
    
    # Calculate the elastic constants at the given temperature and pressure
    # (Assuming linear dependence; if quadratic terms are needed, add them accordingly)
    C11 = C11_0 + dC11_dT * temperature + dC11_dP * pressure
    C12 = C12_0 + dC12_dT * temperature + dC12_dP * pressure
    C13 = C13_0 + dC13_dT * temperature + dC13_dP * pressure
    C33 = C33_0 + dC33_dT * temperature + dC33_dP * pressure
    C44 = C44_0 + dC44_dT * temperature + dC44_dP * pressure
    
    # For hexagonal ice, C66 is given by (C11 - C12)/2.
    C66 = (C11 - C12) / 2.0
    
    # Convert from GPa to Pa (1 GPa = 1e9 Pa)
    conv = 1e9
    C11 *= conv
    C12 *= conv
    C13 *= conv
    C33 *= conv
    C44 *= conv
    C66 *= conv
    
    # Build the 6x6 stiffness tensor in Voigt notation.
    # The tensor for hexagonal symmetry (ice Ih) is:
    # [ [ C11, C12, C13,   0,   0,   0 ],
    #   [ C12, C11, C13,   0,   0,   0 ],
    #   [ C13, C13, C33,   0,   0,   0 ],
    #   [   0,   0,   0, C44,   0,   0 ],
    #   [   0,   0,   0,   0, C44,   0 ],
    #   [   0,   0,   0,   0,   0, C66] ]
    C = np.zeros((6,6))
    C[0,0] = C11; C[0,1] = C12; C[0,2] = C13
    C[1,0] = C12; C[1,1] = C11; C[1,2] = C13
    C[2,0] = C13; C[2,1] = C13; C[2,2] = C33
    C[3,3] = C44
    C[4,4] = C44
    C[5,5] = C66
    
    return C

# -----------------------------------------------------------------------------
def ice_stiffness_gagnon(
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
    
    C[0,0] = 136.813 - 0.28941*temperature - 0.0017827*(temperature**2) \
      + 4.7546*pressure - 0.11307*(pressure**2) 
    C[0,1] = 72.639 - 0.14673*temperature - 0.00090362*(temperature**2) \
      + 7.2109*pressure + 0.37096*(pressure**2)
    C[0,2] = 61.445 - 0.11916*temperature - 0.00073120*(temperature**2) \
      + 2.8083*pressure - .12510*(pressure**2)
    C[2,2] = 157.960 - 0.31129*temperature - 0.0018948*(temperature**2) \
      + 2.8083*pressure - 0.12510*(pressure**2)
    C[3,3] = 29.490 - 0.062874*temperature - 0.00038956*(temperature**2) \
      + 1.5055*pressure + 0.2677*(pressure**2)
    
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
def ice_density(temperature: float, pressure: float = 0.0, method: str = None):
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
    elif method == 'gagnon':
        rho = 922.15 + 10.28e-2 * pressure - 0.250310e-4 * pressure **2 - 0.1445 * temperature + 2.77e-4 * temperature**2
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
            permittivity, conductivity = vrh2(permittivity, conductivity, euler)
            
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
        perm_ast = complex(perm, -complex_perm / 1.2)
        
    permittivity = np.eye(3,3) * perm_ast 
    if method == 'fujita':
        permittivity[2,2] = complex(perm + dP, -complex_perm) 

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
        temperature, density, porosity, 0
    )
    
    rho_water = rho_water_correction(temperature)
    # LWC is a percentage we need it in g/cm3
    lwc = lwc/rho_water 
    # dry snow density is in kg/m3 but we also need it in g/cm3
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
    
    permittivity = np.eye(3,3) * complex(perm, -complex_permittivity)

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
