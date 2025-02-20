import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.animation as anim
import os.path
from typing import Optional
from subprocess import call
from scipy.io import FortranFile
import glob2
import copy
import json

import seidart.routines.materials as mf
from seidart.routines.prjbuild import image2int, update_json, readwrite_json

__all__ = [
    'read_dat',
    'complex2str',
    'str2complex',
    'loadproject',
    'airsurf',
    'cpmlcompute',
    'cpml_parameters',
    'rcxgen',
    'coherdt',
    'coherstf',
    'stfvec2mat',
    'movingsrc',
    'indvar',
    'agc',
    'correct_geometric_spreading',
    'exponential_gain',
    'parameter_profile_1d',
    'CFL', 'clight',
]

# --------------------------------- Globals ------------------------------------
clight = 2.99792458e8  # Speed of light in vacuum (m/s)
eps0 = 8.85418782e-12  # Permittivity of free space
mu0 = 4.0 * np.pi * 1.0e-7  # Permeability of free space
mu_r = 1.0

# These can be changed for control over the cpml parameters
#CPML
#----
sig_opt_scalar = 1.2
alpha_max_scalar = 1.0
NP = 2  # Numerical parameter for CPML
NPA = 2  # Additional numerical parameter for CPML
kappa_max = 5  # Max value for CPML parameter
Rcoef = 0.0010  # Reflection coefficient, used for seismic only

# Courant-Friedrichs-Levy condition
CFL = 1/np.sqrt(3) # 3D CFL but can be changed to 1/np.sqrt(2) for 2D. 

# -------------------------- Function Definitions ------------------------------
def read_dat(
        fn: str, 
        channel: str, 
        domain, 
        single: bool =False
    ) -> np.ndarray:
    """
    Read data from an unformatted Fortran binary file and return it as a numpy 
    array. This function supports reading data in single or double precision, 
    and can handle complex data.
    
    The domain parameter is expected to be an object with attributes 
        `dim`, `nx`, `ny`, and `nz`.
    These are used to determine the shape of the returned numpy array.
    
    :param fn: The filename of the data file to read.
    :type fn: str
    :param channel: Specifies the data channel ('Ex', 'Ey', 'Ez', etc.) 
        to read.
    :type channel: str
    :param domain: The domain object that contains dimensions (dim, nx, ny, nz) 
        of the data.
    :type domain: Domain
    :param single: A flag indicating if the data should be read in 
        single precision. Defaults to False (double precision).
    :type single: bool
    :return: The data read from the file, reshaped according to the 
        domain dimensions.
    :rtype: np.ndarray
    
    """
    
    if single:
        dtype = np.float32 
    else:
        dtype = np.float64 
    
    if domain.dim == 2.5:
        NX = domain.nz
        NY = domain.ny
        NZ = domain.nx
    else:
        NZ = domain.nx #! This should be nx and the next one should be nz
        NX = domain.nz
    #

    with FortranFile(fn, 'r') as f:
        if single:
            dat = f.read_reals(dtype = np.float32)
            # dat = np.fromfile(fn, dtype = np.float32)
        else:
            dat = f.read_reals(dtype = np.float64)
            # dat = np.fromfile(fn, dtype = np.float64)

    if domain.dim == 2.5:
        dat = dat.reshape(NX, NY, NZ)
    else:
        dat = dat.reshape(NX, NZ)
    
    f.close()
    return(dat)

# =============================================================================
# ============================== Useful Functions =============================
# =============================================================================
# ------------------------------------------------------------------------------
# After computing tensor coefficients we want to append them to the given text
# file. In order to do this, we need to create a new file then move it to the
# original file

def complex2str(complex_vector: np.ndarray) -> str:
    """
    Convert a numpy array of complex numbers into a string representation, 
    with the first element representing an identifier and the subsequent 
    elements representing complex numbers.

    This function extracts the real and imaginary parts of each complex number 
    in the input vector, converts them to strings, and concatenates them with 
    'j' as the delimiter. The first real number of the input vector is treated 
    as an identifier (ID) and is converted separately.

    :param complex_vector: A numpy array of complex numbers where the first 
        element is an ID.
    :type complex_vector: np.ndarray
    :return: A string representation of the ID followed by the complex numbers 
        in the vector.
    :rtype: str

    The returned string format is 'ID,real1jimag1,real2jimag2,...' where ID is 
    the integer ID, and realN and imagN are the real and imaginary parts of the
    Nth complex number, respectively.
    """
    reals = complex_vector.real[1:]
    comps = complex_vector.imag[1:]
    
    reals = [format(x, '.2e') for x in reals]
    comps = [format(x, '.2e') for x in comps]
    
    id = complex_vector.real[0].astype(int)
    m = len(reals)
    complex_string_vector = np.zeros([m], dtype = object)
    for ind in range(m):
        complex_string_vector[ind] = str(reals[ind]) + 'j' + str(comps[ind])
    
    return( str(id) + ',' + ','.join(complex_string_vector))

def str2complex(strsplit: list) -> np.ndarray:
    """
    Convert a list of strings into a numpy array of complex numbers, where the
    first string is treated as an identifier (ID) and the subsequent strings 
    represent complex numbers.

    The function splits each string by 'j' to separate the real and imaginary 
    parts of each complex number, then constructs a numpy array of complex 
    numbers from these parts. The first element of the array is the ID.

    :param strsplit: A list of strings, where the first string is an ID and the 
        remaining strings represent complex numbers in 'realjimag' format.
    :type strsplit: list
    :return: A numpy array of complex numbers with the ID as the first element.
    :rtype: np.ndarray

    Note that the input list should follow the format ['ID', 'real1jimag1', 
    'real2jimag2', ...], where 'ID' is an integer identifier, and 'realNjimagN' 
    are the real and imaginary parts of the Nth complex number.
    """
    # strsplit = line.split(',')
    id = int(strsplit[0])
    m = len(strsplit)
    complex_vector = np.zeros([m], dtype = complex)
    complex_vector[0] = id
    for ind in range(1,m):
        realvalue, complexvalue = strsplit[ind].split('j')
        complex_vector[ind] = complex(float(realvalue), float(complexvalue))
    
    return(complex_vector)
        

# ==============================================================================
# ========================= Read/Assign Project File ===========================
# ==============================================================================
# ------------------------------------------------------------------------------
def loadproject(
        project_file: str, 
        domain, material, seismic, electromag
    ):
    """
    Loads project settings from a project file and assigns values to domain, 
    material, seismic, and electromagnetic objects.

    This function parses a specified project file, reading configurations and parameters
    for the domain, materials, seismic modeling, and electromagnetic properties. These
    parameters are then assigned to the provided objects.

    :param project_file: The path to the project configuration file.
    :type project_file: str
    :param domain: An object to hold domain-specific parameters and configurations.
    :param material: An object to hold material-specific parameters and configurations.
    :param seismic: An object to hold seismic modeling parameters and configurations.
    :param electromag: An object to hold electromagnetic modeling parameters and configurations.
    :return: A tuple containing the updated domain, material, seismic, and electromagnetic objects.
    :rtype: tuple

    The function assumes the project file contains specific prefixes for different
    types of configurations (e.g., 'I' for images, 'D' for domain configurations, 'S' for
    seismic modeling parameters, etc.).
    """
    stiffness_columns = np.array([
            "c11", "c12", "c13", "c14", "c15", "c16",
            "c22", "c23", "c24", "c25", "c26",
            "c33", "c34", "c35", "c36",
            "c44", "c45", "c46", 
            "c55", "c56",
            "c66", "rho"
    ])
    attenuation_columns = np.array([
        'gamma_x', 'gamma_y', 'gamma_z', 'gamma_yz', 'gamma_xz', 'gamma_xy', 'reference_frequency'
    ])
    conductivity_columns = np.array([
            "s11", "s12", "s13", 
            "s22", "s23", 
            "s33"
    ])
    permittivity_columns = np.array([
            "e11", "e12", "e13", 
            "e22", "e23", 
            "e33"
    ])
    material_columns = np.array([
        'id', 'Name', 'RGB', 'Temperature', 'Density', 'Porosity', 
        'Water_Content', 'is_anisotropic', 'ANGfiles'
    ])
    # domain, material, seismic, electromag are the objects that we will assign
    # values to
    data = readwrite_json(project_file)
    
    # ----------------------------- Domain Values ------------------------------
    (
        domain.dim, domain.nx, domain.ny, domain.nz, 
        domain.dx, domain.dy, domain.dz, domain.cpml, 
        domain.nmats, domain.imfile 
    ) = list(data['Domain'].values())
    im,__ = image2int(data['Domain']['image_file'])
    domain.geometry = im.transpose().astype(int)
    
    # ---------------------------- Seismic Values ------------------------------
    (
        seismic.dt, seismic.time_steps,
        seismic.x, seismic.y, seismic.z, 
        seismic.xind, seismic.yind, seismic.zind,
        seismic.f0, seismic.theta, seismic.phi, 
        seismic.source_amplitude, seismic.source_type    
    ) = list(data['Seismic']['Source'].values())
    
    (
        electromag.dt, electromag.time_steps,
        electromag.x, electromag.y, electromag.z,
        electromag.xind, electromag.yind, electromag.zind,
        electromag.f0, electromag.theta, electromag.phi,
        electromag.source_amplitude, electromag.source_type    
    ) = list(data['Electromagnetic']['Source'].values())

    # --------------------------- Indexed Values -------------------------------    
    # preallocate
    
    seismic.stiffness_coefficients = pd.DataFrame(
        np.zeros([domain.nmats, 22]), columns = stiffness_columns
    )
    seismic.attenuation_coefficients = pd.DataFrame(
        np.zeros([domain.nmats,7]), columns = attenuation_columns
    )
    electromag.permittivity_coefficients = pd.DataFrame(
        np.zeros([domain.nmats,6]), columns = permittivity_columns
    )
    electromag.conductivity_coefficients = pd.DataFrame(
        np.zeros([domain.nmats,6]), columns = conductivity_columns
    )
    
    seismic.fref = np.zeros([domain.nmats])
    material.material_list =pd.DataFrame(
        np.zeros([domain.nmats, 9], dtype = object), columns = material_columns
    )
    for ind in range(domain.nmats):
        # First do seismic
        coefs = list(data['Seismic']['Stiffness_Coefficients'][ind].values())
        seismic.stiffness_coefficients.loc[ind] = np.array(coefs[1:])
        coefs = list(data['Seismic']['Attenuation'][ind].values())
        seismic.attenuation_coefficients.loc[ind] = np.array(coefs[1:])
        # Now do electromag
        coefs = list(data['Electromagnetic']['Permittivity_Coefficients'][ind].values())
        electromag.permittivity_coefficients.loc[ind] = np.array(coefs[1:])
        coefs = list(data['Electromagnetic']['Conductivity_Coefficients'][ind].values())
        electromag.conductivity_coefficients.loc[ind] = np.array(coefs[1:])
        # Material values
        vals = list(data['Materials'][ind].values())
        material.material_list.loc[ind] = vals
    
    seismic.project_file = project_file 
    electromag.project_file = project_file
    
    seismic.is_seismic = True 
    electromag.is_seismic = False
    material.sort_material_list() 
    return domain, material, seismic, electromag

# ------------------------------------------------------------------------------
def airsurf(material, domain, N: int = 2) -> np.ndarray:
    """
    Generates a gradient matrix to simulate the air surface in a domain, based 
    on material properties.

    This function calculates a gradient matrix representing the transition 
    between air and non-air materials within a domain. It is particularly 
    useful for creating a surface with variable density or other properties 
    near the air-material interface.

    :param material: The material object containing material lists and 
        properties.
    :type material: Material
    :param domain: The domain object containing the geometry of the simulation 
        area.
    :type domain: Domain
    :param N: The number of gradational steps in the air surface simulation, 
        defaults to 2.
    :type N: int
    :return: A gradient matrix representing the air surface within the domain.
    :rtype: np.ndarray

    The function uses the domain's geometry to determine air interfaces and 
    applies a gradational approach to simulate the air surface effect. The 
    gradient matrix can be used to adjust physical properties like density near 
    the surface.
    """
    # This can be generalized a little better, but for now...
    try:
        airnum = np.where(material.material == 'air')[0][0].astype(int)
    except:
        airnum = False
    
    if airnum:
        gradmatrix = (domain.geometry != airnum).astype(int)
        # Take the gradient in both directions
        gradz = np.diff(gradmatrix, axis = 0)
        gradx = np.diff(gradmatrix, axis = 1)

        # For grady we will append a column of zeros at the beginning so that the value
        # 1 is located at the interface but on the non-air side
        gradzpos = np.row_stack([np.zeros([gradz.shape[1] ]),gradz])
        # For gradx we will append a row of zeros at the beginning
        gradxpos = np.column_stack([np.zeros([gradx.shape[0] ]),gradx])
        # -1 also means that there is an air interface. We will need to append to the
        # end of the array, then we can just flip the sign
        gradzneg = (np.row_stack( [gradz, np.zeros([gradz.shape[1]]) ] ) )
        gradxneg = (np.column_stack( [gradx, np.zeros([gradx.shape[0]]) ] ) )

        # At the surface we want to have 15% of density
        grad = np.zeros( [gradx.shape[0], gradz.shape[1], N] )
        grad[:,:,0] = gradzpos + gradxpos - gradzneg - gradxneg
        grad[ grad[:,:,0]>0, 0] = 0.15

        # We will make the change gradational by splitting the difference each step
        # For instance, 1 - 0.15 = 0.85, so the next percentage will be
        # 0.85/2 + 0.15 and so on
        pct = np.zeros([N])
        pct[0] = 0.15

        for ind in range(1, N):
            pct[ind] = pct[ind-1] + (1-pct[ind-1])/2
            gradzpos = np.row_stack( [np.zeros([gradz.shape[1] ]),gradzpos] )[:-1,:]
            gradxpos = np.column_stack( [np.zeros( [ gradx.shape[0] ] ),gradxpos ])[:,:-1]

            gradzneg = (np.row_stack( [ gradzneg, np.zeros( [gradz.shape[1] ]) ] )[1:,:])
            gradxneg = (np.column_stack( [ gradxneg, np.zeros( [gradx.shape[0] ]) ] )[:,1:])
            grad[:,:, ind] = gradzpos + gradxpos - gradzneg - gradxneg
            grad[ grad[:,:,ind] > 0, ind] = pct[ind]

        gradcomp = np.zeros( [grad.shape[0], grad.shape[1] ])
        for ind in range(N-1, -1, -1):
            gradcomp[ grad[:,:,ind] > 0] = pct[ind]

        gradcomp[ gradcomp == 0] = 1
    else:
        gradcomp = np.ones([int(domain.nx), int(domain.nz) ])
    
    f = FortranFile('density_gradient.dat', 'w')
    f.write_record(gradcomp)
    f.close()
        
    return(gradcomp)

def cpmlcompute(
        modelclass, 
        domain, 
        direction: str, 
        half: bool = False, 
    ) -> None:
    """
    Computes CPML parameters for a given direction and updates model/domain.

    :param modelclass: The model class instance to update.
    :param domain: The domain class instance to update.
    :param direction: Direction to compute CPML ('x', 'y', or 'z').
    :param half: Flag to compute half CPML parameters. Defaults to False.
    :type modelclass: Model
    :type domain: Domain
    :type direction: str
    :type half: bool
    """
    global NP, NPA, sig_opt_scalar, kappa_max

    # For 2D models, we don't need to compute the cpml in the y-direction
    if domain.dim == 2 and direction == 'y':
        return 
    
    nx = domain.nx + 2*domain.cpml
    nz = domain.nz + 2*domain.cpml
    if domain.dim == 2.5:
        ny = domain.ny + 2*domain.cpml
        deltamin = np.min([domain.dx, domain.dy, domain.dz]) 
    else:
        deltamin = np.min([domain.dx, domain.dz]) 

    # Allocate space
    if direction == 'x':
        N = int(nx)
        dx = float(domain.dx)
    elif direction == 'y':
        N = int(ny)
        dx = float(domain.dy) 
    else:
        N = int(nz)
        dx = float(domain.dz)
    
    # -----------------------------------------------------------------------------
    # Compute the distance along the absorbing boundary relative to the end of the 
    # original model space. 
    dist = dx * np.arange(0, domain.cpml)
    if half:
        dist = dist + dx/2
    
    dist = dx*domain.cpml - dist
    dist = dist/(dx*domain.cpml)

    # Compute the maximum sigma, and alpha values for the CPML.  
    if modelclass.is_seismic:
        alpha_max = np.pi*modelclass.f0
        quasi_cp_max = 0.7 * deltamin / (2.0 * modelclass.dt)
        sig_max = - np.log(Rcoef) * (NP+1) * quasi_cp_max / (2.0 * domain.cpml )
         # This seems to work well even at higher frequencies
        sigma, kappa, alpha, acoeff, bcoeff = cpml_parameters(
            sig_max, alpha_max, kappa_max, dist, N, modelclass.dt, is_seismic = True
        )
    else:
        # We will use the maximum permittivity coefficient and assume that the 
        # magnetic permeability is 1. We can use a different value
        sig_max = sig_opt_scalar * \
            ((NP + 1) / (deltamin * ((mu0/eps0)**0.5) ) ) #!!! We need to multiply eps0 by the relative permattivity for a better estimate
        alpha_max = alpha_max_scalar * 2 * np.pi * eps0 * modelclass.f0 
        sigma, kappa, alpha, acoeff, bcoeff = cpml_parameters(
            sig_max, alpha_max, kappa_max, dist, N, modelclass.dt
        )

    # Save the results to a fortran binary
    if half:
        sigma.tofile('sigma' + direction + '_half_cpml.dat')
        kappa.tofile('kappa' + direction + '_half_cpml.dat')
        alpha.tofile('alpha' + direction + '_half_cpml.dat')
        acoeff.tofile('acoef' + direction + '_half_cpml.dat')
        bcoeff.tofile('bcoef' + direction + '_half_cpml.dat')
    else:
        sigma.tofile('sigma' + direction + '_cpml.dat')
        kappa.tofile('kappa' + direction + '_cpml.dat')
        alpha.tofile('alpha' + direction + '_cpml.dat')
        acoeff.tofile('acoef' + direction + '_cpml.dat')
        bcoeff.tofile('bcoef' + direction + '_cpml.dat')
        
# -----------------------------------------------------------------------------
def cpml_parameters(
        sig_max: float, 
        alpha_max: float, 
        kappa_max: float, 
        distance, 
        N: int, 
        dt: float,
        is_seismic = False,
    ):
    """
    """
    kappa = np.ones([N])
    alpha = np.zeros([N])
    sigma = np.zeros([N])
    acoeff = np.zeros([N])
    bcoeff = np.zeros([N])

    # Compute in the x, and z directions
    for ind in range(0, len(distance)):
        # From 0
        sigma[ind] = sig_max * (distance[ind]**NP)
        kappa[ind] = 1.0 + (kappa_max - 1.0) * distance[ind]**NP
        alpha[ind] = alpha_max * (1 - distance[ind])**NPA
        # From the end 
        sigma[-(ind+1)] = sig_max*distance[ind]**NP
        kappa[-(ind+1)] = 1 + (kappa_max - 1) * distance[ind]**NP
        alpha[-(ind+1)] = alpha_max * (1 - distance[ind])**NPA
        
        if is_seismic:
            bcoeff[ind] = np.exp( - (sigma[ind] / kappa[ind] + alpha[ind]) * dt)
            bcoeff[-(ind+1)] = np.exp(- (sigma[-(ind+1)] / kappa[-(ind+1)] + alpha[-(ind+1)]) * dt)
        else:
            bcoeff[ind] = np.exp( - (sigma[ind] / kappa[ind] + alpha[ind]) * (dt/eps0) )
            bcoeff[-(ind+1)] = np.exp(- (sigma[-(ind+1)] / kappa[-(ind+1)] + alpha[-(ind+1)]) * (dt/eps0) )
            

    # Compute the a-coefficients 
    alpha[np.where(alpha < 0.0)] = 0.0
    indices = np.where(np.abs(sigma) > 1.0e-6)
    acoeff[indices] = sigma[indices] * (bcoeff[indices] - 1) / \
            (kappa[indices] * sigma[indices] + kappa[indices] * alpha[indices] )

    return sigma, kappa, alpha, acoeff, bcoeff

# -----------------------------------------------------------------------------
def ecpml_parameters2(domain, modelclass, dist, direction):
    """
    Compute the boundary values for variable permittivity and conductivity values at
    the boundary. 
    """
    global mu0, NP, NPA

    nx = domain.nx + 2 * domain.cpml 
    nz = domain.nz + 2 * domain.cpml 
    
    if direction == 'x':
        perm_ind = 1
        cond_ind = 7
        nx -= 1
    if direction == 'z':
        perm_ind = 6
        cond_ind = 12
        nz -= 1
    
    sigma = np.zeros([nx, nz])
    alpha = sigma.copy() 
    kappa = np.ones([nx, nz])
    acoeff = sigma.copy() 
    bcoeff = sigma.copy() 

    m,n = domain.geometry.shape

    for ii in range(0, domain.cpml):
        for jj in range(domain.cpml, nz-domain.cpml):
            perm1 = modelclass.tensor_coefficients_original[
                domain.geometry[0,jj-domain.cpml]
            ][perm_ind]
            cond1 = modelclass.tensor_coefficients_original[domain.geometry[0,jj-domain.cpml]][cond_ind].real
            perm2 = modelclass.tensor_coefficients_original[domain.geometry[-1,jj-domain.cpml]][perm_ind]
            cond2 = modelclass.tensor_coefficients_original[domain.geometry[-1,jj-domain.cpml]][cond_ind].real
            if cond1 > 0:
                sig_max1 = cond1
            elif perm1.imag > 0 and cond1 == 0:
                sig_max1 = eps0 * modelclass.f0 * np.sqrt(perm1).imag
            else:
                sig_max1 =  (NP+1) / (150.0 * np.pi * domain.dx * np.sqrt(perm1.real) )
            
            if cond2 > 0:
                sig_max2 = cond2
            elif perm2.imag > 0 and cond2 == 0:
                sig_max2 = eps0 * modelclass.f0 * np.sqrt(perm2).imag
            else:
                sig_max2 =  (NP+1) / (150.0 * np.pi * domain.dx * np.sqrt(perm2.real) )
            
            # In places where there is no complex permittivity, we need to assign alpha to a small non-zero value
            alpha_max1 = 2 * np.pi * modelclass.f0 * eps0 * perm1.real / 2
            alpha_max2 = 2 * np.pi * modelclass.f0 * eps0 * perm1.real / 2
            kappa_max1 = 1 + (alpha_max1 / modelclass.f0) * np.sqrt(mu0 * perm1.real)
            kappa_max2 = 1 + (alpha_max2 / modelclass.f0) * np.sqrt(mu0 * perm2.real)

            #!!! This is for testing 
            sig_max1 = 0.8 * ( (NP+1) / (domain.dx * (mu0/eps0 )**0.5) )
            sig_max2 = 0.8 * ( (NP+1) / (domain.dx * (mu0/eps0 )**0.5) )
            
            alpha_max1 = 2 * np.pi * eps0 * modelclass.f0
            alpha_max2 = 2 * np.pi * eps0 * modelclass.f0
            kappa_max1 = 5.0
            kappa_max2 = 5.0
            #!!!
            sigma[ii,jj] = sig_opt_scalar * sig_max1 * (dist[ii]**NP)
            sigma[-(ii+1),jj] = sig_opt_scalar * sig_max2 * ( (1-dist[ii] )**NPA)
            alpha[ii,jj] = alpha_max1 * (1 - dist[ii])**NPA
            alpha[-(ii+1),jj] = alpha_max2 * (1 - dist[ii])**NPA
            kappa[ii,jj] = 1.0 + (kappa_max1 - 1.0) * dist[ii]**NP
            kappa[-(ii+1),jj] = 1 + (kappa_max2 - 1) * dist[ii]**NP

            bcoeff[ii,jj] = np.exp( 
                -( (sigma[ii,jj] / kappa[ii,jj]) + alpha[ii,jj] ) * (modelclass.dt/eps0)
            )
            bcoeff[-(ii+1),jj] = np.exp( 
                -( (sigma[-(ii+1),jj] / kappa[-(ii+1),jj] + alpha[-(ii+1),jj]) * (modelclass.dt/eps0) )
            )

        for jj in range(domain.cpml, nx-domain.cpml):
            perm1 = modelclass.tensor_coefficients_original[domain.geometry[jj-domain.cpml,0]][perm_ind]
            cond1 = modelclass.tensor_coefficients_original[domain.geometry[jj-domain.cpml,0]][cond_ind].real
            perm2 = modelclass.tensor_coefficients_original[domain.geometry[jj-domain.cpml,-1]][perm_ind]
            cond2 = modelclass.tensor_coefficients_original[domain.geometry[jj-domain.cpml,-1]][cond_ind].real

            # In places where there is no complex permittivity, we need to assign alpha to a small non-zero value
            if cond1 > 0:
                sig_max1 = cond1
            elif perm1.imag > 0 and cond2 == 0.0:
                sig_max1 = eps0 * modelclass.f0 * np.sqrt(perm1).imag
            else:
                sig_max1 =  (NP+1) / (150.0 * np.pi * domain.dx * np.sqrt(perm1.real) )
            
            if cond2 > 0:
                sig_max2 = cond2
            elif perm2.imag > 0 and cond2 == 0.0:
                sig_max2 = eps0 * modelclass.f0 * np.sqrt(perm2).imag            
            else:
                sig_max2 =  (NP+1) / (150.0 * np.pi * domain.dx * np.sqrt(perm2.real) )
            
            # In places where there is no complex permittivity, we need to assign alpha to a small non-zero value
            if perm1.imag == 0.0:
                alpha_max1 = 0.02
            else:
                alpha_max1 = modelclass.f0 * np.sqrt(perm1).imag / 2

            if perm2.imag == 0.0:
                alpha_max2 = 0.02
            else:
                alpha_max2 = modelclass.f0 * np.sqrt(perm2).imag / 2 

            kappa_max1 = 1 + (alpha_max1 / modelclass.f0) * np.sqrt(mu0 * perm1.real)
            kappa_max2 = 1 + (alpha_max2 / modelclass.f0) * np.sqrt(mu0 * perm2.real) 

            #!!! This is for testing 
            sig_max1 = 0.8 * ( (NP+1) / (domain.dx * (mu0/eps0 )**0.5) )
            sig_max2 = 0.8 * ( (NP+1) / (domain.dx * (mu0/eps0 )**0.5) )
            
            alpha_max1 = 2 * np.pi * eps0 * modelclass.f0
            alpha_max2 = 2 * np.pi * eps0 * modelclass.f0
            kappa_max1 = 5.0
            kappa_max2 = 5.0
            #!!!
            sigma[jj,ii] = sig_opt_scalar * sig_max1 * (dist[ii]**NP)
            sigma[jj,-(ii+1)] = sig_opt_scalar * sig_max2 * ( (1-dist[ii] )**NPA)
            alpha[jj,ii] = alpha_max1 * (1 - dist[ii])**NPA
            alpha[jj,-(ii+1)] = alpha_max2 * (1 - dist[ii])**NPA
            kappa[jj,ii] = 1.0 + (kappa_max - 1.0) * dist[ii]**NP
            kappa[jj,-(ii+1)] = 1 + (kappa_max - 1) * dist[ii]**NP

            bcoeff[jj,ii] = np.exp( 
                -(sigma[jj,ii] / kappa[jj,ii] + alpha[jj,ii]) * modelclass.dt
            )
            bcoeff[jj,-(ii+1)] = np.exp( 
                -(sigma[jj,-(ii+1)] / kappa[jj,-(ii+1)] + alpha[jj,-(ii+1)]) * modelclass.dt
            )

    alpha[np.where(alpha < 0.0)] = 0.0
    indices = np.where(np.abs(sigma) > 1.0e-6)
    acoeff[indices] = sigma[indices] * (bcoeff[indices] - 1.0) / \
        (kappa[indices] * sigma[indices] + kappa[indices]**2 * alpha[indices] )

    return sigma, alpha, kappa, acoeff, bcoeff
 
# -----------------------------------------------------------------------------
def ecpml_parameters3(domain, modelclass, dist, direction):
    """
    """
    global mu0, NP, NPA
    nx = domain.nx + 2 * domain.cpml 
    ny = domain.ny + 2 * domain.cpml
    nz = domain.nz + 2 * domain.cpml
    if direction == 'x':
        perm_ind = 1
        cond_ind = 7
        nx -= 1
    if direction == 'y':
        perm_ind = 4
        cond_ind = 10
        ny -= 1
    if direction == 'z':
        perm_ind = 6
        cond_ind = 12
        nz -= 1
    
    sigma = np.zeros([nx, ny, nz])
    alpha = sigma.copy() 
    kappa = np.ones([nx, ny, nz])
    acoeff = sigma.copy() 
    bcoeff = sigma.copy()

    m,n = domain.geometry.shape

    for ii in range(0, domain.cpml):
        for jj in range(domain.cpml, nz-domain.cpml):
            perm1 = modelclass.tensor_coefficients_original[
                domain.geometry[0,jj-domain.cpml]
            ][perm_ind]
            cond1 = modelclass.tensor_coefficients_original[
                domain.geometry[0,jj-domain.cpml]
            ][cond_ind].real
            perm2 = modelclass.tensor_coefficients_original[
                domain.geometry[-1,jj-domain.cpml]
            ][perm_ind]
            cond2 = modelclass.tensor_coefficients_original[
                domain.geometry[-1,jj-domain.cpml]
            ][cond_ind].real
            if cond1 > 0:
                sig_max1 = cond1
            elif perm1.imag > 0:
                sig_max1 = eps0 * modelclass.f0 * np.sqrt(perm1).imag
            else:
                sig_max1 =  (NP+1) / (150.0 * np.pi * domain.dx * np.sqrt(perm1.real) )
            
            if cond2 > 0:
                sig_max2 = cond2
            elif perm2.imag > 0:
                sig_max2 = eps0 * modelclass.f0 * np.sqrt(perm2).imag
            else:
                sig_max2 =  (NP+1) / (150.0 * np.pi * domain.dx * np.sqrt(perm2.real) )
            
            # In places where there is no complex permittivity, we need to assign alpha to a small non-zero value
            if perm1.imag == 0.0:
                alpha_max1 = 0.02
            else:
                alpha_max1 = modelclass.f0 * np.sqrt(perm1).imag / 2

            if perm2.imag == 0.0:
                alpha_max2 = 0.02
            else:
                alpha_max2 = modelclass.f0 * np.sqrt(perm2).imag / 2 

            kappa_max1 = 1 + (alpha_max1 / modelclass.f0) * np.sqrt(mu0 * perm1.real)
            kappa_max2 = 1 + (alpha_max2 / modelclass.f0) * np.sqrt(mu0 * perm2.real)

            sigma[ii,:,jj] = sig_opt_scalar * sig_max1 * (dist[ii]**NP)
            sigma[-(ii+1),:,jj] = sig_opt_scalar * sig_max2 * ( (1-dist[ii] )**NPA)
            alpha[ii,:,jj] = alpha_max1 * (1 - dist[ii])**NPA
            alpha[-(ii+1),:,jj] = alpha_max2 * (1 - dist[ii])**NPA
            kappa[ii,:,jj] = 1.0 + (kappa_max1 - 1.0) * dist[ii]**NP
            kappa[-(ii+1),:,jj] = 1 + (kappa_max2 - 1) * dist[ii]**NP

            bcoeff[ii,:,jj] = np.exp( 
                -(sigma[ii,:,jj] / kappa[ii,:,jj] + alpha[ii,:,jj]) * modelclass.dt
            )
            bcoeff[-(ii+1),:,jj] = np.exp( 
                -(sigma[-(ii+1),:,jj] / kappa[-(ii+1),:,jj] + alpha[-(ii+1),:,jj]) * modelclass.dt
            )

        for jj in range(domain.cpml, nx-domain.cpml):
            perm1 = modelclass.tensor_coefficients_original[
                domain.geometry[jj-domain.cpml,0]
            ][perm_ind]
            cond1 = modelclass.tensor_coefficients_original[
                domain.geometry[jj-domain.cpml,0]
            ][cond_ind].real
            perm2 = modelclass.tensor_coefficients_original[
                domain.geometry[jj-domain.cpml,-1]
            ][perm_ind]
            cond2 = modelclass.tensor_coefficients_original[
                domain.geometry[jj-domain.cpml,-1]
            ][cond_ind].real
            if cond1 > 0:
                sig_max1 = cond1
            elif perm1.imag > 0:
                sig_max1 = eps0 * modelclass.f0 * np.sqrt(perm1).imag
            else:
                sig_max1 =  (NP+1) / (150.0 * np.pi * dx * np.sqrt(perm1.real) )
            
            if cond2 > 0:
                sig_max2 = cond2
            elif perm2.imag > 0:
                sig_max2 = eps0 * modelclass.f0 * np.sqrt(perm2).imag
            else:
                sig_max2 =  (NP+1) / (150.0 * np.pi * dx * np.sqrt(perm2.real) )

            # In places where there is no complex permittivity, we need to assign alpha to a small non-zero value
            if perm1.imag == 0.0:
                alpha_max1 = 0.02
            else:
                alpha_max1 = modelclass.f0 * np.sqrt(perm1).imag / 2

            if perm2.imag == 0.0:
                alpha_max2 = 0.02
            else:
                alpha_max2 = modelclass.f0 * np.sqrt(perm2).imag / 2 
            
            kappa_max1 = 1 + (alpha_max1 / modelclass.f0) * np.sqrt(mu0 * perm1.real)
            kappa_max2 = 1 + (alpha_max2 / modelclass.f0) * np.sqrt(mu0 * perm2.real) 

            sigma[jj,:,ii] = sig_opt_scalar * sig_max1 * (dist[ii]**NP)
            sigma[jj,:,-(ii+1)] = sig_opt_scalar * sig_max2 * ( (1-dist[ii] )**NPA)
            alpha[jj,:,ii] = alpha_max1 * (1 - dist[ii])**NPA
            alpha[jj,:,-(ii+1)] = alpha_max2 * (1 - dist[ii])**NPA
            kappa[jj,:,ii] = 1.0 + (kappa_max1 - 1.0) * dist[ii]**NP
            kappa[jj,:,-(ii+1)] = 1 + (kappa_max2 - 1) * dist[ii]**NP

            bcoeff[jj,ii] = np.exp( 
                -(sigma[jj,:,ii] / kappa[jj,:,ii] + alpha[jj,:,ii]) * modelclass.dt
            )
            bcoeff[jj,:,-(ii+1)] = np.exp( 
                -(sigma[jj,:,-(ii+1)] / kappa[jj,:,-(ii+1)] + alpha[jj,:,-(ii+1)]) * modelclass.dt
            )

        for jj in range(domain.cpml, nx-domain.cpml):
            for kk in range(domain.cpml, nz-domain.cpml):
                perm = modelclass.tensor_coefficients_original[
                    domain.geometry[jj-domain.cpml,kk-domain.cpml]
                ][perm_ind]
                cond = modelclass.tensor_coefficients_original[
                    domain.geometry[jj-domain.cpml,kk-domain.cpml]
                ][cond_ind].real
                if cond > 0:
                    sig_max = cond
                elif perm.imag > 0:
                    sig_max = eps0 * modelclass.f0 * np.sqrt(perm).imag
                else:
                    sig_max =  (NP+1) / (150.0 * np.pi * dx * np.sqrt(perm.real) )

                # In places where there is no complex permittivity, we need to assign alpha to a small non-zero value
                if perm.imag == 0.0:
                    alpha_max = 0.02
                else:
                    alpha_max = modelclass.f0 * np.sqrt(perm1).imag / 2


                kappa_max = 1 + (alpha_max / modelclass.f0) * np.sqrt(mu0 * perm.real)
                sigma[jj,ii,kk] = sig_opt_scalar * sig_max1 * (dist[ii]**NP)
                sigma[jj,-(ii+1),kk] = sig_opt_scalar * sig_max2 * ( (1-dist[ii] )**NPA)
                alpha[jj,ii,kk] = alpha_max1 * (1 - dist[ii])**NPA
                alpha[jj,-(ii+1),kk] = alpha_max2 * (1 - dist[ii])**NPA
                kappa[jj,ii,kk] = 1.0 + (kappa_max1 - 1.0) * dist[ii]**NP
                kappa[jj,-(ii+1),kk] = 1 + (kappa_max2 - 1) * dist[ii]**NP

                bcoeff[jj,ii,kk] = np.exp( 
                    -(sigma[ii,jj] / kappa[ii,jj] + alpha[ii,jj]) * modelclass.dt
                )
                bcoeff[jj,-(ii+1),kk] = np.exp( 
                    -(sigma[jj,-(ii+1),kk] / kappa[jj,-(ii+1),kk] + alpha[jj,-(ii+1),kk]) * modelclass.dt
                )

    
    alpha[np.where(alpha < 0.0)] = 0.0
    indices = np.where(np.abs(sigma) > 1.0e-6)
    acoeff[indices] = sigma[indices] * (bcoeff[indices] - 1) / \
        (kappa[indices] * sigma[indices] + kappa[indices] * alpha[indices] )

    return sigma, alpha, kappa, acoeff, bcoeff

# ------------------------------------------------------------------------------
def rcxgen(
        rgb: list, 
        domain, material, 
        filename: str = 'receivers.xyz'
    ):
    """
    Creates a receiver list from a specified RGB value found in the model 
    image, outputting to a CSV file.
    
    This function searches the domain's geometry for occurrences of a given RGB 
    value (specified as a list of integers representing red, green, and blue 
    components). It then generates a list of receiver coordinates based on the 
    locations where the RGB value is found. Currently, the function is set up 
    for 2D receivers only, with Y-coordinates set to zero.

    :param rgb: A list of three integers representing the RGB value to search 
        for in the domain's geometry.
    :type rgb: list
    :param domain: An object representing the domain, including its geometry 
        and spatial discretization parameters (dz, dx).
    :param material: An object containing the material list and associated RGB 
        values.
    :param filename: The name of the output CSV file containing the receiver 
        coordinates, defaults to 'receivers.xyz'.
    :type filename: str
    :return: A numpy array of receiver coordinates where each row represents 
        [Z, Y, X] coordinates.
    :rtype: np.ndarray

    The function converts the RGB list into a string format, looks up the 
        corresponding integer value from the
    material list, finds all occurrences in the domain's geometry, and outputs 
        the coordinates to a specified CSV file.
    """
    rgb = '/'.join(np.array(rgb).astype(str))
    rgbint = int(material.material_list[material.rgb == rgb,0])
    z,x = np.where(domain.geometry == rgbint)
    y = np.zeros([len(x)])
    xyz = np.stack([z*domain.dz,y,x*domain.dx], axis = 1)
    df = pd.DataFrame(xyz, columns = ['X', 'Y', 'Z'])
    df.to_csv(filename, index = False)
    return(xyz)


# ============================= Source Functions ===============================
# ------------------------------------------------------------------------------
def coherdt(
        alpha: float, 
        v: float, 
        n: int, 
        dx: float, 
        loc: str = 'left'
    ) -> np.ndarray:
    """
    Calculates a vector of time delays based on the angle of incidence, 
    velocity, and spatial discretization.

    This function determines the time delay for coherent addition of traces 
    based on the propagation angle (alpha), velocity (v), and the 
    discretization in space (dx) along the domain length.

    :param alpha: The angle of incidence in degrees, counter clockwise from the 
        x-plane. Range: (-90, 90).
    :type alpha: float
    :param v: Velocity in meters/second.
    :type v: float
    :param n: The number of spatial points in the domain.
    :type n: int
    :param dx: The spatial discretization in meters.
    :type dx: float
    :param loc: Specifies the reference point's location for time delay 
        calculation, defaults to 'left'.
    :type loc: str
    :return: A numpy array containing the time delays for each spatial point in 
        the domain.
    :rtype: np.ndarray

    The function adjusts angles outside the (-90, 90) range and calculates time 
    delays accordingly.
    """
    # Return a vector of time delays along the length of the domain. The angle
    # alpha, given in degrees between 90 and -90, is counter clockwise from the
    # x plane. Velocity and dx are in meters/second and meters, respectively.
    if alpha < 0:
        self.dz = None
        self.cpml = None
        self.write = None
        self.image_file = None
        self.exit_status = 1

        # Flag to specify whether the all inputs are fulfilled
        self.seismic_model = False
        self.electromag_model = False

    def para_check(self):
        self.exit_s
    if alpha > 90.0 or alpha < -90:
        alpha = 180 - np.sign(alpha) * alpha
    else:
        print('Alpha must be on the interval (-90,90)')
        quit()

    x = np.arange(0, n*dx, dx)
    dt = x * np.sin( np.deg2rad(alpha) ) / v
    return(dt)

# ------------------------------------------------------------------------------
def coherstf(
        timediff: np.ndarray, 
        source_function: np.ndarray, 
        dt: float, 
        m: int, n: int, cpml: int, 
        bottom: bool = False
    ) -> None:
    """
    Applies time shifts to a source time function along the domain and saves it 
    as an m-by-n-by-p matrix.

    :param timediff: A numpy array containing time differences for time 
        shifting the source function.
    :type timediff: np.ndarray
    :param sf: The source time function as a numpy array.
    :type sf: np.ndarray
    :param dt: The time interval between successive samples in the source 
        function.
    :type dt: float
    :param m: Number of indices in the y-direction (height of the domain).
    :type m: int
    :param n: Number of indices in the x-direction (width of the domain).
    :type n: int
    :param cpml: The size of the Convolutional Perfectly Matched Layer (CPML) 
        padding.
    :type cpml: int
    :param bottom: Specifies whether to place the source at the bottom of the 
        domain. Defaults to False.
    :type bottom: bool
    :return: None

    This function time shifts the source time function for each node in the 
    x-direction, considering CPML adjustments.
    
    'topleft' origin is (0,0), 
    'topright' origin is (x_n, 0)
    'bottomleft' origin is (0, y_n), 
    'bottomright' origin is (x_n, y_n)
    """
    p = len(source_function)
    sfmat = np.zeros([m, n, p], order='F') # m indices are in the y-direction
    cpmlmat = np.zeros([2*cpml + m, 2*cpml + n, p])
    sfarray = np.zeros([m, p])  # p indices are in the t-direction
    ndiff = int(np.round(timediff/dt))
    for ind in range(0, m):
        sfarray[ind,:] = np.roll(source_function, timediff[i])
        sfarray[ind,0:timediff[i]] == 0
    if bottom == True:
        sfmat[m,:,:] = sfarray[:]
    else:
        sfmat[0,:,:] = sfarray[:]
    cpmlmat[cpml:,cpml:,:] = sfmat[:,:,:]
    cpmlmat.T.tofile('sourcefunctionmatrix.dat')

# ------------------------------------------------------------------------------
def stfvec2mat(
        source_function: np.ndarray, 
        xind: int, zind: int, 
        m: int, n: int, 
        cpml: int, 
        yind: int = None
    ) -> None:
    """
    Converts a source time function vector to a matrix for 2D or 3D 
    simulations.

    :param sf: The source time function vector.
    :type sf: np.ndarray
    :param xind: The index in the x-direction where the source is located.
    :type xind: int
    :param zind: The index in the z-direction where the source is located.
    :type zind: int
    :param m: The number of indices in the y-direction (height of the domain).
    :type m: int
    :param n: The number of indices in the x-direction (width of the domain).
    :type n: int
    :param cpml: The size of the Convolutional Perfectly Matched Layer (CPML) 
        padding.
    :type cpml: int
    :param yind: The index in the y-direction where the source is located, for 
        3D simulations. Defaults to None for 2D simulations.
    :type yind: Optional[int]
    :return: None

    This function arranges the source time function in the spatial domain and applies CPML padding.
    """
    # If 2D, y = None. So far we only have 2D. 3D will be implemented later
    p = len(st)
    sfmat = np.zeros([m,n,p], order = 'F')
    sfmat[zind, xind,:] = source_function[:]
    cpmlmat = np.zeros([2*cpml + m, 2*cpml + n, p], order = 'F')
    cpmlmat[cpml:,cpml:,:] = sfmat[:,:,:]
    cpmlmat.T.tofile('sourcefunctionmatrix.dat')

# ------------------------------------------------------------------------------
def movingsrc(sf: np.ndarray, txlocs: np.ndarray) -> None:
    """
    Simulates a moving source given a stationary source time function.

    :param sf: The stationary source time function as a numpy array.
    :type sf: np.ndarray
    :param txlocs: A numpy array containing the transmission locations over time.
    :type txlocs: np.ndarray
    :return: None
    
    This function is intended as a placeholder to simulate a moving source by manipulating the stationary source time function according to the provided transmission locations.
    """
    pass  # Implementation yet to be done or specified.

# ----------------------------- Plotting Functions -----------------------------
def indvar(
        model_object, domain
    ) -> tuple[np.ndarray, Optional[np.ndarray], np.ndarray, np.ndarray]:
    """
    Generates spatial and temporal grids based on domain and model object 
    attributes.

    :param model_object: An object containing model-specific parameters, including 
        time steps and time interval (dt).
    :param domain: An object containing domain-specific parameters such as 
        spatial dimensions (nx, ny, nz) and discretizations (dx, dy, dz).
    :return: A tuple of numpy arrays representing the grids in x, y (None if ny 
        is not set), z, and t dimensions.
    :rtype: tuple

    The function calculates the spatial grid points in the x, z, and optionally 
    y dimensions, and temporal grid points based on the provided domain and 
    model parameters. The y-dimension grid is generated only if the domain 
    parameters include 'ny'.
    """
    nx = int(domain.nx[0])
    nz = int(domain.nz[0])
    dx = float(domain.dx[0])
    dz = float(domain.dz[0])
    dt = float(model_object.dt[0])
    nt = int(model_object.time_steps[0])

    x = np.linspace(0, dx * (nx - 1), nx)
    z = np.linspace(0, dz * (nz - 1), nz)
    t = np.linspace(0, dt * (nt - 1), nt)
    try:
        y = np.linspace(
            0, float(domain.dy[0]) * (int(domain.ny) - 1), int(domain.ny)
        )
    except:
        y = None

    return(x,y,z,t)

# ------------------------------------------------------------------------------
def parameter_profile_1d(
        domain, material, model, 
        indice: int, parameter: str = 'velocity'
    ):
    """
    Get the values of a 1D profile for a parameter (i.e. velocity in the x,y,z
    directions with respect to z). The available options are:
        'velocity', 'temperature', 'density', 'lwc', 'conductivity' (em only)

    Defaults to 'velocity' which returns 4 m-by-1 arrays. All other options will 
    return 2 m-by-1 arrays. 

    :param domain: The domain object of the model
    :type domain: Domain
    :param material: The material object of the model
    type material: Material
    :param model: The model object
    :type model: Model
    :param indice: Specify the x-indice to pull the 1D profile
    :type indice: int
    :param parameter: Specify which parameter to use
    :type parameter: str
    """

    profile = domain.geometry[indice,:]
    n = len(profile)
    if material.material is None:
        material.sort_material_list() 
    
    attribute_map = {
        'temperature': 'temp',
        'density': 'rho',
        'lwc': 'lwc',
        'porosity': 'porosity'
    }

     
    if parameter == 'velocity':
        if model.is_seismic:
            output_values = np.empty([n, 4])
            tensor_coefficients = model.tensor_coefficients[:,:-1].astype(float)
            # T = np.zeros([6,6])
            for ind in range(n):
                T = tensor_coefficients[profile[ind],:]
                rho = material.rho[profile[ind]]
                output_values[ind,:] = mf.tensor2velocities(
                    T, rho, seismic = True
                )
        else:
            output_values = np.empty([n, 3])
            if parameter == 'conductivity':
                tensor_coefficients = model.tensor_coefficients[:,7].real.astype(float)
                for ind in range(n):
                    T = tensor_coefficients[profile[ind],:] 
                    output_values[ind,:] = T[np.array([0,3,5])]
            else:
                tensor_coefficients = model.tensor_coefficients[:,1:7].real.astype(float)
                # T = np.zeros([3,3])  
                for ind in range(n):
                    T = tensor_coefficients[profile[ind],:]
                    output_values[ind,:] = mf.tensor2velocities(
                        T, seismic = False
                    )
    else:
        attribute_name = attribute_map.get(parameter)
        vals = getattr(material, attribute_name, None)
        output_values = np.empty([n])
        for ind in range(n):
            output_values[ind] = vals[profile[ind]]

    z = np.arange(0,n) * domain.dz
    return z, output_values 

# ============================ Processing Functions ============================
def agc(ts: np.ndarray, k: int, agctype: str) -> np.ndarray:
    """
    Applies auto-gain control (AGC) normalization to a time series using a 
    specified window length and type. This function normalizes the input time 
    series using a running window approach, with the normalization type 
    determined by `agctype`. It supports standard deviation ("std"), mean 
    ("mean"), or median ("median") based normalization. The function modifies 
    the series to have a uniform amplitude characteristic over time.

    :param ts: The input time series as a numpy array.
    :type ts: np.ndarray
    :param k: The length of the window used for normalization.
    :type k: int
    :param agctype: The type of normalization to apply, choices are "std", 
        "mean", or "median".
    :type agctype: str
    :return: The time series after applying AGC normalization.
    :rtype: np.ndarray
    """
    n = len(ts)

    k2 = int(k/2) # This will floor the number if k is even; (k-1)/2
    if np.mod(k, 2) == 0: # even numbers need will not have a centered window
        k = int( k + 1)

    stat = np.ones([n])
    # Normalize
    if agctype == "std":
        for i in range(0, k2):
            stat[i] = np.std( abs( ts[0:(i+k2)] ) )
            stat[ (n-k2+i) ] = np.std( abs(ts[ (n-2*k2+i):n] ) )
        for i in range(k2,n-k2):
            stat[i] = np.std( abs( ts[ (i-k2):(i+k2) ] ) )
    elif agctype == "mean":
        for i in range(0, k2):
            stat[i] = np.mean( abs( ts[0:(i+k2)] ) )
            stat[ (n-k2+i) ] = np.mean( abs(ts[ (n-2*k2+i):n] ) )
        for i in range(k2,n-k2):
            stat[i] = np.mean( abs( ts[ (i-k2):(i+k2) ] ) )
    else:
        for i in range(0, k2):
            stat[i] = np.std( ts[i:(i+k2)] )
            stat[ (n-k2+i) ] = np.std( ts[ (n-2*k2+i):n] )
        for i in range(k2,n-k2):
            stat[i] = np.std( ts[ (i-k2):(i+k2) ] )

    stat[stat == 0] = 1
    ts = ts/stat
    return ts

# ------------------------------------------------------------------------------
def correct_geometric_spreading(
        time_series: np.ndarray, 
        distances: np.ndarray, 
        exponent: float = 1.0
    ) -> np.ndarray:
    """
    Corrects a seismic time series for geometric spreading.

    :param time_series: The seismic time series data. Should be a numpy array of 
        shape (n_samples, n_traces).
    :type time_series: np.ndarray
    :param distances: The distances from the source to each trace. Should be a 
        numpy array of shape (n_traces,).
    :type distances: np.ndarray
    :param exponent: The exponent used in the geometric spreading correction. 
        Defaults to 1.0.
    :type exponent: float, optional

    :raises ValueError: 
        If the length of distances does not match the number of traces.

    :return: 
        The corrected seismic time series. Will be a numpy array of shape (n_samples, n_traces).
    :rtype: np.ndarray

    **Example usage:**

    .. code-block:: python

        # Example seismic time series (n_samples x n_traces)
        time_series = np.random.rand(1000, 10)  # Replace with your actual data

        # Example distances (n_traces,)
        distances = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])  # Replace with your actual distances

        # Apply geometric spreading correction
        corrected_time_series = correct_geometric_spreading(time_series, distances)

        print(corrected_time_series)
    """
    if len(time_series.shape) == 1:
        time_series = time_series[:, np.newaxis]
        
    n_samples, n_traces = time_series.shape
    if distances.shape[0] != n_traces:
        raise ValueError(
            "The length of distances must match the number of traces."
        )

    # Calculate the geometric spreading correction factor for each trace
    correction_factors = distances ** exponent

    # Apply the correction to each trace
    corrected_time_series = time_series / correction_factors

    return corrected_time_series

# ------------------------------------------------------------------------------
def exponential_gain(time_series: np.ndarray, alpha: float) -> np.ndarray:
    """
    Applies an exponential gain to a time series to correct for attenuation.

    :param time_series: The original time series data. Should be a numpy array 
        of shape (n_samples,).
    :type time_series: np.ndarray
    :param alpha: 
        The attenuation coefficient.
    :type alpha: float

    :return: The corrected time series.
    :rtype: np.ndarray
    """
    n_samples = time_series.shape[0]
    time = np.arange(n_samples)
    
    # Calculate exponential gain
    gain = np.exp(alpha * time)
    
    # Apply the exponential gain to the time series
    corrected_time_series = time_series * gain

    return corrected_time_series
