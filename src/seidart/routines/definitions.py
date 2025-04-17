import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.animation as anim
import os.path
from typing import Optional
from subprocess import call
from scipy.io import FortranFile
from scipy.signal import hilbert, correlate
from scipy.signal.windows import tukey

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
    'rotate_to_zrt',
    'rotate_to_qlt',
    'compute_envelope',
    'polarization_analysis',
    'agc',
    'correct_geometric_spreading',
    'exponential_gain',
    'parameter_profile_1d',
    'plot_3c',
    'plot_hodogram',
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
# sig_opt_scalar = 1.2
# alpha_max_scalar = 1.0
# NP = 2  # Numerical parameter for CPML
# NPA = 2  # Additional numerical parameter for CPML
# kappa_max = 5  # Max value for CPML parameter
# Rcoef = 0.0010  # Reflection coefficient, used for seismic only

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
        domain.nmats, 
        domain.sig_opt_scalar, domain.alpha_max_scalar, domain.kappa_max,
        domain.NP, domain.NPA, domain.Rcoef,
        domain.imfile 
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
        alpha_max = domain.alpha_max_scalar * np.pi*modelclass.f0
        quasi_cp_max = 0.7 * deltamin / (2.0 * modelclass.dt)
        sig_max = - np.log(domain.Rcoef) * (domain.NP+1) * quasi_cp_max / (2.0 * domain.cpml )
         # This seems to work well even at higher frequencies
        sigma, kappa, alpha, acoeff, bcoeff = cpml_parameters(
            sig_max, alpha_max, domain.kappa_max, 
            dist, N, domain.NP, domain.NPA, modelclass.dt, is_seismic = True
        )
    else:
        # We will use the maximum permittivity coefficient and assume that the 
        # magnetic permeability is 1. We can use a different value
        sig_max = domain.sig_opt_scalar * \
            ((domain.NP + 1) / (deltamin * ((mu0/eps0)**0.5) ) ) #!!! We need to multiply eps0 by the relative permattivity for a better estimate
        alpha_max = domain.alpha_max_scalar * 2 * np.pi * eps0 * modelclass.f0 
        sigma, kappa, alpha, acoeff, bcoeff = cpml_parameters(
            sig_max, alpha_max, domain.kappa_max, 
            dist, N, domain.NP, domain.NPA, modelclass.dt
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
        NP: int,
        NPA: int, 
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

# ============================= Plotting Functions =============================
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

# ------------------------------------------------------------------------------
def plot_3c(time, 
            data1, 
            data2=None, 
            data3=None, 
            data_label1='Data1',
            data_label2='Data2',
            data_label3='Data3',
            source_wavelet=None, 
            arrivals1=None, 
            arrivals2=None,
            arrivals3=None,
            component_labels=None,
            arrival_labels1=None, 
            arrival_labels2=None,
            arrival_labels3=None,
            arrival_colors1=None, 
            arrival_colors2=None,
            arrival_colors3=None,
            color_data1='black', 
            color_data2='blue',
            color_data3='green',
            envelope_color1='lightgray', 
            envelope_color2='lightgray',
            envelope_color3='lightgray',
            plot_envelope=True,
            figsize=(10, 8),
            show_legend=True,
            xlims = None,
            yaxis_rotation = 90,
            envelope_pad_width=50):
    """
    Plot 3-component seismograms with vertical dashed lines for expected arrival times.
    
    For each component channel, the function will:
      1. Optionally cross-correlate the channel with a provided source wavelet.
      2. Compute a power envelope (via Hilbert transform with padding to reduce edge effects).
      3. Plot the envelope as a shaded region and overlay the cross-correlated seismogram.
      
    The function supports up to three datasets:
      - data1 is required.
      - data2 and data3 are optional.
      
    Separate arrival times, labels, and line colors can be specified for each dataset.
    
    Parameters
    ----------
    time : numpy.ndarray
        1D array of time values.
    data1 : numpy.ndarray
        2D array (n_samples x 3) for the first set of seismogram components.
    data2 : numpy.ndarray, optional
        2D array (n_samples x 3) for the second set of seismogram components.
    data3 : numpy.ndarray, optional
        2D array (n_samples x 3) for the third set of seismogram components.
    source_wavelet : numpy.ndarray, optional
        1D array containing the source wavelet. If provided, cross correlation is performed.
    arrivals1 : dict or list/tuple, optional
        Arrival times for data1. If a list/tuple, must have 3 elements in order [P, S1, S2].
    arrivals2 : dict or list/tuple, optional
        Arrival times for data2. If a list/tuple, must have 3 elements in order [P, S1, S2].
    arrivals3 : dict or list/tuple, optional
        Arrival times for data3. If a list/tuple, must have 3 elements in order [P, S1, S2].
    component_labels : list of str, optional
        Labels for the 3 components. Default is ['Vertical', 'Radial', 'Transverse'].
    arrival_labels1 : dict, optional
        Arrival labels for data1. Default is {'P': 'P1', 'S1': 'S1_1', 'S2': 'S2_1'}.
    arrival_labels2 : dict, optional
        Arrival labels for data2. Default is {'P': 'P2', 'S1': 'S1_2', 'S2': 'S2_2'}.
    arrival_labels3 : dict, optional
        Arrival labels for data3. Default is {'P': 'P3', 'S1': 'S1_3', 'S2': 'S2_3'}.
    arrival_colors1 : dict, optional
        Colors for vertical dashed lines for data1 arrivals. Default is {'P': 'red', 'S1': 'red', 'S2': 'red'}.
    arrival_colors2 : dict, optional
        Colors for vertical dashed lines for data2 arrivals. Default is {'P': 'blue', 'S1': 'blue', 'S2': 'blue'}.
    arrival_colors3 : dict, optional
        Colors for vertical dashed lines for data3 arrivals. Default is {'P': 'green', 'S1': 'green', 'S2': 'green'}.
    color_data1 : str, optional
        Color for the data1 trace (default 'black').
    color_data2 : str, optional
        Color for the data2 trace (default 'blue').
    color_data3 : str, optional
        Color for the data3 trace (default 'green').
    envelope_color1 : str, optional
        Fill color for the envelope of data1 (default 'lightgray').
    envelope_color2 : str, optional
        Fill color for the envelope of data2 (default 'lightgray').
    envelope_color3 : str, optional
        Fill color for the envelope of data3 (default 'lightgray').
    figsize : tuple, optional
        Figure size.
    envelope_pad_width : int, optional
        Padding width (in samples) for the Hilbert transform envelope computation.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure handle.
    """
    # Set default component labels if not provided.
    if component_labels is None:
        component_labels = [r'$\mathrm{Vertical}$', r'$\mathrm{Radial}$', r'$\mathrm{Transverse}$']
    
    # Process arrivals for data1.
    if arrivals1 is not None:
        if isinstance(arrivals1, (list, tuple, np.ndarray)):
            if len(arrivals1) != 3:
                raise ValueError("arrivals1 must have three elements: [P, S1, S2].")
            arrivals1 = {'P': arrivals1[0], 'S1': arrivals1[1], 'S2': arrivals1[2]}
    # Process arrivals for data2.
    if arrivals2 is not None:
        if isinstance(arrivals2, (list, tuple, np.ndarray)):
            if len(arrivals2) != 3:
                raise ValueError("arrivals2 must have three elements: [P, S1, S2].")
            arrivals2 = {'P': arrivals2[0], 'S1': arrivals2[1], 'S2': arrivals2[2]}
    # Process arrivals for data3.
    if arrivals3 is not None:
        if isinstance(arrivals3, (list, tuple, np.ndarray)):
            if len(arrivals3) != 3:
                raise ValueError("arrivals3 must have three elements: [P, S1, S2].")
            arrivals3 = {'P': arrivals3[0], 'S1': arrivals3[1], 'S2': arrivals3[2]}
    
     # Set default arrival labels.
    if arrival_labels1 is None:
        arrival_labels1 = {'P': r'$P_{1}$', 'S1': r'$S_{1,1}$', 'S2': r'$S_{2,1}$'}
    if arrival_labels2 is None:
        arrival_labels2 = {'P': r'$P_{2}$', 'S1': r'$S_{1,2}$', 'S2': r'$S_{2,2}$'}
    if arrival_labels3 is None:
        arrival_labels3 = {'P': r'$P_{3}$', 'S1': r'$S_{1,3}$', 'S2': r'$S_{2,3}$'}
        
    # Set default arrival colors.
    if arrival_colors1 is None:
        arrival_colors1 = {'P': 'red', 'S1': 'red', 'S2': 'red'}
    if arrival_colors2 is None:
        arrival_colors2 = {'P': 'blue', 'S1': 'blue', 'S2': 'blue'}
    if arrival_colors3 is None:
        arrival_colors3 = {'P': 'green', 'S1': 'green', 'S2': 'green'}
    
    if not plot_envelope:
        envelope_alpha = 0.0
    else: 
        envelope_alpha = 0.4
    
    # Create figure with 3 subplots (one per component).
    fig, axs = plt.subplots(3, 1, figsize=figsize, sharex=True)
    
    # Loop over the three components.
    for i in range(3):
        # Process data1.
        if source_wavelet is not None:
            cc1 = correlate(data1[:, i], source_wavelet, mode='same')
            cc1 = cc1 / np.max(np.abs(cc1))
        else:
            cc1 = data1[:, i]
        try:
            env1 = compute_envelope(cc1, pad_width=envelope_pad_width)
        except ValueError:
            env1 = np.abs(hilbert(cc1))
        axs[i].fill_between(time, -env1, env1, color=envelope_color1, alpha=envelope_alpha)
        axs[i].plot(time, cc1, color=color_data1, label=data_label1 if i == 0 else None, lw = 2, alpha = 0.75)
        
        # Process data2 if provided.
        if data2 is not None:
            if source_wavelet is not None:
                cc2 = correlate(data2[:, i], source_wavelet, mode='same')
                cc2 = cc2 / np.max(np.abs(cc2))
            else:
                cc2 = data2[:, i]
            try:
                env2 = compute_envelope(cc2, pad_width=envelope_pad_width)
            except ValueError:
                env2 = np.abs(hilbert(cc2))
            axs[i].fill_between(time, -env2, env2, color=envelope_color2, alpha=envelope_alpha)
            axs[i].plot(time, cc2, color=color_data2, label=data_label2 if i == 0 else None, lw = 2, alpha = 0.75)
        
        # Process data3 if provided.
        if data3 is not None:
            if source_wavelet is not None:
                cc3 = correlate(data3[:, i], source_wavelet, mode='same')
                cc3 = cc3 / np.max(np.abs(cc3))
            else:
                cc3 = data3[:, i]
            try:
                env3 = compute_envelope(cc3, pad_width=envelope_pad_width)
            except ValueError:
                env3 = np.abs(hilbert(cc3))
            axs[i].fill_between(time, -env3, env3, color=envelope_color3, alpha=envelope_alpha)
            axs[i].plot(time, cc3, color=color_data3, label=data_label3 if i == 0 else None, lw = 2, alpha = 0.75)
        
        axs[i].set_ylabel(component_labels[i], fontsize=12, rotation = yaxis_rotation)
        axs[i].grid(True, linestyle='--', alpha=0.5)
        
        # Plot vertical dashed lines for arrivals in data1.
        if arrivals1 is not None:
            for key, atime in arrivals1.items():
                lbl = arrival_labels1[key] if i == 0 else None
                axs[i].axvline(atime, color=arrival_colors1.get(key, 'red'), linestyle='--', label=lbl)
        
        # Plot vertical dashed lines for arrivals in data2.
        if data2 is not None and arrivals2 is not None:
            for key, atime in arrivals2.items():
                lbl = arrival_labels2[key] if i == 0 else None
                axs[i].axvline(atime, color=arrival_colors2.get(key, 'blue'), linestyle=':', label=lbl)
        
        # Plot vertical dashed lines for arrivals in data3.
        if data3 is not None and arrivals3 is not None:
            for key, atime in arrivals3.items():
                lbl = arrival_labels3[key] if i == 0 else None
                axs[i].axvline(atime, color=arrival_colors3.get(key, 'green'), linestyle='-.', label=lbl)
        
        # Add legend only on the first subplot.
        if i == 0:
            if show_legend:
                axs[i].legend(loc='upper right')
    
    if xlims is not None:
        for ax in axs:
            ax.set_xlim(xlims)
    
    axs[-1].set_xlabel("Time (s)", fontsize = 12)
    fig.tight_layout()
    return fig, axs

# ------------------------------------------------------------------------------
def plot_hodogram(data, dt, window=None, components=('L', 'Q', 'T'), title='Hodogram'):
    """
    Plots hodograms (particle motion plots) from given seismic data and returns the figure and axes objects.
    
    Parameters
    ----------
    data : ndarray (N x 3)
        Seismic data array with columns corresponding to [L, Q, T] components.
    dt : float
        Time interval between samples (seconds).
    window : tuple (start, end), optional
        Time window (in seconds) to select data. If None, uses full data.
    components : tuple of str
        Labels for components, default ('L', 'Q', 'T').
    title : str
        Plot title.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    axs : ndarray of matplotlib.axes.Axes
        Array containing the axes objects.
    """
    time = np.arange(data.shape[0]) * dt
    
    if window:
        idx = np.where((time >= window[0]) & (time <= window[1]))[0]
        data = data[idx]
        time = time[idx]
    
    L, Q, T = data[:,0], data[:,1], data[:,2]
    
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    
    # L vs. Q
    axs[0].plot(L, Q, color='k')
    axs[0].set_xlabel(f'{components[0]}')
    axs[0].set_ylabel(f'{components[1]}')
    axs[0].set_title(f'{components[0]} vs. {components[1]}')
    axs[0].grid(True)
    
    # L vs. T
    axs[1].plot(L, T, color='k')
    axs[1].set_xlabel(f'{components[0]}')
    axs[1].set_ylabel(f'{components[2]}')
    axs[1].set_title(f'{components[0]} vs. {components[2]}')
    axs[1].grid(True)
    
    # Q vs. T
    axs[2].plot(Q, T, color='k')
    axs[2].set_xlabel(f'{components[1]}')
    axs[2].set_ylabel(f'{components[2]}')
    axs[2].set_title(f'{components[1]} vs. {components[2]}')
    axs[2].grid(True)
    
    plt.tight_layout()
    
    return fig, axs



# ============================ Processing Functions ============================
def compute_envelope(signal, pad_width=50):
    """
    Compute the amplitude envelope of a signal using the Hilbert transform.
    Pads the signal (using reflection) to reduce edge effects, then removes the pad.
    """
    if len(signal) < 2*pad_width:
        raise ValueError("Signal too short for the given pad_width.")
    padded = np.pad(signal, pad_width, mode='reflect')
    analytic = hilbert(padded)
    envelope = np.abs(analytic)
    # Remove the padded regions.
    return envelope[pad_width:-pad_width]


# ------------------------------------------------------------------------------
def polarization_analysis(data, dt, M, alpha=0.1):
    """
    Compute rectilinearity, backazimuth, and incidence angle using a sliding Tukey window.
    Handles signal edges with internal padding and shorter windows near boundaries.

    Parameters
    ----------
    data : ndarray (N x 3)
        3-component motion data (e.g., LQT, ZRT, or XYZ).
    dt : float
        Time step in seconds.
    M : int
        Window length (samples).
    alpha : float
        Tukey window shape parameter (default = 0.1).
    
    Returns
    -------
    times : ndarray
        Time vector (seconds).
    rectilinearity : ndarray
        Rectilinearity values.
    backazimuth : ndarray
        Backazimuth angles (degrees from y-axis).
    incidence : ndarray
        Incidence angles (degrees from vertical).
    cone_mask : ndarray
        Boolean array, True where full-length window is used (i.e., outside the cone of influence).
    """
    N = data.shape[0]
    times = np.arange(N) * dt
    
    rectilinearity = np.full(N, np.nan)
    backazimuth = np.full(N, np.nan)
    incidence = np.full(N, np.nan)
    cone_mask = np.zeros(N, dtype=bool)
    
    half_M = M // 2
    
    for i in range(N):
        # Determine actual window range (shorter near edges)
        start = max(0, i - half_M)
        end = min(N, i + half_M)
        
        segment = data[start:end]
        win_len = end - start
        
        if win_len < 10:
            continue  # Not enough data to compute polarization
        
        # Apply Tukey window of current length
        taper = tukey(win_len, alpha)
        tapered = segment * taper[:, None]
        
        # Covariance and eigen-decomposition
        cov = np.cov(tapered.T)
        eigvals, eigvecs = np.linalg.eigh(cov)
        principal = eigvecs[:, 2]
        
        # Flip to ensure "upward" or consistent direction
        if principal[2] < 0:
            principal *= -1
        
        # Compute metrics
        r = 1.0 - eigvals[1] / eigvals[2]
        inc = np.degrees(np.arccos(principal[2]))
        ba = np.degrees(np.arctan2(principal[0], principal[1]))
        if ba < 0:
            ba += 360
        
        # Save results
        rectilinearity[i] = r
        backazimuth[i] = ba
        incidence[i] = inc
        
        # Flag if this is a full window
        cone_mask[i] = (win_len == M)
        
    return times, rectilinearity, backazimuth, incidence, cone_mask

# --------------------------------------------------------------------------
def rotate_to_zrt(data, source=None, receiver=None, direction=None):
    """
    Rotate seismogram data from Cartesian (vx, vy, vz) to Z, R, T components.
    
    Parameters
    ----------
    data : numpy.ndarray
        Array of shape (n_samples, 3) containing the seismogram data in (vx, vy, vz).
    source : array_like, optional
        Source coordinates as a 3-element array [x, y, z].
    receiver : array_like, optional
        Receiver coordinates as a 3-element array [x, y, z].
    direction : array_like, optional
        Propagation direction vector (3 elements). If provided, its horizontal 
        projection is used to define the radial direction.
        
    Returns
    -------
    rotated_data : numpy.ndarray
        Array of shape (n_samples, 3) containing the rotated seismogram data in (Z, R, T).
    
    Notes
    -----
    - If 'direction' is provided, it is normalized and its horizontal (x,y) components
    are used for the radial direction.
    - If 'direction' is not provided, both 'source' and 'receiver' must be given. The 
    horizontal difference (receiver - source) is used to compute the azimuth.
    - The vertical unit vector is assumed to be [0, 0, 1].
    - The transverse component is computed as the cross product of the vertical and radial directions.
    """
    # Determine the radial unit vector
    if direction is not None:
        direction = np.asarray(direction, dtype=float)
        norm_dir = np.linalg.norm(direction)
        if norm_dir == 0:
            raise ValueError("The provided direction vector has zero length.")
        direction_normalized = direction / norm_dir
        
        # Use the horizontal projection of the direction vector for the radial component
        radial_h = np.array([direction_normalized[0], direction_normalized[1], 0.0])
        norm_radial_h = np.linalg.norm(radial_h)
        if norm_radial_h < 1e-8:
            # If the horizontal part is nearly zero, choose an arbitrary horizontal direction.
            radial_unit = np.array([1.0, 0.0, 0.0])
        else:
            radial_unit = radial_h / norm_radial_h
    elif source is not None and receiver is not None:
        source = np.asarray(source, dtype=float)
        receiver = np.asarray(receiver, dtype=float)
        horizontal_diff = receiver[:2] - source[:2]  # Only x and y components
        norm_horizontal = np.linalg.norm(horizontal_diff)
        if norm_horizontal < 1e-8:
            # If source and receiver are nearly vertically aligned, choose an arbitrary horizontal direction.
            radial_unit = np.array([1.0, 0.0, 0.0])
        else:
            radial_unit = np.array([horizontal_diff[0] / norm_horizontal,
                                    horizontal_diff[1] / norm_horizontal, 0.0])
    else:
        raise ValueError("Either 'direction' or both 'source' and 'receiver' must be provided.")
    
    # Define the vertical unit vector (assumed to be [0, 0, 1])
    Z_unit = np.array([0.0, 0.0, 1.0])
    # Compute the transverse unit vector as the cross product of Z and radial
    T_unit = np.cross(Z_unit, radial_unit)
    norm_T = np.linalg.norm(T_unit)
    if norm_T < 1e-8:
        # If T is degenerate (should rarely happen), use an alternative horizontal direction.
        T_unit = np.array([0.0, 1.0, 0.0])
    else:
        T_unit = T_unit / norm_T
    
    # Build the rotation matrix.
    # Each row is one of the new unit vectors, so that:
    #   rotated_data = [Z, R, T] = data @ rotation_matrix.T
    rotation_matrix = np.vstack([Z_unit, radial_unit, T_unit])
    
    # Apply the rotation matrix to each sample in the seismogram.
    rotated_data = np.dot(data, rotation_matrix.T)
    
    return rotated_data

def rotate_to_qlt(data, source_location, receiver_location, backazimuth = None, incidence = None):
    """
    Rotate data from (x, y, z) to (L, Q, T) components.

    Parameters
    ----------
    data : ndarray (N x 3)
        Input data with columns [vx, vy, vz].
    source_location : ndarray (3,)
        Source coordinates [x, y, z].
    receiver_location : ndarray (3,)
        Receiver coordinates [x, y, z].

    Returns
    -------
    ba : float
        Backazimuth angle in radians.
    inc : float
        Inclination angle in radians.
    lqt : ndarray (N x 3)
        Rotated data [L, Q, T].
    """ 
    
    # Compute displacement vector (receiver - source)
    d = receiver_location - source_location
    dx, dy, dz = d
    
    # Horizontal distance
    if incidence is None:
        H = np.sqrt(dx**2 + dy**2)
        
        # Inclination angle (positive downward)
        inc = np.arctan2(dz, H)
    else:
        inc = incidence 
    
    # Backazimuth angle (clockwise from +y)
    if backazimuth is None:
        ba = np.arctan2(dx, dy)
        if ba < 0:
            ba += 2 * np.pi
    else: 
        ba = backazimuth
    
    # Trigonometric components
    sin_inc, cos_inc = np.sin(inc), np.cos(inc)
    sin_ba, cos_ba = np.sin(ba), np.cos(ba)
    
    # Rotation matrix from XYZ to LQT
    # L aligned along propagation direction
    R = np.array([
        [sin_inc * sin_ba,  sin_inc * cos_ba, cos_inc],     # L-direction (along propagation)
        [cos_inc * sin_ba,  cos_inc * cos_ba, -sin_inc],    # Q-direction (SV-plane, perpendicular to propagation)
        [cos_ba,           -sin_ba,           0.0]          # T-direction (SH, perpendicular horizontally)
    ])
    
    # Apply rotation
    lqt = data @ R.T
    
    return ba, inc, lqt

# ------------------------------------------------------------------------------
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
