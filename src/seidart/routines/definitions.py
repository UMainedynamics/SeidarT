import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.animation as anim
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import os.path
from typing import Tuple, Optional, Dict, List, Union
from subprocess import call
from scipy.io import FortranFile
from scipy.signal import hilbert, correlate, butter, filtfilt, sosfiltfilt, firwin, minimum_phase, savgol_filter

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
    'compute_polar_energy',
    'compute_phase_gain',
    'agc',
    'correct_geometric_spreading',
    'exponential_gain',
    'plot_3c',
    'plot_hodogram',
    'compute_dispersion',
    'compute_dispersion_image',
    'compute_fk_spectrum', 
    'plot_fk_spectrum',
    'plot_dispersion',
    'plot_dispersion_image',
    'compute_fk_spectrum',
    'plot_fk_spectrum',
    'CFL', 
    'clight',
]

# --------------------------------- Globals ------------------------------------
clight = 2.99792458e8  # Speed of light in vacuum (m/s)
eps0 = 8.85418782e-12  # Permittivity of free space
mu0 = 4.0 * np.pi * 1.0e-7  # Permeability of free space
mu_r = 1.0


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
        domain.alpha_max_scalar, domain.kappa_max,
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
        half: bool = False,
    ) -> None:
    """
    Computes CPML parameters for a given direction and updates model/domain.

    :param modelclass: The model class instance to update.
    :param domain: The domain class instance to update.
    # :param direction: Direction to compute CPML ('x', 'y', or 'z').
    :param half: Flag to compute half CPML parameters. Defaults to False.
    :type modelclass: Model
    :type domain: Domain
    # :type direction: str
    :type half: bool
    """
    nx = domain.nx + 2*domain.cpml
    nz = domain.nz + 2*domain.cpml
    if domain.dim == 2.5:
        # ny = domain.ny + 2*domain.cpml
        deltamin = np.min([domain.dx, domain.dy, domain.dz]) 
    else:
        deltamin = np.min([domain.dx, domain.dz]) 
    
    # -----------------------------------------------------------------------------
    # Compute the distance along the absorbing boundary relative to the end of the 
    # original model space. 
    distx = domain.dx * np.arange(0, domain.cpml)
    distz = domain.dz * np.arange(0, domain.cpml)
    if half:
        distx = distx + domain.dx/2
        distz = distz + domain.dz/2
    
    distx = domain.dx*domain.cpml - distx
    distz = domain.dz*domain.cpml - distz
    distx = distx/(domain.dx*domain.cpml)
    distz = distz/(domain.dz*domain.cpml)
    
    # Create lookup array indexed by material ID
    lookup_array = np.zeros(domain.nmats)
    for mat_id, v in modelclass.max_velocity_per_material.items():
        lookup_array[mat_id] = v
    
    # Now map geometry to velocities
    velocity_map = lookup_array[domain.geometry]
    
    # Compute the maximum sigma, and alpha values for the CPML.  
    if modelclass.is_seismic:
        alpha_max = domain.alpha_max_scalar * np.pi*modelclass.f0
        quasi_cp_max = 0.7 * velocity_map / 2.0 
        sig_max = - np.log(domain.Rcoef) * (domain.NP+1) * quasi_cp_max / (2.0 * domain.cpml )
        # This seems to work well even at higher frequencies
        sigma, kappa, alpha, acoef, bcoef = cpml_parameters(
            sig_max, alpha_max, domain.kappa_max, nx, nz,
            distx, distz, domain.NP, domain.NPA, modelclass.dt, is_seismic = True
        )
    else:
        # We will use the maximum permittivity coefficient and assume that the 
        # magnetic permeability is 1. We can use a different value
        # sig_max = domain.sig_opt_scalar * ((domain.NP + 1) / (deltamin * ((mu0/eps0)**0.5) ) ) #!!! We need to multiply eps0 by the relative permattivity for a better estimate
        # eps_r_min = modelclass.permittivity_coefficients[['e11', 'e22', 'e33']].min().min()
        # c_max = clight * np.sqrt(1/eps_r_min) 
        # sig_max = -(domain.NP + 1) * np.log(domain.Rcoef) * eps0 * c_max / (2 * domain.cpml * deltamin)
        c_max = clight / np.sqrt(velocity_map)
        sig_max = -(domain.NP + 1) * np.log(domain.Rcoef) * eps0 * c_max / (2 * domain.cpml * deltamin)
        alpha_max = domain.alpha_max_scalar * np.pi * eps0 * modelclass.f0 
        sigma, kappa, alpha, acoef, bcoef = cpml_parameters(
            sig_max, alpha_max, domain.kappa_max, nx, nz,
            distx, distz, domain.NP, domain.NPA, modelclass.dt
        )
    
    # Save the results to a fortran binary
    if half:
        sigma_fn = 'sigma_half_cpml.dat'
        kappa_fn = 'kappa_half_cpml.dat'
        alpha_fn = 'alpha_half_cpml.dat'
        acoef_fn = 'acoef_half_cpml.dat'
        bcoef_fn = 'bcoef_half_cpml.dat'
    else:
        sigma_fn = 'sigma_cpml.dat'
        kappa_fn = 'kappa_cpml.dat'
        alpha_fn = 'alpha_cpml.dat'
        acoef_fn = 'acoef_cpml.dat'
        bcoef_fn = 'bcoef_cpml.dat'
    
    f = FortranFile(sigma_fn, 'w')
    f.write_record(sigma.T)
    f.close()
    f = FortranFile(kappa_fn, 'w')
    f.write_record(kappa.T)
    f.close()
    f = FortranFile(alpha_fn, 'w')
    f.write_record(alpha.T)
    f.close()
        
    f = FortranFile(acoef_fn, 'w')
    f.write_record(acoef.T)
    f.close()
    f = FortranFile(bcoef_fn, 'w')
    f.write_record(bcoef.T)
    f.close()
        
    return sigma, kappa, alpha, acoef, bcoef
        
# -----------------------------------------------------------------------------
def cpml_parameters(
        sig_max: float, 
        alpha_max: float, 
        kappa_max: float, 
        nx: int,
        nz: int,
        distancex,
        distancez, 
        NP: int,
        NPA: int, 
        dt: float,
        is_seismic = False,
    ):
    """
    """
    kappa = np.ones([nx,nz])
    alpha = np.zeros([nx,nz])
    sigma = np.zeros([nx,nz])
    acoeff = np.zeros([nx,nz])
    bcoeff = np.zeros([nx,nz])
    
    distx_max = distancex.max() 
    distz_max = distancez.max()
    
    m = len(distancex)
    n = len(distancez)
    # Compute in the x, and z direction
    for ind in range(0, m):
        sigma[ind,m:-m] = sig_max[0,:] * (distancex[ind]**NP)
        kappa[ind,:] = 1.0 + (kappa_max - 1.0) * (distancex[ind]/distx_max)**NP
        alpha[ind,:] = alpha_max * (1 - distancex[ind])**NPA
        
        # From the end 
        sigma[-(ind+1),m:-m] = sig_max[-1,:]*distancex[ind]**NP
        kappa[-(ind+1),:] = 1 + (kappa_max - 1) * (distancex[ind]/distx_max)**NP
        alpha[-(ind+1),:] = alpha_max * (1 - distancex[ind])**NPA
    
    for ind in range(0, n):        
        sigma[n:-n,ind] = sig_max[:,0] * (distancez[ind]**NP)
        kappa[:,ind] = 1.0 + (kappa_max - 1.0) * (distancez[ind]/distz_max)**NP
        alpha[:,ind] = alpha_max * (1 - distancez[ind])**NPA
                    
        sigma[n:-n,-(ind+1)] = sig_max[:,-1]*distancez[ind]**NP
        kappa[:,-(ind+1)] = 1 + (kappa_max - 1) * (distancez[ind]/distz_max)**NP
        alpha[:,-(ind+1)] = alpha_max * (1 - distancez[ind])**NPA
    
    for ii in range(m):
        for jj in range(n):
            r = np.sqrt(distancex[ii]**2 + distancez[jj]**2)
            sigma[ii, jj]           = sig_max[0, 0]   * (r ** NP)
            sigma[-(ii+1), jj]      = sig_max[-1, 0]  * (r ** NP)
            sigma[ii, -(jj+1)]      = sig_max[0, -1]  * (r ** NP)
            sigma[-(ii+1), -(jj+1)] = sig_max[-1, -1] * (r ** NP)
            
            kappa[ii, jj]           = 1.0 + (kappa_max - 1.0) * (r ** NP)
            kappa[-(ii+1), jj]      = 1.0 + (kappa_max - 1.0) * (r ** NP)
            kappa[ii, -(jj+1)]      = 1.0 + (kappa_max - 1.0) * (r ** NP)
            kappa[-(ii+1), -(jj+1)] = 1.0 + (kappa_max - 1.0) * (r ** NP)
            
            alpha[ii, jj]           = alpha_max * ((1.0 - r) ** NPA)
            alpha[-(ii+1), jj]      = alpha_max * ((1.0 - r) ** NPA)
            alpha[ii, -(jj+1)]      = alpha_max * ((1.0 - r) ** NPA)
            alpha[-(ii+1), -(jj+1)] = alpha_max * ((1.0 - r) ** NPA)
    
    if is_seismic:
        bcoeff = np.exp( - (sigma / kappa + alpha) * dt)
    else:
        bcoeff = np.exp( - (sigma / kappa + alpha) * (dt/eps0) )
    
    # Compute the a-coefficients 
    alpha[np.where(alpha < 0.0)] = 0.0
    indices = np.where(np.abs(sigma) > 1.0e-6)
    acoeff[indices] = sigma[indices] * (bcoeff[indices] - 1) / \
            (kappa[indices] * sigma[indices] + kappa[indices] * alpha[indices] )
    
    return sigma, kappa, alpha, acoeff, bcoeff


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

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
def compute_polar_energy(E_x, E_y, dt, angles, to_db=False, corrected=False):
    """
    Compute the integrated energy of the horizontal field projection 
    onto multiple polarization angles for a 3C time series.
    
    Parameters
    ----------
    E_x : array_like
        1D array of the x-component time series.
    E_y : array_like
        1D array of the y-component time series.
    dt : float
        Time-step interval (seconds).
    angles : array_like
        1D array of angles (radians) at which to project the field.
    to_db : bool, optional
        If True, convert energies to dB scale: 10*log10(E/E.max()).
    
    Returns
    -------
    energies : np.ndarray
        Array of shape (len(angles),) containing the energy at each angle,
        or in dB if `to_db=True`.
    """
    E_x = np.asarray(E_x)
    E_y = np.asarray(E_y)
    energies = np.empty(len(angles))
    
    for i, theta in enumerate(angles):
        # Project onto the polarization direction theta
        E_theta = E_x * np.cos(theta) + E_y * np.sin(theta)
        # Compute energy = integral of |E|^2 dt
        energies[i] = np.sum(np.abs(E_theta)**2) * dt
    
    if to_db:
        # Normalize to max and convert to dB
        energies = 10 * np.log10(energies / np.max(energies))
    
    if corrected:
        baseline = np.cos(angles)**2
        corr = energies / baseline 
        return 10*np.log10(corr / corr.max() )
    else:
        return energies

# ------------------------------------------------------------------------------
def compute_phase_gain(ts_ref, ts, f0, dt):
    """
    Compute phase shift and gain (in dB) between two time series at a specific frequency.
    
    Parameters
    ----------
    ts_ref : array_like
        Reference time series (1D array).
    ts : array_like
        Comparison time series (same shape as ts_ref).
    f0 : float
        Center frequency in Hz at which to compute phase and gain.
    dt : float
        Sampling interval in seconds.
    
    Returns
    -------
    phase_deg : float
        Phase shift of `ts` relative to `ts_ref` at f0, in degrees.
    gain_db : float
        Gain in decibels: 20 * log10(|H|), where H is the complex ratio of spectral amplitudes.
    """
    ts_ref = np.asarray(ts_ref)
    ts = np.asarray(ts)
    if ts_ref.shape != ts.shape:
        raise ValueError("Input time series must have the same shape")
    N = ts_ref.size
    
    # Compute real FFTs
    Y_ref = np.fft.rfft(ts_ref)
    Y = np.fft.rfft(ts)
    freqs = np.fft.rfftfreq(N, d=dt)
    
    # Find the index of the bin closest to f0
    idx = np.argmin(np.abs(freqs - f0))
    
    # Complex transfer at f0
    H = Y[idx] / Y_ref[idx]
    
    gain_db = 20 * np.log10(np.abs(H))
    phase_deg = np.angle(H, deg=True)
    
    return phase_deg, gain_db

# ------------------------------------------------------------------------------
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
def exponential_gain(time_series: np.ndarray, dt: float, alpha: float) -> np.ndarray:
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


# ------------------------------------------------------------------------------
def compute_fk_spectrum(
        data: np.ndarray,
        dt: float,
        dx: float,
        taper: bool = True,
        ntfft: Optional[int] = None,
        nxfft: Optional[int] = None,
        fmin: float = 0.0,
        fmax: Optional[float] = None,
        kmin: Optional[float] = None,
        kmax: Optional[float] = None,
        wavenumber_units: str = 'cycles'  # 'cycles' or 'radians'
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the slant-stack f–k power spectrum of a 2D shot gather.
    
    Returns
    -------
    freqs : 1D array of frequencies (Hz)
    ks     : 1D array of wavenumbers in either cycles/m or rad/m
    P      : 2D power spectrum array shape (len(freqs), len(ks))
    """
    nt, nx = data.shape
    ntfft = ntfft or nt
    nxfft = nxfft or nx
    
    # optional taper to suppress leakage
    dat = data.copy()
    if taper:
        dat *= np.hanning(nt)[:, None] * np.hanning(nx)[None, :]
    
    # 2D FFT and power
    F2 = np.fft.fft2(dat, s=(ntfft, nxfft))
    Pfull = np.abs(F2)**2
    
    # frequency axis
    freqs_full = np.fft.rfftfreq(ntfft, dt)
    if fmax is None:
        fmask = freqs_full >= fmin
    else:
        fmask = (freqs_full >= fmin) & (freqs_full <= fmax)
    freqs = freqs_full[fmask]
    
    # wavenumber axis (cycles/m)
    ks_full = np.fft.fftshift(np.fft.fftfreq(nxfft, dx))
    # convert to rad/m if requested
    if wavenumber_units == 'radians':
        ks_full = ks_full * 2 * np.pi
    
    # apply k bounds
    if kmin is None: kmin = 0.0
    if kmax is None: kmax = np.abs(ks_full).max()
    kmask = (np.abs(ks_full) >= kmin) & (np.abs(ks_full) <= kmax)
    ks = ks_full[kmask]
    
    # extract and shift power matrix
    P = Pfull[: len(freqs_full), :]
    P = P[fmask, :]
    P = np.fft.fftshift(P, axes=1)
    P = P[:, kmask]
    
    return freqs, ks, P

# -----------------------------------------------------------------------------

def plot_fk_spectrum(
        freqs: np.ndarray,
        ks: np.ndarray,
        P: np.ndarray,
        fig: Optional[Figure] = None,
        ax: Optional[Axes] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        log_scale: bool = False,
        smooth: bool = False,
        smooth_sigma: Tuple[float, float] = (1.0, 1.0),
        mode_lines: Optional[Dict[str, Union[float, Tuple[float,str]]]] = None,
        wavenumber_units: str = 'cycles'  # match compute_fk_spectrum
    ) -> Tuple[Figure, Axes]:
    """
    Display an f–k power spectrum, with optional log-scale, smoothing,
    and overlaid theoretical mode-velocity lines in either cycles/m or rad/m.
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    
    # optional smoothing
    Pplot = P.copy()
    if smooth:
        Pplot = gaussian_filter(Pplot, sigma=smooth_sigma)
    
    eps = np.finfo(float).eps
    if log_scale:
        Pplot = 10 * np.log10(Pplot + eps)
        cbar_label = 'Power (dB)'
    else:
        cbar_label = 'Power'
    
    pcm = ax.pcolormesh(ks, freqs, Pplot, shading='auto', vmin=vmin, vmax=vmax)
    fig.colorbar(pcm, ax=ax, label=cbar_label)
    
    # overlay mode lines
    if mode_lines:
        for label, spec in mode_lines.items():
            if isinstance(spec, (tuple, list)):
                v, col = spec
            else:
                v, col = spec, 'white'
            # compute k_line in correct units
            if wavenumber_units == 'radians':
                kline = 2 * np.pi * freqs / v
            else:
                kline = freqs / v
            ax.plot(kline, freqs, '--', color=col, lw=1, label=label)
            ax.plot(-kline, freqs, '--', color=col, lw=1)
        # ax.legend(loc='upper right')
    
    xlabel = 'Wavenumber (rad/m)' if wavenumber_units=='radians' else 'Wavenumber (1/m)'
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Frequency (Hz)')
    ax.set_title('f–k Spectrum')
    return fig, ax

# ------------------------------------------------------------------------------
# SURFACE WAVES
def compute_dispersion(
        seismograms: np.ndarray, 
        dx: float, 
        dt: float,
        fmin: float, 
        fmax: float, 
        nfreq: int=50
    ) -> Tuple[np.ndarray,np.ndarray]:
    """Linear-phase fit with 2π-branch correction and physical clamp."""
    nt, nx = seismograms.shape
    spec = np.fft.rfft(seismograms, axis=0)
    freq_axis = np.fft.rfftfreq(nt, dt)
    freqs = np.linspace(fmin, fmax, nfreq)
    vs = np.full(nfreq, np.nan)
    x = np.arange(nx)*dx
    L = x[-1]-x[0]
    # max physically resolvable velocity (Nyquist spatial): dx/dt
    max_vel = dx / dt
    for i,f in enumerate(freqs):
        idx = np.argmin(np.abs(freq_axis-f))
        ph = np.unwrap(np.angle(spec[idx,:]))
        raw = np.polyfit(x,ph,1)[0]
        impl = raw*L
        branch = np.round(impl/(2*np.pi))
        slope = raw - branch*(2*np.pi/L)
        k = -slope
        v = 2*np.pi*f/ k
        # clamp to positive and below max_vel
        if v > 0 and v <= max_vel:
            vs[i] = v
        else:
            vs[i] = np.nan
    return freqs, vs

# ============= Slant-stack energy image and smooth picks =============
def compute_dispersion_image(
        seismograms: np.ndarray,
        dx: float,
        dt: float,
        fmin: float,
        fmax: float,
        nv: int = 100,
        nfreq: int = 100,
        vmin: float = None,
        vmax: float = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """MASW slant-stack with taper, smoothing, and sub-bin peak picking."""
    # taper in time to reduce spectral leakage
    nt, nx = seismograms.shape
    win = np.hanning(nt)
    tapered = seismograms * win[:, None]
    # FFT in time
    spec = np.fft.rfft(tapered, axis=0)
    freq_axis = np.fft.rfftfreq(nt, dt)
    
    # initial phase-fit picks for auto bounding
    freqs, phase_picks = compute_dispersion(seismograms, dx, dt, fmin, fmax, nfreq)
    # determine vmin/vmax: user override or auto
    if vmin is not None and vmax is not None:
        vmin, vmax = vmin, vmax
    else:
        pos = phase_picks[np.isfinite(phase_picks) & (phase_picks > 0)]
        if pos.size > 0:
            vmin = pos.min() * 0.8
            vmax = pos.max() * 1.2
        else:
            vmin = dx/dt * 0.1
            vmax = dx/dt * 2.0
    if vmin >= vmax:
        raise ValueError(f"Invalid velocity bounds: vmin={vmin}, vmax={vmax}")
    
    vs = np.linspace(vmin, vmax, nv)
    image = np.zeros((len(freqs), len(vs)))
    x = np.arange(nx) * dx
    
    # build slant-stack energy image
    for i, f in enumerate(freqs):
        idx = np.argmin(np.abs(freq_axis - f))
        A = spec[idx, :]
        omega = 2 * np.pi * f
        for j, v in enumerate(vs):
            image[i, j] = np.abs(np.sum(A * np.exp(1j * omega * x / v)))**2
    
    # normalize per-frequency
    image /= (image.max(axis=1, keepdims=True) + np.finfo(float).eps)
    # smooth image to reduce striping
    from scipy.ndimage import gaussian_filter
    image = gaussian_filter(image, sigma=(1, 2))
    
    # sub-bin (parabolic) peak picking
    picks = np.zeros_like(freqs)
    dv = vs[1] - vs[0]
    for i in range(len(freqs)):
        row = image[i]
        imax = np.nanargmax(row)
        if 0 < imax < len(vs) - 1:
            y0, y1, y2 = row[imax-1], row[imax], row[imax+1]
            dp = 0.5 * (y0 - y2) / (y0 - 2*y1 + y2)
            picks[i] = vs[imax] + dp * dv
        else:
            picks[i] = vs[imax]
    # mask outside bounds
    picks = np.where((picks >= vmin) & (picks <= vmax), picks, np.nan)
    
    return freqs, vs, image, picks

# ----------------------------------------------------------------------------

def plot_dispersion(freqs: np.ndarray, velocities: np.ndarray) -> None:
    """
    Plot the dispersion curve freq vs phase velocity.
    """
    plt.figure()
    plt.plot(freqs, velocities, '-o')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Phase velocity (m/s)')
    plt.title('Dispersion curve')
    plt.grid(True)
    plt.show()


def plot_dispersion_image(
        freqs: np.ndarray,
        vs: np.ndarray,
        image: np.ndarray,
        picks_freqs: Optional[np.ndarray] = None,
        picks_vs: Optional[np.ndarray] = None,
        colormap = 'inferno'
    ) -> Tuple[Figure, Axes]:
    """
    Plot a MASW-style dispersion image and optional picked dispersion curve.
    
    Parameters
    ----------
    freqs : (Nf,) array
        Frequencies (Hz).
    vs : (Nv,) array
        Phase velocities (m/s).
    image : (Nf, Nv) array
        Normalized energy map (frequency × velocity).
    picks_freqs : (Np,) array, optional
        Frequencies of picked dispersion curve.
    picks_vs : (Np,) array, optional
        Phase velocities of picked dispersion curve.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # plot the energy map
    pcm = ax.pcolormesh(freqs, vs, image.T, shading='auto', cmap = colormap)
    cbar = fig.colorbar(pcm, ax=ax, label='Normalized Energy')
    
    # labels and title
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Phase Velocity (m/s)')
    ax.set_title('Dispersion Image')
    
    # overlay picks if provided
    if picks_freqs is not None and picks_vs is not None:
        ax.scatter(
            picks_freqs, picks_vs,
            c='r', marker='.', s=10,
            label='Picked Curve'
        )
        ax.legend(loc='upper right')
    
    return fig, ax

# ----------------------------------------------------------------------------
def apply_time_taper(data, dt, t0=0.18, t1=0.25):
    """
    Apply a cosine taper to suppress data after t0–t1 (in seconds).
    """
    nt, nx = data.shape
    time = np.arange(nt) * dt
    taper = np.ones_like(time)

    # Apply taper between t0 and t1
    mask = (time >= t0) & (time <= t1)
    taper[mask] = 0.5 * (1 + np.cos(np.pi * (time[mask] - t0) / (t1 - t0)))
    taper[time > t1] = 0.0

    return data * taper[:, None]