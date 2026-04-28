import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.animation as anim
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import os.path
from typing import Tuple, Optional, Dict, List, Union
from numpy.typing import ArrayLike, NDArray
from subprocess import call
from scipy.io import FortranFile
from scipy.signal import hilbert, correlate, butter, filtfilt, sosfiltfilt, firwin, minimum_phase, savgol_filter
from scipy.ndimage import gaussian_filter
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
    'sponge_boundary',
    'rotate_sdr',
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
    domain.dim = data['Domain']['dim']
    domain.nx = data['Domain']['nx']
    domain.ny = data['Domain']['ny']
    domain.nz = data['Domain']['nz']
    domain.dx = data['Domain']['dx']
    domain.dy = data['Domain']['dy']
    domain.dz = data['Domain']['dz']
    domain.cpml = data['Domain']['cpml']
    domain.nmats = data['Domain']['nmats']
    domain.alpha_max_scalar = data['Domain']['alpha_max_scalar']
    domain.kappa_max = data['Domain']['kappa_max']
    domain.NP = data['Domain']['NP']
    domain.NPA = data['Domain']['NPA']
    domain.Rcoef = data['Domain']['Rcoef']
    domain.imfile = data['Domain']['image_file']
    # (
    #     domain.dim, domain.nx, domain.ny, domain.nz, 
    #     domain.dx, domain.dy, domain.dz, domain.cpml, 
    #     domain.nmats, 
    #     domain.alpha_max_scalar, domain.kappa_max,
    #     domain.NP, domain.NPA, domain.Rcoef,
    #     domain.imfile 
    # ) = list(data['Domain'].values())
    im,__ = image2int(data['Domain']['image_file'])
    domain.geometry = im.transpose().astype(int)
    
    # ---------------------------- Seismic Values ------------------------------
    
    seismic.dt = data['Seismic']['Source']['dt']
    seismic.time_steps = data['Seismic']['Source']['time_steps'] 
    seismic.x = data['Seismic']['Source']['x']
    seismic.y = data['Seismic']['Source']['y']
    seismic.z = data['Seismic']['Source']['z']
    seismic.f0 = data['Seismic']['Source']['source_frequency']
    seismic.theta = data['Seismic']['Source']['x-z_rotation']
    seismic.phi = data['Seismic']['Source']['x-y_rotation']
    seismic.psi = data['Seismic']['Source']['y-z_rotation']
    seismic.source_amplitude = data['Seismic']['Source']['amplitude']
    seismic.source_wavelet = data['Seismic']['Source']['source_wavelet']
    seismic.source_type = data['Seismic']['Source']['source_type']
    
    electromag.dt = data['Electromagnetic']['Source']['dt']
    electromag.time_steps = data['Electromagnetic']['Source']['time_steps'] 
    electromag.x = data['Electromagnetic']['Source']['x']
    electromag.y = data['Electromagnetic']['Source']['y']
    electromag.z = data['Electromagnetic']['Source']['z']
    electromag.f0 = data['Electromagnetic']['Source']['source_frequency']
    electromag.theta = data['Electromagnetic']['Source']['x-z_rotation']
    electromag.phi = data['Electromagnetic']['Source']['x-y_rotation']
    electromag.psi = data['Electromagnetic']['Source']['y-z_rotation']
    electromag.source_amplitude = data['Electromagnetic']['Source']['amplitude']
    electromag.source_wavelet = data['Electromagnetic']['Source']['source_wavelet']
    electromag.source_type = data['Electromagnetic']['Source']['source_type']
    
    # (
    #     electromag.dt, electromag.time_steps,
    #     electromag.x, electromag.y, electromag.z,
    #     electromag.xind, electromag.yind, electromag.zind,
    #     electromag.f0, electromag.theta, electromag.phi, electromag.psi,
    #     electromag.source_amplitude, electromag.source_wavelet    
    # ) = list(data['Electromagnetic']['Source'].values())

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

# ------------------------------------------------------------------------------
def cpmlcompute(
        modelclass, 
        domain, 
        half: bool = False,
        velocity_scaling_factor = 1.0,
        pml_smoothing: bool = "global",
        smoothing_window: int = 7, 
        smoothing_passes: int = 1
    ) -> None:
    """
    Computes CPML parameters for a given direction and updates model/domain.

    :param modelclass: The model class instance to update.
    :param domain: The domain class instance to update.
    # :param direction: Direction to compute CPML ('x', 'y', or 'z').
    :param half: Flag to compute half CPML parameters. Defaults to False.
    :param pml_smoothing: Specify the type of smoothing to apply to the PML
        layer to avoid impedence contrasts and instability within the region. 
        Options are: 
            None - no smoothing is applied and cpml parameters will reflect the inner modeling domain boundary values
            global (Default) - set the max values from the max velocity in the domain 
            smoothed - computes the max from the boundary value and smooths within the PML using the arithmetic mean
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
    
    modelclass.compute_max_velocities(dim = domain.dim) 
    
    # build separate max velocities for x and z directions
    lookup_array_x = np.zeros(domain.nmats)
    lookup_array_z = np.zeros(domain.nmats)
    
    modelclass.compute_max_velocities(dim = domain.dim, directions='x')
    for mat_id, v in modelclass.vmax_x.items():
        lookup_array_x[mat_id] = v
    
    modelclass.compute_max_velocities(dim = domain.dim, directions='z')
    for mat_id, v in modelclass.vmax_z.items():
        lookup_array_z[mat_id] = v
    
    velocity_map_x = lookup_array_x[domain.geometry]
    velocity_map_z = lookup_array_z[domain.geometry]
    
    # Compute the maximum sigma, and alpha values for the CPML.  
    if modelclass.is_seismic:    
        alpha_max = domain.alpha_max_scalar * np.pi*modelclass.f0
        quasi_cp_max_x = velocity_scaling_factor * 0.7 * velocity_map_x / 2.0
        quasi_cp_max_z = velocity_scaling_factor *  0.7 * velocity_map_z / 2.0
        sig_max_x = - np.log(domain.Rcoef) * (domain.NP+1) * quasi_cp_max_x / (2.0 * domain.cpml)
        sig_max_z = - np.log(domain.Rcoef) * (domain.NP+1) * quasi_cp_max_z / (2.0 * domain.cpml)
    else:
        c_max_x = clight / np.sqrt(velocity_map_x)
        c_max_z = clight / np.sqrt(velocity_map_z) 
        sig_max_x = -(domain.NP + 1) * np.log(domain.Rcoef) * eps0 * c_max_x / (2 * domain.cpml * deltamin)
        sig_max_z = -(domain.NP + 1) * np.log(domain.Rcoef) * eps0 * c_max_z / (2 * domain.cpml * deltamin)
        alpha_max = domain.alpha_max_scalar * np.pi * eps0 * modelclass.f0 
    
    # -------------------- apply pml_smoothing strategy --------------------
    # We only need the tangential boundary lines because cpml_parameters
    # consumes sig_max_x[0, :], sig_max_x[-1, :], sig_max_z[:, 0], sig_max_z[:, -1].
    mode = (pml_smoothing or "global").lower()
    
    if mode == "global":
        if modelclass.is_seismic:
            vgx = float(np.nanmax(velocity_map_x))
            vgz = float(np.nanmax(velocity_map_z))
            qc_x = velocity_scaling_factor * 0.7 * vgx / 2.0
            qc_z = velocity_scaling_factor * 0.7 * vgz / 2.0
            s_x = - np.log(domain.Rcoef) * (domain.NP + 1) * qc_x / (2.0 * domain.cpml)
            s_z = - np.log(domain.Rcoef) * (domain.NP + 1) * qc_z / (2.0 * domain.cpml)
        else:
            cmax_x = float(np.nanmax(clight / np.sqrt(np.maximum(velocity_map_x, 1.0))))
            cmax_z = float(np.nanmax(clight / np.sqrt(np.maximum(velocity_map_z, 1.0))))
            s_x = -(domain.NP + 1) * np.log(domain.Rcoef) * eps0 * cmax_x / (2.0 * domain.cpml * deltamin)
            s_z = -(domain.NP + 1) * np.log(domain.Rcoef) * eps0 * cmax_z / (2.0 * domain.cpml * deltamin)
        
        # overwrite with constants (most robust)
        sig_max_x = np.full_like(sig_max_x, s_x, dtype=float)
        sig_max_z = np.full_like(sig_max_z, s_z, dtype=float)
    elif mode == "smoothed":
        # Smooth tangential lines only (avoid jagged contrasts along the PML):
        # - for x-sides, smooth along z: rows 0 and -1
        # - for z-sides, smooth along x: cols 0 and -1
        # (The variation across thickness comes from the rx/rz profiles inside cpml_parameters.)
        sig_max_x[0,  :]  = _smooth1d(sig_max_x[0,  :], k=7, passes=2)
        sig_max_x[-1, :]  = _smooth1d(sig_max_x[-1, :], k=7, passes=2)
        sig_max_z[:, 0]   = _smooth1d(sig_max_z[:, 0],  k=7, passes=2)
        sig_max_z[:, -1]  = _smooth1d(sig_max_z[:, -1], k=7, passes=2)
        
        # Ensure non-negativity
        np.maximum(sig_max_x, 0.0, out=sig_max_x)
        np.maximum(sig_max_z, 0.0, out=sig_max_z)
    else:
        # Do nothing
        print("No cpml smoothing specified.")
    
    
    sigma, kappa, alpha, acoef, bcoef = cpml_parameters(
        sig_max_x, sig_max_z, alpha_max, domain.kappa_max, nx, nz,
        distx, distz, domain.NP, domain.NPA, modelclass.dt, is_seismic=False
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

# ------------------------------------------------------------------------------
def smooth1d(x: ArrayLike, k: int = 7, passes: int = 1) -> NDArray[np.floating]:
    """
    Smooth a 1D array using an arithmetic moving average.
    
    The function applies a length-``k`` moving-average filter to the input
    one or more times. The window size is coerced to an odd integer of at
    least 3 to keep the kernel symmetric.
    
    :param x: Input 1D data. Any array-like object convertible to a NumPy array.
    :type x: ArrayLike
    :param k: Nominal window length for the moving average. Must be positive;
              values are coerced to the nearest odd integer greater than or
              equal to 3.
    :type k: int
    :param passes: Number of times to apply the moving-average filter.
    :type passes: int
    :returns: Smoothed data with the same shape as the input.
    :rtype: numpy.typing.NDArray[numpy.floating]
    """
    # arithmetic moving average smoother
    k = max(3, int(k) | 1)  # ensure odd >=3
    kernel = np.ones(k, dtype=float) / k
    y = x.astype(float)
    for _ in range(int(passes)):
        y = np.convolve(y, kernel, mode="same")
    return y

# ------------------------------------------------------------------------------
def cpml_parameters(
        sig_max_x: NDArray[np.floating],
        sig_max_z: NDArray[np.floating],
        alpha_max: float,
        kappa_max: float,
        nx: int,
        nz: int,
        distancex: NDArray[np.floating],
        distancez: NDArray[np.floating],
        NP: float,
        NPA: float,
        dt: float,
        is_seismic: bool = False,
    ):
    """
    Compute CPML damping parameters on a 2D grid.

    This routine constructs the conductivity, stretching, and auxiliary
    coefficients for a convolutional perfectly matched layer (CPML) in
    two dimensions, using user-specified maximum values and polynomial
    grading in the damping profiles.

    :param sig_max_x: Maximum conductivity values along the x boundaries.
                      Typically a 2-by-``nz`` array, with first and last
                      rows corresponding to the low and high x boundaries.
    :type sig_max_x: numpy.typing.NDArray[numpy.floating]
    :param sig_max_z: Maximum conductivity values along the z boundaries.
                      Typically an ``nx``-by-2 array, with first and last
                      columns corresponding to the low and high z boundaries.
    :type sig_max_z: numpy.typing.NDArray[numpy.floating]
    :param alpha_max: Maximum value of the CPML :math:`\\alpha` parameter.
    :type alpha_max: float
    :param kappa_max: Maximum value of the CPML :math:`\\kappa` stretching
                      parameter.
    :type kappa_max: float
    :param nx: Number of grid points in the x direction.
    :type nx: int
    :param nz: Number of grid points in the z direction.
    :type nz: int
    :param distancex: Distance profile in the x direction, used to grade
                      the CPML profiles. Expected to be normalized or
                      normalizable to :math:`[0, 1]`.
    :type distancex: numpy.typing.NDArray[numpy.floating]
    :param distancez: Distance profile in the z direction, used to grade
                      the CPML profiles. Expected to be normalized or
                      normalizable to :math:`[0, 1]`.
    :type distancez: numpy.typing.NDArray[numpy.floating]
    :param NP: Polynomial order for the :math:`\\sigma` and :math:`\\kappa`
               grading.
    :type NP: float
    :param NPA: Polynomial order for the :math:`\\alpha` grading.
    :type NPA: float
    :param dt: Time step size.
    :type dt: float
    :param is_seismic: If ``True``, use a seismic formulation for the
                       damping coefficients. If ``False``, use an
                       electromagnetic formulation that depends on the
                       permittivity constant ``eps0``.
    :type is_seismic: bool

    :returns: Tuple ``(sigma, kappa, alpha, acoeff, bcoeff)``, where each
              element is a 2D array of shape ``(nx, nz)``:
              ``sigma`` is the CPML conductivity profile,
              ``kappa`` is the CPML stretching profile,
              ``alpha`` is the CPML attenuation profile,
              ``acoeff`` and ``bcoeff`` are auxiliary CPML coefficients.
    :rtype: tuple[numpy.typing.NDArray[numpy.floating], \
    numpy.typing.NDArray[numpy.floating], \
    numpy.typing.NDArray[numpy.floating], \
    numpy.typing.NDArray[numpy.floating], \
    numpy.typing.NDArray[numpy.floating]]
    """
    kappa = np.ones([nx, nz])
    alpha = np.zeros([nx, nz])
    sigma = np.zeros([nx, nz])
    acoeff = np.zeros([nx, nz])
    bcoeff = np.zeros([nx, nz])
    
    distx = distancex / distancex.max()
    distz = distancez / distancez.max()

    rx = distx ** NP
    rz = distz ** NP
    m = len(distx)
    n = len(distz)

    for ii in range(m):
        for jj in range(n):
            # r = np.sqrt(distx[ii]**2 + distz[jj]**2)
            sigma[ii, jj] = sig_max_x[0, 0] * (rx[ii]) + sig_max_z[0, 0] * (rz[jj])
            sigma[-(ii + 1), jj] = sig_max_x[-1, 0] * (rx[ii]) + sig_max_z[-1, 0] * (rz[jj])
            sigma[ii, -(jj + 1)] = sig_max_x[0, -1] * (rx[ii]) + sig_max_z[0, -1] * (rz[jj])
            sigma[-(ii + 1), -(jj + 1)] = sig_max_x[-1, -1] * (rx[ii]) + sig_max_z[-1, -1] * (rz[jj])

            kappa[ii, jj] = (1.0 + (kappa_max - 1.0) * (rx[ii])) * (1.0 + (kappa_max - 1.0) * (rz[jj]))
            kappa[-(ii + 1), jj] = (1.0 + (kappa_max - 1.0) * (rx[ii])) * (1.0 + (kappa_max - 1.0) * (rz[jj]))
            kappa[ii, -(jj + 1)] = (1.0 + (kappa_max - 1.0) * (rx[ii])) * (1.0 + (kappa_max - 1.0) * (rz[jj]))
            kappa[-(ii + 1), -(jj + 1)] = (1.0 + (kappa_max - 1.0) * (rx[ii])) * (1.0 + (kappa_max - 1.0) * (rz[jj]))

            alpha[ii, jj] = alpha_max * ((1.0 - distx[ii]) ** NPA) * ((1.0 - distz[jj]) ** NPA)
            alpha[-(ii + 1), jj] = alpha_max * ((1.0 - distx[ii]) ** NPA) * ((1.0 - distz[jj]) ** NPA)
            alpha[ii, -(jj + 1)] = alpha_max * ((1.0 - distx[ii]) ** NPA) * ((1.0 - distz[jj]) ** NPA)
            alpha[-(ii + 1), -(jj + 1)] = alpha_max * ((1.0 - distx[ii]) ** NPA) * ((1.0 - distz[jj]) ** NPA)

    # Compute in the x, and z direction
    for ind in range(0, m):
        sigma[ind, m:-m] = sig_max_x[0, :] * (rx[ind])
        kappa[ind, :] = 1.0 + (kappa_max - 1.0) * (rx[ind])
        alpha[ind, :] = alpha_max * (1 - distx[ind]) ** NPA

        # From the end
        sigma[-(ind + 1), m:-m] = sig_max_x[-1, :] * rx[ind]
        kappa[-(ind + 1), :] = 1 + (kappa_max - 1) * (rx[ind])
        alpha[-(ind + 1), :] = alpha_max * (1 - distx[ind]) ** NPA

    for ind in range(0, n):
        sigma[n:-n, ind] = sig_max_z[:, 0] * (rz[ind])
        kappa[:, ind] = 1.0 + (kappa_max - 1.0) * (rz[ind])
        alpha[:, ind] = alpha_max * (1 - distz[ind]) ** NPA

        sigma[n:-n, -(ind + 1)] = sig_max_z[:, -1] * rz[ind]
        kappa[:, -(ind + 1)] = 1 + (kappa_max - 1) * (rz[ind])
        alpha[:, -(ind + 1)] = alpha_max * (1 - distz[ind]) ** NPA

    if is_seismic:
        bcoeff = np.exp(-(sigma / kappa + alpha) * dt)
    else:
        bcoeff = np.exp(-(sigma / kappa + alpha) * (dt / eps0))

    # Compute the a-coefficients
    alpha[np.where(alpha < 0.0)] = 0.0
    indices = np.where(np.abs(sigma) > 1.0e-6)
    acoeff[indices] = sigma[indices] * (bcoeff[indices] - 1) / (
        kappa[indices] * sigma[indices] + kappa[indices] * alpha[indices]
    )

    return sigma, kappa, alpha, acoeff, bcoeff

# ------------------------------------------------------------------------------
def sponge_boundary(
        domain, eta_max: float, eta: NDArray[np.floating]
    ) -> NDArray[np.floating]:
    """
    Apply a polynomial sponge boundary to a 2D damping field.
    
    A polynomial sponge profile is added near the domain boundaries to
    gradually damp outgoing waves. The strength and thickness of the
    sponge are controlled by attributes of ``domain``.
    
    :param domain: Simulation domain object providing sponge parameters.
                   Must define ``cpml`` (sponge thickness in grid points)
                   and ``NPS`` (polynomial order for the sponge profile).
    :type domain: Any
    :param eta_max: Maximum sponge damping coefficient applied at the
                    boundaries.
    :type eta_max: float
    :param eta: Existing 2D damping field to be modified in place by
                adding the sponge boundary contribution.
    :type eta: numpy.typing.NDArray[numpy.floating]
    :returns: The updated damping field with sponge boundary applied.
    :rtype: numpy.typing.NDArray[numpy.floating]
    """
    distx = np.arange(0, domain.cpml) / (domain.cpml - 1)
    distz = np.arange(0, domain.cpml) / (domain.cpml - 1)
    
    m = len(distx)
    n = len(distz)
    
    for ii in range(m):
        for jj in range(n):
            eta[ii, jj]           += eta_max * ((1.0 - distx[ii]) ** domain.NPS) * ((1.0 - distz[jj]) ** domain.NPS)
            eta[-(ii + 1), jj]    += eta_max * ((1.0 - distx[ii]) ** domain.NPS) * ((1.0 - distz[jj]) ** domain.NPS)
            eta[ii, -(jj + 1)]    += eta_max * ((1.0 - distx[ii]) ** domain.NPS) * ((1.0 - distz[jj]) ** domain.NPS)
            eta[-(ii + 1), -(jj + 1)] += eta_max * ((1.0 - distx[ii]) ** domain.NPS) * ((1.0 - distz[jj]) ** domain.NPS)
    
    # Left/right
    for ii in range(domain.cpml):
        eta[ii, :]      += eta_max * (1.0 - distx[ii]) ** domain.NPS
        eta[-(ii + 1), :] += eta_max * (1.0 - distx[ii]) ** domain.NPS
    
    # Top/bottom
    for jj in range(domain.cpml):
        eta[:, jj]      += eta_max * (1.0 - distz[jj]) ** domain.NPS
        eta[:, -(jj + 1)] += eta_max * (1.0 - distz[jj]) ** domain.NPS
    
    return eta

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
def compute_envelope(signal: ArrayLike, pad_width: int = 50) -> NDArray[np.floating]:
    """
    Compute the amplitude envelope of a 1D signal using the Hilbert transform.
    
    The signal is first padded by reflection at both ends to reduce edge
    effects in the Hilbert transform. The envelope is then computed on
    the padded signal and the padding is removed before returning.

    :param signal: Real-valued input signal.
    :type signal: ArrayLike
    :param pad_width: Number of samples to pad on each side using
                      reflection. The input length must be at least
                      ``2 * pad_width``.
    :type pad_width: int
    :returns: Amplitude envelope of the input signal with the same length
              as the original (unpadded) signal.
    :rtype: numpy.typing.NDArray[numpy.floating]
    :raises ValueError: If the input signal is shorter than
                        ``2 * pad_width``.
    
    .. note::
    
       Reflection padding helps mitigate boundary artifacts introduced
       by the Hilbert transform. If the signal is very short or strongly
       non-stationary near the edges, you may need to adjust
       ``pad_width`` to obtain a stable envelope.

    """
    if len(signal) < 2 * pad_width:
        raise ValueError("Signal too short for the given pad_width.")
    padded = np.pad(signal, pad_width, mode="reflect")
    analytic = hilbert(padded)
    envelope = np.abs(analytic)
    # Remove the padded regions.
    return envelope[pad_width:-pad_width]

# ------------------------------------------------------------------------------
def polarization_analysis(
        data: NDArray[np.floating],
        dt: float,
        M: int,
        alpha: float = 0.1,
    ) -> tuple[
        NDArray[np.floating],
        NDArray[np.floating],
        NDArray[np.floating],
        NDArray[np.floating],
        NDArray[np.bool_],
    ]:
    """
    Compute rectilinearity, backazimuth, and incidence angle using a sliding Tukey window.

    The analysis uses a sliding window over 3-component motion data to
    estimate polarization attributes via covariance eigen-decomposition.
    Near the signal edges, shorter windows are used, and a cone-of-influence
    mask is returned to indicate where the full window length is available.

    :param data: 3-component motion data with shape ``(N, 3)`` (e.g., LQT,
                 ZRT, or XYZ).
    :type data: numpy.typing.NDArray[numpy.floating]
    :param dt: Time step in seconds.
    :type dt: float
    :param M: Window length in samples.
    :type M: int
    :param alpha: Tukey window shape parameter in ``[0, 1]``. Default is
                  ``0.1``.
    :type alpha: float
    :returns: A tuple ``(times, rectilinearity, backazimuth, incidence,
              cone_mask)`` where:
              ``times`` is the time vector in seconds,
              ``rectilinearity`` contains rectilinearity values,
              ``backazimuth`` contains backazimuth angles in degrees from
              the y-axis,
              ``incidence`` contains incidence angles in degrees from the
              vertical,
              ``cone_mask`` is a boolean array that is ``True`` where a
              full-length window is used.
    :rtype: tuple[
        numpy.typing.NDArray[numpy.floating],
        numpy.typing.NDArray[numpy.floating],
        numpy.typing.NDArray[numpy.floating],
        numpy.typing.NDArray[numpy.floating],
        numpy.typing.NDArray[numpy.bool_]
    ]

    .. note::

       The polarization attributes are computed from the eigenvectors and
       eigenvalues of the windowed covariance matrix. The principal
       eigenvector is forced to point in a consistent (upward) direction
       by flipping it when its vertical component is negative.

       Near the series boundaries, the effective window is shorter than
       ``M`` samples. In these regions, the estimates may be less stable;
       the ``cone_mask`` output can be used to restrict interpretation to
       times where the full window length is available.

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
def compute_polar_energy(
        E_x: ArrayLike,
        E_y: ArrayLike,
        dt: float,
        angles: ArrayLike,
        to_db: bool = False,
        corrected: bool = False,
    ) -> NDArray[np.floating]:
    """
    Compute the integrated energy of the horizontal field projection.

    The horizontal field is projected onto multiple polarization angles
    for a 3-component time series (using the x and y components), and the
    energy is integrated over time for each projection angle.

    :param E_x: 1D array of the x-component time series.
    :type E_x: ArrayLike
    :param E_y: 1D array of the y-component time series.
    :type E_y: ArrayLike
    :param dt: Time-step interval in seconds.
    :type dt: float
    :param angles: 1D array of polarization angles in radians at which to
                   project the horizontal field.
    :type angles: ArrayLike
    :param to_db: If ``True``, convert energies to dB scale using
                  ``10 * log10(E / E.max())``. Default is ``False``.
    :type to_db: bool
    :param corrected: If ``True``, apply a cosine-squared baseline
                      correction based on the projection geometry and
                      return the corrected energies in dB. Default is
                      ``False``.
    :type corrected: bool
    :returns: Array of shape ``(len(angles),)`` containing the energy
              at each projection angle. If ``to_db`` is ``True`` and
              ``corrected`` is ``False``, the values are in dB relative
              to the maximum energy. If ``corrected`` is ``True``, the
              returned values are baseline-corrected and expressed in dB.
    :rtype: numpy.typing.NDArray[numpy.floating]

    .. note::
    
       For each angle :math:`\\theta`, the horizontal field is projected
       as :math:`E_\\theta = E_x \\cos(\\theta) + E_y \\sin(\\theta)`,
       and the energy is computed as the discrete integral
       :math:`\\sum |E_\\theta|^2 \\Delta t`.

       When ``corrected`` is ``True``, the raw energy pattern is divided
       by a cosine-squared baseline (``cos(angles)**2``) to account for
       purely geometric projection effects before converting to dB. This
       can be useful when comparing polarization strength across angles.

    """
    E_x = np.asarray(E_x)
    E_y = np.asarray(E_y)
    angles = np.asarray(angles)
    energies = np.empty(len(angles))
    
    for i, theta in enumerate(angles):
        # Project onto the polarization direction theta
        E_theta = E_x * np.cos(theta) + E_y * np.sin(theta)
        # Compute energy = integral of |E|^2 dt
        energies[i] = np.sum(np.abs(E_theta) ** 2) * dt
    
    if to_db:
        # Normalize to max and convert to dB
        energies = 10 * np.log10(energies / np.max(energies))
    
    if corrected:
        baseline = np.cos(angles) ** 2
        corr = energies / baseline
        return 10 * np.log10(corr / corr.max())
    else:
        return energies

# ------------------------------------------------------------------------------
def compute_phase_gain(
        ts_ref: ArrayLike,
        ts: ArrayLike,
        f0: float,
        dt: float,
    ) -> tuple[float, float]:
    """
    Compute phase shift and gain between two time series at a specific frequency.
    
    The function computes the real FFT of both input time series, locates
    the frequency bin closest to ``f0``, and evaluates the complex transfer
    function between the spectra. From this transfer function, it returns
    the phase shift (in degrees) and gain (in decibels).
    
    :param ts_ref: Reference time series (1D array).
    :type ts_ref: ArrayLike
    :param ts: Comparison time series, with the same shape as ``ts_ref``.
    :type ts: ArrayLike
    :param f0: Center frequency in Hz at which to compute phase and gain.
    :type f0: float
    :param dt: Sampling interval in seconds.
    :type dt: float
    :returns: A tuple ``(phase_deg, gain_db)`` where:
              ``phase_deg`` is the phase shift of ``ts`` relative to
              ``ts_ref`` at ``f0``, in degrees, and
              ``gain_db`` is the gain in decibels, computed as
              ``20 * log10(|H|)`` with ``H`` the complex spectral ratio.
    :rtype: tuple[float, float]
    :raises ValueError: If the input time series do not have the same shape
                        or are empty.
    
    .. note::

       The frequency bin used for the calculation is chosen as the one
       whose center frequency is closest to ``f0`` in the real FFT
       frequency grid (as returned by ``numpy.fft.rfftfreq``). For more
       accurate estimates at arbitrary frequencies, a longer time series
       or additional spectral interpolation may be required.

       If the reference spectrum is (near) zero at the selected bin,
       the transfer function may be ill-conditioned and the resulting
       gain or phase can be unstable. It is recommended to check that
       the reference signal has sufficient energy near ``f0``.

    """
    ts_ref = np.asarray(ts_ref)
    ts = np.asarray(ts)
    
    if ts_ref.shape != ts.shape:
        raise ValueError("Input time series must have the same shape")
    if ts_ref.size == 0:
        raise ValueError("Input time series must be non-empty")
    
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
def rotate_to_zrt(
        data: NDArray[np.floating],
        source: ArrayLike | None = None,
        receiver: ArrayLike | None = None,
        direction: ArrayLike | None = None,
    ) -> NDArray[np.floating]:
    """
    Rotate seismogram data from Cartesian components to Z, R, T components.
    
    The input 3-component seismogram data (``vx, vy, vz``) are rotated into
    a local coordinate system defined by the vertical, radial, and transverse
    directions. The radial direction is determined either from a propagation
    direction vector or from the source–receiver geometry.
    
    :param data: Seismogram data in Cartesian components, with shape
                 ``(n_samples, 3)`` and components ordered as ``(vx, vy, vz)``.
    :type data: numpy.typing.NDArray[numpy.floating]
    :param source: Source coordinates as a 3-element array ``[x, y, z]``.
                   Used together with ``receiver`` when ``direction`` is not
                   provided.
    :type source: array_like or None
    :param receiver: Receiver coordinates as a 3-element array ``[x, y, z]``.
                     Used together with ``source`` when ``direction`` is not
                     provided.
    :type receiver: array_like or None
    :param direction: Propagation direction vector (3 elements). If provided,
                      its horizontal projection defines the radial direction.
    :type direction: array_like or None
    :returns: Rotated seismogram data with shape ``(n_samples, 3)`` and
              components ordered as ``(Z, R, T)``.
    :rtype: numpy.typing.NDArray[numpy.floating]
    :raises ValueError: If neither ``direction`` nor a valid ``source`` and
                        ``receiver`` pair is provided, if the direction
                        vector has zero length, or if the source–receiver
                        horizontal separation is effectively zero.
    
    .. note::

       If ``direction`` is provided, it is normalized and its horizontal
       :math:`(x, y)` components are used to define the radial direction.

       If ``direction`` is not provided, both ``source`` and ``receiver``
       must be given. The horizontal difference ``receiver - source`` is
       used to determine the azimuth and radial direction.

       The vertical unit vector is assumed to be :math:`[0, 0, 1]`, and the
       transverse direction is computed as the cross product of the vertical
       and radial directions, ensuring an approximately right-handed
       coordinate system.
    
    .. code-block:: python
    
       # Example: rotate from vx, vy, vz to Z, R, T using source/receiver
       data = np.random.randn(1000, 3)  # (vx, vy, vz)
       source = [0.0, 0.0, 0.0]
       receiver = [10.0, 5.0, 0.0]
        
       zrt = rotate_to_zrt(data, source=source, receiver=receiver)

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
            radial_unit = np.array(
                [horizontal_diff[0] / norm_horizontal, horizontal_diff[1] / norm_horizontal, 0.0]
            )
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

# ------------------------------------------------------------------------------
def rotate_to_qlt(
        data: NDArray[np.floating],
        source_location: ArrayLike,
        receiver_location: ArrayLike,
        backazimuth: float | None = None,
        incidence: float | None = None,
    ) -> tuple[float, float, NDArray[np.floating]]:
    """
    Rotate 3-component data from Cartesian (x, y, z) to (L, Q, T) components.

    The input 3-component data are rotated into the longitudinal (L),
    quasi-vertical (Q), and transverse (T) directions using geometry
    defined by the source and receiver locations and, optionally, given
    backazimuth and incidence angles.

    :param data: Input data with shape ``(N, 3)`` and columns ``[vx, vy, vz]``.
    :type data: numpy.typing.NDArray[numpy.floating]
    :param source_location: Source coordinates ``[x, y, z]``.
    :type source_location: array_like
    :param receiver_location: Receiver coordinates ``[x, y, z]``.
    :type receiver_location: array_like
    :param backazimuth: Backazimuth angle in radians (clockwise from the
                        positive y-axis). If ``None``, it is computed
                        from the source–receiver geometry.
    :type backazimuth: float or None
    :param incidence: Incidence (inclination) angle in radians, positive
                      downward from horizontal. If ``None``, it is computed
                      from the source–receiver geometry.
    :type incidence: float or None
    :returns: A tuple ``(ba, inc, lqt)`` where:
              ``ba`` is the backazimuth angle in radians,
              ``inc`` is the inclination angle in radians,
              ``lqt`` is the rotated data with shape ``(N, 3)`` and
              columns ``[L, Q, T]``.
    :rtype: tuple[float, float, numpy.typing.NDArray[numpy.floating]]

    .. note::

       When ``backazimuth`` or ``incidence`` is not provided, both are
       derived from the displacement vector ``receiver_location -
       source_location``. The inclination is taken as the arctangent of
       the vertical offset over the horizontal distance, and the
       backazimuth is measured clockwise from the positive y-axis.
        
       The rotation matrix is constructed so that the L direction is
       aligned with the propagation direction, Q lies in the vertical
       plane containing the propagation direction, and T is the
       horizontal transverse (SH) direction.

    .. code-block:: python
    
       # Example: rotate from vx, vy, vz to L, Q, T using geometry
       data = np.random.randn(1000, 3)  # [vx, vy, vz]
       src = np.array([0.0, 0.0, 0.0])
       rcv = np.array([10.0, 5.0, -2.0])
        
       ba, inc, lqt = rotate_to_qlt(data, src, rcv)
        
    """
    # Compute displacement vector (receiver - source)
    d = receiver_location - source_location
    dx, dy, dz = d
    
    # Horizontal distance and inclination
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
    R = np.array(
        [
            [sin_inc * sin_ba, sin_inc * cos_ba, cos_inc],      # L-direction
            [cos_inc * sin_ba, cos_inc * cos_ba, -sin_inc],     # Q-direction
            [cos_ba, -sin_ba, 0.0],                             # T-direction
        ]
    )
    
    # Apply rotation
    lqt = data @ R.T
    
    return ba, inc, lqt

# ------------------------------------------------------------------------------
def rotate_sdr(
        strike: float,
        dip: float,
        rake: float,
        is_degree: bool = True,
    ) -> NDArray[np.floating]:
    """
    Create a 3D rotation matrix from strike, dip, and rake angles.
    
    The rotation matrix is constructed from three sequential rotations
    about the global axes corresponding to strike, dip, and rake angles.
    This is commonly used to transform between fault-plane and global
    coordinate systems in seismological applications.

    :param strike: Strike angle. Interpreted as degrees if
                   ``is_degree`` is ``True``, otherwise radians.
    :type strike: float
    :param dip: Dip angle. Interpreted as degrees if ``is_degree`` is
                ``True``, otherwise radians.
    :type dip: float
    :param rake: Rake (slip) angle. Interpreted as degrees if
                 ``is_degree`` is ``True``, otherwise radians.
    :type rake: float
    :param is_degree: If ``True``, input angles are treated as degrees
                      and converted to radians internally. If
                      ``False``, angles are assumed to be in radians.
    :type is_degree: bool
    :returns: Rotation matrix of shape ``(3, 3)`` formed as
              ``Rz(strike) @ Ry(dip) @ Rx(rake)``.
    :rtype: numpy.typing.NDArray[numpy.floating]
    
    .. note::

       The rotation is applied as a product of three basic rotations:
       about the z-axis by ``strike``, about the y-axis by ``dip``, and
       about the x-axis by ``rake``. The precise physical interpretation
       depends on your coordinate and fault-plane conventions.

       When using this matrix to rotate vectors ``v`` in the original
       coordinate system, the transformed vectors can be obtained as
       ``v_rot = R @ v``.
    
    .. code-block:: python
    
       # Example: build a rotation matrix for a given focal mechanism
       strike, dip, rake = 210.0, 30.0, 90.0  # degrees
       R = rotate_sdr(strike, dip, rake)
    
       v = np.array([1.0, 0.0, 0.0])
       v_rot = R @ v  # rotated vector
    
    """
    if is_degree:
        s, d, r = np.deg2rad([strike, dip, rake])
    else:
        s, d, r = strike, dip, rake
    
    Rz = np.array(
        [
            [np.cos(s), -np.sin(s), 0.0],
            [np.sin(s), np.cos(s), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    
    Ry = np.array(
        [
            [np.cos(d), 0.0, np.sin(d)],
            [0.0, 1.0, 0.0],
            [-np.sin(d), 0.0, np.cos(d)],
        ]
    )
    
    Rx = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, np.cos(r), -np.sin(r)],
            [0.0, np.sin(r), np.cos(r)],
        ]
    )
    
    return Rz @ Ry @ Rx

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
        data: NDArray[np.floating],
        dt: float,
        dx: float,
        taper: bool = True,
        ntfft: Optional[int] = None,
        nxfft: Optional[int] = None,
        fmin: float = 0.0,
        fmax: Optional[float] = None,
        kmin: Optional[float] = None,
        kmax: Optional[float] = None,
        wavenumber_units: str = "cycles",
    ) -> Tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """
    Compute the slant-stack f–k power spectrum of a 2D shot gather.

    The input gather is assumed to have time along the first axis and
    space (offset) along the second axis. A 2D FFT is computed, and the
    power spectrum is returned as a function of frequency and spatial
    wavenumber with optional band-limiting in both dimensions.

    :param data: Input shot gather with shape ``(nt, nx)`` where ``nt`` is
                 the number of time samples and ``nx`` is the number of
                 spatial traces.
    :type data: numpy.typing.NDArray[numpy.floating]
    :param dt: Time sampling interval in seconds.
    :type dt: float
    :param dx: Spatial sampling interval in meters.
    :type dx: float
    :param taper: If ``True``, apply a separable Hann taper in time and
                  space to reduce spectral leakage. Default is ``True``.
    :type taper: bool
    :param ntfft: Number of time samples for the FFT. If ``None``, use
                  ``nt`` (no zero-padding).
    :type ntfft: int or None
    :param nxfft: Number of spatial samples for the FFT. If ``None``, use
                  ``nx`` (no zero-padding).
    :type nxfft: int or None
    :param fmin: Minimum frequency in Hz to include in the output.
    :type fmin: float
    :param fmax: Maximum frequency in Hz to include in the output. If
                 ``None``, include up to the Nyquist frequency.
    :type fmax: float or None
    :param kmin: Minimum absolute wavenumber to include in the output,
                 expressed in the units specified by ``wavenumber_units``.
                 If ``None``, defaults to ``0.0``.
    :type kmin: float or None
    :param kmax: Maximum absolute wavenumber to include in the output,
                 expressed in the units specified by ``wavenumber_units``.
                 If ``None``, defaults to the maximum available wavenumber.
    :type kmax: float or None
    :param wavenumber_units: Units of the output wavenumber axis:
                             ``"cycles"`` for cycles per meter (default)
                             or ``"radians"`` for radians per meter.
    :type wavenumber_units: str
    :returns: A tuple ``(freqs, ks, P)`` where:
              ``freqs`` is a 1D array of frequencies in Hz,
              ``ks`` is a 1D array of spatial wavenumbers in the requested
              units, and
              ``P`` is a 2D power spectrum array with shape
              ``(len(freqs), len(ks))``.
    :rtype: tuple[
        numpy.typing.NDArray[numpy.floating],
        numpy.typing.NDArray[numpy.floating],
        numpy.typing.NDArray[numpy.floating]
    ]
    
    .. note::
    
       The temporal frequency axis is constructed using
       ``numpy.fft.rfftfreq`` and restricted to the range
       ``[fmin, fmax]``. The spatial wavenumber axis is built from
       ``numpy.fft.fftfreq`` and ``fftshift``, and optionally converted
       to radians per meter.

       Only wavenumbers whose absolute value lies within
       ``[kmin, kmax]`` are retained, and the corresponding columns are
       selected from the shifted power spectrum.
    
    .. code-block:: python
    
       # Example: compute an f–k spectrum from synthetic data
       nt, nx = 1024, 64
       dt, dx = 0.004, 10.0  # 4 ms, 10 m
       data = np.random.randn(nt, nx)
    
       freqs, ks, P = compute_fk_spectrum(
           data,
           dt=dt,
           dx=dx,
           fmin=2.0,
           fmax=60.0,
           kmin=0.0,
           kmax=0.2,
           wavenumber_units="cycles",
       )
    
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
    Pfull = np.abs(F2) ** 2
    
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
    if wavenumber_units == "radians":
        ks_full = ks_full * 2 * np.pi
    
    # apply k bounds
    if kmin is None:
        kmin = 0.0
    if kmax is None:
        kmax = np.abs(ks_full).max()
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
        freqs: NDArray[np.floating],
        ks: NDArray[np.floating],
        P: NDArray[np.floating],
        fig: Optional[Figure] = None,
        ax: Optional[Axes] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        log_scale: bool = False,
        smooth: bool = False,
        smooth_sigma: Tuple[float, float] = (1.0, 1.0),
        mode_lines: Optional[Dict[str, Union[float, Tuple[float, str]]]] = None,
        wavenumber_units: str = "cycles",
    ) -> Tuple[Figure, Axes]:
    """
    Display an f–k power spectrum.
    
    The f–k power spectrum is shown as a 2D color mesh of power versus
    frequency and wavenumber, with options for logarithmic scaling,
    Gaussian smoothing, and overlaying theoretical mode-velocity lines
    in either cycles/m or rad/m.
    
    :param freqs: 1D array of frequencies in Hz, typically the output
                  from :func:`compute_fk_spectrum`.
    :type freqs: numpy.typing.NDArray[numpy.floating]
    :param ks: 1D array of spatial wavenumbers in either cycles/m or
               rad/m, matching the ``wavenumber_units`` setting.
    :type ks: numpy.typing.NDArray[numpy.floating]
    :param P: 2D power spectrum array with shape
              ``(len(freqs), len(ks))``.
    :type P: numpy.typing.NDArray[numpy.floating]
    :param fig: Existing Matplotlib :class:`matplotlib.figure.Figure`
                to draw into. If ``None``, a new figure is created.
    :type fig: matplotlib.figure.Figure or None
    :param ax: Existing Matplotlib :class:`matplotlib.axes.Axes` to draw
               into. If ``None``, a new axes is created.
    :type ax: matplotlib.axes.Axes or None
    :param vmin: Minimum value for the color scale. If ``None``, it is
                 determined automatically by Matplotlib.
    :type vmin: float or None
    :param vmax: Maximum value for the color scale. If ``None``, it is
                 determined automatically by Matplotlib.
    :type vmax: float or None
    :param log_scale: If ``True``, display power in decibels using
                      ``10 * log10(P)``. Default is ``False``.
    :type log_scale: bool
    :param smooth: If ``True``, apply a Gaussian filter to ``P`` before
                   plotting. Default is ``False``.
    :type smooth: bool
    :param smooth_sigma: Standard deviations of the Gaussian smoothing
                         kernel in the (frequency, wavenumber) directions.
    :type smooth_sigma: tuple[float, float]
    :param mode_lines: Optional dictionary defining theoretical mode
                       lines to overlay. Keys are labels; values are
                       either velocities (float) or ``(velocity, color)``
                       tuples. Velocities are in m/s.
    :type mode_lines: dict[str, float or tuple[float, str]] or None
    :param wavenumber_units: Units of ``ks`` and the horizontal axis:
                             ``"cycles"`` for cycles/m (default) or
                             ``"radians"`` for rad/m.
    :type wavenumber_units: str
    :returns: The figure and axes used for the plot.
    :rtype: tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
    
    .. note::
    
       When ``log_scale`` is ``True``, a small epsilon is added to the
       power spectrum before converting to dB to avoid taking the
       logarithm of zero. Smoothing can help visualize coherent
       energy trends but may blur fine spectral features.
    
       Mode-velocity lines are drawn as pairs of symmetric curves
       :math:`\\pm k(f)` for each velocity, computed using
       :math:`k = 2 \\pi f / v` for wavenumbers in radians/m or
       :math:`k = f / v` for cycles/m, depending on
       ``wavenumber_units``.
    
    .. code-block:: python
    
       # Example: plot an f–k spectrum with a fundamental mode line
       freqs, ks, P = compute_fk_spectrum(data, dt=0.004, dx=10.0)
    
       fig, ax = plot_fk_spectrum(
           freqs,
           ks,
           P,
           log_scale=True,
           smooth=True,
           smooth_sigma=(1.0, 2.0),
           mode_lines={"fundamental": 300.0},  # velocity in m/s
       )
    
       plt.show()
    
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
        cbar_label = "Power (dB)"
    else:
        cbar_label = "Power"
    
    pcm = ax.pcolormesh(ks, freqs, Pplot, shading="auto", vmin=vmin, vmax=vmax)
    fig.colorbar(pcm, ax=ax, label=cbar_label)
    
    # overlay mode lines
    if mode_lines:
        for label, spec in mode_lines.items():
            if isinstance(spec, (tuple, list)):
                v, col = spec
            else:
                v, col = spec, "white"
            # compute k_line in correct units
            if wavenumber_units == "radians":
                kline = 2 * np.pi * freqs / v
            else:
                kline = freqs / v
            ax.plot(kline, freqs, "--", color=col, lw=1, label=label)
            ax.plot(-kline, freqs, "--", color=col, lw=1)
        # ax.legend(loc='upper right')
    
    xlabel = "Wavenumber (rad/m)" if wavenumber_units == "radians" else "Wavenumber (1/m)"
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title("f–k Spectrum")
    
    return fig, ax

# ------------------------------------------------------------------------------
# SURFACE WAVES
def compute_dispersion(
        seismograms: NDArray[np.floating],
        dx: float,
        dt: float,
        fmin: float,
        fmax: float,
        nfreq: int = 50,
    ) -> Tuple[NDArray[np.floating], NDArray[np.floating]]:
    """
    Estimate phase-velocity dispersion from a 2D seismogram gather.
    
    A linear phase-vs-distance trend is fitted for each selected frequency
    bin, with :math:`2\\pi` branch correction applied to unwrap the overall
    phase slope. The resulting wavenumber is converted to phase velocity,
    which is then clamped to a physically plausible range.

    :param seismograms: Input seismogram gather with shape ``(nt, nx)``,
                        where ``nt`` is the number of time samples and
                        ``nx`` is the number of spatial traces.
    :type seismograms: numpy.typing.NDArray[numpy.floating]
    :param dx: Spatial sampling interval in meters.
    :type dx: float
    :param dt: Time sampling interval in seconds.
    :type dt: float
    :param fmin: Minimum frequency in Hz for dispersion estimation.
    :type fmin: float
    :param fmax: Maximum frequency in Hz for dispersion estimation.
    :type fmax: float
    :param nfreq: Number of frequency samples between ``fmin`` and
                  ``fmax`` to evaluate. Default is ``50``.
    :type nfreq: int
    :returns: A tuple ``(freqs, vs)`` where:
              ``freqs`` is a 1D array of frequencies in Hz, and
              ``vs`` is a 1D array of phase velocities in m/s for each
              frequency, with invalid or out-of-range estimates set to
              ``NaN``.
    :rtype: tuple[
        numpy.typing.NDArray[numpy.floating],
        numpy.typing.NDArray[numpy.floating]
    ]
    
    .. note::

       For each target frequency, the complex spectrum is extracted
       across all traces, the phase is unwrapped, and a linear fit
       of phase versus distance is computed. A :math:`2\\pi` branch
       correction ensures that the fitted slope is consistent over
       the aperture.
    
       The maximum physically resolvable velocity is taken as
       ``dx / dt`` based on the spatial and temporal sampling.
       Phase-velocity estimates outside ``(0, max_vel]`` are rejected
       and stored as ``NaN`` in the output.
    
    .. code-block:: python

       # Example: compute a dispersion curve from a synthetic gather
       nt, nx = 2048, 48
       dt, dx = 0.004, 5.0
       seismograms = np.random.randn(nt, nx)
    
       freqs, vs = compute_dispersion(
           seismograms,
           dx=dx,
           dt=dt,
           fmin=2.0,
           fmax=40.0,
           nfreq=60,
       )

    """
    nt, nx = seismograms.shape
    spec = np.fft.rfft(seismograms, axis=0)
    freq_axis = np.fft.rfftfreq(nt, dt)
    freqs = np.linspace(fmin, fmax, nfreq)
    vs = np.full(nfreq, np.nan)
    
    x = np.arange(nx) * dx
    L = x[-1] - x[0]
    # max physically resolvable velocity (Nyquist spatial): dx/dt
    max_vel = dx / dt
    
    for i, f in enumerate(freqs):
        idx = np.argmin(np.abs(freq_axis - f))
        ph = np.unwrap(np.angle(spec[idx, :]))
        raw = np.polyfit(x, ph, 1)[0]
        impl = raw * L
        branch = np.round(impl / (2 * np.pi))
        slope = raw - branch * (2 * np.pi / L)
        k = -slope
        v = 2 * np.pi * f / k
        # clamp to positive and below max_vel
        if v > 0 and v <= max_vel:
            vs[i] = v
        else:
            vs[i] = np.nan
    
    return freqs, vs

# ============= Slant-stack energy image and smooth picks =============
def compute_dispersion_image(
        seismograms: NDArray[np.floating],
        dx: float,
        dt: float,
        fmin: float,
        fmax: float,
        nv: int = 100,
        nfreq: int = 100,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
    ) -> Tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """
    Compute a MASW-style dispersion image and phase-velocity picks.
    
    A multi-channel analysis of surface waves (MASW) slant-stack is
    formed over a range of trial phase velocities and frequencies. The
    resulting f–v image is normalized, smoothed, and used to perform
    sub-bin parabolic peak picking of the dominant dispersion curve.

    :param seismograms: Input seismogram gather with shape ``(nt, nx)``,
                        where ``nt`` is the number of time samples and
                        ``nx`` is the number of spatial traces.
    :type seismograms: numpy.typing.NDArray[numpy.floating]
    :param dx: Spatial sampling interval in meters.
    :type dx: float
    :param dt: Time sampling interval in seconds.
    :type dt: float
    :param fmin: Minimum frequency in Hz for the dispersion image.
    :type fmin: float
    :param fmax: Maximum frequency in Hz for the dispersion image.
    :type fmax: float
    :param nv: Number of phase-velocity samples between ``vmin`` and
               ``vmax``. Default is ``100``.
    :type nv: int
    :param nfreq: Number of frequency samples between ``fmin`` and
                  ``fmax`` for the image and picks. Default is ``100``.
    :type nfreq: int
    :param vmin: Minimum phase velocity in m/s for the image. If
                 ``None``, it is estimated from initial phase picks.
    :type vmin: float or None
    :param vmax: Maximum phase velocity in m/s for the image. If
                 ``None``, it is estimated from initial phase picks.
    :type vmax: float or None
    :returns: A tuple ``(freqs, vs, image, picks)`` where:
              ``freqs`` is a 1D array of frequencies in Hz,
              ``vs`` is a 1D array of trial phase velocities in m/s,
              ``image`` is a 2D normalized dispersion image with shape
              ``(len(freqs), len(vs))``, and
              ``picks`` is a 1D array of sub-bin phase-velocity picks
              (m/s) at each frequency, with out-of-bounds values set
              to ``NaN``.
    :rtype: tuple[
        numpy.typing.NDArray[numpy.floating],
        numpy.typing.NDArray[numpy.floating],
        numpy.typing.NDArray[numpy.floating],
        numpy.typing.NDArray[numpy.floating]
    ]
    :raises ValueError: If the automatically or manually determined
                        velocity bounds satisfy ``vmin >= vmax``.
    
    .. note::
    
       Initial phase-velocity estimates are obtained via
       :func:`compute_dispersion` and used to set automatic bounds when
       ``vmin`` and ``vmax`` are not provided. If these picks are
       unavailable, fallback bounds are based on the sampling ratio
       ``dx / dt``.
    
       For each frequency, complex amplitudes are slant-stacked over
       trial velocities using an exponential steering term. The image
       is normalized per frequency, smoothed with a Gaussian filter,
       and a parabolic fit around the maximum value in each row is used
       to obtain sub-bin phase-velocity picks.
    
    .. code-block:: python
    
       # Example: build a dispersion image and picks from a gather
       nt, nx = 2048, 48
       dt, dx = 0.004, 5.0
       seismograms = np.random.randn(nt, nx)
    
       freqs, vs, image, picks = compute_dispersion_image(
           seismograms,
           dx=dx,
           dt=dt,
           fmin=2.0,
           fmax=40.0,
           nv=120,
           nfreq=80,
       )
    
    """
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
            vmin = dx / dt * 0.1
            vmax = dx / dt * 2.0
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
            image[i, j] = np.abs(np.sum(A * np.exp(1j * omega * x / v))) ** 2
    
    # normalize per-frequency
    image /= image.max(axis=1, keepdims=True) + np.finfo(float).eps
    # smooth image to reduce striping
    image = gaussian_filter(image, sigma=(1, 2))
    
    # sub-bin (parabolic) peak picking
    picks = np.zeros_like(freqs)
    dv = vs[1] - vs[0]
    for i in range(len(freqs)):
        row = image[i]
        imax = np.nanargmax(row)
        if 0 < imax < len(vs) - 1:
            y0, y1, y2 = row[imax - 1], row[imax], row[imax + 1]
            dp = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
            picks[i] = vs[imax] + dp * dv
        else:
            picks[i] = vs[imax]
    # mask outside bounds
    picks = np.where((picks >= vmin) & (picks <= vmax), picks, np.nan)
    
    return freqs, vs, image, picks

# ----------------------------------------------------------------------------
def plot_dispersion(
        freqs: NDArray[np.floating], velocities: NDArray[np.floating]
    ) -> None:
    """
    Plot a dispersion curve of frequency versus phase velocity.
    
    A simple line plot is created showing phase velocity as a function
    of frequency, suitable for visualizing dispersion picks or curves.

    :param freqs: 1D array of frequencies in Hz.
    :type freqs: numpy.typing.NDArray[numpy.floating]
    :param velocities: 1D array of phase velocities in m/s, with the
                       same shape as ``freqs``.
    :type velocities: numpy.typing.NDArray[numpy.floating]
    :returns: This function does not return anything; it creates a
              Matplotlib figure and displays it.
    :rtype: None
    
    .. note::
    
       This function calls :func:`matplotlib.pyplot.show`, which will
       block execution in some environments. If you need more control
       over the figure (e.g., saving to file or embedding in a GUI),
       consider modifying the code to return the figure and axes
       instead of calling ``plt.show()`` directly.
    
    .. code-block:: python
    
       # Example: plot dispersion picks
       freqs = np.linspace(2.0, 40.0, 50)
       velocities = 300.0 + 5.0 * freqs  # synthetic trend
    
       plot_dispersion(freqs, velocities)
    
    """
    plt.figure()
    plt.plot(freqs, velocities, "-o")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Phase velocity (m/s)")
    plt.title("Dispersion curve")
    plt.grid(True)
    plt.show()

# ----------------------------------------------------------------------------
def plot_dispersion_image(
        freqs: NDArray[np.floating],
        vs: NDArray[np.floating],
        image: NDArray[np.floating],
        picks_freqs: Optional[NDArray[np.floating]] = None,
        picks_vs: Optional[NDArray[np.floating]] = None,
        colormap: str = "inferno",
    ) -> Tuple[Figure, Axes]:
    """
    Plot a MASW-style dispersion image with an optional picked curve.
    
    The function displays a frequency–velocity (f–v) image of normalized
    energy and, if provided, overlays a picked dispersion curve on top
    of the image.
    
    :param freqs: 1D array of frequencies in Hz with shape ``(Nf,)``.
    :type freqs: numpy.typing.NDArray[numpy.floating]
    :param vs: 1D array of phase velocities in m/s with shape ``(Nv,)``.
    :type vs: numpy.typing.NDArray[numpy.floating]
    :param image: 2D normalized energy map with shape ``(Nf, Nv)``,
                  typically the output of :func:`compute_dispersion_image`.
    :type image: numpy.typing.NDArray[numpy.floating]
    :param picks_freqs: Optional 1D array of frequencies in Hz for a
                        picked dispersion curve.
    :type picks_freqs: numpy.typing.NDArray[numpy.floating] or None
    :param picks_vs: Optional 1D array of phase velocities in m/s for a
                     picked dispersion curve, matching ``picks_freqs``.
    :type picks_vs: numpy.typing.NDArray[numpy.floating] or None
    :param colormap: Matplotlib colormap name used for the image
                     (default ``"inferno"``).
    :type colormap: str
    :returns: The figure and axes containing the dispersion image plot.
    :rtype: tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
    
    .. note::
    
       The energy map is plotted using
       :func:`matplotlib.axes.Axes.pcolormesh` with transposed data
       so that frequency appears on the horizontal axis and phase
       velocity on the vertical axis.
    
       When both ``picks_freqs`` and ``picks_vs`` are provided, the
       picked curve is overlaid as red scatter points and a legend is
       added in the upper-right corner.
    
    .. code-block:: python
    
       # Example: visualize a dispersion image and picked curve
       freqs, vs, image, picks = compute_dispersion_image(
           seismograms,
           dx=5.0,
           dt=0.004,
           fmin=2.0,
           fmax=40.0,
       )
    
       fig, ax = plot_dispersion_image(
           freqs,
           vs,
           image,
           picks_freqs=freqs,
           picks_vs=picks,
           colormap="inferno",
       )
        
       plt.show()
    
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # plot the energy map
    pcm = ax.pcolormesh(freqs, vs, image.T, shading="auto", cmap=colormap)
    fig.colorbar(pcm, ax=ax, label="Normalized Energy")
    
    # labels and title
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Phase Velocity (m/s)")
    ax.set_title("Dispersion Image")
    
    # overlay picks if provided
    if picks_freqs is not None and picks_vs is not None:
        ax.scatter(
            picks_freqs,
            picks_vs,
            c="r",
            marker=".",
            s=10,
            label="Picked Curve",
        )
        ax.legend(loc="upper right")
    
    return fig, ax

# ----------------------------------------------------------------------------
def apply_time_taper(
        data: NDArray[np.floating],
        dt: float,
        t0: float = 0.18,
        t1: float = 0.25,
    ) -> NDArray[np.floating]:
    """
    Apply a cosine taper to suppress data after a given time window.
    
    A one-sided cosine taper is applied between ``t0`` and ``t1`` (in
    seconds), smoothly reducing the amplitudes to zero. All samples
    after ``t1`` are fully muted.
    
    :param data: Input data array with shape ``(nt, nx)``, where ``nt``
                 is the number of time samples and ``nx`` the number of
                 traces or components.
    :type data: numpy.typing.NDArray[numpy.floating]
    :param dt: Time sampling interval in seconds.
    :type dt: float
    :param t0: Start time of the taper in seconds. Default is ``0.18``.
    :type t0: float
    :param t1: End time of the taper in seconds. Default is ``0.25``.
               Samples after ``t1`` are set to zero.
    :type t1: float
    :returns: Tapered data array with the same shape as the input.
    :rtype: numpy.typing.NDArray[numpy.floating]
    :raises ValueError: If ``t1`` is less than or equal to ``t0``.
    
    .. note::
    
       The taper is equal to 1.0 before ``t0``, follows a cosine shape
       between ``t0`` and ``t1``, and is 0.0 after ``t1``. This is
       useful for suppressing late-arriving noise or boundary
       reflections without introducing sharp truncation artifacts.
    
    .. code-block:: python
    
       # Example: mute late-time energy in a gather
       nt, nx = 2000, 24
       dt = 0.004  # 4 ms
       data = np.random.randn(nt, nx)
    
       tapered = apply_time_taper(data, dt, t0=0.5, t1=0.8)
    
    """
    if t1 <= t0:
        raise ValueError(f"t1 must be greater than t0 (got t0={t0}, t1={t1})")
    
    nt, nx = data.shape
    time = np.arange(nt) * dt
    taper = np.ones_like(time)
    
    # Apply taper between t0 and t1
    mask = (time >= t0) & (time <= t1)
    taper[mask] = 0.5 * (1 + np.cos(np.pi * (time[mask] - t0) / (t1 - t0)))
    taper[time > t1] = 0.0
    
    return data * taper[:, None]
