import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg 
import matplotlib.animation as anim
import os 
import os.path 
from scipy.io import FortranFile
from subprocess import call
import json 
from typing import Optional
import trimesh

import seidart.routines.materials as mf
from seidart.routines import sourcefunction
from seidart.routines.definitions import *
from seidart.routines.prjbuild import readwrite_json
import seidart.routines.materials as mf

__all__ = [
    'Domain',
    'Material',
    'Model',
    'AnimatedGif'
]

# ------------------------------------------------------------------------------
class Domain:
    def __init__(self):
        super().__init__()
        self.build()

    def build(self):
        # Initialize variables
        self.geometry = None
        self.dim = None
        self.nx = None
        self.ny = None
        self.nz = None
        self.dx = None
        self.dy = None
        self.dz = None
        self.cpml = None
        self.cpml_attenuation = 0.0 # Default attenuation in the cpml region. Higher is more attenuation.
        # CPML parameters for the domain
        self.sig_opt_scalar = None
        self.alpha_opt_scalar = None 
        self.NP = None 
        self.NPA = None 
        self.NPS = 1.0 # Sponge exponent
        self.kappa_max = 5
        self.Rcoef = None # This is the reflection coefficient and has replaced sig_opt_scalar 
        # Some more values that might be of interest 
        self.write = None
        self.image_file = None
        self.exit_status = 1

        # Flag to specify whether the all inputs are fulfilled
        self.seismic_model = False
        self.electromag_model = False

    def parameter_check(self):
        self.exit_status = 0
        # make sure there's a geometry. This implies whether nx, and nz exist
        if self.geometry is None:
            self.exit_status = 1
            print('No geometry loaded\n')
        if not self.dx or self.dx is None:
            self.exit_status = 1
            print('No step size in the x-direction')
        if not self.dz or self.dz is None:
            self.exit_status = 1
            print('No step size in the z-direction')
        if not self.cpml or self.cpml is None:
            self.exit_status = 1
            print('No cpml thickness is given')

        # Check 2.5D
        if self.dim == '2.5' and exit_status == 0:
            if self.ny is None or self.ny == 'n/a':
                print('No dimension given for y-direction. Assigning default ny=3.')
                self.ny = 3
            if self.dy is None or self.dy == 'n/a':
                self.dy = np.min( [int(self.dx), int(self.dz)])
                print('No step size given for y-direction. Assigning min(dx,dz).')

        # Convert variables to required types. Ny and Dy are already converted if given
        if self.exit_status == 0:
            self.dim = float(self.dim)
            self.nx = int(self.nx)
            self.nz = int(self.nz)
            self.dx = float(self.dx)
            self.dz = float(self.dz)
            self.cpml = int(self.cpml)
        else:
            print('\n Domain inputs are not satisfied. I can"t go on anymore. \n')
            # quit()

# ------------------------------------------------------------------------------
class Material:
    """
    A class to manage materials for simulation purposes.

    Attributes
    ----------
    material_list : pd.DataFrame
        An array to store the list of materials.
        A flag indicating whether the materials were read in successfully.
    material : numpy.ndarray or None
        Stores material information.
    rgb : numpy.ndarray or None
        Stores RGB values for materials (not used currently).
    temp : numpy.ndarray or None
        Stores temperature values for materials.
    rho : numpy.ndarray or None
        Stores density values for materials.
    porosity : numpy.ndarray or None
        Stores porosity values for materials.
    lwc : numpy.ndarray or None
        Stores liquid water content values for materials.
    is_anisotropic : numpy.ndarray or None
        Indicates if the material is anisotropic.
    angfile : numpy.ndarray or None
        Stores ANG file paths for anisotropic materials.
    functions : Any
        Stores material processing functions.
    """
    # initialize the class
    def __init__(self) -> None:
        """
        Initializes the Material class by building the initial structure.
        """
        super().__init__()
        self.build()

    def build(self) -> None:
        """
        Initializes the material attributes with default values.
        """
        self.material_list = None

        # We will assign each of the list variables
        self.material = None
        self.rgb = None
        self.temp = None
        # self.attenuation = None
        self.rho = None
        self.porosity = None
        self.lwc = None
        self.is_anisotropic = None
        self.angfile = None
        # self.permittivity_coefficients = None 
        # self.conductivity_coefficients = None

        # The processing functions
        self.functions = mf

    def sort_material_list(self) -> None:
        """
        Sorts the material list based on the material properties.
        """
        self.material = self.material_list['Name'].to_numpy()
        self.rgb = self.material_list['RGB'].to_numpy()
        self.temp = self.material_list['Temperature'].to_numpy()
        self.rho = self.material_list['Density'].to_numpy()
        self.porosity = self.material_list['Porosity'].to_numpy()
        self.lwc = self.material_list['Water_Content'].to_numpy()
        self.is_anisotropic = self.material_list['is_anisotropic'].to_numpy()
        self.angfile = self.material_list['ANGfiles'].to_numpy()
        

    def parameter_check(self) -> None:
        """
        Checks the parameters of the material list for completeness.

        It ensures that necessary fields are provided and checks for the 
        presence of .ANG files for anisotropic materials.
        """
        # The fields for the materials in the input are defined as:
        # 'id, R/G/B, Temp., Dens., Por., WC, Anis, ANG_File'
        # but the R/G/B column is deleted

        if len(self.material_list) > 0:
            # Check to make sure the necessary fields are provided
            check = 0

            for row in self.material_list[:,0:6]:
                for val in row:
                    if not val:
                        check = check + 1

        if check == 0:
            file_check = 0
            for ind in range(0, self.material_list.shape[0]):
                if self.material_list['is_anisotropic'][ind] == 'True' and \
                    not self.material_list[ind,7] or \
                        self.material_list[ind,7] == 'n/a':
                    file_check = file_check + 1
        else:
            print('Material inputs aren"t satisfied.')

        if check == 0:
            if file_check != 0:
                print('No .ANG file specified for anisotropic material')
    
# ------------------------------------------------------------------------------
class Model:
    """
    A class to manage the simulation model configuration.

    Attributes
    ----------
    dt : float or None
        The time step size.
    time_steps : int or None
        The total number of time steps in the simulation.
    x : float or None
        The x-coordinate of the source location.
    y : float or None
        The y-coordinate of the source location. Optional, defaults to None.
    z : float or None
        The z-coordinate of the source location.
    f0 : float or None
        The source frequency.
    theta : float or None
        The angle of incidence in the xz-plane. Optional, defaults to 0 if 
        unspecified.
    phi : float or None
        The angle of incidence in the xy-plane. Optional, defaults to 0 if `y` 
        is specified and `phi` is unspecified.
    src : Any
        The source information. Type is unspecified.
    tensor_coefficients : numpy.ndarray or None
        The tensor coefficients for the simulation. Optional, but required for 
        tensor-based simulations.
    compute_coefficients : bool
        A flag indicating whether to compute coefficients. Defaults to True.
    attenuation_coefficients : numpy.ndarray or None
        The attenuation coefficients for the simulation. Optional.
    fref : float or None
        The reference frequency for attenuation. Optional.
    attenuation_fadjust : bool or None
        Flag to adjust frequency for attenuation. Optional.
    exit_status : int
        Status code to indicate the success or failure of parameter checks.
    """
    def __init__(self) -> None:
        """
        Initializes the Model class by building the initial configuration.
        """
        super().__init__()
        self.project_file = None
        self.dt = None
        self.time_steps = None
        self.x = None
        self.xind = None
        self.y = None
        self.yind = None
        self.z = None
        self.zind = None
        self.f0 = None
        self.theta = None
        self.phi = None
        self.source_amplitude = None 
        self.source_wavelet = None
        self.source_type = None
        self.src = None
        self.moment_tensor = None
        self.attenuation_fadjust = None
        self.exit_status = 0
        self.is_seismic = None
        self.step_limit = None
        self.sourcefunction_x = None 
        self.sourcefunction_y = None 
        self.sourcefunction_z = None
        self.initialcondition_x = None 
        self.initialcondition_y = None 
        self.initialcondition_z = None 
        self.compute_coefficients = True
        self.domain_density = None
        self.air_gradient_integer = 2
        self.use_multimodal = False
        self.center = False
        self.use_broadband = False 
        self.broadband_fmin = None
        self.broadband_fmax = None
        self.source_n_octaves = 2
        self.sourcefunction = sourcefunction.pointsource # This will change
        self.get_seismic = mf.get_seismic 
        self.get_perm = mf.get_perm
        self.CFL = 1/np.sqrt(2) # Default minimum. If the dimension is 2.5, then this will change automatically
        self.density_method = 'geometric' # options are 'none', 'harmonic', 'geometric', and 'arithmetic'
        self.vmax_x = None
        self.vmax_z = None
        self.velocity_scaling_factor = 1.0 
        self.voigt_to_full_tensor = mf.voigt_to_full_tensor
        
    def kband_check(self, domain):
        '''
        Check to make sure the discretization values of the domain satisfy the 
        wavenumber bandlimiter as described in Volvcan et al. (2024).

        :param material: 
        :type material: Material 
        :param domain:
        :type domain: Domain
        '''
        
        if self.is_seismic:
            velocities = np.zeros([domain.nmats, 4])
            for ind in range(domain.nmats):
                row = self.stiffness_coefficients.loc[ind]
                rho = row['rho']
                tensor = row.drop('rho').to_numpy()
                velocities[ind,:] = mf.tensor2velocities(
                    tensor, rho,
                    seismic = True
                )
        else: 
            velocities = np.zeros([domain.nmats, 3])
            for ind in range(domain.nmats):
                velocities[ind,:] = mf.tensor2velocities(
                    self.permittivity_coefficients.loc[ind].to_numpy().real, seismic = False
                )
        
        lambda_values = velocities/self.f0
        # For things such as air and water, there will be zero valued velocities
        # which will be taken as the step_limit. Remove the zero values.
        lambda_values = lambda_values[lambda_values > 0]
        step_limit = lambda_values.min()/4 
        self.step_limit = step_limit 
        if domain.dim == 2.5:
            step_max = np.array([domain.dx, domain.dy, domain.dz]).max() 
        else: 
            step_max = np.array([domain.dx, domain.dz]).max()
        
        if step_max > step_limit:
            print(
                f'Wavenumber bandlimit is not met. Reduce maximum spatial step to {step_limit}'
            )
        else:
            print(
                f'Wavenumber bandlimit is satisfied for step limit = {step_limit}'
            )
    
    # --------------------------------------------------------------------------
    def parameter_check(self) -> None:
        """
        Performs parameter checks for essential simulation settings and updates 
        the exit status accordingly.
        """
        if not self.time_steps:
            self.exit_status = 1
            print('Number of time steps aren"t satisfied.')
        
        if not self.x or not self.z:
            self.exit_status = 1
            print('No source location is specified.')
        
        if not self.f0:
            self.exit_status = 1
            print('No source frequency is specified.')
        
        # in case theta or phi aren't specified we can assign defaults
        if not self.theta:
            self.theta = 0
        
        # if y is specified but not phi
        if self.y and not self.phi:
            self.phi = 0
    
    # --------------------------------------------------------------------------
    def build(
            self, 
            material: Material, 
            domain: Domain, 
            recompute_tensors: bool = True,
            write_tensor: bool = True,
        ) -> None:
        """
        Checks the status of the modeling classes and appends coefficients to the
        project file if necessary.

        :param material: The material class instance.
        :type material: Material
        :param domain: The domain class instance.
        :type domain: Domain
        :param append_to_json: Flag to append coefficients to the project file. 
                                Defaults to True.
        :type append_to_json: bool
        """
        global clight
        
        self.xind = int(self.x / domain.dx )
        self.yind = int(self.y / domain.dy )
        self.zind = int(self.z / domain.dz )
                
        # Set the CFL to its max value if it exceeds the 
        if domain.dim == 2.5:
            self.CFL = np.min([self.CFL, 1/np.sqrt(3)])
        
        # Make sure all of the initial conditions are input 
        if self.initialcondition_x  is None:
            if domain.dim == 2.5:
                self.initialcondition_x = np.zeros([
                    domain.nx+2*domain.cpml, 
                    domain.ny+2*domain.cpml, 
                    domain.nz+2*domain.cpml
                ])
            else:
                self.initialcondition_x = np.zeros([
                    domain.nx+2*domain.cpml, domain.nz+2*domain.cpml
                ])
        if self.initialcondition_z  is None:
            if domain.dim == 2.5:
                self.initialcondition_z = np.zeros([
                    domain.nx+2*domain.cpml, 
                    domain.ny+2*domain.cpml, 
                    domain.nz+2*domain.cpml
                ])
            else:
                self.initialcondition_z = np.zeros([
                    domain.nx+2*domain.cpml, domain.nz+2*domain.cpml
                ])
        if self.initialcondition_y  is None:
            if domain.dim == 2.5:
                self.initialcondition_y = np.zeros([
                    domain.nx+2*domain.cpml, 
                    domain.ny+2*domain.cpml, 
                    domain.nz+2*domain.cpml
                ])
            else:
                self.initialcondition_y = np.zeros([
                    domain.nx+2*domain.cpml, domain.nz+2*domain.cpml
                ])
        
        # ----------------------
        if self.exit_status == 0 and recompute_tensors:
            # The coefficients aren't provided but the materials are so we can 
            # compute them
            # Assign the materials to their respective corners
            if self.is_seismic:
                print('Computing the stiffness coefficients.')
                self.get_seismic(self, material)
            else:
                print('Computing the permittivity and conductivity coefficients.')    
                self.get_perm(self, material)
        
        # Always recalculate dt
        self.compute_max_velocities(dim = domain.dim) 
        max_vel = max(self.max_velocity_per_material.values() ) 
        
        if not self.is_seismic:
            max_vel = clight/ \
                np.sqrt(self.permittivity_coefficients[['e11', 'e22', 'e33']].min().min() ) 
        
        if domain.dim == 2.5 or domain.dim == 3:
            denom  = np.sqrt((1/domain.dx)**2 + (1/domain.dy)**2 + (1/domain.dz)**2) 
        else:
            denom = np.sqrt((1/domain.dx)**2 + (1/domain.dz)**2) 
        
        self.dt = self.CFL / ( max_vel * denom )  
                  
        src_nyquist = 1/(2*self.f0)
        if src_nyquist < self.dt:
            print(
                '''Nyquist is not small enough for the source frequency. Change
                the source frequency or decrease the spatial step size'''
            )
        
        # ----------------------
        # We need to set the 
        print("Creating the source function")
        (
            __, self.sourcefunction_x, self.sourcefunction_y, 
            self.sourcefunction_z, __ 
        ) = self.sourcefunction(
            self, 
            broadband = self.use_broadband,
            multimodal = self.use_multimodal,
            center = self.center,
            fmin = self.broadband_fmin,
            fmax = self.broadband_fmax,
            num_octaves = self.source_n_octaves
        )
                
        direction = ['x', 'y', 'z']
        # Compute CPML
        print('Computing CPML boundary values and writing outputs to Fortran files.')
        for d in direction:
            __ = cpmlcompute(
                self, domain, 
                velocity_scaling_factor = self.velocity_scaling_factor, 
            )
            __ = cpmlcompute(self, domain, half = True,
                velocity_scaling_factor = self.velocity_scaling_factor, 
            )
        
        # ----------------------
        # Write out the tensor components to file
        if write_tensor:
            print('Writing tensor components to individual .dat files.')
            # Compute the density gradient at boundaries with air. If the model is a
            # a single material, this will return an error. Booooo, errors!
            if domain.nmats > 1:
                self.domain_density = airsurf(material, domain, self.air_gradient_integer)
                self.domain_density = np.ones([domain.nx, domain.nz]) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            else:
                self.domain_density = np.ones([domain.nx, domain.nz])
            
            # Write out the arrays of tensor coefficients
            if self.is_seismic:
                self.tensor2dat(self.stiffness_coefficients, domain)
                # gamma/attenuation coefficients need to be scaled by the 
                # dominant frequency of the source
                self.atten2dat(self.attenuation_coefficients * self.f0, domain, self.gamma_max)
            else:
                self.tensor2dat(self.permittivity_coefficients, domain) 
                self.tensor2dat(self.conductivity_coefficients, domain)
        else:
            print('Tensor components have already been written to .dat files.')
            return
        
    def writemask(self, material: Material, domain: Domain):
        """
        
        """
        # Define any materials that need to be masked
        masked_materials = ['air', 'water', 'oil']
        
        masked_array = np.zeros([domain.nx, domain.nz])
        extended_masked_array = np.zeros(
            [domain.nx + 2*domain.cpml, domain.nz + 2*domain.cpml]
        )
        for mm in masked_materials:
            material_id = material.material_list['id'].loc[material.material_list['Name'] == 'air'].to_numpy()
            masked_array = masked_array + (domain.geometry == material_id).astype(int)
        
        masked_array[masked_array > 0] = 1
        extended_masked_array[
                domain.cpml:domain.nx+domain.cpml,domain.cpml:domain.nz+domain.cpml
        ] = masked_array
        extended_masked_array[0:domain.cpml,:] = extended_masked_array[domain.cpml+1,:]
        extended_masked_array[domain.nx+domain.cpml:,:] = extended_masked_array[domain.nx+domain.cpml-1,:]
        extended_masked_array[:,0:domain.cpml] = extended_masked_array[:,domain.cpml+1].reshape(-1,1)
        extended_masked_array[:,domain.nz+domain.cpml:] = extended_masked_array[:,domain.nz+domain.cpml-1].reshape(-1,1)
        f = FortranFile('geometry_mask.dat', 'w')
        f.write_record(extended_masked_array.T)
        f.close()
        self.geometry_mask = extended_masked_array
        return 
        
    def tensor2dat(self, tensor, domain):
        """

        """
        columns = tensor.columns
        for col in columns:
            fn = col + '.dat'
            
            if domain.dim == 3.0:
                coef_array = np.zeros([domain.nx, domain.ny, domain.nz])
                extended_array = np.zeros([domain.nx+2*domain.cpml, domain.nz+2*domain.cpml])
                for ii in range(domain.nx):
                    for jj in range(domain.ny):
                        for kk in range(domain.nz):
                            coef_array[ii,jj,kk] = tensor[col][domain.geometry[ii,jj,kk]]
                
                # Extend values into the pml
                extended_array[
                        domain.cpml:domain.nx+domain.cpml,
                        domain.cpml:domain.ny+domain.cpml,
                        domain.cpml:domain.nz+domain.cpml
                    ] = coef_array
                extended_array[0:domain.cpml,:,:] = extended_array[domain.cpml+1,:,:]
                extended_array[domain.nx+domain.cpml:,:,:] = extended_array[domain.nx+domain.cpml-1,:,:]
                extended_array[:,0:domain.cpml,:] = extended_array[:,domain.cpml+1,:].reshape(-1,1)
                extended_array[:,domain.ny+domain.cpml:,:] = extended_array[:,domain.ny+domain.cpml-1,:].reshape(-1,1)
                extended_array[:,:,0:domain.cpml] = extended_array[:,:,domain.cpml+1].reshape(-1,1)
                extended_array[:,:,domain.nz+domain.cpml:] = extended_array[:,:,domain.nz+domain.cpml-1].reshape(-1,1)
            
            else:
                coef_array = np.zeros([domain.nx, domain.nz])
                extended_array = np.zeros([domain.nx+2*domain.cpml, domain.nz+2*domain.cpml])
                for ii in range(domain.nx):
                    for jj in range(domain.nz):
                        coef_array[ii,jj] = tensor[col][domain.geometry[ii,jj]]
                # Extend values into the pml
                extended_array[
                        domain.cpml:domain.nx+domain.cpml,domain.cpml:domain.nz+domain.cpml
                    ] = coef_array
                extended_array[0:domain.cpml,:] = extended_array[domain.cpml+1,:]
                extended_array[domain.nx+domain.cpml:,:] = extended_array[domain.nx+domain.cpml-1,:]
                extended_array[:,0:domain.cpml] = extended_array[:,domain.cpml+1].reshape(-1,1)
                extended_array[:,domain.nz+domain.cpml:] = extended_array[:,domain.nz+domain.cpml-1].reshape(-1,1)
            
            f = FortranFile(fn, 'w')
            f.write_record(extended_array.T)
            f.close()
        
    def atten2dat(self, tensor, domain, eta_max = None):
        """

        """
        columns = tensor.columns
        coef_array = np.zeros([domain.nx, domain.nz])
        extended_array = np.zeros([domain.nx+2*domain.cpml, domain.nz+2*domain.cpml])
        for col in columns[:-1]:
            fn = col + '.dat' 
            for ii in range(domain.nx):
                for jj in range(domain.nz):
                    coef_array[ii,jj] = tensor[col][domain.geometry[ii,jj]]
            
            # Extend values into the pml
            extended_array[
                    domain.cpml:domain.nx+domain.cpml,domain.cpml:domain.nz+domain.cpml
                ] = coef_array
            extended_array[0:domain.cpml,:] = extended_array[domain.cpml+1,:]
            extended_array[domain.nx+domain.cpml:,:] = extended_array[domain.nx+domain.cpml-1,:]
            extended_array[:,0:domain.cpml] = extended_array[:,domain.cpml+1].reshape(-1,1)
            extended_array[:,domain.nz+domain.cpml:] = extended_array[:,domain.nz+domain.cpml-1].reshape(-1,1)
            if eta_max:
                extended_array = sponge_boundary(domain, self.gamma_max, extended_array)
                        
            f = FortranFile(fn, 'w')
            f.write_record(extended_array.T)
            f.close()
        
        return(extended_array)

    # --------------------------------------------------------------------------
    def add_noise(self, domain, scalar_amplitude: float = 1e0):
        """
        
        :param domain: The domain instance  
        :type domain: Domain
        :param scalar_amplitude: The scalar value to multiply the noise by. The
            range of noise values are between 0 and 1. 
        :type scalar_amplitude: float
        
        """
        # Create a gaussian noise model and superpose it onto the initial 
        # conditions
        if domain.dim == 2.5:
            self.initialcondition_x = self.initialcondition_x + \
                scalar_amplitude * \
                    np.random.normal(0, 1, [domain.nx, domain.ny, domain.nz])
            self.initialcondition_y = self.initialcondition_y + \
                scalar_amplitude * \
                    np.random.normal(0, 1, [domain.nx, domain.ny, domain.nz])
            self.initialcondition_z = self.initialcondition_z + \
                scalar_amplitude * \
                    np.random.normal(0, 1, [domain.nx, domain.ny, domain.nz])
        else:
            self.initialcondition_x = self.initialcondition_x + \
                scalar_amplitude * \
                    np.random.normal(0, 1, [domain.nx, domain.nz])
            self.initialcondition_y = self.initialcondition_y + \
                scalar_amplitude * \
                    np.random.normal(0, 1, [domain.nx, domain.nz])
            self.initialcondition_z = self.initialcondition_z + \
                scalar_amplitude * \
                    np.random.normal(0, 1, [domain.nx, domain.nz])
    
    # --------------------------------------------------------------------------
    def save_to_json(self, jsonfile = None):
        # First we can append the tensor_coefficients
        if self.is_seismic:
            section = "Seismic"
            self.append_coefficients(self.stiffness_coefficients, "Stiffness")
        else:
            section = "Electromagnetic"
            self.append_coefficients(
                self.permittivity_coefficients, "Permittivity"
            )
            self.append_coefficients(
                self.conductivity_coefficients, "Conductivity"
            )
        
        # Load the json file as a dictionary so we can edit the fields
        # if jsonfile:
        #     json_dictionary = readwrite_json(jsonfile)
        # else:
        json_dictionary = readwrite_json(self.project_file)
        
        # Add the time step 
        json_dictionary[section]['Source']['dt'] = self.dt
        json_dictionary[section]['Source']['x'] = self.x
        json_dictionary[section]['Source']['y'] = self.y
        json_dictionary[section]['Source']['z'] = self.z
        json_dictionary[section]['Source']['xind'] = self.xind
        json_dictionary[section]['Source']['yind'] = self.yind
        json_dictionary[section]['Source']['zind'] = self.zind
        json_dictionary[section]['Source']['source_frequency'] = self.f0
        json_dictionary[section]['Source']['x-z_rotation'] = self.theta
        json_dictionary[section]['Source']['x-y_rotation'] = self.phi
        json_dictionary[section]['Source']['amplitude'] = self.source_amplitude
        json_dictionary[section]['Source']['source_wavelet'] = self.source_wavelet
        
        if self.source_type:
            json_dictionary[section]['Source']['source_type'] = self.source_type
        
        if jsonfile:
            readwrite_json(jsonfile, json_dictionary)
        else:
            readwrite_json(self.project_file, json_dictionary)

    # --------------------------------------------------------------------------
    def append_coefficients(
            self, tensor: np.ndarray, tensor_name: str = None 
        ):
        """
        Appends coefficients to the JSON dictionary based on the provided tensor. 
        
        :param tensor: A numpy array containing the tensor coefficients to append.
        :type tensor: np.ndarray
        :param tensor_name: The parameter that the tensor describes. Available 
            options are Stiffness, Permittivity, Conductivity
        :type tensor_name: 
        
        :return: None
        """
        if tensor_name not in ['Stiffness', 'Permittivity', 'Conductivity']:
            raise ValueError(
                f"""You have input {tensor_name} as the input. The tensor_name must be 
                either 'Stiffness', 'Permittivity', 'Conductivity'. This is case 
                invariant, but spelling must be accurate."""
            )
        
        if tensor_name == "Stiffness":
            section = "Seismic"
        else: 
            section = "Electromagnetic"
        
        field = f"{tensor_name}_Coefficients"
        sub_field = tensor.columns.values
        material_id = np.arange(tensor.shape[0])

        json_dictionary = readwrite_json(self.project_file)     
        for ii in material_id:
            for component in sub_field:
                json_dictionary[section][field][ii][component] = tensor[component][ii]
        
        # Make sure the the indice values for the x, y, and z are added
        readwrite_json(self.project_file, json_dictionary)
    
    # --------------------------------------------------------------------------
    def run(self, num_threads = 1, jsonfile = None):
        '''
        
        '''
        # Write the initial conditions to their respective .dat files 
        if self.is_seismic:
            M = 'V'
        else:
            M = 'E'
        
        f = FortranFile(f'initialcondition{M}x.dat', 'w')
        f.write_record(self.initialcondition_x.T)
        f.close()
        
        f = FortranFile(f'initialcondition{M}y.dat', 'w')
        f.write_record(self.initialcondition_y.T)
        f.close()
        
        f = FortranFile(f'initialcondition{M}z.dat', 'w')
        f.write_record(self.initialcondition_z.T)
        f.close()
        
        # Update the JSON
        if not jsonfile:
            jsonfile = self.project_file 
        
        self.save_to_json(jsonfile = jsonfile)
        
        # Run it
        if num_threads > 1:
            env = os.environ.copy() 
            env['OMP_NUM_THREADS'] = str(num_threads)
            call([
                    'seidartfdtd', 
                    jsonfile, 
                    f'seismic={str(self.is_seismic).lower()}',
                    f'density_method={str(self.density_method).lower()}'
                ],
                env=env 
            )
        else:
            call([
                'seidartfdtd', 
                jsonfile, 
                f'seismic={str(self.is_seismic).lower()}',
                f'density_method={str(self.density_method).lower()}'
            ])
    
    # --------------------------------------------------------------------------
    def plotsource(self):
        """
        """
        if self.is_seismic:
            y_units = '(m/s)'
        else:
            y_units = ''
        time_vector = self.dt * np.arange(len(self.sourcefunction_x))
        fig, (axx, axy, axz) = plt.subplots(nrows = 3)
        axx.plot(time_vector, self.sourcefunction_x) 
        axx.set_ylabel('X ' + y_units)
        axx.set_xlabel('time (s)')
        axy.plot(time_vector, self.sourcefunction_y)
        axy.set_ylabel('Y ' + y_units)
        axy.set_xlabel('time (s)')
        axz.plot(time_vector, self.sourcefunction_z)
        axz.set_ylabel('Z ' + y_units)
        axz.set_xlabel('time (s)')
        plt.show()
    
    # --------------------------------------------------------------------------
    def tensor_test(self):
        '''
        Test for positive definiteness and ellipticity so that we can ensure that
        our domain is physically stable in the continuum mechanics sense
        '''
        if self.is_seismic:
            nmats = self.stiffness_coefficients.shape[0]
            pd_ell= np.zeros([nmats, 2], dtype = bool)
            
            for ii in range(nmats):
                row = self.stiffness_coefficients.loc[ii]
                C = np.array([
                    [row['c11'], row['c12'], row['c13'], row['c14'], row['c15'], row['c16']],
                    [row['c12'], row['c22'], row['c23'], row['c24'], row['c25'], row['c26']],
                    [row['c13'], row['c23'], row['c33'], row['c34'], row['c35'], row['c36']],
                    [row['c14'], row['c24'], row['c34'], row['c44'], row['c45'], row['c46']],
                    [row['c15'], row['c25'], row['c35'], row['c45'], row['c55'], row['c56']],
                    [row['c16'], row['c26'], row['c36'], row['c46'], row['c56'], row['c66']]
                ])
                eigenvalues = np.linalg.eigvalsh(C) 
                min_eigenvalue = np.min(eigenvalues) 
                pd_ell[ii,0] = np.all(eigenvalues > 0) 
                pd_ell[ii,1] = min_eigenvalue > 0  
            
            self.positive_definite_and_elliptical = pd_ell
        else: 
            nmats = self.permittivity_coefficients.shape[0] 
            pd_pd = np.zeros([nmats,2], dtype = bool) 
            for ii in range(nmats):
                rowp = self.permittivity_coefficients.loc[ii] 
                rowc = self.conductivity_coefficients.loc[ii]
                P = np.array([
                    [row['e11'], row['e12'], row['e13']],
                    [row['e12'], row['e22'], row['e23']],
                    [row['e13'], row['e23'], row['e33']]
                ])
                C = np.array([
                    [row['s11'], row['s12'], row['s13']],
                    [row['s12'], row['s22'], row['s23']],
                    [row['s13'], row['s23'], row['s33']]
                ])
                eigenvalues = np.linalg.eigvalsh(P) 
                pd_pd[ii,0] = np.all(eigenvalues > 0) 
                eigenvalues = np.linalg.eigvalsh(C) 
                pd_pd[ii,1] = np.all(eigenvalues > 0) 
            
            self.positive_definite = pd_pd
    # --------------------------------------------------------------------------
    
    # --------------------------------------------------------------------------
    def get_christoffel_matrix(self, material_indice, direction):
        if self.is_seismic:
            row = self.stiffness_coefficients.loc[material_indice]
        else:
            row_epsilon = self.permittivity_coefficients.loc[material_indice]
            row_sigma = self.conductivity_coefficients.loc[material_indice]
            row = pd.concat([row_epsilon, row_sigma])
        
        return mf.get_christoffel_matrix(row, self.f0, self.is_seismic, direction)
    
    # --------------------------------------------------------------------------
    def compute_christoffel_directions(self, material_index, n_theta=30, n_phi=60):
        """
        Compute eigenvectors from the Christoffel matrix over a grid of directions.

        Parameters
        ----------
        material_index : int
            Index of the material from stiffness_coefficients.
        n_theta : int
            Number of points in inclination (from 0 to pi/2, lower hemisphere).
        n_phi : int
            Number of azimuthal angles (from 0 to 2*pi).

        Stores
        ------
        self.christoffel_solutions : list of dict
            Each dict contains direction vector, eigenvectors, eigenvalues, velocities.
        """
        self.christoffel_solutions = []

        for i in range(n_theta):
            theta = (i + 0.5) * (np.pi / 2) / n_theta  # 0 to pi/2
            for j in range(n_phi):
                phi = j * 2 * np.pi / n_phi  # 0 to 2pi

                # Convert to Cartesian
                n = np.array([
                    np.sin(theta) * np.cos(phi),
                    np.sin(theta) * np.sin(phi),
                    np.cos(theta)
                ])

                Gamma, eigvals, eigvecs, velocities = self.get_christoffel_matrix(material_index, n)
                self.christoffel_solutions.append({
                    'direction': n,
                    'eigenvectors': eigvecs,
                    'eigenvalues': eigvals,
                    'velocities': velocities
                })

    # --------------------------------------------------------------------------
    def compute_max_velocities(self, dim = 2, directions=None):
        """
        Loop over all materials in the model and compute the maximum phase velocity
        for each material using the Christoffel matrix.

        :param directions: List of direction keys to use. Defaults to ['x', 'y', 'z', 'xy', 'xz', 'yz'].
        :return: Dictionary {material_id: max_velocity}
        """
        
        if directions is None:
            if dim == 2:
                directions = ['x', 'z', 'xz']
            else:
                directions = ['x', 'y', 'z', 'xy', 'xz', 'yz']
        
        max_velocity_per_material = {}
        
        # Loop over all material IDs present in the stiffness/permittivity tables
        if self.is_seismic:
            material_ids = self.stiffness_coefficients.index 
        else: 
            material_ids = self.permittivity_coefficients.index
        
        for mat_id in material_ids:
            local_max = 0.0
            for direction in directions:
                _, _, _, velocities = self.get_christoffel_matrix(mat_id, direction)
                local_max = max(local_max, np.max(velocities))
            max_velocity_per_material[mat_id] = local_max
            
        if not self.is_seismic:
            for key, value in max_velocity_per_material.items():
                max_velocity_per_material[key] = 1/value
        
        if direction == 'x':
            self.vmax_x = max_velocity_per_material 
        elif direction == 'z': 
            self.vmax_z = max_velocity_per_material 
        else:
            self.max_velocity_per_material = max_velocity_per_material
    
    # --------------------------------------------------------------------------
    def plot_lower_hemisphere_polarizations(self):
        """
        Plot lower hemisphere projections of P, SV, and SH polarization directions
        using previously computed Christoffel solutions.

        Returns
        -------
        fig, axs : Matplotlib figure and axes objects.
        """
        def project(v):
            x, y, z = v
            # Flip to lower hemisphere if needed
            if z < 0:
                v = -v
            r = np.sqrt(1 - v[2]**2)
            return v[0] * r, v[1] * r

        fig, axs = plt.subplots(1, 3, figsize=(12, 4), subplot_kw=dict(aspect='equal'))

        titles = ['P-wave', 'SV-wave', 'SH-wave']
        colors = ['r', 'b', 'g']

        for i, ax in enumerate(axs):
            ax.set_title(titles[i])
            ax.set_xlim(-1.05, 1.05)
            ax.set_ylim(-1.05, 1.05)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.axhline(0, color='k', linewidth=0.5)
            ax.axvline(0, color='k', linewidth=0.5)
            ax.add_patch(plt.Circle((0, 0), 1, color='gray', fill=False, linewidth=0.5))

            for sol in self.christoffel_solutions:
                vec = sol['eigenvectors'][:, i]  # i = 0: P, 1: SV, 2: SH
                x, y = project(vec)
                ax.plot(x, y, marker='o', color=colors[i], markersize=2, alpha=0.7)

        fig.tight_layout()
        return fig, axs

    # ------------------------------------------------------------------------------
    def parameter_profile_1d(
            self, domain, material, 
            indice: int, parameter: str = 'velocity', direction = 'z'
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
            
            output_values = np.empty([n, 3])
            if self.is_seismic:
                for ind in range(n):
                    __, __, __, output_values[ind,:] = self.get_christoffel_matrix(profile[ind], direction)
            else:                      
                for ind in range(n):
                    T = self.permittivity_coefficients[profile[ind],:]
                    output_values[ind,:] = mf.tensor2velocities(
                        T, seismic = False
                    )
        elif parameter == 'conductivity':
            for ind in range(n):
                T = em.conductivity_coefficients.loc[profile[ind],:] 
                output_values[ind,:] = T[np.array([0,3,5])]
        else:
            attribute_name = attribute_map.get(parameter)
            vals = getattr(material, attribute_name, None)
            output_values = np.empty([n])
            for ind in range(n):
                output_values[ind] = vals[profile[ind]]

        z = np.arange(0,n) * domain.dz
        self.parameter_z = z.copy()
        self.parameter_output_values = output_values.copy()
        return  z, output_values
    
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------    
    def plot_unit_sphere_wave_speeds(self, material_indice, num_samples=100):
        """
        Plot wave speeds on a unit sphere for a given material.

        Parameters:
            self: The class containing material properties and methods.
            material_indice (int): Index of the material in the DataFrame.
            num_samples (int): Number of unit sphere sample points.

        Returns:
            None (displays figure).
        """
        
        # Generate spherical coordinates (θ = azimuthal, φ = polar)
        theta = np.linspace(0, 2 * np.pi, num_samples)  # 0 to 360 degrees
        phi = np.linspace(0, np.pi, num_samples // 2)   # 0 to 180 degrees
        
        theta, phi = np.meshgrid(theta, phi)  # Create a meshgrid
        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)
        
        # Initialize velocity storage
        Vp_values = np.zeros_like(x)
        Vs1_values = np.zeros_like(x)
        Vs2_values = np.zeros_like(x)
        
        # Loop over each (θ, φ) direction
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                # Propagation direction (normalized)
                n = np.array([x[i, j], y[i, j], z[i, j]])
                n /= np.linalg.norm(n)  # Ensure unit vector
                
                # Compute Christoffel matrix for this direction
                Gamma = self.get_christoffel_matrix(material_indice, n)
                eigenvalues, _ = np.linalg.eigh(Gamma)  # Solve eigenvalue problem
                
                # Compute wave speeds (v = sqrt(eigenvalue))
                velocities = np.sqrt(np.abs(eigenvalues))
                velocities.sort()  # Sort eigenvalues into (S1, S2, P)
                
                Vs1_values[i, j] = velocities[0]  # Slowest S-wave
                Vs2_values[i, j] = velocities[1]  # Second S-wave
                Vp_values[i, j] = velocities[2]   # Fastest P-wave
        
        # Find max velocity for annotation
        max_vp = np.max(Vp_values)
        max_vs1 = np.max(Vs1_values)
        max_vs2 = np.max(Vs2_values)
        
        # Create figure with subplots
        fig = plt.figure(figsize=(12, 8))
        
        # Plot P-wave velocity
        ax1 = fig.add_subplot(131, projection='3d')
        ax1.plot_surface(x, y, z, facecolors=plt.cm.viridis(Vp_values / max_vp), rstride=1, cstride=1, alpha=0.8)
        ax1.set_title(f'P-Wave Velocity\nMax: {max_vp:.2f} m/s')
        
        # Plot S-wave 1 velocity
        ax2 = fig.add_subplot(132, projection='3d')
        ax2.plot_surface(x, y, z, facecolors=plt.cm.viridis(Vs1_values / max_vs1), rstride=1, cstride=1, alpha=0.8)
        ax2.set_title(f'S-Wave 1 Velocity\nMax: {max_vs1:.2f} m/s')
        
        # Plot S-wave 2 velocity
        ax3 = fig.add_subplot(133, projection='3d')
        ax3.plot_surface(x, y, z, facecolors=plt.cm.viridis(Vs2_values / max_vs2), rstride=1, cstride=1, alpha=0.8)
        ax3.set_title(f'S-Wave 2 Velocity\nMax: {max_vs2:.2f} m/s')
        
        # Adjust layout
        plt.tight_layout()
        plt.show()
        
        return fig, ax1, ax2, ax3

# ------------------------------------------------------------------------------
class AnimatedGif:
    """
    A class to create an animated GIF from a series of images.

    :param size: The size of the animated GIF in pixels (width, height), 
        defaults to (640, 480)
    :type size: Tuple[int, int], optional

    Attributes:
        fig (matplotlib.figure.Figure): The figure object for the animation.
        images (List): A list of image frames to be included in the GIF.
        background (List): Background image data.
        source_location (List[int, int]): The location of the source point in 
        the plot.
        nx (int): Width of the plot.
        nz (int): Height of the plot.
        output_format (int): Format of the output file. 0 for GIF, 1 for other 
        formats.
    """
    def __init__(self, size=(640,480) ):
        """
        Initialize the AnimatedGif class with a figure size.
        """
        self.fig = plt.figure()
        self.fig.set_size_inches(size[0]/100, size[1]/100)
        ax = self.fig.add_axes([0, 0, 1, 1], frameon=False, aspect=1)
        ax.set_xticks([])
        ax.set_yticks([])
        self.images = []
        self.background = []
        self.source_location = []
        self.nx = size[0]
        self.nz = size[1]
        self.output_format = 0

    def add(self, image, label='', extent=None ):
        """
        Adds an image to the animation sequence.

        :param image: The image data to add.
        :type image: np.ndarray
        :param label: The label to add to the image, defaults to ''.
        :type label: str, optional
        :param extent: The extent of the image, defaults to None.
        :type extent: Tuple[int, int, int, int], optional
        """
        bound = np.max([abs(np.min(image)),abs(np.max(image))])
        plt_im = plt.imshow(
            image,cmap='seismic',
            animated=True,
            extent=(0, (self.nx), (self.nz), 0),
            vmin=-bound,vmax=bound
        )
        plt_bg = plt.imshow(
            self.background,
            alpha = 0.3,
            extent=extent,
            animated = True
        )
        plt.scatter(
            self.source_location[0],
            self.source_location[1],
            marker = '*',
            s = 30,
            linewidths = 1,
            edgecolor = (0.2, 0.2, 0.2, 1 )
        )
        plt_txt = plt.text(
            extent[0] + 20,
            extent[2] + 20,
            label,
            color='red'
        ) # Lower left corner
        self.images.append([plt_im, plt_bg, plt_txt])

    def save(self, filename, frame_rate = 50):
        """
        Save the animated GIF.

        :param filename: The name of the output file.
        :type filename: str
        :param frame_rate: The frame rate of the animation, defaults to 50.
        :type frame_rate: int, optional
        """
        animation = anim.ArtistAnimation(self.fig,
                                        self.images,
                                        interval = frame_rate,
                                        blit = True
                                        )
        if self.output_format == 1:
            animation.save(filename, dpi = 300)
        else:
            animation.save(filename, dpi = 300, writer = 'imagemagick')


# ------------------------------------------------------------------------------
class VolumeBuilder:
    """
    Discretize an OBJ scene onto a regular FD grid with region labels.
    
    Parameters
    ----------
    obj_path : str
        Path to the OBJ file.
    priority : dict
        Dict mapping group names to priority (higher = later overwrite).
        Example: {"ice": 0, "air": 1, "base": 2, "heterogeneity": 3}
    x_min, x_max, dx : float
    y_min, y_max, dy : float
    z_min, z_max, dz : float
        Grid extents and spacings in OBJ coordinates.
    """
    
    def __init__(
            self,
            obj_path,
            priority,
            x_min, x_max, dx,
            y_min, y_max, dy,
            z_min, z_max, dz,
            background_label=0,
            special_geometries = None
        ):
        self.obj_path = obj_path
        base, _ = obj_path.rsplit(".", 1)
        self.mtl_path = base + ".mtl"
        self.priority = priority
        self.background_label = background_label
        
        self._parse_mtl_kd() 
        
        # Load scene
        scene = trimesh.load(self.obj_path, process=False)
        if isinstance(scene, trimesh.Scene):
            self.geoms = scene.geometry  # dict: name -> Trimesh [web:70][web:75]
        else:
            self.geoms = {"mesh": scene}
        
        # Build grid
        self.xs = np.arange(x_min, x_max + 0.5 * dx, dx)
        self.ys = np.arange(y_min, y_max + 0.5 * dy, dy)
        self.zs = np.arange(z_min, z_max + 0.5 * dz, dz)
        self.dx, self.dy, self.dz = dx, dy, dz 
        self.nx, self.ny, self.nz = len(self.xs), len(self.ys), len(self.zs)
        
        X, Y, Z = np.meshgrid(self.xs, self.ys, self.zs, indexing="ij")
        self.points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
        self.labels_1d = np.full(self.points.shape[0],
                                 self.background_label,
                                 dtype=np.int32)
        self.label_grid = None
        
        if special_geometries is None:
            self.special_geometries = ()
        elif isinstance(special_geometries, str):
            self.special_geometries = (special_geometries,)   # single-element tuple
        elif isinstance(special_geometries, np.ndarray):
            self.special_geometries = tuple( str(x) for x in special_geometries.ravel() )
        else:
            self.special_geometries = tuple(special_geometries)
        
    
    # -------------------------
    # internal helpers
    # -------------------------
    def _parse_mtl_kd(self):
        """
        Parse a .mtl file and return {material_name: (r, g, b)} from Kd lines.
        r,g,b are floats in [0,1].
        """
        kd_map = {}
        current = None
        
        with open(self.mtl_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if parts[0] == "newmtl" and len(parts) > 1:
                    current = parts[1]
                elif parts[0] == "Kd" and current is not None and len(parts) >= 4:
                    r, g, b = map(float, parts[1:4])
                    kd_map[current] = (r, g, b)
        
        # Build label -> RGB (float triple) by matching tag substrings to material names
        label_to_rgb = {}
        for tag, label in self.priority.items():
            rgb = None
            for mtl_name, kd in kd_map.items():
                if tag in mtl_name.lower():
                    rgb = kd
                    break
            if rgb is None:
                # Default color if no material matches (e.g. gray)
                rgb = (0.5, 0.5, 0.5)
            
            rgb_f = np.array(rgb) 
            rgb_arr = (rgb_f * 255 ).astype(int)
            label_to_rgb[label] = rgb_arr    
            
        self.label_to_rgb = label_to_rgb
        self.kd_map = kd_map
    
    def _match_tag(self, name: str):
        """
        Find which priority tag applies to this geometry name.
        
        Returns
        -------
        tag : str or None
            The matched tag (key in self.priority) or None if no match.
        """
        lower = name.lower()
        for tag in self.priority.keys():
            if tag in lower:
                return tag
        return None
    
    def _get_priority_for_name(self, name: str) -> int:
        """
        Get the numeric priority for a geometry name.
        """
        tag = self._match_tag(name)
        if tag is None:
            return -1  # lowest
        return self.priority[tag]
    
    def _get_label_for_name(self, name: str) -> int:
        """
        Label ID = priority value associated with matched tag.
        """
        tag = self._match_tag(name)
        if tag is None:
            return self.background_label
        return self.priority[tag]
    
    # -------------------------
    # main API
    # -------------------------
    
    def label_domain(self):
        """
        Fill self.labels_1d and self.label_grid using priority-based overwrite.
        
        Returns
        -------
        label_grid : np.ndarray, shape (nx, ny, nz)
        """
        names_sorted = sorted(self.geoms.keys(), key=self._get_priority_for_name)
        
        for name in names_sorted:
            mesh = self.geoms[name]
            tag = self._match_tag(name)
            prio = self._get_priority_for_name(name)
            label = self._get_label_for_name(name)
            
            if tag is None:
                print(f"Skipping '{name}' (no priority tag match).")
                continue
            
            print(f"Processing '{name}' tag='{tag}' label={label} priority={prio}")
            
            bounds_min, bounds_max = mesh.bounds
            
            coarse = (
                (self.points[:, 0] >= bounds_min[0]) & (self.points[:, 0] <= bounds_max[0]) &
                (self.points[:, 1] >= bounds_min[1]) & (self.points[:, 1] <= bounds_max[1]) &
                (self.points[:, 2] >= bounds_min[2]) & (self.points[:, 2] <= bounds_max[2])
            )
            
            candidate_idx = np.where(coarse)[0]
            if candidate_idx.size == 0:
                print("  No points in bounding box; skipping.")
                continue

            pts_candidate = self.points[candidate_idx]

            # name-based special-geometry detection
            lower = name.lower()
            is_special = False
            for geometry in self.special_geometries:
                if geometry and geometry in lower:
                    is_special = True
                    break

            if is_special:
                inside_local = mesh.contains(pts_candidate)
                self.labels_1d[candidate_idx[inside_local]] = label
                print(
                    f"  {lower}: coarse {candidate_idx.size}, "
                    f"inside {inside_local.sum()}"
                )
            else:
                # For axis-aligned cubes (ice, air, base), AABB is enough
                self.labels_1d[candidate_idx] = label
                print(f"  Bulk region: labeled {candidate_idx.size} cells (AABB)")

        label_grid = self.labels_1d.reshape((self.nx, self.ny, self.nz))
        # This all needs to be transposed to match the correct coordinates in Fortran
        # label_grid = np.flip(np.transpose(label_grid, (1, 2, 0)), axis = 0)
        
        self.label_grid = label_grid
        
        f = FortranFile('geometry.dat', 'w')
        f.write_record(np.asfortranarray(self.label_grid).astype(np.int32))
        f.close()
        
        return self.label_grid
    
    # -------------------------
    # visualization helper
    # -------------------------
    
    def plot_slice(self, index, plane="xy"):
        """
        Show a 2D slice of a 3D array along one of the principal planes.
        plane: 'xy', 'xz', or 'yz'
        """
        Nx, Ny, Nz = self.label_grid.shape
        if plane == "xy":
            assert 0 <= index < Nz
            img = self.label_grid[index, :, :]
            xlabel, ylabel = "ix", "iy"
            title = f"Slice plane=xy, iz={index}"
        elif plane == "xz":
            assert 0 <= index < Ny
            img = self.label_grid[:, index, :]
            xlabel, ylabel = "ix", "iz"
            title = f"Slice plane=xz, iy={index}"
        elif plane == "yz":
            assert 0 <= index < Nx
            img = self.label_grid[:, :, index]
            xlabel, ylabel = "iy", "iz"
            title = f"Slice plane=yz, ix={index}"
        else:
            raise ValueError("plane must be one of 'xy', 'xz', 'yz'")
        
        plt.figure(figsize=(6, 5))
        plt.imshow(img)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()
        plt.show()
    
    # ----------------------------------
    # Create the project file
    # ----------------------------------
    def prjbuild(self, output_json):
        """
        Generate the project file like the 2/2.5D uses
        """
        self.material_ids = np.unique(self.label_grid)
        self.nmats = len(self.material_ids)
        template = generate_template(
            self.nmats, 
            Domain_nx = self.nx, Domain_nz = self.nz, Domain_ny = self.ny,
            Domain_dx = self.dx, Domain_dz = self.dz, Domain_dy = self.dy,
        )
        updates = {}
        
        for i, mid in enumerate(self.material_ids):
            rgb = self.label_to_rgb.get(mid, np.array([128, 128, 128]))
            updates[("Materials", mid, "rgb")]  = '/'.join(rgb.astype(str))
        
        updates[("Domain", "image_file")] = self.obj_path 
        updated_template = update_json(template, updates) 
        
        with open(output_json, 'w') as f: 
            json.dump(updated_template, f, indent = 4) 
