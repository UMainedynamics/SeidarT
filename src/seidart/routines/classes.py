import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg 
import matplotlib.animation as anim
import os 
import os.path 
from scipy.io import FortranFile
from subprocess import call

from typing import Optional
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
        self.kappa_max = 5 
        self.Rcoef = None # Seismic only parameter
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
        self.source_type = None
        self.src = None
        self.attenuation_fadjust = None
        self.exit_status = 0
        self.is_seismic = None
        self.kappa_max = np.zeros([3])
        self.sig_max = np.zeros([3])
        self.alpha_max = np.zeros([3]) 
        self.kappa_max_half = np.zeros([3])
        self.sig_max_half = np.zeros([3])
        self.alpha_max_half = np.zeros([3])
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
        self.sourcefunction = sourcefunction.pointsource # This will change
        self.get_seismic = mf.get_seismic 
        self.get_perm = mf.get_perm
        self.CFL = 1/np.sqrt(2) # Default minimum. If the dimension is 2.5, then this will change automatically
        
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
                velocities[ind,:] = mf.tensor2velocities(
                    self.stiffness_coefficients.loc[ind][:-1].to_numpy(), 
                    self.stiffness_coefficients['rho'][ind], 
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
        
        if self.exit_status == 0 and recompute_tensors:
            # The coefficients aren't provided but the materials are so we can 
            # compute them
            # Assign the materials to their respective corners
            max_vels = np.zeros([domain.nmats])
            
            if self.is_seismic:
                print('Computing the stiffness coefficients.')
                self.get_seismic(self, material)
                
                # We need to compute dt from the Courant number. We can use the 
                # maximum tensor value, and the maximum density even if they don't 
                # correspond to the same material.
                
                for ind in range(domain.nmats):
                    max_vels[ind] = np.sqrt( (self.stiffness_coefficients.loc[ind]/ \
                        self.stiffness_coefficients.rho.loc[ind]).max() )
                    
                max_rho = self.stiffness_coefficients["rho"].max()
                if domain.dim == 2:
                    self.dt = self.CFL / np.sqrt(1/domain.dx**2 + 1/domain.dz**2) / max_vels.max() / 2
                else:
                    self.dt = self.CFL / np.sqrt(1/domain.dx**2 + 1/domain.dy**2 + 1/domain.dz**2) / max_vels.max() / 2
            else:
                print('Computing the permittivity and conductivity coefficients.')
                
                self.get_perm(self, material)
                
                for ind in range(domain.nmats):
                    max_vels[ind] = clight / \
                        np.sqrt(self.permittivity_coefficients[['e11', 'e22', 'e33']].loc[ind].min().real)
                        
                self.dt = CFL * np.min([domain.dx, domain.dz]) / max_vels.max()
            
            # The time step needs to satisfy the Courant number and also have a nyquist
            # that will resolve the source frequency
            src_nyquist = 1/(2*self.f0)
            if src_nyquist < self.dt:
                print(
                    '''Nyquist is not small enough for the source frequency. Change
                    the source frequency or decrease the spatial step size'''
                )
        
        # We need to set the 
        print("Creating the source function")
        (
            __, self.sourcefunction_x, self.sourcefunction_y, 
            self.sourcefunction_z, __ 
        ) = self.sourcefunction(self)
        
        direction = ['x', 'y', 'z']
        # Compute CPML
        print('Computing CPML boundary values and writing outputs to Fortran files.')
        for d in direction:
            cpmlcompute(self, domain, d)
            cpmlcompute(self, domain, d, half = True)
        
        # Write out the tensor components to file
        if write_tensor:
            print('Writing tensor components to individual .dat files.')
            # Compute the density gradient at boundaries with air. If the model is a
            # a single material, this will return an error. Booooo, errors!
            if domain.nmats > 1:
                self.domain_density = airsurf(material, domain, self.air_gradient_integer)
            else:
                self.domain_density = np.ones([domain.nx, domain.nz])
            
            # Write out the arrays of tensor coefficients
            if self.is_seismic:
                self.tensor2dat(self.stiffness_coefficients, domain)
                # gamma/attenuation coefficients need to be scaled by the 
                # dominant frequency of the source
                self.tensor2dat(self.attenuation_coefficients * self.f0, domain)
            else:
                self.tensor2dat(self.permittivity_coefficients, domain) 
                self.tensor2dat(self.conductivity_coefficients, domain)
        else:
            print('Tensor components have already been written to .dat files.')
            return
        
    def tensor2dat(self, tensor, domain):
        """

        """
        columns = tensor.columns
        coef_array = np.zeros([domain.nx, domain.nz])
        extended_array = np.zeros([domain.nx+2*domain.cpml, domain.nz+2*domain.cpml])
        for col in columns:
            for ii in range(domain.nx):
                for jj in range(domain.nz):
                    coef_array[ii,jj] = tensor[col][domain.geometry[ii,jj]]
                    fn = col + '.dat' 
            
            if col == 'rho':
                coef_array = coef_array*self.domain_density
            
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
        
        # # For em models, we need to put in the initial conditions for the magnetic field
        # if not self.is_seismic:
        #     for direction in ['x', 'y', 'z']:
        #         f = FortranFile(f'initialconditionH{direction}.dat', 'w')
        #         f.write_record(self.initialcondition_)
        #         f.close() 
        
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
                f'seismic={str(self.is_seismic).lower()}'],
                env=env 
            )
        else:
            call([
                'seidartfdtd', 
                jsonfile, 
                f'seismic={str(self.is_seismic).lower()}'
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
    def voigt_to_full_tensor(self, C_voigt):
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
    
    # --------------------------------------------------------------------------
    def get_christoffel_matrix(self, material_indice, direction):
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
        
        # Convert 6x6 Voigt tensor to 3x3x3x3 full tensor
        voigt_to_tensor = {
            0: (0, 0), 1: (1, 1), 2: (2, 2),  # Normal components
            3: (1, 2), 4: (0, 2), 5: (0, 1)   # Shear components
        }
        # Define a factor for each Voigt index:
        voigt_factor = {0: 1.0, 1: 1.0, 2: 1.0, 3: 1/np.sqrt(2), 4: 1/np.sqrt(2), 5: 1/np.sqrt(2)}
        
        if self.is_seismic:
            C_full = np.zeros((3, 3, 3, 3))
            # Extract stiffness coefficients from the specified row
            row = self.stiffness_coefficients.iloc[material_indice]
            C = np.array([
                [row['c11'], row['c12'], row['c13'], row['c14'], row['c15'], row['c16']],
                [row['c12'], row['c22'], row['c23'], row['c24'], row['c25'], row['c26']],
                [row['c13'], row['c23'], row['c33'], row['c34'], row['c35'], row['c36']],
                [row['c14'], row['c24'], row['c34'], row['c44'], row['c45'], row['c46']],
                [row['c15'], row['c25'], row['c35'], row['c45'], row['c55'], row['c56']],
                [row['c16'], row['c26'], row['c36'], row['c46'], row['c56'], row['c66']]
            ])
            
            rho = row['rho']  # Extract density
            C_full = self.voigt_to_full_tensor(C)
            Gamma = np.zeros((3, 3))
            for i in range(3):
                for j in range(3):
                    Gamma[i, j] = np.sum(C_full[i, :, j, :] * np.outer(n, n))/rho
        else:
             # Extract permittivity values from the specified row
            row_epsilon = self.permittivity_coefficients.iloc[material_indice]
            epsilon = np.array([
                [row_epsilon['e11'], row_epsilon['e12'], row_epsilon['e13']],
                [row_epsilon['e12'], row_epsilon['e22'], row_epsilon['e23']],
                [row_epsilon['e13'], row_epsilon['e23'], row_epsilon['e33']]
            ])

            # Extract conductivity values from the specified row
            row_sigma = self.conductivity_coefficients.iloc[material_indice]
            sigma = np.array([
                [row_sigma['s11'], row_sigma['s12'], row_sigma['s13']],
                [row_sigma['s12'], row_sigma['s22'], row_sigma['s23']],
                [row_sigma['s13'], row_sigma['s23'], row_sigma['s33']]
            ])

            # Compute the effective permittivity tensor
            omega = 2 * np.pi * self.f0  # Angular frequency
            epsilon_eff = epsilon + 1j * sigma / omega  # Complex permittivity

            # Compute inverse of effective permittivity tensor
            epsilon_eff_inv = np.linalg.inv(epsilon_eff)
            
            # Construct the electromagnetic Christoffel matrix
            Gamma = np.zeros((3, 3), dtype=complex)
            for i in range(3):
                for j in range(3):
                    Gamma[i, j] = np.sum(epsilon_eff_inv[i, :] * n[j] * n)
        
        # Compute eigenvalues (which correspond to v^2 for wave propagation)
        eigenvalues, eigenvectors = np.linalg.eigh(Gamma)
        
        # Convert to phase velocities (v = sqrt(v^2)) ensuring non-negative values
        velocities = np.sqrt(np.abs(eigenvalues))  # Take absolute value before sqrt
        
        # Sort eigenvalues to match P-wave (fastest), S1, and S2 waves
        velocities.sort()
        
        return Gamma, eigenvalues, eigenvectors, velocities

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
