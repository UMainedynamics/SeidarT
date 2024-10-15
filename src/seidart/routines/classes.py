import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg 
import matplotlib.animation as anim 
import os.path 

from typing import Optional
import seidart.routines.materials as mf
from seidart.routines import sourcefunction
from seidart.routines.definitions import *

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
    material_list : numpy.ndarray
        An array to store the list of materials.
    material_flag : bool
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
        self.material_list = np.array([]) # initialize
        self.material_flag = False # Whether the materials were read in

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
        m = len(self.material_list)
        self.material = np.zeros([m], dtype = 'U12')
        self.temp = np.zeros([m], dtype = float)
        self.rho = np.zeros([m], dtype = float)
        self.porosity = np.zeros([m], dtype = float)
        self.lwc = np.zeros([m], dtype = float)
        self.is_anisotropic = np.zeros([m], dtype = bool)
        self.angfile = np.zeros([m], dtype = object)
        
        for ind in range(m):
            (
                __, self.material[ind], self.temp[ind], self.rho[ind], 
                self.porosity[ind] , self.lwc[ind] , self.is_anisotropic[ind], 
                self.angfile[ind] 
            ) = self.material_list[ind]

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
                if self.material_list[ind,6] == 'True' and \
                    not self.material_list[ind,7] or \
                        self.material_list[ind,7] == 'n/a':
                    file_check = file_check + 1
        else:
            print('Material inputs aren"t satisfied.')

        if check == 0:
            if file_check == 0:
                self.material_flag = True
            else:
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
        self.dt = None
        self.time_steps = None
        self.x = None
        self.y = None
        self.z = None
        self.f0 = None
        self.theta = None
        self.phi = None
        self.src = None
        self.fref = None
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
        self.build()
    
    def build(self) -> None:
        """
        Initializes the simulation model attributes with default values.
        """
        self.sourcefunction = sourcefunction.pointsource # This will change
        self.get_seismic = mf.get_seismic 
        self.get_perm = mf.get_perm
            
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
    def status_check(
            self, 
            material: Material, 
            domain: Domain, 
            projectfile,
            append_to_json: bool = True
        ) -> None:
        """
        Checks the status of the modeling classes and appends coefficients to the
        project file if necessary.

        :param material: The material class instance.
        :type material: Material
        :param domain: The domain class instance.
        :type domain: Domain
        :param projectfile: Path to the project file.
        :type projectfile: str
        :param append_to_json: Flag to append coefficients to the project file. 
                                Defaults to True.
        :type append_to_json: bool
        """
        global CFL
        
        # Set the CFL to its max value if it exceeds the 
        if domain.dim == 2.5:
            CFL = np.min([CFL, 1/np.sqrt(3)])
        else:
            CFL = np.min([CFL, 1/np.sqrt(2)])
        
        # Make sure all of the initial conditions are input 
        if self.initialcondition_x  is None:
            if domain.dim == 2.5:
                self.initialcondition_x = np.zeros(domain.nx, domain.ny, domain.nz)
            else:
                self.initialcondition_x = np.zeros([domain.nx, domain.nz])
        if self.initialcondition_z  is None:
            if domain.dim == 2.5:
                self.initialcondition_z = np.zeros([domain.nx, domain.ny, domain.nz])
            else:
                self.initialcondition_z = np.zeros([domain.nx, domain.nz])
        if self.initialcondition_y  is None:
            if domain.dim == 2.5:
                self.initialcondition_y = np.zeros([domain.nx, domain.ny, domain.nz])
            else:
                self.initialcondition_y = np.zeros([domain.nx, domain.nz])
                
        
        if self.exit_status == 0 and \
            material.material_flag and \
                append_to_json:
            # The coefficients aren't provided but the materials are so we can 
            # compute them
            # Assign the materials to their respective corners
            if self.is_seismic:
                print('Computing the stiffness coefficients.')
                self.get_seismic(self, material)
                
                # We need to compute dt from the Courant number. We can use the 
                # maximum tensor value, and the maximum density even if they don't 
                # correspond to the same material.
                ind = np.where(self.stiffness_coefficients.max() == self.stiffness_coefficients)
                max_rho = self.stiffness_coefficients[ ind[0][0], -1]
                self.dt = CFL * np.min([domain.dx, domain.dz]) / np.sqrt(3.0 * self.stiffness_coefficients.max()/max_rho )
                append_coefficients(
                    projectfile, 
                    self.stiffness_coefficients, tensor_name = 'Stiffness', 
                    dt = self.dt
                )
            else:
                print('Computing the permittivity and conductivity coefficients.')
                
                self.get_perm(self, material)

                self.dt = CFL * np.min([domain.dx, domain.dz]) / \
                    ( clight/ \
                        np.sqrt(np.min(
                            [
                                self.permittivity_coefficients[:,1].real.astype(float).min(), 
                                self.permittivity_coefficients[:,4].real.astype(float).min()
                            ]
                        )) 
                    )

                append_coefficients(
                    projectfile, 
                    self.permittivity_coefficients, 
                    tensor_name = 'Permittivity', 
                    dt = self.dt
                )
                append_coefficients(
                    projectfile, 
                    self.conductivity_coefficients, 
                    tensor_name = 'Conductivity', 
                    dt = self.dt
                )
                
            
            # The time step needs to satisfy the Courant number and also have a nyquist
            # that will resolve the source frequency
            src_nyquist = 1/(2*self.f0)
            if src_nyquist < self.dt:
                print(
                    '''Nyquist is not small enough for the source frequency. Change
                    the source frequency or decrease the spatial step size'''
                )

            print("Finished. Appending to project file.\n")
    
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
    def run(self, project_file, material, domain):
        direction = ['x', 'y', 'z']
        # Compute CPML
        print('Computing CPML boundary values and writing outputs to Fortran files.')
        for d in direction:
            cpmlcompute(self, domain, d)
            cpmlcompute(self, domain, d, half = True)
        
        # Compute the density gradient at boundaries with air. 
        self.domain_density = airsurf(material, domain, self.air_gradient_integer)
        
        # Write the initial conditions to their respective .dat files 
        if self.is_seismic:
            M = 'V'
        else:
            M = 'E'
        
        f = FortranFile(f'initialcondition{M}y.dat', 'w')
        f.write_record(self.initialcondition_y)
        f.close()
        
        f = FortranFile(f'initialcondition{M}y.dat', 'w')
        f.write_record(self.initialcondition_y)
        f.close()
        
        f = FortranFile(f'initialcondition{M}z.dat', 'w')
        f.write_record(self.initialcondition_z)
        f.close()
        
        call(['seidartfdtd', project_file, f'seismic={self.is_seismic}'])

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
