import numpy as np
import pandas as pd
import seidart.routines.materials as mf
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.animation as anim
import os.path
from typing import Optional
from subprocess import call
from scipy.io import FortranFile
import glob2
import copy

__all__ = [
    'read_dat',
    'image2int',
    'complex2str',
    'str2complex',
    'append_coefficients',
    'loadproject',
    'prepme',
    'coefs2prj',
    'airsurf',
    'rcxgen',
    'coherdt',
    'coherstf',
    'stfvec2mat',
    'movingsrc',
    'indvar',
    'agc',
    'correct_geometric_spreading',
    'exponential_gain',
    'Domain',
    'Material',
    'Model',
    'AnimatedGif',
    'kband_check',
    'parameter_profile_1d',
]

# =============================================================================
# =========================== Define Class Variables ==========================
# =============================================================================

# We need to define some class variables
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
        self.imfile = None
        self.exit_status = 1

        # Flag to specify whether the all inputs are fulfilled
        self.seismic_model = False
        self.electromag_model = False

    def para_check(self):
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
            print('No step size in the y-direction')
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
    angfiles : numpy.ndarray or None
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
        self.angfiles = None

        # The processing functions
        self.functions = mf

    def sort_material_list(self) -> None:
        """
        Sorts the material list based on the material properties.
        """
        self.material = self.material_list[:,1]
        self.temp = self.material_list[:,2].astype(float)
        self.rho = self.material_list[:,3].astype(float)
        self.porosity = self.material_list[:,4].astype(float)
        self.lwc = self.material_list[:,5].astype(float)
        self.is_anisotropic = self.material_list[:,6] == 'True'
        self.angfiles = self.material_list[:,7]

    def para_check(self) -> None:
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
        self.build()
    
    def build(self) -> None:
        """
        Initializes the simulation model attributes with default values.
        """
        self.dt = None
        self.time_steps = None
        self.x = None
        self.y = None
        self.z = None
        self.f0 = None
        self.theta = None
        self.phi = None
        self.src = None
        self.tensor_coefficients = None
        self.compute_coefficients = True
        self.attenuation_coefficients = None
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
        
    
    def tensor_check(self) -> None:
        """
        Checks if tensor coefficients are specified and valid. Disables 
        coefficient computation if valid.
        """
        # If the tensors are there
        check = 0
        # if not self.tensor_coefficients:
        # 	print('ldkjf')
        for row in self.tensor_coefficients:
            for val in row:
                if not val:
                    check = check + 1
        
        if check == 0:
            self.compute_coefficients = False
    
    def para_check(self) -> None:
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
def kband_check(modelclass, domain):
    '''
    Check to make sure the discretization values of the domain satisfy the 
    wavenumber bandlimiter as described in Volvcan et al. (2024).

    :param modelclass:
    :type modelclass: Model 
    :param material: 
    :type material: Material 
    :param domain:
    :type domain: Domain
    '''
    
    if modelclass.is_seismic:
        velocities = np.zeros([domain.nmats, 4])
        for ind in range(domain.nmats):
            velocities[ind,:] = mf.tensor2velocities(
                modelclass.tensor_coefficients[ind,1:-1].astype(float), 
                modelclass.tensor_coefficients[ind,-1].astype(float), 
                seismic = True
            )
    else: 
        velocities = np.zeros([domain.nmats, 3])
        for ind in range(domain.nmats):
            velocities[ind,:] = mf.tensor2velocities(
                modelclass.tensor_coefficients[ind,1:7].real, seismic = False
            )
    
    lambda_values = velocities/modelclass.f0
    # For things such as air and water, there will be zero valued velocities
    # which will be taken as the step_limit. Remove the zero values.
    lambda_values = lambda_values[lambda_values > 0]
    step_limit = lambda_values.min()/4 
    modelclass.step_limit = step_limit 
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


# -------------------------- Function Definitions ------------------------------
def read_dat(
        fn: str, 
        channel: str, 
        domain: Domain, 
        is_complex: bool, 
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
    :param is_complex: A flag indicating if the data is complex.
    :type is_complex: bool
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
        NX = domain.nz
        NZ = domain.nx
    #

    with FortranFile(fn, 'r') as f:
        if is_complex:
            if single:
                dat = f.read_reals(dtype=np.float32)#.reshape(nx, nz)
                dat = dat[:(NX*NZ)] + 1j * dat[(NX*NZ):]
            else:
                # Read complex data directly as double precision
                dat = f.read_complex(dtype=np.complex128)#.reshape(nx, nz)
        else:
            if single:
                dat = f.read_reals(dtype = np.float32)
            else:
                dat = f.read_reals(dtype = np.float64)

    if domain.dim == 2.5:
        dat = dat.reshape(NX, NY, NZ)
    else:
        dat = dat.reshape(NX, NZ)

    f.close()
    return(dat)

# =============================================================================
# ============================== Useful Functions =============================
# =============================================================================

def image2int(image_filename: str) -> np.ndarray:
    """
    Convert an image file into an integer array where each unique RGB color 
    is mapped to a unique integer.

    This function reads an image from a file, converts RGB colors to a single 
    integer value by combining the RGB channels, then maps each unique color in
    the image to a unique integer. This can be useful for image analysis tasks 
    where colors represent different categories or labels.

    :param image_filename: The path to the image file.
    :type image_filename: str
    :return: A 2D array where each element represents a unique integer mapping 
        of the original RGB colors.
    :rtype: np.ndarray

    The conversion of RGB to a single value is performed using the formula: 
    65536*red + 255*green + blue. Note that this function assumes the input 
    image is in a format readable by matplotlib.image.imread (e.g., PNG, JPG), 
    and that it has three color channels (RGB).
    """
    # read the image
    img = mpimg.imread(image_filename)

    # Convert RGB to a single value
    rgb_int = np.array(65536*img[:,:,0] +  255*img[:,:,1] + img[:,:,2])

    # Get the unique values of the image
    rgb_uni = np.unique(rgb_int)

    # mat_id = np.array( range(0, len(rgb_uni) ) )
    for ind in range(0, len(rgb_uni) ):
        rgb_int[ rgb_int == rgb_uni[ind] ] = ind

    return rgb_int.astype(int)

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
        
def append_coefficients(
        prjfile: str, tensor: np.ndarray, CP: str = None, dt : float = 1.0 
    ):
    """
    Appends coefficients to a project file based on the provided tensor, 
    optionally using a line identifier.
    
    This function reads a project file, appending new coefficients derived from 
    a tensor to it. The coefficients are appended based on a line identifier 
    (`CP`), which specifies the type of modeling (e.g., 'C' for stiffness, 'P' 
    for permittivity). The function also appends the time step (`dt`) for each 
    modeling type specified by `CP`.
    
    Note: The project file is updated in place, with the original file replaced by 
    the modified one. The `tensor` can be either a complex or a real numpy array. The function 
    handles these differently when appending to the file. This function relies on an external command (`mv`) to replace the 
    original file, which requires that the script has permissions to execute 
    system commands.
    
    :param prjfile: The path to the project file to be modified.
    :type prjfile: str
    :param tensor: A numpy array containing the tensor coefficients to append.
    :type tensor: np.ndarray
    :param CP: The line identifier indicating the type of modeling ('C' or 'P'), optional.
    :type CP: str, optional
    :param dt: The time step to append for the modeling type, defaults to 1.
    :type dt: float
    :return: None
    """
    newfile = 'newfile.txt'
    
    with open(prjfile, 'r') as ogprj, open(newfile, 'a') as temp:
        # Determine modeling type based on CP
        mt = 'S' if CP == 'C' else 'E'
        
        for line in ogprj.readlines():
            if line.startswith(CP):
                line = line.strip().split(',')
                if np.iscomplexobj(tensor):
                    temp.write(
                        CP + ',' + \
                            complex2str(tensor[int(float(line[1])), :]) + \
                                '\n'
                    )
                else:
                    temp.write(
                        CP + ',' + \
                            ','.join(map(str, tensor[int(float(line[1])), :])) + \
                                '\n'
                    )
            elif line.startswith(mt) and line[2:4] == 'dt':
                temp.write(mt + ',dt,' + str(dt) + '\n')
            else:
                temp.write(line)
    
    # Replace the original project file with the new file
    call('mv ' + newfile + ' ' + prjfile, shell=True)

# ==============================================================================
# ========================= Read/Assign Project File ===========================
# ==============================================================================
# ------------------------------------------------------------------------------
def loadproject(
        project_file: str, 
        domain: Domain, material: Material, seismic: Model, electromag: Model
    ) -> tuple[Domain, Material, Model, Model]:
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
    # domain, material, seismic, electromag are the objects that we will assign
    # values to
    f = open(project_file)

    # Let the if train begin
    for line in f:
        if line[0] == 'I':
            # There's a trailing new line value
            im = image2int(line[2:-1])
            domain.geometry = im.transpose().astype(int)
            # Get the image file
            domain.imfile = line[2:-1]
        
        # All domain inputs must be input except for ny and dy
        if line[0] == 'D':
            temp = line.split(',')
            
            if temp[1] == 'dim':
                domain.dim = float( (temp[2].rsplit())[0])
            if temp[1] == 'nx':
                domain.nx = int( (temp[2].rsplit())[0])
            if temp[1] == 'ny':
                try:
                    domain.ny = int( (temp[2].rsplit())[0])
                except:
                    pass
            if temp[1] == 'nz':
                domain.nz = int( (temp[2].rsplit())[0])
            if temp[1] == 'dx':
                domain.dx = float( (temp[2].rsplit())[0])
            if temp[1] == 'dy':
                try:
                    domain.dy = float( (temp[2].rsplit())[0])
                except:
                    pass
            if temp[1] == 'dz':
                domain.dz = float( (temp[2].rsplit())[0])
            if temp[1] == 'cpml':
                domain.cpml = int( (temp[2].rsplit())[0])
            # if temp[1] == 'write':
            #     domain.write = temp[2].rsplit()
            if temp[1] == 'nmats':
                domain.nmats = int( (temp[2].rsplit())[0])
            if temp[1] == 'cpml_alpha':
                domain.cpml_attenuation = float( (temp[2].rsplit())[0])
        
        if line[0] == 'S':
            temp = line.split(',')
            if temp[1] == 'dt':
                try :
                    seismic.dt = float( (temp[2].rsplit())[0])
                except:
                    pass
            if temp[1] == 'time_steps':
                seismic.time_steps = float( (temp[2].rsplit())[0])
            if temp[1] == 'x':
                seismic.x = float( (temp[2].rsplit())[0])
            if temp[1] == 'y':
                seismic.y = float( (temp[2].rsplit())[0])
            if temp[1] == 'z':
                seismic.z = float( (temp[2].rsplit())[0])
            if temp[1] == 'f0':
                seismic.f0 = float( (temp[2].rsplit())[0])
            if temp[1] == 'theta':
                seismic.theta = float( (temp[2].rsplit())[0])
            if temp[1] == 'phi':
                seismic.phi = float( (temp[2].rsplit())[0])

        if line[0] == 'E':
            temp = line.split(',')
            if temp[1] == 'dt':
                try:
                    electromag.dt = float( (temp[2].rsplit())[0])
                except:
                    pass
            if temp[1] == 'time_steps':
                electromag.time_steps = float( (temp[2].rsplit())[0])
            if temp[1] == 'x':
                electromag.x = float( (temp[2].rsplit())[0])
            if temp[1] == 'y':
                electromag.y = float( (temp[2].rsplit())[0])
            if temp[1] == 'z':
                electromag.z = float( (temp[2].rsplit())[0])
            if temp[1] == 'f0':
                electromag.f0 = float( (temp[2].rsplit())[0])
            if temp[1] == 'theta':
                electromag.theta = float( (temp[2].rsplit())[0])
            if temp[1] == 'phi':
                electromag.phi = float( (temp[2].rsplit())[0])

        if line[0] == 'M':
            line = line[0:-1]
            temp = line.split(',')
            if temp[1] == '0':
                material.material_list = temp[1:]
                material.rgb = temp[3]
            else:
                material.material_list = np.row_stack( (material.material_list, temp[1:]))
                material.rgb = np.append( material.rgb, temp[3] )

        if line[0] == 'C':
            temp = line.split(',')
            # We need to remove the '\n' at the end. Whether the coefficients are
            # given results in a different string
            try:
                temp[-1] = temp[-1].rsplit()[0]
            except:
                temp[-1] = ''

            if temp[1] == '0' or temp[1] == '0.0':
                seismic.tensor_coefficients = temp[1:]
            else:
                seismic.tensor_coefficients = np.row_stack((seismic.tensor_coefficients, temp[1:]))

        # Permittivity coefficients
        if line[0] == 'P':
            if 'j' in line:
                is_complex = True 
            else: 
                is_complex = False 
            
            temp = line.split(',')
            try:
                temp[-1] = temp[-1].rsplit()[0] # An index error will be given if coefficients are provided
            except:
                temp[-1] = ''

            if temp[1] == '0' or temp[1] == '0.0':
                if is_complex:
                    electromag.tensor_coefficients = str2complex(temp[1:])
                else:
                    electromag.tensor_coefficients = temp[1:]
                
            else:
                if is_complex:
                    electromag.tensor_coefficients = np.row_stack( 
                        (electromag.tensor_coefficients, str2complex(temp[1:]) )
                    )
                else:
                    electromag.tensor_coefficients = np.row_stack( 
                        (electromag.tensor_coefficients, temp[1:])
                    )
                    
        # Attenuation coefficients
        if line[0] == 'A':
            temp = line.split(',')
            try: 
                seismic.attenuation_coefficients = np.row_stack(
                    (
                        seismic.attenuation_coefficients,
                        np.array(temp[2:5], dtype = float)
                    )
                )
            except:
                seismic.attenuation_coefficients = np.array(
                    temp[2:5], dtype = float
                )
                
            # electromag.fref = float(temp[8])
            seismic.fref = float(temp[5])

    f.close()
    seismic.is_seismic = True 
    electromag.is_seismic = False
    return domain, material, seismic, electromag


# ------------------------------------------------------------------------------
# Make sure variables are in the correct type for Fortran
def prepme(
        model_object: Model, 
        domain: Domain, 
        complex_tensor: bool = True
    ) -> tuple:
    """
    Prepare modeling objects and domain configurations for simulation, ensuring
    that all parameters are in the correct format for Fortran processing.

    This function adjusts the type of various parameters in the modeling and domain objects,
    such as converting time steps to integers or ensuring that tensor coefficients are in the
    correct numerical format (complex or float). Additionally, it arranges source and domain
    parameters in the correct order based on the dimensionality of the domain.

    :param model_object: The modeling object containing simulation parameters and configurations.
    :param domain: The domain object containing spatial configurations and parameters.
    :param complex_tensor: A flag indicating whether the tensor coefficients should be treated as complex numbers.
    :type complex_tensor: bool
    :return: A tuple containing the updated modeling object and domain object.
    :rtype: tuple

    The function specifically adjusts parameters for compatibility with Fortran-based
    simulation engines, accounting for differences in array indexing and data types.
    """
    # Check if there are no errors and the coefficients have been computed
    model_object.time_steps = int(model_object.time_steps)
    model_object.f0 = float(model_object.f0)
    model_object.theta = float(model_object.theta)
    model_object.x = float(model_object.x)
    model_object.z = float(model_object.z)
    # We don't want to get rid of the complex coefficients if they are there, but we need 
    # make sure that we are only using the real components if we aren't solving for the complex
    # equations
    model_object.tensor_coefficients_original = model_object.tensor_coefficients.copy()
    if not complex_tensor:
        model_object.tensor_coefficients = model_object.tensor_coefficients.astype(float)
    
    # Put source and domain parameters in correct order
    if domain.dim == 2.5:
        # There are additional values we need to assign
        domain.ny = int(domain.ny)
        domain.dy = float(domain.dy)
        model_object.y = float(model_object.y)
        model_object.phi = float(model_object.phi)

        model_object.src = np.array(
            [
                model_object.x/domain.dx,
                model_object.y/domain.dy,
                model_object.z/domain.dz
            ]
        ).astype(int)
    else:
        model_object.src = np.array(
            [
                model_object.x/domain.dx,
                model_object.z/domain.dz
            ]
        ).astype(int)

    return(model_object, domain)

# ------------------------------------------------------------------------------
# Append coefficients
def coefs2prj(model_object, material, domain, modeltype):
    pass


# ------------------------------------------------------------------------------
def airsurf(material: Material, domain: Domain, N: int = 2) -> np.ndarray:
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
    airnum = material.material_list[material.material_list[:,1] == 'air', 0]

    if airnum:
        airnum = int(airnum[0])
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

    return(gradcomp)

# ------------------------------------------------------------------------------
def rcxgen(
        rgb: list, 
        domain: Domain, material: Material, 
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
        self.imfile = None
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
        model_object: Model, domain: Domain
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
        domain: Domain, material: Material, model: Model, 
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
