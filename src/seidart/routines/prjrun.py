import argparse
import os.path
import os
import numpy as np
import matplotlib.image as mpimg
from subprocess import call
from typing import Tuple
# import seidart.routines.materials as mf
from seidart.routines.definitions import *

# Modeling modules
from seidart.fortran.cpmlfdtd import cpmlfdtd

# Global Constants
clight = 2.99792458e8  # Speed of light in vacuum (m/s)
NP = 2  # Numerical parameter for CPML
NPA = 2  # Additional numerical parameter for CPML
k_max = 1.1e1  # Max value for CPML parameter
eps0 = 8.85418782e-12  # Permittivity of free space
mu0 = 4.0 * np.pi * 1.0e-7  # Permeability of free space
Rcoef = 0.0010  # Reflection coefficient, used for seismic only

# ============================ Create the objects =============================
# Let's initiate the domain and check to make sure everything is good to go
def domain_initialization(
        prjfile: str
    ) -> Tuple[Domain, Material, Model, Model]:
    """
    Initializes the simulation domain and materials, and prepares seismic and 
    electromagnetic models based on a project file.

    This function reads simulation parameters from a specified project file and 
    uses them to initialize the domain, material properties, and two models 
    (seismic and electromagnetic). It performs checks to ensure all necessary 
    parameters are specified and modifies the material list by removing RGB 
    values. 

    :param prjfile: The path to the project file containing initialization 
        parameters for the domain, materials, and models.
    :type prjfile: str
    :return: A tuple containing initialized Domain, Material, and two Model 
        instances (seismic and electromagnetic).
    :rtype: Tuple[Domain, Material, Model, Model]

    Each Model instance is checked to ensure required parameters are provided, 
    and the material list is pruned to exclude RGB values before final checks 
    are performed.
    """

    domain, material, seismic, electromag = loadproject(
        prjfile,
        Domain(), 
        Material(),
        Model(),
        Model()
    )

    # Model initialization and checks
    domain.para_check()
    seismic.para_check()
    electromag.para_check()
    seismic.tensor_check()
    electromag.tensor_check()

    # Modify material list and perform final checks
    material.material_list = np.delete(material.material_list, 2, axis=1)
    material.para_check()

    return(domain, material, seismic, electromag)

# -----------------------------------------------------------------------------
def status_check(
        modelclass: Model, 
        material: Material, 
        domain: Domain, 
        prjfile: str, 
        seismic: bool = True, 
        append_to_prjfile: bool = True
    ) -> None:
    """
    Checks the status of the modeling classes and appends coefficients to the
    project file if necessary.

    :param modelclass: The model class to check.
    :param material: The material class instance.
    :param domain: The domain class instance.
    :param prjfile: Path to the project file.
    :param seismic: Flag indicating if the model is seismic. Defaults to True.
    :param append_to_prjfile: Flag to append coefficients to the project file. 
                              Defaults to True.
    :type modelclass: Union[Model, None]
    :type material: Material
    :type domain: Domain
    :type prjfile: str
    :type seismic: bool
    :type append_to_prjfile: bool
    """
    if modelclass.exit_status == 0 and \
        material.material_flag and \
            append_to_prjfile:
        # The coefficients aren't provided but the materials are so we can 
        # compute them
        # Assign the materials to their respective corners
        material.sort_material_list()
        
        if seismic:
            print('Computing the stiffness coefficients.')
            tensor = material.functions.get_seismic(
                temp = material.temp, 
                rho = material.rho,
                porosity = material.porosity,
                lwc = material.lwc, 
                anisotropic = material.is_anisotropic,
                angfile = material.angfiles, 
                material_name = material.material
            )
            modelclass.tensor_coefficients = tensor
        else:
            print('Computing the permittivity and conductivity coefficients.')
            
            tensor = material.functions.get_perm(
                material, modelclass
            )
            modelclass.tensor_coefficients = tensor
        
        # Before we append the coefficients to the text file let's round to the second decimal
        tensor = np.round(tensor, 2)
        if seismic:
            ind = np.where(tensor.max() == tensor)
            max_rho = tensor[ ind[0][0], -1]
        
        # We're going to find the lines marked 'C' or 'P' and input the values there

        if seismic:
            modelclass.dt = np.min([domain.dx, domain.dz]) / np.sqrt(3.0 * tensor.max()/max_rho )
            append_coefficients(prjfile, tensor, CP = 'C', dt = modelclass.dt)
        else:
            modelclass.dt = np.min([domain.dx, domain.dz]) / \
                (2.0 * clight/ \
                    np.sqrt(np.min(
                        [
                            tensor[:,1].astype(float).min(), 
                            tensor[:,4].astype(float).min()
                        ]
                    )) 
                )
            append_coefficients(prjfile, tensor, CP = 'P', dt = modelclass.dt)
        
        # The time step needs to satisfy the Courant number and also have a nyquist
        # that will resolve the source frequency
        src_nyquist = 1/(2*modelclass.f0)
        if src_nyquist < modelclass.dt:
            print(
                '''Nyquist is not small enough for the source frequency. Change
                the source frequency or decrease the spatial step size'''
            )

        print("Finished. Appending to project file.\n")

# -----------------------------------------------------------------------------
def cpmlcompute(
        modelclass: Model, 
        domain: Domain, 
        direction: str, 
        half: bool = False, 
        seismic: bool = True
    ) -> None:
    """
    Computes CPML parameters for a given direction and updates model/domain.

    :param modelclass: The model class instance to update.
    :param domain: The domain class instance to update.
    :param direction: Direction to compute CPML ('x', 'y', or 'z').
    :param half: Flag to compute half CPML parameters. Defaults to False.
    :param seismic: Flag indicating if the model is seismic. Defaults to True.
    :type modelclass: Model
    :type domain: Domain
    :type direction: str
    :type half: bool
    :type seismic: bool
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

    quasi_cp_max = 0.7* deltamin / (2.0 * modelclass.dt)
    alpha_max = np.pi*modelclass.f0
    if seismic:
        sig_max = - np.log(Rcoef) * (NP+1) * quasi_cp_max / (2.0 * domain.cpml )
    else:
        sig_max = 0.7 * (NP+1) / (dx * np.sqrt(mu0/eps0) )

    kappa = np.ones([N])
    alpha = np.zeros([N])
    sigma = np.zeros([N])
    acoeff = np.zeros([N])
    bcoeff = np.zeros([N])

    # Compute in the x, and z directions
    for ind in range(0, domain.cpml):
        # From 0
        sigma[ind] = sig_max*dist[ind]**NP
        kappa[ind] = 1.0 + (k_max - 1.0) * dist[ind]**NP
        alpha[ind] = alpha_max * (1 - dist[ind])**NPA
        sigma[-(ind+1)] = sig_max*dist[ind]**NP
        kappa[-(ind+1)] = 1 + (k_max - 1) * dist[ind]**NP
        alpha[-(ind+1)] = alpha_max * (1 - dist[ind])**NPA
        bcoeff[-(ind+1)] = np.exp(- (sigma[-(ind+1)] / kappa[-(ind+1)] + alpha[-(ind+1)]) * modelclass.dt)
        bcoeff[ind] = np.exp( - (sigma[ind] / kappa[ind] + alpha[ind]) * modelclass.dt)

    # Compute the a-coefficients 
    alpha[np.where(alpha < 0.0)] = 0.0
    indices = np.where(np.abs(sigma) > 1.0e-6)
    acoeff[indices] = sigma[indices] * (bcoeff[indices] - 1) / \
            (kappa[indices] * sigma[indices] + kappa[indices] * alpha[indices] )

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

# ================================== SEISMIC ==================================
def runseismic(
        modelclass: Model, 
        material: Material, 
        domain: Domain, 
        single_precision: bool = True
    ) -> None:
    """
    Runs the seismic model using the initialized modelclass, material, and domain.

    :param modelclass: The seismic model class instance.
    :param material: The material class instance.
    :param domain: The domain class instance.
    :param single_precision: Flag for single precision data processing.
    :type modelclass: Model
    :type material: Material
    :type domain: Domain
    :type single_precision: bool
    """

    modelclass, domain = prepme(modelclass, domain, complex_tensor = False)
    direction = ['x', 'y', 'z']
    # Compute CPML
    print(direction)
    print('computing cpml')
    for d in direction:
        cpmlcompute(modelclass, domain, d, seismic = True)
        cpmlcompute(modelclass, domain, d, half = True, seismic = True)
    
    # We need to set a density gradient at air interfaces because high
    # density gradients lead to numerical instability
    rhograd = airsurf(material, domain, 2)
    # Write the coefficient images to a fortran file
    
    # Create the stiffness .dat files
    cpmlfdtd.stiffness_write(
        domain.geometry + 1,
        modelclass.tensor_coefficients,
        domain.cpml,
        rhograd,
        domain.nx,
        domain.nz
    )
    # attenuation_coefficients = modelclass.attenuation_coefficients * \
    #     modelclass.dt / (np.ones([domain.nmats, 3]) * \
    #         modelclass.attenuation_coefficients.astype(float) ).T
    # Create the attenuation .dat files
    cpmlfdtd.attenuation_write(
        domain.geometry + 1,
        modelclass.attenuation_coefficients,
        domain.cpml,
        # domain.nx, 
        # domain.nz,
        domain.cpml_attenuation
    )
    
    if domain.dim == 2.5:
        print('Running 2.5D model')
        # Run the FDTD
        cpmlfdtd.seismic25(
            domain.nx + 2*domain.cpml, 
            domain.ny + 2*domain.cpml, 
            domain.nz + 2*domain.cpml,
            domain.dx, domain.dy, domain.dz,
            domain.cpml,
            modelclass.src,
            modelclass.time_steps,
            single_precision
        )
    else:
            print('Running 2D model')
            cpmlfdtd.seismic2(
                domain.nx + 2*domain.cpml, 
                domain.nz + 2*domain.cpml,
                domain.dx, domain.dz,
                domain.cpml,
                modelclass.src,
                modelclass.time_steps,
                single_precision
            )     

# =============================================================================
def runelectromag(
        modelclass: Model, 
        material: Material, 
        domain: Domain, 
        use_complex_equations: bool = False,
        single_precision: bool = True
    ) -> None:
    """
    Runs the electromagnetic model with options for complex equations and precision.

    :param modelclass: The electromagnetic model class instance.
    :param material: The material class instance.
    :param domain: The domain class instance.
    :param use_complex_equations: Flag to use complex equations in modeling.
    :param single_precision: Flag for single precision data processing.
    :type modelclass: Model
    :type material: Material
    :type domain: Domain
    :type use_complex_equations: bool
    :type single_precision: bool
    """

    modelclass, domain = prepme(
        modelclass, domain, complex_tensor = use_complex_equations
    )
    direction = ['x', 'y', 'z']
    # Compute CPML
    print(direction)
    print('computing cpml')
    for d in direction:
        cpmlcompute(modelclass, domain, d, seismic = False)
        cpmlcompute(modelclass, domain, d, half = True, seismic = False)
    
    if use_complex_equations:
        cpmlfdtd.permittivity_write_c(
            domain.geometry+1,
            modelclass.tensor_coefficients,
            domain.cpml,
            domain.nx, 
            domain.nz
        )
        if domain.dim == 2.5:
            print('Running complex 3D model.')
            cpmlfdtd.electromag25c(
                domain.nx + 2*domain.cpml, 
                domain.ny + 2*domain.cpml, 
                domain.nz + 2*domain.cpml,
                domain.dx, domain.dy, domain.dz,
                domain.cpml,
                modelclass.src,
                modelclass.time_steps,
                single_precision
            )
        else:
            print('Running complex 2D model')
            cpmlfdtd.electromag2c(
                domain.nx + 2*domain.cpml,
                domain.nz + 2*domain.cpml,
                domain.dx, domain.dz,
                domain.cpml,
                modelclass.src,
                modelclass.time_steps,
                single_precision
            )
    else:  
        cpmlfdtd.permittivity_write(
                domain.geometry+1,
                modelclass.tensor_coefficients.real,
                domain.cpml,
                domain.nx, 
                domain.nz
            )
        if domain.dim == 2.5:
            print('Running 2.5D model')
            cpmlfdtd.electromag25(
                domain.nx + 2*domain.cpml, 
                domain.ny + 2*domain.cpml, 
                domain.nz + 2*domain.cpml,
                domain.dx, domain.dy, domain.dz,
                domain.cpml,
                modelclass.src,
                modelclass.time_steps,
                single_precision
            )
        else:
            print('Running 2D model')
            cpmlfdtd.electromag2(
                domain.nx + 2*domain.cpml,
                domain.nz + 2*domain.cpml,
                domain.dx, domain.dz,
                domain.cpml,
                modelclass.src,
                modelclass.time_steps,
                single_precision
            )
        
def main(
        prjfile: str, 
        model_type: str, 
        append_to_prjfile: bool, 
        pwd: str, 
        double_precision: bool, 
        use_complex_equations: bool
    ) -> None:
    """
    Main function orchestrating domain initialization, status checking, and 
    model execution based on provided command line arguments.

    Initializes the simulation domain and material properties from a project file,
    checks the status of seismic and electromagnetic models, and runs the specified
    model(s) with configured settings. The function allows for re-computation and 
    appending of material coefficients to the project file and supports both
    seismic and electromagnetic model types with options for precision and complex
    equations.

    :param prjfile: The path to the project file containing simulation settings.
    :param model_type: Specifies the type of model to run ('s' for seismic, 'e' 
                       for electromagnetic).
    :param append_to_prjfile: Flag indicating whether to append computed coefficients 
                              to the project file.
    :param pwd: The working directory where the project file is located.
    :param double_precision: Flag indicating whether to use double precision for 
                             simulation outputs.
    :param use_complex_equations: Flag indicating whether to use complex equations 
                                  for electromagnetic simulations.
    :type prjfile: str
    :type model_type: str
    :type append_to_prjfile: bool
    :type pwd: str
    :type double_precision: bool
    :type use_complex_equations: bool
    """
    domain, material, seismic, electromag = domain_initialization(prjfile)
    status_check(
        seismic, 
        material,
        domain,
        prjfile, 
        seismic=True, 
        append_to_prjfile = append_to_prjfile
    )
    status_check(
        electromag, 
        material, 
        domain,
        prjfile, 
        seismic=False, 
        append_to_prjfile = append_to_prjfile,
        use_complex_equations = use_complex_equations
    )
    
    
    if model_type == 's':
        runseismic(
            seismic, material, domain, 
            double_precision
        )
    if model_type == 'e':
        runelectromag(
            electromag, material, domain, 
            use_complex_equations = use_complex_equations, 
            double_precision = double_precision
        )

# -------------------------- Command Line Arguments ---------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        This script will read the parameters given from a project file and run 
        the specified models. The SeidarT software requires a .PNG image that is used to construct the model domain for seismic and
        electromagnetic wave propagation. Given the image file, a project file
        will be constructed which contains all the necessary parameters to be
        read in to the finite differences time domain modeling schemes."""
    )

    parser.add_argument(
        '-p', '--prjfile', nargs=1, type=str, required = True,
        help='the full file path for the project file', default=None
    )

    parser.add_argument(
        '-m', '--model', nargs = 1, type = str, required = False,
        help = """Specify whether to run the seismic (s), or electromagnetic (e), 
        or none (default = n)""",
        default = 'n'
    )

    parser.add_argument(
        '-a', '--append',
        action='store_true', required = False,
        help = """Append/recompute the coefficients to the permittivity and
        stiffness matrices; DEFAULT is True"""
    )
    
    parser.add_argument(
        '-d', '--double_precision', 
        action='store_false', required = False,
        help = """Specify double precision output of the simulation. If 
        complex, the outputs are real valued of each of the components of the 
        complex value. The default is True."""
    )
    
    parser.add_argument(
        '-c', '--use_complex_equations', action = 'store_true', required = False,
        help = """Flag whether to use the complex permittivity in the model
        simulation. The complex permittivity will be computed for ice and snow
        in all situations, but unless specified here, the complex permittivity 
        will only be used to compute the conductivity."""
    )
    
    parser.add_argument(
        '-o', '--output_directory_path', required = False, 
        nargs = 1, type = str, 
        help = '''Specify the output directory folder path to write the Fortran
        outputs to. '''
    )
    # parser.add_argument(
    #     '-A', '--attenuation', action='store_true', required = False,
    #     help = """Specify whether to include attenuation in the model. The 
    #     default is False."""
    # )

    # Get the arguments
    args = parser.parse_args()
    prjfile = ''.join(args.prjfile)
    model_type = ''.join(args.model)
    append_to_prjfile = args.append
    pwd = os.path.dirname(prjfile)
    double_precision = args.double_precision
    use_complex_equations = args.use_complex_equations
     
    # attenuation = args.attenuation
    #
    main(
        prjfile, 
        model_type, 
        append_to_prjfile, 
        pwd, 
        double_precision, 
        use_complex_equations
    )