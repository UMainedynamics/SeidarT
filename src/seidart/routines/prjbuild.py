#!/usr/bin/env python3
#
# This script will read an image and build the template project file template
# to be used in the seisarT program
#
# -----------------------------------------------------------------------------

import argparse
import numpy as np
import matplotlib.image as mpimg
from typing import Tuple
from scipy.io import FortranFile
import json 

__all__ = [
    'prjbuild',
    'generate_template',
    'update_json',
    'image2int',
]

# ------------------------------------------------------------------------------
def generate_template(nmats, **kwargs):
    """
    """
    template = {
        "Domain": {
            "dim": 2,
            "nx": None,
            "ny": 1,
            "nz": None,
            "dx": None,
            "dy": 1,
            "dz": None,
            "cpml": 10,
            "nmats": nmats,
            "alpha_max_scalar": 1.0,
            "kappa_max": 5,
            "NP": 2,
            "NPA": 2,
            "Rcoef": 0.0010,
            "image_file": None
        },
        "Materials": [
            {
                "id": i,
                "name": None,
                "rgb": f"{-i}/{-i}/{-i}",
                "temperature": None,
                "density": None,
                "porosity": None,
                "water_content": None,
                "is_anisotropic": None,
                "euler_angles": None
            } for i in range(nmats)
        ],
        "Seismic": {
            "Source": {
                "dt": 1e-12,
                "time_steps": 0,
                "x": 1.0,
                "y": 1.0,
                "z": 1.0,
                "xind": 0,
                "yind": 0,
                "zind": 0,
                "source_frequency": 1.0,
                "x-z_rotation": 0,
                "x-y_rotation": 0,
                "amplitude": 1.0,
                "source_type": "gaus1"
            },
            "Attenuation": [
                {
                    "id": i,
                    "gamma_x": 0.0,
                    "gamma_y": 0.0,
                    "gamma_z": 0.0,
                    "gamma_yz": 0.0,
                    "gamma_xz": 0.0,
                    "gamma_xy": 0.0,
                    "reference_frequency": 1.0
                } for i in range(nmats)
            ],
            "Stiffness_Coefficients": [
                {
                    "id": i,
                    "c11": 0.0,
                    "c12": 0.0,
                    "c13": 0.0,
                    "c14": 0.0,
                    "c15": 0.0,
                    "c16": 0.0,
                    "c22": 0.0,
                    "c23": 0.0,
                    "c24": 0.0,
                    "c25": 0.0,
                    "c26": 0.0,
                    "c33": 0.0,
                    "c34": 0.0,
                    "c35": 0.0,
                    "c36": 0.0,
                    "c44": 0.0,
                    "c45": 0.0,
                    "c46": 0.0,
                    "c55": 0.0,
                    "c56": 0.0,
                    "c66": 0.0,
                    "rho": 0.0
                } for i in range(nmats)
            ]
        },
        "Electromagnetic": {
            "Source": {
                "dt": 1e-12,
                "time_steps": 0,
                "x": 1.0,
                "y": 1.0,
                "z": 1.0,
                "xind": 0,
                "yind": 0,
                "zind": 0,
                "source_frequency": 1.0,
                "x-z_rotation": 0,
                "x-y_rotation": 0,
                "amplitude": 1.0,
                "source_type": "gaus1"
            },
            "Permittivity_Coefficients": [
                {
                    "id": i,
                    "e11": 0.0,
                    "e12": 0.0,
                    "e13": 0.0,
                    "e22": 0.0,
                    "e23": 0.0,
                    "e33": 0.0
                } for i in range(nmats)
            ],
            "Conductivity_Coefficients": [
                {
                    "id": i,
                    "s11": 0.0,
                    "s12": 0.0,
                    "s13": 0.0,
                    "s22": 0.0,
                    "s23": 0.0,
                    "s33": 0.0
                } for i in range(nmats)
            ]
        }
    }
    
    # Apply updates from **kwargs
    for kwarg_key, value in kwargs.items():
        path = parse_kwargs_to_path(kwarg_key)
        set_value_in_dict(template, path, value)
    
    return template

# ------------------------------------------------------------------------------
def parse_kwargs_to_path(kwarg_key):
    """
    Converts a kwarg like 'Domain_nx' or 'Electromagnetic.Conductivity.s11' into 'Domain.nx' 
    to handle nested dictionary paths.
    
    It also handles array indexing by converting 'Electromagnetic_0_Conductivity_0_s11' into 
    'Electromagnetic[0].Conductivity[0].s11'.
    """
    if '.' in kwarg_key:
        # If dot notation is used, return it as is
        return kwarg_key
    
    # Handle underscore notation by converting it to dot notation with array indexing
    keys = kwarg_key.split('_')
    path = []
    for key in keys:
        if key.isdigit():
            path[-1] = f'{path[-1]}[{key}]'  # Handles array index conversion
        else:
            path.append(key)
    return '.'.join(path)

# ------------------------------------------------------------------------------
# Helper function to set a value in a nested dictionary using dot notation
def set_value_in_dict(data, path, value):
    keys = path.split('.')
    d = data
    for key in keys[:-1]:
        # Handle list indices if they appear in the path
        if '[' in key and ']' in key:
            key, idx = key.split('[')
            idx = int(idx[:-1])  # Extract index and cast to integer
            d = d[key][idx]
        else:
            d = d[key]
    # Final key to set the value
    if '[' in keys[-1] and ']' in keys[-1]:
        key, idx = keys[-1].split('[')
        idx = int(idx[:-1])
        d[key][idx] = value
    else:
        d[keys[-1]] = value

# ------------------------------------------------------------------------------
def set_value_by_id(data, section, id_value, field, value):
    """
    Sets a value in a list of dictionaries where the dictionary has a specific `id`.
    
    Parameters:
    - data (dict): The dictionary to update.
    - section (str): The top-level section to look for the id (e.g., 'Materials' or 'Seismic').
    - id_value (int): The id to look for.
    - field (str): The field to update (e.g., 'rgb' or 'Attenuation.gamma_x').
    - value: The value to set.
    """
    # Find the dictionary in the section where 'id' matches id_value
    for item in data[section]:
        if item['id'] == id_value:
            # Split the field if it's a nested field (like 'Attenuation.gamma_x')
            keys = field.split('.')
            d = item
            for key in keys[:-1]:
                d = d[key]
            d[keys[-1]] = value
            break

# ------------------------------------------------------------------------------
def update_json(template, updates):
    """
    Updates the template JSON based on a dictionary of updates where the keys are a mixture of:
    - dot notation paths
    - underscore notation paths
    - tuple notation, where the index or `id` can appear at different positions in the tuple.
    
    Parameters:
    - template (dict): The JSON template to update.
    - updates (dict): A dictionary where keys are paths (dot/underscore notation) or tuples, 
                      and values are the new values to set.
                      
    # Example updates with tuple, dot, and underscore notation
    updates = {
        # Tuple with id at second index
        ("Materials", 0, "rgb"): "255/0/0",   
        # Tuple with id at third index
        ("Seismic", "Attenuation", 1, "gamma_x"): 0.1,  
        # Dot notation
        "Domain.nx": 800,  
        # Underscore notation
        "Electromagnetic_0_Conductivity_0_s11": 0.002
    }
    
    # Apply the updates
    updated_template = update_json(template, updates)
    """
    for key, value in updates.items():
        if isinstance(key, tuple):
            # If it's a tuple, determine the structure of the tuple and handle it accordingly
            if len(key) == 2:
                path = '.'.join(key)
                set_value_in_dict(template, path, value)
            elif len(key) == 3:
                # Case 1: (section, id, field) e.g. ("Materials", 0, "rgb")
                section, id_value, field = key
                set_value_by_id(template, section, id_value, field, value)
            elif len(key) == 4:
                # Case 2: (section, field, id_position, value) where id_position is dynamic
                section, field, id_value, sub_field = key
                set_value_by_id(template, section, id_value, f"{field}.{sub_field}", value)
            else:
                raise ValueError(f"Unsupported tuple format: {key}")
        else:
            # If it's not a tuple, it's a regular dot or underscore notation update
            path = parse_kwargs_to_path(key)
            set_value_in_dict(template, path, value)
    
    print("JSON file updated with the provided values.")
    return template

# ------------------------ Some Necessary Definitions -------------------------
def image2int(imfilename: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Converts an image file to a 2D array of integer values representing unique
    RGB combinations and returns the unique RGB values.

    :param imfilename: The path to the image file.
    :type imfilename: str
    :return: A tuple containing the 2D array of integer values and the array of unique RGB values.
    :rtype: Tuple[np.ndarray, np.ndarray]
    """
    # read the image
    img = mpimg.imread(imfilename)
    # Convert RGB to a single value
    rgb_int = np.array(65536*img[:,:,0] +  255*img[:,:,1] + img[:,:,2])
    # Get the unique values of the image
    rgb_uni = np.unique(rgb_int)
    # We want the unique rgb values too
    rgb = np.zeros( [len(rgb_uni), 3] )
    # reshape the image. We know it's three channels
    img_vect = np.zeros( [np.prod(rgb_int.shape), 3] )
    img_vect[:,0] = np.reshape(img[:, :, 0], np.prod(np.shape(img[:, :, 0]) ) )
    img_vect[:,1] =	np.reshape(img[:, :, 1], np.prod(np.shape(img[:, :, 1]) ) )
    img_vect[:,2] =	np.reshape(img[:, :, 2], np.prod(np.shape(img[:, :, 2]) ) )
    
    for ind in range(0, len(rgb_uni) ):
        rgb_ind = np.reshape(rgb_int == rgb_uni[ind], [np.prod(rgb_int.shape)])
        rgb[ind,:] = (img_vect[rgb_ind,:])[0,:]
        rgb_int[ rgb_int == rgb_uni[ind] ] = ind
    
    if np.max(rgb) <= 1.0:
        rgb = rgb * 255
        rgb = rgb.astype(int)
    
    rgb_int = rgb_int.astype(int)
    f = FortranFile('geometry.dat', 'w')
    f.write_record(rgb_int)
    f.close()
    
    return( rgb_int, rgb)

# ------------------------------------------------------------------------------
def readwrite_json(jsonfile: str, data: dict = None):
    """
    
    """
    if data:
        with open(jsonfile, 'w') as file:
            json.dump(data, file, indent = 4)
    else:
        with open(jsonfile, 'r') as file:
            data = json.load(file)
        return(data)
    
    file.close()
# ------------------------------------------------------------------------------

def prjbuild(
        image_file: str, 
        outputjson: str
    ) -> None:
    """
    Generates a project file (.prj) based on an input image file. This file
    contains domain parameters, material parameters, attenuation parameters,
    seismic parameters, and electromagnetic parameters derived from the image.

    :param image_file: The path to the input image file.
    :param outputjson: The path where the project file is to be saved.
    :type image_file: str
    :type outputjson: str
    
    
    """
    # ----- Read the image file
    im, rgb = image2int(image_file)
    im = im.transpose()
    material_id = np.unique(im)
    
    nx, nz = im.shape
    nmats = len(material_id)
    template = generate_template(nmats, Domain_nx = nx, Domain_nz = nz)
        
        # -------------------------------------------------------------------------
    updates = {}
    for i in range(nmats):
        updates[("Materials", material_id[i], "rgb")] = '/'.join(rgb[i].astype(str))
    
    updates[("Domain","image_file")] = image_file
    updated_template = update_json(template, updates)
    # -------------------------------------------------------------------------
    
    # Save the generated template JSON file
    with open(outputjson, 'w') as f:
        json.dump(updated_template, f, indent=4)
    

def buildwizard(jsonfile):
    pass 

def main():
    # -------------------------- Command Line Arguments ------------------------
    parser = argparse.ArgumentParser(description="""The SeidarT software 
        requires a .PNG image that is used to construct the model domain for 
        seismic and electromagnetic wave propagation. Given the image file, a 
        project file will be constructed which contains all the necessary 
        parameters to be read in to the finite differences time domain modeling 
        schemes.""" )
    
    parser.add_argument(
        '-i','--imagefile', 
        nargs=1, type=str, required = True,
        help='the full file path for the image', default=None
    )

    parser.add_argument(
        '-o', '--outputjson',
        nargs=1, type=str, required = False, default = 'steven_bernsen_rules.prj',
        help = """Name of output file path with extension .json and excluding
        the full path directory."""
    )

    # Get the arguments
    args = parser.parse_args()
    image_file = ''.join(args.imagefile)
    outputjson = ''.join(args.outputjson)
    prjbuild(image_file, outputjson)


if __name__ == "__main__":
    main()