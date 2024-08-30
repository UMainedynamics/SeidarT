#!/usr/bin/env python

import os
import shutil

def main():
    print()
    print('Running Post Install')
    # Use CONDA_PREFIX if it's set; otherwise, fall back to MESON_INSTALL_PREFIX or a default
    install_prefix =  os.getenv('MESON_INSTALL_PREFIX', '') or os.getenv('CONDA_PREFIX')
    
    # Define paths
    install_path = os.path.join(install_prefix)
    target_path = os.path.join('seidart', 'fortran')

    print(f'Install path is set to {install_path}')
    print(f'Target path is set to {target_path}')
    # Ensure the target directory exists
    os.makedirs(target_path, exist_ok=True)

    # Copy the shared library to the target path
    for file in os.listdir(install_path):
        print(file)
        if file.endswith('.so'):
            shutil.copy(os.path.join(install_path, file), target_path)

if __name__ == "__main__":
    main()
