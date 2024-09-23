from skbuild import setup
from setuptools import find_packages

setup(
    name='seidart',
    version='2.3.2',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    package_data={
        'seidart': ['fortran/*.so'],  # Include shared libraries
    },
)
