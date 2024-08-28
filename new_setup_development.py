from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import numpy.f2py 

# class CustomBuildExt(build_ext):
#     def build_extensions(self):
#         # Add the Fortran compiler
#         for ext in self.extensions:
#             if ext.name == 'seidart.fortran.cpmlfdtd':
#                 ext.extra_f90_compile_args = ['-std=f95']
#         super().build_extensions()

setup(
    name='seidart',
    version='2.2.3',
    packages=[
        'seidart',
        'seidart.fortran',
        'seidart.routines',
        'seidart.simulations',
        'seidart.visualization'
    ],
    ext_modules=[
        Extension(
            name='seidart.fortran.cpmlfdtd',
            sources=['src/seidart/fortran/cpmlfdtd.f95'],
            include_dirs=[numpy.get_include()],
        ),
    ],
    # cmdclass={'build_ext': CustomBuildExt},
    entry_points={
        'console_scripts': [
            'prjbuild=seidart.routines.prjbuild:main',
            'prjrun=seidart.routines.prjrun:main',
            'arraybuild=seidart.routines.arraybuild:main',
            'sourcefunction=seidart.routines.sourcefunction:main',
            'rcxdisplay=seidart.visualization.rcxdisplay:main',
            'im2anim=seidart.visualiztion.im2anim:build_animation',
            # other entry points...
        ]
    },
    install_requires=[
        'numpy',
        'setuptools',
        'wheel',
        'matplotlib',
        'seaborn'
    ]
)
