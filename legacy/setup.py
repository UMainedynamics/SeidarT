from numpy.distutils.core import setup, Extension

fortran_sources = [
    'src/seidart/fortran/cpmlfdtd.f95'
]

fortran_extension = Extension(name='seidart.fortran.cpmlfdtd',
                              sources=fortran_sources)

setup(
    name='seidart',
    version='0.0.6',
    packages=[
        'seidart', 
        'seidart.fortran', 
        'seidart.routines', 
        'seidart.simulations', 
        'seidart.visualization'
    ],
    ext_modules=[fortran_extension]
)
