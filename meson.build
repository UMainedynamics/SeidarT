project('seidart', 'fortran', version : '2.2.3')

python = import('python').find_installation('python3')

# Define the Fortran source files
seidart_sources = files('src/seidart/fortran/cpmlfdtd.f95')

# Build the Python extension module directly, without creating a separate static library
python.extension_module('seidart.fortran.cpmlfdtd',
    sources: seidart_sources,
    install: true,
)
