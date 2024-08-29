# Project definition
project('seidart', 'fortran', 'c',
  version: '2.2.3',
  default_options: ['warning_level=2'])

# Define the Fortran sources
seidart_sources = files('src/seidart/fortran/cpmlfdtd.f95')

# Compile to a shared library
seidart_lib = shared_library('cpmlfdtd', seidart_sources,
                             install: true,
                             install_dir: meson.current_source_dir() / 'src' / 'seidart' / 'fortran')

# Install the shared library directly in the package directory
install_data(seidart_lib, install_dir: meson.current_source_dir() / 'src' / 'seidart' / 'fortran')
