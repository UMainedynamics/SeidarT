cmake_minimum_required(VERSION 3.13)

project(seidart LANGUAGES Fortran)

# Set policy for deprecated Python modules
if(POLICY CMP0148)
  cmake_policy(SET CMP0148 NEW)
endif()

# Find Python and NumPy packages
find_package(PythonExtensions REQUIRED)
find_package(NumPy REQUIRED)

# Add the Fortran module
add_library(cpmlfdtd MODULE cpmlfdtd.f95)

# Include directories for NumPy
target_include_directories(cpmlfdtd PRIVATE ${NumPy_INCLUDE_DIRS})

# Mark the library as a Python extension module
python_extension_module(cpmlfdtd)

# Specify the installation directory
install(TARGETS cpmlfdtd LIBRARY DESTINATION .)
