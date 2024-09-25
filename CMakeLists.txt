cmake_minimum_required(VERSION 3.15)

project(seidart LANGUAGES Fortran)

# Module name and Fortran source file
set(module_name "cpmlfdtd")
set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/src/seidart/fortran/cpmlfdtd.f95")

# Define the correct output file name based on OS and architecture
set(osname "${CMAKE_SYSTEM_NAME}")
string(TOLOWER ${osname} osname_lower)
set(arch ${CMAKE_SYSTEM_PROCESSOR})

# Base name follows Python conventions
set(basefile "${module_name}.cpython-311")

# Filename generation logic for Linux and macOS
if(${osname_lower} STREQUAL "linux")
    set(filename "${basefile}-${arch}-${osname_lower}-gnu.so")
else()
    set(filename "${basefile}-${osname_lower}.so")
endif()

# -----------------------------------
# Call f2py to compile the Fortran code into a Python extension module
add_custom_command(
    OUTPUT ${filename}
    COMMAND ${PYTHON_EXECUTABLE} -m numpy.f2py -m ${module_name} -c ${fortran_src_file} -o ${filename}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${fortran_src_file}
)

add_custom_target(${module_name}_build ALL
    DEPENDS ${filename}
)

# -----------------------------------
# Install the generated shared library to the Python package's Fortran folder
install(FILES ${filename} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.11/site-packages/seidart/fortran)
