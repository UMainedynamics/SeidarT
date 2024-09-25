cmake_minimum_required(VERSION 3.15)

project(seidart LANGUAGES Fortran)

set(module_name "cpmlfdtd")
set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/src/seidart/fortran/cpmlfdtd.f95")

# -----------------------------------
# Define the correct output file name
set(osname "${CMAKE_SYSTEM_NAME}") 
string(TOLOWER ${osname} osname_lower)
set(arch ${CMAKE_SYSTEM_PROCESSOR})

set(basefile "${module_name}.cpython-311")

# We need support for Linux, and MacOS from Monterey to Sequoia
if(${osname_lower} STREQUAL "linux")
    set(filename "${basefile}-${arch}-${osname_lower}-gnu.so")
else()
    set(filename "${basefile}-${osname_lower}.so")
endif()

#set(generated_module_file "${CMAKE_CURRENT_BINARY_DIR}/${filename}")
#set(generated_module_file "${CMAKE_CURRENT_BINARY_DIR}/${module_name}.cpython-311-darwin.so")

# -----------------------------------
# Create a shared library from the Fortran source file
add_library(${module_name} MODULE ${fortran_src_file})

# Set the output file name for the shared library
set_target_properties(${module_name} PROPERTIES 
    PREFIX ""               # No 'lib' prefix
    SUFFIX ""            # Ensure '.so' extension
    OUTPUT_NAME ${filename} # Set correct output name
)
# -----------------------------------
# Install the generated shared library to the Python package's Fortran folder
install(TARGETS ${module_name}
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.11/site-packages/seidart/fortran)

