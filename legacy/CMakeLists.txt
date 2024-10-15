cmake_minimum_required(VERSION 3.15)

project(seidart LANGUAGES Fortran C)

# Module name and Fortran source file
set(module_name "cpmlfdtd")
set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/src/seidart/fortran/cpmlfdtd.f95")
set(c_wrapper_file "${CMAKE_CURRENT_SOURCE_DIR}/src/seidart/fortran/cpmlfdtd_wrapper.c")

# Find Python development files
find_package(Python3 REQUIRED COMPONENTS Development NumPy)

# Include the Python and NumPy headers
include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})

# ------------------------------------------------------------------------------
# Define the correct output file name based on OS and architecture
set(osname "${CMAKE_SYSTEM_NAME}")
string(TOLOWER ${osname} osname_lower)
set(arch ${CMAKE_SYSTEM_PROCESSOR})

# Base name follows Python conventions
set(basefile "${module_name}.cpython-311")

# Filename generation logic for Linux and macOS
if(${osname_lower} STREQUAL "linux")
    set(filename "${basefile}-${arch}-${osname_lower}-gnu.so")
elseif(${osname_lower} STREQUAL "darwin")
    set(filename "${basefile}-${osname_lower}")
else()
    message(FATAL_ERROR "Unsupported OS: ${osname_lower}.so")
endif()

# ------------------------------------------------------------------------------
# Add the Fortran source file to the target
add_library(${module_name} MODULE ${fortran_src_file} ${c_wrapper_file})

# Remove the 'lib' prefix from the shared object name
set_target_properties(${module_name} PROPERTIES
    PREFIX ""  # This removes the 'lib' prefix
    SUFFIX ""
    OUTPUT_NAME ${filename}  # Filename without the .so extension
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
# Full path to the generated shared object
set(full_output_path "${CMAKE_CURRENT_BINARY_DIR}/${filename}")

# Make sure the correct compiler flags are set for position-independent code
target_compile_options(${module_name} PRIVATE -fPIC)


# set(CMAKE_INSTALL_PREFIX "${Python3_SITELIB}" CACHE PATH "Installation Directory" FORCE)
message(STATUS "CMAKE_INSTALL_PREFIX is set to: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "CMAKE_SEIDART_INSTALL_PREFIX is set to: $ENV{CMAKE_SEIDART_INSTALL_PREFIX}")
message(STATUS "FULL_OUTPUT_PATH for generated so is set to: ${full_output_path}")

# Install the shared object to the desired Python package location
install(FILES ${full_output_path} DESTINATION ${CMAKE_INSTALL_PREFIX})

# execute_process(
#     COMMAND ${Python3_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_path('purelib'))"
#     OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
#     OUTPUT_STRIP_TRAILING_WHITESPACE
# )

# # Install the shared object to the Python site-packages directory
# install(FILES ${full_output_path} DESTINATION ${PYTHON_SITE_PACKAGES}/seidart/fortran)
