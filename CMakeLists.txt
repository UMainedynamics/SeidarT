cmake_minimum_required(VERSION 3.15)

project(seidart LANGUAGES Fortran)

set(F2PY_EXECUTABLE ${PYTHON_EXECUTABLE} -m numpy.f2py)
set(f2py_module_name "cpmlfdtd")

set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/src/seidart/fortran/cpmlfdtd.f95")

# -----------------------------------
# Define the correct output file name
set(osname "${CMAKE_SYSTEM_NAME}") 
string(TOLOWER ${osname} osname_lower)
set(arch ${CMAKE_SYSTEM_PROCESSOR})
set(basefile "${f2py_module_name}.cpython-311-${arch}-${osname_lower}")
# We need support for Linux, and MacOS from Monterey to Sequoia
if(${osname_lower} STREQUAL "linux")
    set(filename "${basefile}-gnu.so")
else()
    set(filename "${basefile}.so")
endif()

set(generated_module_file "${CMAKE_CURRENT_BINARY_DIR}/${filename}")
#set(generated_module_file "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}.cpython-311-darwin.so")

add_custom_command(
    OUTPUT ${generated_module_file}
    COMMAND ${F2PY_EXECUTABLE} -m ${f2py_module_name} -c ${fortran_src_file}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${fortran_src_file}
)

add_custom_target(${f2py_module_name}_build ALL
    DEPENDS ${generated_module_file}
)

# Install the generated module file directly
install(FILES ${generated_module_file} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.11/site-packages/seidart/fortran)
