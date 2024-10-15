project('seidart', ['fortran','c'],
  version : '2.3.2',
  license : 'GPL-3.0-or-later'
)

# Locate Python interpreter
python_mod = import('python')
python = python_mod.find_installation(pure: false)
python_dep = python.dependency() 

# ------------------------------------------------------------------------------
# build_dir = meson.current_build_dir()

# # Include the build directory where .mod files are generated
# inc_dirs = include_directories(build_dir)

# ------------------------------------------------------------------------------
incdir_numpy = run_command(python, 
    ['-c', 'import numpy; print(numpy.get_include())'], #'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
    check : true
).stdout().strip() 

incdir_f2py = run_command(python,
    ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

# ------------------------------------------------------------------------------
# Compile the Fortran file to generate .mod files
# fortran_module_compilation = custom_target(
#   'compile_cpmlfdtd_mod',
#   input: 'src/seidart/fortran/cpmlfdtd.f95',
#   output: 'cpmlfdtd.mod',  # Expected .mod file
#   command: ['gfortran', '-c', '@INPUT@', '-J' + meson.build_root()],  # Compile with -J to specify the .mod directory
#   build_by_default: true
# )

# Run the shell script to compile modules
module_compilation = run_command('sh', 'mod_build.sh', meson.current_build_dir(), check: true)
if module_compilation.returncode() != 0
    error('Failed to compile Fortran modules')
endif

# Specify the Fortran source file
fortran_sources = custom_target('generate_cpmlfdtd_c_and_f90', #'cpmlfdtd.c',
  input : ['src/seidart/fortran/cpmlfdtd.f95'],
  output : ['cpmlfdtdmodule.c', 'cpmlfdtd-f2pywrappers2.f90'],
  command : [python, '-m', 'numpy.f2py', '@INPUT@', '-m', 'cpmlfdtd', '--build-dir', '@OUTDIR@', '--lower'],
  build_by_default : true
)

python.extension_module('cpmlfdtd',
    [fortran_sources, incdir_f2py / 'fortranobject.c'],
    include_directories: inc_np,
    dependencies : python_dep,
    install : true
)

# ------------------------------------------------------------------------------
install_dir = python.get_install_dir()
message('Python install directory: ' + install_dir)

# Install the Fortran module within the 'seidart.fortran' package
install_subdir('src/seidart/routines', install_dir: install_dir / 'seidart')
install_subdir('src/seidart/fortran', install_dir: install_dir / 'seidart')
install_subdir('src/seidart/visualization', install_dir: install_dir / 'seidart')
install_subdir('src/seidart/simulations', install_dir: install_dir / 'seidart')

# Install the license file (GNU GPL)
install_data('LICENSE.txt', install_dir: get_option('datadir') / 'licenses')
