import ctypes
import importlib.resources as resources

# Use `files()` to locate the shared object file
lib_path = resources.files('seidart.fortran').joinpath('cpmlfdtd.so')

# Convert to a string and load the shared object using ctypes
_cpmlfdtd = ctypes.CDLL(str(lib_path))

# Define Python wrappers for the Fortran subroutines
# You can define argument types if necessary, like this:
# _cpmlfdtd.permittivity_write.argtypes = [ctypes.c_double, ctypes.c_int]
# _cpmlfdtd.permittivity_write.restype = None  # Assuming Fortran subroutine has no return

def permittivity_write(*args):
    return _cpmlfdtd.permittivity_write(*args)

def permittivity_write_c(*args):
    return _cpmlfdtd.permittivity_write_c(*args)

def attenuation_write(*args):
    return _cpmlfdtd.attenuation_write(*args)

def stiffness_write(*args):
    return _cpmlfdtd.stiffness_write(*args)

def seismic2(*args):
    return _cpmlfdtd.seismic2(*args)

def seismic25(*args):
    return _cpmlfdtd.seismic25(*args)

def electromag2(*args):
    return _cpmlfdtd.electromag2(*args)

def electromag25(*args):
    return _cpmlfdtd.electromag25(*args)

def electromag2c(*args):
    return _cpmlfdtd.electromag2c(*args)

def electromag25c(*args):
    return _cpmlfdtd.electromag25c(*args)
