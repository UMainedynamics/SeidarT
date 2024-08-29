project('seidart', 'fortran', 
        version : '2.2.3',
        default_options : ['warning_level=3'])

py3_mod = import('python3')
py3_dep = dependency('python3')

py3_mod.extension_module('seidart',
    'src/seidart/fortran/cpmlfdtd.f95',
    dependencies : py3_dep
)