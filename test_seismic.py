import numpy as np
import pandas as pd
from seidart.routines.domain import Domain
from seidart.routines.material import Material
from seidart.routines.classes import Biot

# Initialize mock domain and material
domain = Domain()
domain.dim = 2.0
domain.nx, domain.nz = 100, 100
domain.dx, domain.dz = 1.0, 1.0
domain.cpml = 10
domain.geometry = np.zeros((100, 100), dtype=int)

material = Material()
material.material = ['sand']
material.temp = [20.0]
material.rho = [2000.0]
material.porosity = [0.3]
material.lwc = [0.0]

# Run Biot coefficient builder to ensure no crashes
b = Biot()
b.build(material, domain, write_tensor=False)
print("CPML usage skipped successfully. Stiffness and attenuation coefficients generated.")
