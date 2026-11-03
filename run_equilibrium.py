"""p-GaN gate HEMT at thermal equilibrium with DEVSIM.

    python geometry.py          # once: creates the mesh pGaN_HEMT_simple.msh
    python run_equilibrium.py   # solves Poisson, writes output/ (ParaView / Tecplot files)

Parameters: hemt/config.py
"""
import os
import devsim as ds

from hemt.config import DEVICE, OUT_DIR
from hemt import mesh, equilibrium

os.makedirs(OUT_DIR, exist_ok=True)

equilibrium.set_global_options()
mesh.build()
equilibrium.assign_parameters()
equilibrium.initial_guess()
equilibrium.setup_poisson()
equilibrium.setup_interfaces()
equilibrium.setup_contacts()
equilibrium.solve()

# files for ParaView (.vtm) and Tecplot (.dat): Potential, Ec, Ev, Electrons, Holes, ...
ds.write_devices(file=os.path.join(OUT_DIR, "hemt_equilibrium"),     type="vtk")
ds.write_devices(file=os.path.join(OUT_DIR, "hemt_equilibrium.dat"), type="tecplot")
print(f"results in {OUT_DIR}")
