"""Gmsh mesh import: regions, interfaces and contacts (names defined in geometry.py)."""
import devsim as ds

from hemt.config import DEVICE, MESH, MESH_FILE, REGION_MATERIAL

INTERFACES = {
    # name : (region0, region1)
    "pGaN_AlGaN": ("pGaN",  "AlGaN"),
    "AlGaN_GaN":  ("AlGaN", "GaN"),
}

CONTACTS = {
    # name : region
    "gate":   "pGaN",
    "source": "GaN",
    "drain":  "GaN",
    "body":   "GaN",
}


def build():
    ds.create_gmsh_mesh(mesh=MESH, file=MESH_FILE)
    for reg, mat in REGION_MATERIAL.items():
        ds.add_gmsh_region(gmsh_name=reg, mesh=MESH, region=reg, material=mat)
    for con, reg in CONTACTS.items():
        ds.add_gmsh_contact(gmsh_name=con, mesh=MESH, name=con, region=reg, material="metal")
    for itf, (r0, r1) in INTERFACES.items():
        ds.add_gmsh_interface(gmsh_name=itf, mesh=MESH, name=itf, region0=r0, region1=r1)
    ds.finalize_mesh(mesh=MESH)
    ds.create_device(mesh=MESH, device=DEVICE)

    print("Regions:    ", ds.get_region_list(device=DEVICE))
    print("Contacts:   ", ds.get_contact_list(device=DEVICE))
    print("Interfaces: ", ds.get_interface_list(device=DEVICE))
