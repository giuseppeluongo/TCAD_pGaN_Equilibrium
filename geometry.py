"""Geometry and mesh of a simple p-GaN gate HEMT (Gmsh), for the DEVSIM equilibrium model.

Cross section (y grows downwards, units um, saved in cm):

         gate
       +------+
       | pGaN |
  +----+------+-------------------+----+
  |S   |        AlGaN barrier     |   D|      S, D: ohmic contacts on the GaN
  +----+--------------------------+----+  <- 2DEG (AlGaN/GaN interface)
  |                 GaN                |
  +------------------------------------+
                  body

Run:  python geometry.py   ->  pGaN_HEMT_simple.msh
"""
import gmsh

UM_TO_CM = 1.0e-4       # Gmsh works in um, DEVSIM in cm

# ---------------- dimensions (um) ----------------
L_C     = 0.5           # ohmic contact length
L_SG    = 1.0           # source-gate distance
L_G     = 1.0           # gate length
L_GD    = 3.0           # gate-drain distance
T_PGAN  = 0.070         # p-GaN thickness
T_ALGAN = 0.015         # AlGaN barrier thickness
T_GAN   = 1.0           # GaN thickness

x0 = 0.0
x1 = x0 + L_C           # end of source contact / start of the barrier
x2 = x1 + L_SG          # gate start
x3 = x2 + L_G           # gate end
x4 = x3 + L_GD          # end of the barrier / start of drain contact
x5 = x4 + L_C           # device end

y0 = 0.0                # pGaN / AlGaN interface
y1 = y0 + T_ALGAN       # AlGaN / GaN interface (2DEG)
y2 = y1 + T_GAN         # bottom (body contact)

gmsh.initialize(["-nopopup"])
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("pGaN_HEMT_simple")
occ = gmsh.model.occ

# ---------------- 1. regions ----------------
pGaN  = occ.addRectangle(x2, y0 - T_PGAN, 0, L_G, T_PGAN)
AlGaN = occ.addRectangle(x1, y0, 0, x4 - x1, T_ALGAN)
GaN   = occ.addRectangle(x0, y1, 0, x5 - x0, T_GAN)
occ.fragment([(2, GaN), (2, AlGaN), (2, pGaN)], [])     # shared nodes on the interfaces
occ.synchronize()

# ---------------- 2. physical groups ----------------
eps = 1e-4
def curves_in(xmin, ymin, xmax, ymax):
    return [t for d, t in gmsh.model.getEntitiesInBoundingBox(
        xmin - eps, ymin - eps, -eps, xmax + eps, ymax + eps, eps, 1)]
def surfaces_in(xmin, ymin, xmax, ymax):
    return [t for d, t in gmsh.model.getEntitiesInBoundingBox(
        xmin - eps, ymin - eps, -eps, xmax + eps, ymax + eps, eps, 2)]

# regions
gmsh.model.addPhysicalGroup(2, surfaces_in(x2, y0 - T_PGAN, x3, y0), name="pGaN")
gmsh.model.addPhysicalGroup(2, surfaces_in(x1, y0, x4, y1),          name="AlGaN")
gmsh.model.addPhysicalGroup(2, surfaces_in(x0, y1, x5, y2),          name="GaN")
# interfaces
gmsh.model.addPhysicalGroup(1, curves_in(x2, y0, x3, y0), name="pGaN_AlGaN")
gmsh.model.addPhysicalGroup(1, curves_in(x1, y1, x4, y1), name="AlGaN_GaN")
# contacts
gmsh.model.addPhysicalGroup(1, curves_in(x2, y0 - T_PGAN, x3, y0 - T_PGAN), name="gate")
gmsh.model.addPhysicalGroup(1, curves_in(x0, y1, x1, y1), name="source")
gmsh.model.addPhysicalGroup(1, curves_in(x4, y1, x5, y1), name="drain")
gmsh.model.addPhysicalGroup(1, curves_in(x0, y2, x5, y2), name="body")

# an empty physical group gives no Gmsh error, but DEVSIM refuses the file
for dim, tag in gmsh.model.getPhysicalGroups():
    assert gmsh.model.getEntitiesForPhysicalGroup(dim, tag).size > 0, \
        f"empty physical group: {gmsh.model.getPhysicalName(dim, tag)}"

# ---------------- 3. mesh size ----------------
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
F = gmsh.model.mesh.field

def near_curves(curves, size_min, dist_min, dist_max, size_max=0.25):
    """Fine mesh close to some curves, coarser further away."""
    d = F.add("Distance")
    F.setNumbers(d, "CurvesList", curves)
    F.setNumber(d, "Sampling", 2000)
    t = F.add("Threshold")
    F.setNumber(t, "InField", d)
    F.setNumber(t, "SizeMin", size_min); F.setNumber(t, "SizeMax", size_max)
    F.setNumber(t, "DistMin", dist_min); F.setNumber(t, "DistMax", dist_max)
    return t

def inside_box(size, xmin, xmax, ymin, ymax, thickness=0.1, size_out=0.25):
    """Uniform mesh size inside a rectangle."""
    b = F.add("Box")
    F.setNumber(b, "VIn", size); F.setNumber(b, "VOut", size_out)
    F.setNumber(b, "XMin", xmin); F.setNumber(b, "XMax", xmax)
    F.setNumber(b, "YMin", ymin); F.setNumber(b, "YMax", ymax)
    F.setNumber(b, "ZMin", -1);   F.setNumber(b, "ZMax", 1)
    F.setNumber(b, "Thickness", thickness)
    return b

fields = [
    near_curves(curves_in(x0, y1, x5, y1), 0.002, 0.010, 0.300),   # 2DEG
    near_curves(curves_in(x2, y0, x3, y0), 0.002, 0.010, 0.200),   # pGaN / AlGaN
    near_curves(curves_in(x2, y0 - T_PGAN, x3, y0 - T_PGAN), 0.001, 0.010, 0.200),   # gate
    inside_box(0.003, x1, x4, y0, y1),                              # AlGaN barrier
    inside_box(0.004, x2, x3, y0 - T_PGAN, y0),                     # pGaN
]
f_min = F.add("Min")
F.setNumbers(f_min, "FieldsList", fields)
F.setAsBackgroundMesh(f_min)

# ---------------- 4. mesh + save ----------------
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)       # DEVSIM reads MSH 2.2 only
gmsh.option.setNumber("Mesh.ScalingFactor", UM_TO_CM)   # um -> cm in the saved file
gmsh.model.mesh.generate(2)
gmsh.write("pGaN_HEMT_simple.msh")
print("nodes:", len(gmsh.model.mesh.getNodes()[0]))
gmsh.finalize()
