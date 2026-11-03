"""Model parameters. Textbook / literature values for a generic p-GaN gate HEMT."""
import os

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
q    = 1.602176634e-19      # C
kB   = 1.380649e-23         # J/K
eps0 = 8.854187817e-14      # F/cm
T    = 300.0                # K
Vt   = kB * T / q           # V, thermal voltage

# ---------------------------------------------------------------------------
# Device, mesh and output
# ---------------------------------------------------------------------------
DEVICE    = "hemt"
MESH      = "hemt_mesh"
MESH_FILE = "pGaN_HEMT_simple.msh"     # created by geometry.py (coordinates in cm)
OUT_DIR   = os.path.abspath("output")

# ---------------------------------------------------------------------------
# Materials
# ---------------------------------------------------------------------------
#   Eg: band gap [eV]; dChi: electron affinity relative to GaN [eV]
#   Nc, Nv: effective densities of states [cm^-3]
MATERIALS = {
    "GaN":   dict(Permittivity=8.9  * eps0, Eg=3.39, dChi=0.00,  Nc=2.24e18, Nv=2.51e19),
    "AlGaN": dict(Permittivity=8.82 * eps0, Eg=3.82, dChi=-0.30, Nc=2.58e18, Nv=2.60e19),   # x_Al ~ 0.2
}

# Region -> material
REGION_MATERIAL = {"pGaN": "GaN", "AlGaN": "AlGaN", "GaN": "GaN"}

# ---------------------------------------------------------------------------
# Doping, with incomplete ionization
# ---------------------------------------------------------------------------
#   ND, ED, gD: donors (ED below Ec)      NA, EA, gA: acceptors (EA above Ev)
#   Mg in GaN is a DEEP acceptor (EA ~ 0.17 eV): only a few % is ionized at 300 K
DOPING = {
    "pGaN":  dict(ND=0.0,  ED=0.02, gD=2.0, NA=3e19, EA=0.17, gA=4.0),
    "AlGaN": dict(ND=1e15, ED=0.02, gD=2.0, NA=0.0,  EA=0.17, gA=4.0),
    "GaN":   dict(ND=1e15, ED=0.02, gD=2.0, NA=0.0,  EA=0.17, gA=4.0),
}

# ---------------------------------------------------------------------------
# Polarization
# ---------------------------------------------------------------------------
# Net polarization sheet charge (spontaneous + piezoelectric) of the strained AlGaN
# barrier, ~1e13 cm^-2 for x_Al ~ 0.2 (Ambacher et al., J. Appl. Phys. 85, 3222, 1999).
#   bottom of the barrier (AlGaN/GaN):  +sigma  -> attracts the 2DEG
#   top of the barrier (pGaN/AlGaN):    -sigma
SIGMA = 1.0e13              # cm^-2
SHEET_CHARGE = {
    # interface : (region where the charge is placed, sigma/q)
    "AlGaN_GaN":  ("GaN",  +SIGMA),
    "pGaN_AlGaN": ("pGaN", -SIGMA),
}

# ---------------------------------------------------------------------------
# Gate (Schottky on p-GaN)
# ---------------------------------------------------------------------------
WF_GATE = 5.1       # eV, metal work function (Ni-like)
CHI_GAN = 4.1       # eV, GaN electron affinity
