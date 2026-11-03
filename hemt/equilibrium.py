"""Thermal equilibrium: Poisson equation with polarization charges.

At equilibrium the Fermi level is flat and taken as the energy reference (E_F = 0).
The only unknown is the electrostatic potential psi [V]:
    Ec = -psi - dChi          Ev = Ec - Eg                        [eV]
    n  = Nc*exp((Ec_F - Ec)/kT) = Nc*exp(( psi + dChi)/Vt)        (Boltzmann)
    p  = Nv*exp((Ev - E_F)/kT)  = Nv*exp((-psi - dChi - Eg)/Vt)
    div(eps*grad psi) = -q*(p - n + ND+ - NA-) - rho_pol
"""
import math
import devsim as ds

from hemt.config import DEVICE, Vt, q, MATERIALS, REGION_MATERIAL, DOPING, SHEET_CHARGE, WF_GATE, CHI_GAN


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def get_par(region, name):
    return ds.get_parameter(device=DEVICE, region=region, name=name)


def node_model(region, name, eq):
    """Node model + its derivative with respect to the potential (for Newton)."""
    ds.node_model(device=DEVICE, region=region, name=name, equation=eq)
    ds.node_model(device=DEVICE, region=region, name=f"{name}:Potential",
                  equation=f"simplify(diff({eq}, Potential))")


def edge_model(region, name, eq):
    ds.edge_model(device=DEVICE, region=region, name=name, equation=eq)
    for n in ("n0", "n1"):
        ds.edge_model(device=DEVICE, region=region, name=f"{name}:Potential@{n}",
                      equation=f"simplify(diff({eq}, Potential@{n}))")


# ---------------------------------------------------------------------------
# parameters and initial guess
# ---------------------------------------------------------------------------
def set_global_options():
    ds.set_parameter(name="ElectronCharge", value=q)
    ds.set_parameter(name="Vt", value=Vt)


def assign_parameters():
    """Material and doping parameters of each region (config.py)."""
    for reg in ds.get_region_list(device=DEVICE):
        for name, val in {**MATERIALS[REGION_MATERIAL[reg]], **DOPING[reg]}.items():
            ds.set_parameter(device=DEVICE, region=reg, name=name, value=val)


def neutral_potential(region):
    """Potential that makes the region charge-neutral (bisection on the net charge)."""
    ND, ED, gD = (get_par(region, k) for k in ("ND", "ED", "gD"))
    NA, EA, gA = (get_par(region, k) for k in ("NA", "EA", "gA"))
    Eg, dChi, Nc, Nv = (get_par(region, k) for k in ("Eg", "dChi", "Nc", "Nv"))

    def net_charge(psi):                        # rho/q [cm^-3]
        eta_n = (psi + dChi) / Vt
        eta_p = (-psi - dChi - Eg) / Vt
        n, p = Nc * math.exp(eta_n), Nv * math.exp(eta_p)
        nd = ND / (1 + gD * math.exp(eta_n + ED / Vt))
        na = NA / (1 + gA * math.exp(eta_p + EA / Vt))
        return p - n + nd - na

    lo, hi = -dChi - Eg - 1.0, -dChi + 1.0     # from below Ev to above Ec
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if net_charge(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi)


def initial_guess():
    for reg in ds.get_region_list(device=DEVICE):
        ds.node_solution(device=DEVICE, region=reg, name="Potential")
        ds.set_node_value(device=DEVICE, region=reg, name="Potential", value=neutral_potential(reg))
        ds.edge_from_node_model(device=DEVICE, region=reg, node_model="Potential")


# ---------------------------------------------------------------------------
# Poisson equation in the bulk
# ---------------------------------------------------------------------------
def setup_poisson():
    for reg in ds.get_region_list(device=DEVICE):
        node_model(reg, "Ec", "-Potential - dChi")
        node_model(reg, "Ev", "Ec - Eg")
        node_model(reg, "eta_n", "(Potential + dChi)/Vt")
        node_model(reg, "eta_p", "(-Potential - dChi - Eg)/Vt")
        node_model(reg, "Electrons", "Nc*exp(eta_n)")
        node_model(reg, "Holes",     "Nv*exp(eta_p)")
        # incomplete ionization of donors and acceptors
        node_model(reg, "NDion", "ND/(1 + gD*exp(eta_n + ED/Vt))")
        node_model(reg, "NAion", "NA/(1 + gA*exp(eta_p + EA/Vt))")

        # polarization charge density, filled at the interfaces by add_sheet_charge
        ds.node_solution(device=DEVICE, region=reg, name="PolRho")
        ds.set_node_value(device=DEVICE, region=reg, name="PolRho", value=0.0)

        # DEVSIM form:  div(eps*E) + Z = 0  with  Z = -rho
        node_model(reg, "PotentialNodeCharge",
                   "-ElectronCharge*(Holes - Electrons + NDion - NAion) - PolRho")
        edge_model(reg, "ElectricField", "(Potential@n0 - Potential@n1)*EdgeInverseLength")
        edge_model(reg, "PotentialEdgeFlux", "Permittivity*ElectricField")
        ds.equation(device=DEVICE, region=reg, name="PotentialEquation",
                    variable_name="Potential",
                    node_model="PotentialNodeCharge",
                    edge_model="PotentialEdgeFlux",
                    variable_update="log_damp")


# ---------------------------------------------------------------------------
# Interfaces: continuous potential + polarization sheet charges
# ---------------------------------------------------------------------------
def add_sheet_charge(interface, region, sigma):
    """Sheet charge sigma [cm^-2] spread on the interface nodes of a region:
    rho = q*sigma*SurfaceArea/NodeVolume (its integral is exactly q*sigma*area)."""
    elems = ds.get_element_node_list(device=DEVICE, region=region, interface=interface)
    nodes = sorted({n for e in elems for n in e})
    SA = ds.get_node_model_values(device=DEVICE, region=region, name="SurfaceArea")
    NV = ds.get_node_model_values(device=DEVICE, region=region, name="NodeVolume")
    rho = list(ds.get_node_model_values(device=DEVICE, region=region, name="PolRho"))
    for i in nodes:
        rho[i] += q * sigma * SA[i] / NV[i]
    ds.set_node_values(device=DEVICE, region=region, name="PolRho", values=rho)


def setup_interfaces():
    for itf in ds.get_interface_list(device=DEVICE):
        ds.interface_model(device=DEVICE, interface=itf, name="continuousPotential",
                           equation="Potential@r0 - Potential@r1")
        ds.interface_model(device=DEVICE, interface=itf,
                           name="continuousPotential:Potential@r0", equation="1")
        ds.interface_model(device=DEVICE, interface=itf,
                           name="continuousPotential:Potential@r1", equation="-1")
        ds.interface_equation(device=DEVICE, interface=itf, name="PotentialEquation",
                              interface_model="continuousPotential", type="continuous")
    for itf, (reg, sigma) in SHEET_CHARGE.items():
        add_sheet_charge(itf, reg, sigma)


# ---------------------------------------------------------------------------
# Contacts: fixed potential (all biases = 0 at equilibrium)
#   ohmic (source, drain, body):  psi = charge-neutral potential of the region
#   gate (ideal Schottky):        Ec = E_F + phi_B  ->  psi = -phi_B - dChi
# ---------------------------------------------------------------------------
def contact_potential(contact):
    reg = ds.get_region_list(device=DEVICE, contact=contact)[0]
    if contact == "gate":
        phi_B = WF_GATE - CHI_GAN - get_par(reg, "dChi")     # Schottky barrier for electrons
        return -phi_B - get_par(reg, "dChi")
    return neutral_potential(reg)


def setup_contacts():
    for c in ds.get_contact_list(device=DEVICE):
        ds.set_parameter(device=DEVICE, name=f"{c}_psi", value=contact_potential(c))
        ds.contact_node_model(device=DEVICE, contact=c, name=f"{c}_bc",
                              equation=f"Potential - {c}_psi")
        ds.contact_node_model(device=DEVICE, contact=c, name=f"{c}_bc:Potential",
                              equation="1")
        ds.contact_equation(device=DEVICE, contact=c, name="PotentialEquation",
                            node_model=f"{c}_bc", edge_charge_model="PotentialEdgeFlux")


def solve():
    ds.solve(type="dc", absolute_error=1.0, relative_error=1e-10, maximum_iterations=100)
