import re
from nugridpy.mesa import history_data
import numpy as np

ISOTOPES = ['h1', 'h2', 'he3', 'he4', 'li7', 'be7', 'b8', 'c12', 'c13', 'n13', 'n14', 'n15', 'o16', 'o17', 'o18', 'f19', 'ne22']

def get_element(iso):
    """Returns the element part of an isotope of the form e.g. `c13`"""
    return re.match(r"[a-zA-Z]+", iso)[0]

def get_A(isotope):
     return int(re.match(r"\d+", iso)[0])

    
    
def get_X(star: history_data, iso):
    col = f"surface_{iso}"
    return star.get(col)
    

def get_X0(star: history_data, iso: str):
    """Returns the initial abundance of the given history_data object for the isotope iso"""
    return get_X(star, iso)[0]


def mass_yield(star: history_data, iso: str):
    """Given a history_data object and an isotope, calculates
    the mass yield of the star
    """
    dm = get_dm(star)
    Xs = get_X(star, iso)    
    y = np.sum(Xs * dm)
    return y


def mass_net_yield(star, iso):
    dm = get_dm(star)
    Xs = get_X(star, iso) - get_X0(star, iso)
    y = np.sum(Xs * dm)
    return y

def mass_net_yield_midpoint(star, iso):
    ms = star.get("star_mass")
    dm = -(ms[1:] - ms[:-1])

    Xs = get_X(star, iso) - get_X0(star, iso)

    dX = (Xs[1:] + Xs[:-1])/2 
    y = np.sum(dX * dm)
    return y


def is_primary(star):
    ϵ = 0.01
    if ϵ > np.abs(star.get("star_1_mass")[0] - star.get("star_mass")[0]):
        return True
    if ϵ > np.abs(star.get("star_2_mass")[0] - star.get("star_mass")[0]):
        return False
    raise ValueError("could not determine which star is primary")
    
def is_binary(star):
    return "star_2_mass" in star.cols

def get_dm(star):
    if is_binary(star):
        return get_dm_out_binary(star)
    else:
        return get_dm_single(star)

def get_dm_out_binary(star):
    # is the star donating
    # which star is primary
    s1 = "1"
    s2 = "2"
    if not is_primary(star):
        s1, s2 = s2, s1

    
    donating = (star.get("lg_mtransfer_rate") > -50) # good enough for here, is -99 if no transfer
    donating &= star.get(f"rl_relative_overflow_{s1}") > star.get(f"rl_relative_overflow_{s2}")
    dt = 10**star.get("log_dt")
    mdot_wind = 10**star.get(f"lg_wind_mdot_{s1}")
    mdot_mt = 10**star.get(f"lg_system_mdot_{s1}") + 10**star.get(f'lg_system_mdot_{s2}')
    mdot = mdot_wind  + mdot_mt * donating
    return mdot * dt

def get_dm_single(star):
    # is the star donating
    mdot = 10**star.get(f"log_abs_mdot")
    dt = 10**star.get("log_dt")

    return mdot * dt



def fractional_yield(star, iso, midpoint=False):
    if midpoint:
        return mass_net_yield_midpoint(star, iso) / star.get("star_mass")[0] 
    return  mass_net_yield(star, iso) / star.get("star_mass")[0] 


def isotopic_yields(star, midpoint=False):
    yields = {}
    for iso in ISOTOPES:
        y = fractional_yield(star, iso, midpoint=midpoint)
        yields[iso] = y

    return yields

def elemental_yields(star, midpoint=False):
    yields = {}
    for iso in ISOTOPES:
        ele = get_element(iso)
        y = fractional_yield(star, iso, midpoint=midpoint)
        if ele in yields.keys():
            yields[ele] += y
        else:
            yields[ele] = y

    return yields