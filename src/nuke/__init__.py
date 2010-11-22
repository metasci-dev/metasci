"""MetaSci Nuclear Basics"""

# Prefer importing the system version of isoname
# over the package one.  This is because the system
# version is probably Python bidings of C++ code, 
# and thus faster.
try:
    import isoname
except:
    import IsoName as isoname


def cell_power(Material, SpecificPower, CellVolume, Density):
    """Calculates the power from a given unit cell.

    Args:
        * `Material` (dict): Normalized material dictionary; e.g. UOX, Metal (weight fraction, NOT atom fraction).
        * `SpecificPower` (float): Specific power of material [MW / kgIHM].
        * `CellVolume` (float): Unit cell volume that material occupies [cm^3].
        * `Density` (float): Density of material [g / cm^3].

    Returns:
        * `Power` (float): The power coming from this fuel cell [MW].
    """

    WeightFracIHM = 0.0
    for iso in Material.keys():
        if not ( isoname.zzLL[ isoname.mixed_2_zzaaam(iso)/10000 ] in isoname.ACT ):
            continue
        WeightFracIHM = WeightFracIHM + Material[iso]
    return SpecificPower * (10.0**-3) * Density * WeightFracIHM * CellVolume


def group_collapse(const_g, flux_g, flux = 0.0):
    """Performs a very simple group collapse.

    Args:
        * `const_g` (sequence): Constant as a function of energy group to collapse.
        * `flux_g` (sequence): Flux as a function of energy group.

    Keyword Args:
        * `flux` (float): Total flux, if known.  If not known, this will be calculated from flux_g.

    Returns:
        * `total` (float): One-group value of const.  Calculated via sum(flux_g[g] * const_g[g] / flux).
    """
    if flux == 0.0:
        for f_g in flux_g:
            flux = flux + f_g
    total = 0.0
    for g in range(len(const_g)):
        total = total + (flux_g[g] * const_g[g] / flux  )
    return total


def convolve_initial_fuel_form(stream, chem_form):
    """Convolves a heavy metal mass stream with a chemical form.  The chemical 
    form is a dictionary that should have an "IHM" key, whose value gets replaced
    with the fuel mass stream in the resultant isovec.

    Args:
        * `stream` (MassStream [see BriPy]): An initial heavy metal mass stream.
        * `chem_form` (dict): The chemical form of the fuel.  For instance UOX 
          would be represented by ``{80160: 2.0, "IHM": 1.0}`` while a metal fuel
          would simply be ``{"IHM": 1.0}``.

    Returns:
        * `isovec` (dict): Normalized, mass-weighted isotopic vector of the initial fuel form.
        * `AW` (float): Atomic weight of the fuel mass stream (ie, just the IHM).
        * `MW` (float): Molecular weight of the chemical form after the fuel stream has been 
          inserted into it.  This should always be greater than or equal to AW.
    """

    AW = stream.atomic_weight()

    # Calculate the MW
    if "IHM" not in chem_form:
        raise KeyError("The fuel chemical form must contain the key 'IHM'.")

    MW = 0.0
    for key in chem_form:
        if key == "IHM":
            MW += chem_form[key] * AW
        else:
            key_aw = isoname.nuc_weight(key)
            MW += chem_form[key] * key_aw

    AW_over_MW = AW / MW

    # Calculate the isotopic vector
    isovec = {}
    for key in chem_form:
        if key == "IHM":
            comp = stream.comp
            for iso in comp:
                isovec[iso] = comp[iso] * AW_over_MW
        else:
            key_zz = isoname.mixed_2_zzaaam(key)
            key_aw = isoname.nuc_weight_zzaaam(key_zz)
            isovec[key_zz] = chem_form[key] * key_aw / MW

    return isovec, AW, MW


######################
### MCNP Functions ###
######################

class FailedToParseMCNP(Exception):
    """MCNP parsing error."""
    def __init__(self, line):
        """Initialized with the line at which the  parse attempt failed."""
        self.line = line
        return

    def __str__(self):
        """Useful exception output."""
        return  "Failed to parse the following MCNP line:\n>>> " + self.line

def Line2MCNP(s):
    """Takes a long string and returns an MCNP-valid input card. 
    For cards that are over 80 characters in length, 
    this function returns adds a newline character at most every 80 characters.  
    It also ensures that the card is continued on the next line by using proper indentation.  
    Do not use for multiple cards at once.

    Args:
        * `s` (str): An MCNP card that may be too long.

    Returns:
        * `news` (str): A new MCNP card that has been appropriately broken up into 
          separate lines.
    """

    dist2newline = len(s)		#Distance to the last newline in the string from the back
    while 80 < dist2newline:
        n = 80 - dist2newline
        NotBeenSplit = True
        while NotBeenSplit:
            if s[n] in [" ", ","]:
                s = s[:n] + "\n     " + s[n+1:]
                NotBeenSplit = False
            elif s[n] == "\n":
                print "Reached a newline character without finding a splittable character."
                raise FailedToParseMCNP(s)
            else:
                n = n - 1
        dist2newline = len(s) - s.rfind("\n")
    return s

