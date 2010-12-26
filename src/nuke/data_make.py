"""Functions to make a nuclear data hdf5 file from the raw libraries."""
import os
import re
import math

import numpy as np
import tables as tb

import isoname

########################
### Helper Functions ###
########################

def convert_to_barns(xs, unit):
    unit = unit.lower()
	
    if unit == "b":
        return xs
    elif unit == "mb":
        return xs * (10.0**-3)
    elif unit == "microbarn":
        return xs * (10.0**-6)
    else:
        print("Units {0} could not be converted".format(unit))
        return

def get_xs_from_html_file(nucname, data_dir, XS_Type_Flag, XS_Energy_Flag):
    with open("{0}{1}.html".format(data_dir, nucname), 'r') as f:
        inType = False
        for line in f:
            if XS_Type_Flag in line:
                inType = True

            if inType and ("<li>"+XS_Energy_Flag in line):
                du = line.partition("=")[2].split()
                data = float(du.pop(0))
                unit = ""
                for u in du:
                    unit = unit + u
                unit = unit.partition("\\")[0]
                data = convert_to_barns(data, unit)
                return data

            elif inType and ("</ul>" in line):
                # XS not defined for this energy, returning zero
                return 0.0

    # If the specific XS was not found in trhis file, return zero
    return 0.0


##################
### Decay Data ###
##################

# Decay isotopic description
decay_iso_desc = {
    'from_iso_LL': tb.StringCol(6, pos=0),
    'from_iso_zz': tb.Int32Col(pos=1),

    'half_life':   tb.Float64Col(pos=2),
    'decay_const': tb.Float64Col(pos=3),

    'to_iso_LL': tb.StringCol(6, pos=4),
    'to_iso_zz': tb.Int32Col(pos=5),

    'branch_ratio': tb.Float64Col(pos=6),
    }

def make_decay(h5_file='nuc_data.h5', decay_file='decay.txt'):
    """Makes a decay table and adds it to the hdf5 library.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * decay_file (str): path to the decay text file to load data from.
    """
    # Open the hdf5 file
    decayfile = tb.openFile(h5_file, "a")

    # Get the HDF5 root group
    root = decayfile.root

    # Initialize decay table
    decaytbl = decayfile.createTable(root, "decay", decay_iso_desc, "Isotopic Deacy Information")

    # Now, fill the table:
    row = decaytbl.row

    # Read the text file
    with open(decay_file, 'r') as lib:
        for line in lib:
            ls = line.split()

            from_iso_LL = isoname.mixed_2_LLAAAM(ls[0])
            from_iso_zz = isoname.LLAAAM_2_zzaaam(from_iso_LL)            
       	    row['from_iso_LL'] = from_iso_LL
       	    row['from_iso_zz'] = from_iso_zz

	        # Set halflife
            if ls[1] == "inf":
    	        hl = 10.0**300
            else:
                hl = float(ls[1])
            row['half_life'] = hl

        	# Set decay constants
            row['decay_const'] = math.log(2.0) / hl

            # Set to iso
            if ls[2] == "None":
                to_iso_LL = from_iso_LL
                to_iso_zz = from_iso_zz
    	    else:
                to_iso_LL = isoname.mixed_2_LLAAAM(ls[2])
                to_iso_zz = isoname.LLAAAM_2_zzaaam(to_iso_LL)

            row['to_iso_LL'] = to_iso_LL
            row['to_iso_zz'] = to_iso_zz

            # Set branch ratio
            row['branch_ratio'] = float(ls[3])

            # This injects the Record values
            row.append()

        # Flush the table buffer
        decaytbl.flush()

    # Finally, close the HDF5 file 
    decayfile.close()


############################
### Next, Atomic Weights ###
############################

atomic_weight_desc = {
    'iso_LL': tb.StringCol(itemsize=6, pos=0),
    'iso_zz': tb.IntCol(pos=1),
    'value':  tb.FloatCol(pos=2),
    'error':  tb.FloatCol(pos=3),
    'abund':  tb.FloatCol(pos=4),
    }

def make_atomic_weight(h5_file='nuc_data.h5', data_file='atomic_weight.txt'):
    """Makes an atomic weight table and adds it to the hdf5 library.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_file (str): path to the atomic weight text file to load data from.
    """
    # Open the HDF5 File
    kdb = tb.openFile(h5_file, 'a')

    # Make a new the table
    Atable = kdb.createTable("/", "A", atomic_weight_desc, "Atomic Weight Data [amu]")
    nuc = Atable.row

    with open(data_file, 'r') as f:
        for line in f:
            ls = line.split()
            iso_LL = isoname.mixed_2_LLAAAM(ls[0])
            iso_zz = isoname.LLAAAM_2_zzaaam(iso_LL)

            nuc['iso_LL'] = iso_LL
            nuc['iso_zz'] = iso_zz
            nuc['value'] = float(ls[1])
            nuc['error'] = float(ls[2])
            nuc['abund'] = float(ls[3])

            # Insert nuclide to table
            nuc.append()

    # Ensure that data was written to table
    Atable.flush()

    # Close the hdf5 file
    kdb.close()


#################################
### Now for One Group XS Data ###
#################################

xs_1g_type_flags = {
    "sigma_t":      "Total Cross Section",
    "sigma_e":      "Elastic Scattering Cross Section",
    "sigma_i":      "Total Inelastic Cross Section",
    "sigma_2n":     "(n,2n) Cross Section",
    "sigma_3n":     "(n,3n) Cross Section",
    "sigma_4n":     "(n,4n) Cross Section",
    "sigma_f":      "Total Fission Cross Section",
    "sigma_gamma":  "Radiative Capture Cross Section",
    "sigma_alpha":  "(n,alpha) Cross Section",
    "sigma_proton": "(n,p) Cross Section",
    "sigma_deut":   "(n,d) Cross Section",
    "sigma_trit":   "(n,t) Cross Section",
    }

xs_1g_energy_flags = {
    "Thermal":            "at 0.0253 eV", 
    "ThermalMaxwellAve":  "Maxwell avg. at 0.0253 eV",
    "ResonanceIntegral":  "Resonance integral",
    "FourteenMeV":        "at 14 MeV",
    "FissionSpectrumAve": "Fission spectrum avg.",
    }

xs_1g_desc = {
    'iso_LL':       tb.StringCol(itemsize=6, pos=1),
    'iso_zz':       tb.IntCol(pos=2),
    'sigma_t':      tb.FloatCol(pos=3),
    'sigma_s':      tb.FloatCol(pos=4),
    'sigma_e':      tb.FloatCol(pos=5),
    'sigma_i':      tb.FloatCol(pos=6),
    'sigma_a':      tb.FloatCol(pos=7),
    'sigma_gamma':  tb.FloatCol(pos=8),
    'sigma_f':      tb.FloatCol(pos=9),
    'sigma_alpha':  tb.FloatCol(pos=10),
    'sigma_proton': tb.FloatCol(pos=11),
    'sigma_deut':   tb.FloatCol(pos=12),
    'sigma_trit':   tb.FloatCol(pos=13),
    'sigma_2n':     tb.FloatCol(pos=14),
    'sigma_3n':     tb.FloatCol(pos=15),
    'sigma_4n':     tb.FloatCol(pos=16),
    }

def make_xs_1g(h5_file='nuc_data.h5', data_dir='xs_html/'):
    """Makes an atomic weight table and adds it to the hdf5 library.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_dir (str): path to that holds nuclide html one group neutron 
          cross section files.
    """
    # Get nulcide list from directory
    nuclist = [nuc.partition('.html')[0] for nuc in os.listdir(data_dir)]
    nuclist_zz = isoname.LLAAAM_2_zzaaam_List(nuclist)
    nuclist_zz.sort()
    nuclist = isoname.zzaaam_2_LLAAAM_List(nuclist_zz)

    # Open the HDF5 File
    kdb = tb.openFile(h5_file, 'a')

    # Create neutron group
    if not hasattr(kdb.root, 'neutron'):
        neutron_group = kdb.createGroup('/', 'neutron', 'Neutron Cross Sections')

    # Create xs_1g Group
    xs_1g_group = kdb.createGroup("/neutron", "xs_1g", "One Group Neutron Cross Section Data")

    # Loop through all energy types
    for xsef in xs_1g_energy_flags:
        # Create table for this energy 
        xs_1g_table = kdb.createTable(xs_1g_group, xsef, xs_1g_desc, "({0}) [barns]".format(xs_1g_energy_flags[xsef]))
        nucrow = xs_1g_table.row

        for nuc in nuclist:
            iso_LL = isoname.mixed_2_LLAAAM(nuc)
            iso_zz = isoname.LLAAAM_2_zzaaam(iso_LL)

            nucrow['iso_LL'] = iso_LL
            nucrow['iso_zz'] = iso_zz

            for xstf in xs_1g_type_flags:
                nucrow[xstf] = get_xs_from_html_file(nuc, data_dir, xs_1g_type_flags[xstf], xs_1g_energy_flags[xsef])

            nucrow['sigma_s'] = nucrow['sigma_e'] + nucrow['sigma_i']

            nucrow['sigma_a'] = nucrow['sigma_gamma'] + nucrow['sigma_f'] + nucrow['sigma_alpha'] + \
                                nucrow['sigma_proton'] + nucrow['sigma_deut'] + nucrow['sigma_trit'] + \
                                nucrow['sigma_2n'] + nucrow['sigma_3n'] + nucrow['sigma_4n']

            # Write out this row 
            nucrow.append()

        # Writr out this table
        xs_1g_table.flush()

    # Close the hdf5 file
    kdb.close()


###################################
### Now for Multi-Group XS Data ###
###################################

# These read in cinder.dat
cinder_float = "[\d.+-Ee]+"

def _init_multigroup(kdb):
    """Initializes a multigroup cross-section part of the database.

    Keyword Args:
        * kdb (tables.File): a nuclear data hdf5 file.
    """

    # Create neutron and photon groups
    if not hasattr(kdb.root, 'neutron'):
        neutron_group = kdb.createGroup('/', 'neutron', 'Neutron Cross Sections')

    if not hasattr(kdb.root, 'photon'):
        photon_group = kdb.createGroup('/', 'photon', 'Photon Cross Sections')

    # Create xs_mg groups
    if not hasattr(kdb.root.neutron, 'xs_mg'):
        nxs_mg_group = kdb.createGroup("/neutron", "xs_mg", "Multi-Group Neutron Cross Section Data")

    # Create source groups
    if not hasattr(kdb.root.photon, 'source'):
        gxs_mg_group = kdb.createGroup("/photon", "source", "Multi-Group Photon Source Data")

    # Create fission_yield groups
    if not hasattr(kdb.root.neutron, 'fission_products'):
        nxs_mg_group = kdb.createGroup("/neutron", "fission_products", "Neutron Fission Product Yield Data")

    if not hasattr(kdb.root.photon, 'fission_products'):
        nxs_mg_group = kdb.createGroup("/photon", "fission_products", "Photofission Product Yield Data")



def _get_groups_sizes(raw_data):
    """Gets the number of nuclides and groups in this data file.

    Args:
        * data (str): Input cinder.dat data file as a string.

    Returns:
        * nuclides (int): the number of nuclides in the dataset
        * G_n (int): the number of neutron energy groups in the dataset
        * G_p (int): the number of proton energy groups in the dataset
        * G_g (int): the number of photon energy groups in the dataset
    """
    # Search for the group pattern
    G_pattern = "(\d+) nuclides,\s+(\d+) neutron groups,\s+(\d+) proton groups,\s+(\d+) photon groups"
    m = re.search(G_pattern, raw_data)
    g = m.groups()

    # Convert to ints
    nuclides = int(g[0])
    G_n = int(g[1])
    G_p = int(g[2])
    G_g = int(g[3])

    return nuclides, G_n, G_p, G_g


def cinder_2_zzaaam(iso_c):
    """Converts an isotope from cinder form to zzaaam form.

    Args:
        * iso_c (str): isotope in cinder form.

    Returns:
        * iso_zz (int): isotope in zzaaam form.
    """
    m = int(iso_c[-1])
    zz = int(iso_c[-4:-1])
    aaa = int(iso_c[:-4])

    iso_zz = zz*10000 + aaa*10 + m
    return iso_zz


def make_mg_group_structure(h5_file='nuc_data.h5', data_file='cinder.dat'):
    """Add the energy group bounds arrays to the hdf5 library.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_file (str): path to the cinder.dat data file.
    """
    # Open the HDF5 File
    kdb = tb.openFile(h5_file, 'a')

    # Ensure that the appropriate file structure is present
    _init_multigroup(kdb)

    # Read in cinder data file
    with open(data_file, 'r') as f:
        raw_data = f.read()

    # Get group sizes
    nuclides, G_n, G_p, G_g = _get_groups_sizes(raw_data)

    # Find & write neutron group structure
    n_E_g_pattern = "Neutron group .*, MeV" + ("\s+("+cinder_float+")")*(G_n + 1)
    m = re.search(n_E_g_pattern, raw_data)
    g = m.groups()
    n_E_g = np.array(g, dtype=float)
    kdb.createArray('/neutron/xs_mg', 'E_g', n_E_g, 'Neutron energy group bounds [MeV]')

    # Find & write photon group structure
    g_E_g_pattern = "Gamma structure, MeV" + ("\s+("+cinder_float+")")*(G_g + 1)
    m = re.search(g_E_g_pattern, raw_data)
    g = m.groups()
    g_E_g = np.array(g, dtype=float)
    kdb.createArray('/photon/source', 'E_g', g_E_g, 'Photon energy group bounds [MeV]')

    # Close the hdf5 file
    kdb.close()

# Helpful patterns
from_iso_pattern = "\n#[\s\d]{4}:\s+(\d+).*?\n(_______________________| [\w/-]+ Fission Yield Data)"
to_iso_base = "#[\s\d]{4}:\s+(\d+) produced by the following C-X  \((.{4})\) REF:.*?\n"

absorption_desc = {
    'from_iso_LL': tb.StringCol(6, pos=0),
    'from_iso_zz': tb.Int32Col(pos=1),

    'to_iso_LL': tb.StringCol(6, pos=2),
    'to_iso_zz': tb.Int32Col(pos=3),

    'reaction_type': tb.StringCol(4, pos=4),
    'xs': None, # Should be replaced with tb.Float64Col(shape=(G_n, ), pos=5),    
    }

def make_mg_absorption(h5_file='nuc_data.h5', data_file='cinder.dat'):
    """Adds the absorption reaction rate cross sections to the hdf5 library.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_file (str): path to the cinder.dat data file.
    """
    # Open the HDF5 File
    kdb = tb.openFile(h5_file, 'a')

    # Ensure that the appropriate file structure is present
    _init_multigroup(kdb)

    # Read in cinder data file
    with open(data_file, 'r') as f:
        raw_data = f.read()

    # Get group sizes
    nuclides, G_n, G_p, G_g = _get_groups_sizes(raw_data)

    # Init the neutron absorption table
    absorption_desc['xs'] = tb.Float64Col(shape=(G_n, ), pos=5)
    absorption_table = kdb.createTable('/neutron/xs_mg/', 'absorption', absorption_desc, 
                                       'Neutron absorption reaction rate cross sections [barns]')
    abrow = absorption_table.row

    # Init to_iso_pattern
    to_iso_pattern = to_iso_base + ("\s+("+cinder_float+")")*G_n

    # Iterate through all from isotopes.
    for m_from in re.finditer(from_iso_pattern, raw_data, re.DOTALL):
        from_iso_zz = cinder_2_zzaaam(m_from.group(1))

        # Check matestable state
        if 1 < from_iso_zz%10:
            # Metastable state too high!
            continue
        from_iso_LL = isoname.zzaaam_2_LLAAAM(from_iso_zz)

        # Grab the string for this from_iso in order to get all of the to_isos
        from_iso_part = m_from.group(0)

        # Iterate over all to_isos
        for m_to in re.finditer(to_iso_pattern, from_iso_part):
            to_iso_zz = cinder_2_zzaaam(m_to.group(1))

            # Check matestable state
            if 1 < to_iso_zz%10:
                # Metastable state too high!
                continue
            to_iso_LL = isoname.zzaaam_2_LLAAAM(to_iso_zz)

            # Munge reaction type
            rx_type = m_to.group(2)
            rx_type = rx_type.strip()

            # Setup XS array
            xs = np.array(m_to.groups()[2:], dtype=float)
            assert xs.shape == (G_n, )

            # Write this row to the absorption table
            abrow['from_iso_LL'] = from_iso_LL
            abrow['from_iso_zz'] = from_iso_zz

            abrow['to_iso_LL'] = to_iso_LL
            abrow['to_iso_zz'] = to_iso_zz

            abrow['reaction_type'] = rx_type
            abrow['xs'] = xs

            abrow.append()

        # Flush this from iso
        absorption_table.flush()

    # Close the hdf5 file
    kdb.close()


fission_base = "\n([\s\d]+),([\s\d]+),([\s\d]+)=\(n,f\) yield sets\. If > 0,[\s\d]{4}-gp fisn CX follows\..*?\n"

fission_desc = {
    'iso_LL': tb.StringCol(6, pos=0),
    'iso_zz': tb.Int32Col(pos=1),

    'thermal_yield':     tb.Int8Col(pos=2),
    'fast_yield':        tb.Int8Col(pos=3),
    'high_energy_yield': tb.Int8Col(pos=4),

    'xs': None, # Should be replaced with tb.Float64Col(shape=(G_n, ), pos=5),    
    }

def make_mg_fission(h5_file='nuc_data.h5', data_file='cinder.dat'):
    """Adds the fission reaction rate cross sections to the hdf5 library.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_file (str): path to the cinder.dat data file.
    """
    # Open the HDF5 File
    kdb = tb.openFile(h5_file, 'a')

    # Ensure that the appropriate file structure is present
    _init_multigroup(kdb)

    # Read in cinder data file
    with open(data_file, 'r') as f:
        raw_data = f.read()

    # Get group sizes
    nuclides, G_n, G_p, G_g = _get_groups_sizes(raw_data)

    # Init the neutron absorption table
    fission_desc['xs'] = tb.Float64Col(shape=(G_n, ), pos=5)
    fission_table = kdb.createTable('/neutron/xs_mg/', 'fission', fission_desc, 
                                    'Neutron fission reaction rate cross sections [barns]')
    frow = fission_table.row

    # Init to_iso_pattern
    fission_pattern = fission_base + ("\s+("+cinder_float+")")*G_n

    # Iterate through all from isotopes.
    for m_from in re.finditer(from_iso_pattern, raw_data, re.DOTALL):
        from_iso_zz = cinder_2_zzaaam(m_from.group(1))

        # Check matestable state
        if 1 < from_iso_zz%10:
            # Metastable state too high!
            continue
        from_iso_LL = isoname.zzaaam_2_LLAAAM(from_iso_zz)

        # Grab the string for this from_iso in order to get all of the to_isos
        from_iso_part = m_from.group(0)

        # Grab the fission part
        m_fission = re.search(fission_pattern, from_iso_part)
        if m_fission is None:
            continue

        # Grab yield indexes
        yield_t = int(m_fission.group(1))
        yield_f = int(m_fission.group(2))
        yield_h = int(m_fission.group(3))

        # Grab XS array
        xs = np.array(m_fission.groups()[3:], dtype=float)
        assert xs.shape == (G_n, )

        # Write fission table row
        frow['iso_LL'] = from_iso_LL
        frow['iso_zz'] = from_iso_zz

        frow['thermal_yield']     = yield_t
        frow['fast_yield']        = yield_f
        frow['high_energy_yield'] = yield_h

        frow['xs'] = xs

        # Write out this row
        frow.append()
        fission_table.flush()

    # Close the hdf5 file
    kdb.close()


gamma_decay_base = "gamma spectra from .{3} multiplied by\s+(" + cinder_float + ")\s+to agree w/ Eg=\s*(" +\
                   cinder_float + ")\n"

gamma_decay_desc = {
    'iso_LL': tb.StringCol(6, pos=0),
    'iso_zz': tb.Int32Col(pos=1),

    'energy': tb.Float64Col(pos=2),
    'scaling_factor': tb.Float64Col(pos=3),

    'spectrum': None, # Should be replaced with tb.Float64Col(shape=(G_g, ), pos=4),
    }

def make_mg_gamma_decay(h5_file='nuc_data.h5', data_file='cinder.dat'):
    """Adds the gamma decay spectrum information to the hdf5 library.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_file (str): path to the cinder.dat data file.
    """
    # Open the HDF5 File
    kdb = tb.openFile(h5_file, 'a')

    # Ensure that the appropriate file structure is present
    _init_multigroup(kdb)

    # Read in cinder data file
    with open(data_file, 'r') as f:
        raw_data = f.read()

    # Get group sizes
    nuclides, G_n, G_p, G_g = _get_groups_sizes(raw_data)

    # Init the gamma absorption table
    gamma_decay_desc['spectrum'] = tb.Float64Col(shape=(G_g, ), pos=4)
    gamma_decay_table = kdb.createTable('/photon/source/', 'decay_spectra', gamma_decay_desc, 
                                        'Gamma decay spectrum [MeV]')
    gdrow = gamma_decay_table.row

    # Init to_iso_pattern
    gamma_decay_pattern = gamma_decay_base + ("\s+("+cinder_float+")")*G_g

    # Iterate through all from isotopes.
    for m_from in re.finditer(from_iso_pattern, raw_data, re.DOTALL):
        from_iso_zz = cinder_2_zzaaam(m_from.group(1))

        # Check matestable state
        if 1 < from_iso_zz%10:
            # Metastable state too high!
            continue
        from_iso_LL = isoname.zzaaam_2_LLAAAM(from_iso_zz)

        # Grab the string for this from_iso in order to get all of the to_isos
        from_iso_part = m_from.group(0)

        # Grab the fission part
        m_gd = re.search(gamma_decay_pattern, from_iso_part)
        if m_gd is None:
            continue

        # Grab base data
        scale = float(m_gd.group(1))
        energy = float(m_gd.group(2))

        # Grab spectrum
        spectrum = np.array(m_gd.groups()[2:], dtype=float)
        assert spectrum.shape == (G_g, )

        # Prepare the row
        gdrow['iso_LL'] = from_iso_LL
        gdrow['iso_zz'] = from_iso_zz

        gdrow['energy'] = energy
        gdrow['scaling_factor'] = scale

        gdrow['spectrum'] = spectrum

        # Write out the row
        gdrow.append()
        gamma_decay_table.flush()

    # Close the hdf5 file
    kdb.close()


def make_xs_mg(h5_file='nuc_data.h5', data_file='cinder.dat'):
    """Adds multi-group cross-section data to the hdf5 library from CINDER.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_file (str): path to the cinder.dat data file.
    """

    # Add energy groups to file
    make_mg_group_structure(h5_file='nuc_data.h5', data_file='cinder.dat')

    # Add neutron absorption to file
    make_mg_absorption(h5_file='nuc_data.h5', data_file='cinder.dat')

    # Add fission to file
    make_mg_fission(h5_file='nuc_data.h5', data_file='cinder.dat')

    # Add gamma decay spectrum to file
    make_mg_gamma_decay(h5_file='nuc_data.h5', data_file='cinder.dat')



#######################################
### Make Fission Product Yield Data ###
#######################################


def _init_fission_products(kdb):
    """Initializes the fission product part of the database.

    Keyword Args:
        * kdb (tables.File): a nuclear data hdf5 file.
    """

    # Create neutron and photon groups
    if not hasattr(kdb.root, 'neutron'):
        neutron_group = kdb.createGroup('/', 'neutron', 'Neutron Cross Sections')

    if not hasattr(kdb.root, 'photon'):
        photon_group = kdb.createGroup('/', 'photon', 'Photon Cross Sections')

    # Create fission_yield groups
    if not hasattr(kdb.root.neutron, 'fission_products'):
        nfp_group = kdb.createGroup("/neutron", "fission_products", "Neutron Fission Product Yield Data")

    if not hasattr(kdb.root.photon, 'fission_products'):
        gfp_group = kdb.createGroup("/photon", "fission_products", "Photofission Product Yield Data")


def _get_fp_sizes(raw_data):
    """Gets the number of fission product yield data sets in this file.

    Args:
        * data (str): Input cinder.dat data file as a string.

    Returns:
        * N_n (int): the number of neutron fission product yield datasets in the file
        * N_g (int): the number of photon fission product yield datasets in the file
    """
    # Search for the neutron pattern
    N_n_pattern = "Fission Yield Data.*?\n\s*(\d+)\s+yield sets"
    m_n = re.search(N_n_pattern, raw_data)
    N_n = int(m_n.group(1))

    # Search for the photon pattern
    N_g_pattern = "Photofission Yield Data.*?\n\s*(\d+)\s+yield sets"
    m_g = re.search(N_g_pattern, raw_data)
    N_g = int(m_g.group(1))

    return N_n, N_g


fp_info_desc = {
    'index': tb.Int16Col(pos=0),

    'iso_LL': tb.StringCol(6, pos=1),
    'iso_zz': tb.Int32Col(pos=2),

    'type': tb.StringCol(11, pos=3),
    'mass': tb.Float64Col(pos=4),
    }

fp_type_flag = {
    't': 'thermal', 
    'f': 'fast',
    'h': 'high_energy', 
    's': 'spontaneous', 
    }

iit_pattern = "(\d{1,3})\s+\d{2,3}-\s?([A-Z]{1,2}-[ Mm\d]{3})([tfhs])"
mass_pattern = "\d{1,3}\.\d{1,4}"

nfp_info_pattern = "Fission Yield Data.*?fission products"

def make_neutron_fp_info(h5_file='nuc_data.h5', data_file='cinder.dat'):
    """Adds the neutron fission product yiled info to the hdf5 library.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_file (str): path to the cinder.dat data file.
    """
    # Open the HDF5 File
    kdb = tb.openFile(h5_file, 'a')

    # Ensure that the appropriate file structure is present
    _init_fission_products(kdb)

    # Read in cinder data file
    with open(data_file, 'r') as f:
        raw_data = f.read()

    # Get group sizes
    N_n, N_g = _get_fp_sizes(raw_data)

    # Grab the part of the file that is a neutron fission product yiled info
    m_info = re.search(nfp_info_pattern, raw_data, re.DOTALL)
    nfp_info_raw = m_info.group(0)

    # Grab the index, isotope, and type
    iits = re.findall(iit_pattern, nfp_info_raw)

    # Grab the masses 
    masses = re.findall(mass_pattern, nfp_info_raw)

    # Make sure data is the right size
    assert N_n == len(iits) 
    assert N_n == len(masses)

    # Init the neutron fission product info table
    nfp_table = kdb.createTable('/neutron/fission_products/', 'info', fp_info_desc, 
                                'Neutron Fission Product Yield Information')
    nfprow = nfp_table.row

    # Write out info table rows 
    for m in range(N_n):
        iit = iits[m]
        index = int(iit[0])

        iso_zz = isoname.LLAAAM_2_zzaaam(iit[1])
        # Correct for metastable flag
        if 0 != iso_zz%10:
            iso_zz = iso_zz + 2000

        iso_LL = isoname.zzaaam_2_LLAAAM(iso_zz)
        type = fp_type_flag[iit[2]]
        mass = float(masses[m])

        # Prep row
        nfprow['index'] = index
        nfprow['iso_LL'] = iso_LL
        nfprow['iso_zz'] = iso_zz
        nfprow['type'] = type
        nfprow['mass'] = mass

        # Write out row
        nfprow.append()
        nfp_table.flush()

    # Close the hdf5 file
    kdb.close()


gfp_info_pattern = "Photofission Yield Data.*?fission products"

def make_photon_fp_info(h5_file='nuc_data.h5', data_file='cinder.dat'):
    """Adds the photofission product yiled info to the hdf5 library.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_file (str): path to the cinder.dat data file.
    """
    # Open the HDF5 File
    kdb = tb.openFile(h5_file, 'a')

    # Ensure that the appropriate file structure is present
    _init_fission_products(kdb)

    # Read in cinder data file
    with open(data_file, 'r') as f:
        raw_data = f.read()

    # Get group sizes
    N_n, N_g = _get_fp_sizes(raw_data)

    # Grab the part of the file that is a neutron fission product yiled info
    m_info = re.search(gfp_info_pattern, raw_data, re.DOTALL)
    gfp_info_raw = m_info.group(0)

    # Grab the index, isotope, and type
    iits = re.findall(iit_pattern, gfp_info_raw)

    # Grab the masses 
    masses = re.findall(mass_pattern, gfp_info_raw)

    # Make sure data is the right size
    assert N_g == len(iits) 
    assert N_g == len(masses)

    # Init the neutron fission product info table
    gfp_table = kdb.createTable('/photon/fission_products/', 'info', fp_info_desc, 
                                'Photofission Product Yield Information')
    gfprow = gfp_table.row

    # Write out info table rows 
    for m in range(N_g):
        iit = iits[m]
        index = int(iit[0])

        iso_zz = isoname.LLAAAM_2_zzaaam(iit[1])
        # Correct for metastable flag
        if 0 != iso_zz%10:
            iso_zz = iso_zz + 2000

        iso_LL = isoname.zzaaam_2_LLAAAM(iso_zz)
        type = fp_type_flag[iit[2]]
        mass = float(masses[m])

        # Prep row
        gfprow['index'] = index
        gfprow['iso_LL'] = iso_LL
        gfprow['iso_zz'] = iso_zz
        gfprow['type'] = type
        gfprow['mass'] = mass

        # Write out row
        gfprow.append()
        gfp_table.flush()

    # Close the hdf5 file
    kdb.close()


def make_fission_products(h5_file='nuc_data.h5', data_file='cinder.dat'):
    """Adds fission-product data to the hdf5 library from CINDER.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_file (str): path to the cinder.dat data file.
    """
    # Add neutro info table
    make_neutron_fp_info(h5_file='nuc_data.h5', data_file='cinder.dat')

    # Add neutro info table
    make_photon_fp_info(h5_file='nuc_data.h5', data_file='cinder.dat')



######################################
### Make the data base as a script ###
######################################

if __name__ == "__main__":
    # Clean existing file
    if 'nuc_data.h5' in os.listdir('.'):
        os.remove('nuc_data.h5')

    # Make atomic weights
    make_atomic_weight(h5_file='nuc_data.h5', data_file='atomic_weight.txt')

    # Make decay table
    make_decay(h5_file='nuc_data.h5', decay_file='decay.txt')

    # Make one group xs library
#    make_xs_1g(h5_file='nuc_data.h5', data_dir='xs_html/')

    # Make multi-group xs library
#    make_xs_mg(h5_file='nuc_data.h5', data_file='cinder.dat')

    # Make fission product yield library
    make_fission_products(h5_file='nuc_data.h5', data_file='cinder.dat')
