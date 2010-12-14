"""Functions to make a nuclear data hdf5 file from the raw libraries."""
import os
import math

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

def get_xs_from_html_file(nucname, XS_Type_Flag, XS_Energy_Flag):
	with open("XShtml/{0}.html".format(nucname), 'r') as f:
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
				data = Convert2Barns(data, unit)
				return data
			elif inType and ("</ul>" in line):
				#XS not defined for this energy, returning zero
				return 0.0
	#If the specific XS was not found in trhis file, return zero
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

def make_decay(h5_file="nuc_data.h5", decay_file='decay.txt'):
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

def make_atomic_weight(h5_file="nuc_data.h5", data_file='atomic_weight.txt'):
    """Makes an atomic weight table and adds it to the hdf5 library.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_file (str): path to the atomic weight text file to load data from.
    """
    # Open the HDF5 File
    kdb = tb.openFile(h5_file, 'a')

    # Make a new the table
    Atable = kdb.createTable("/", "A", AtomicWeightDescription, "Atomic Weight Data [amu]")
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

def make_xs_1g(h5_file="nuc_data.h5", data_dir='xs_html/'):
    """Makes an atomic weight table and adds it to the hdf5 library.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_dir (str): path to that holds nuclide html one group neutron 
          cross section files.
    """
    # Get nulcide list from directory
    nuclist = [nuc.partition('.html')[0] for nuc in os.listdir(data_dir)]

    # Open the HDF5 File
    kdb = tb.openFile(h5_file, 'a')

    # Create Group
    xs_1g_group = kdb.createGroup("/", "xs_1g", "One Group Neutron Cross Section Data")

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
                nucrow[xstf] = get_xs_from_html_file(nuc, xs_1g_type_flags[xstf], xs_1g_energy_flags[xsef])

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

