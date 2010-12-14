"""Functions to make a nuclear data hdf5 file from the raw libraries."""

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

def get_xs_from_file(nucname, XS_Type_Flag, XS_Energy_Flag):
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
    'from_iso_LL': tb.StringCol(6, pos=0)
    'from_iso_zz': tb.Int32Col(pos=1)

    'half_life':   tb.Float64Col(pos=2)
    'decay_const': tb.Float64Col(pos=3)

    'to_iso_LL': tb.StringCol(6, pos=4)
    'to_iso_zz': tb.Int32Col(pos=5)

    'branch_ratio': tb.Float64Col(pos=6)	
    }

def make_decay(h5_file="decay.h5", decay_file='decay.txt'):
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
    'iso_LL': tb.StringCol(itemsize=6, pos=0)
    'iso_zz': tb.IntCol(pos=1)
    'value':  tb.FloatCol(pos=2)
    'error':  tb.FloatCol(pos=3)
    'abund':  tb.FloatCol(pos=4)
    }

def make_atomic_weight(h5_file="decay.h5", data_file='atomic_weight.txt'):
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

#######################
### Now for XS Data ###
#######################
XSGroup = kdb.createGroup("/", "XS", "Neutron Cross Section Data")

XS_Type_Flags = {"sigma_t": "Total Cross Section",
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
	"sigma_trit":   "(n,t) Cross Section"}

XS_Energy_Flags = {"Thermal": "at 0.0253 eV", 
	"ThermalMaxwellAve":  "Maxwell avg. at 0.0253 eV",
	"ResonanceIntegral":  "Resonance integral",
	"FourteenMeV":        "at 14 MeV",
	"FissionSpectrumAve": "Fission spectrum avg."}

class XSDescription(tables.IsDescription):
	isoLL        = tables.StringCol(itemsize=6, pos=1)
	isozz        = tables.IntCol(pos=2)
	sigma_t      = tables.FloatCol(pos=3)
	sigma_s      = tables.FloatCol(pos=4)
	sigma_e      = tables.FloatCol(pos=5)
	sigma_i      = tables.FloatCol(pos=6)
	sigma_a      = tables.FloatCol(pos=7)
	sigma_gamma  = tables.FloatCol(pos=8)
	sigma_f      = tables.FloatCol(pos=9)
	sigma_alpha  = tables.FloatCol(pos=10)
	sigma_proton = tables.FloatCol(pos=11)
	sigma_deut   = tables.FloatCol(pos=12)
	sigma_trit   = tables.FloatCol(pos=13)
	sigma_2n     = tables.FloatCol(pos=14)
	sigma_3n     = tables.FloatCol(pos=15)
	sigma_4n     = tables.FloatCol(pos=16)

for xsef in XS_Energy_Flags.keys():
	XStable = kdb.createTable(XSGroup, xsef, XSDescription, "({0}) [barns]".format(XS_Energy_Flags[xsef]))

	nucrow = XStable.row
	for  nuc in stablelist:
		nucrow['isoLL'] = nuc
		nucrow['isozz'] = isoname.LLAAAM_2_zzaaam(nuc)

		for xstf in XS_Type_Flags.keys():
			nucrow[xstf] = GetXSfromFile(nuc, XS_Type_Flags[xstf], XS_Energy_Flags[xsef])

		nucrow['sigma_s'] = nucrow['sigma_e'] + nucrow['sigma_i']
		nucrow['sigma_a'] = nucrow['sigma_gamma'] + nucrow['sigma_f'] + nucrow['sigma_alpha'] + nucrow['sigma_proton'] + \
			nucrow['sigma_deut'] + nucrow['sigma_trit'] + nucrow['sigma_2n'] + nucrow['sigma_3n'] + nucrow['sigma_4n']

		nucrow.append()

	XStable.flush()

#Clean-up
kdb.close()

