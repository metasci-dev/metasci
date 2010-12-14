"""Functions to make a nuclear data hdf5 file from the raw libraries."""

import math

import tables as tb

import isoname

# Decay isotopic description
decay_iso_desc = {
    'from_iso_LL': StringCol(6, pos=0)
    'from_iso_zz': Int32Col(pos=1)

    'half_life':    Float64Col(pos=2)
    'decay_const': Float64Col(pos=3)

    'to_iso_LL': StringCol(6, pos=4)
    'to_iso_zz': Int32Col(pos=5)

    'branch_ratio': Float64Col(pos=6)	
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
    decaytbl = decayfile.createTable(root, "decay", decay_iso_desc)

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



