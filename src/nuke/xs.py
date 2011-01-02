"""This module provides functionality to automatically extract cross-sections
from the nulear database."""

import numpy as np
import tables as tb

import isoname

from . import nuc_data

###############################################################################
### Set up a cross-section cache so the same data isn't loaded repetitively ###
###############################################################################
    
class XSCache(dict):

    def __getitem__(self, key):
        """Overwrites the dict's key lookup by prodiving some custom loading
        from the nuc_data database file."""

        # Grab the energy group structure of the file
        if (key == 'E_n') and ('E_n' not in self):
            with tb.openFile(nuc_data, 'r') as f:
                self['E_n'] = np.array(f.root.neutron.xs_mg.E_g)

        # Return the value requestion
        return super(XSCache, self).__getitem__(key)

xs_cache = XSCache()

def sigma_f(iso, E_g, phi_n):
    """Calculates the neutron fission cross-section for an isotope for a new, lower resolution
    group structure using a higher fidelity flux.  Note that g indexes G, n indexes N, and G < N.

    Args:
        * iso (int or str): An isotope to calculate the fission cross-section for.
        * E_g (sequence of floats): New, lower fidelity energy group structure [MeV]
          that is of length G+1. Ordered from lowest-to-highest energy.
        * phi_n (sequence of floats): The high-fidelity flux [n/cm^2/s] to collapse the fission 
          cross-section over.  Length N.  Ordered from lowest-to-highest energy.

    Returns:
        * sigma_f_g (numpy array): A numpy array of the collapsed fission cross-section.
    """
    pass
