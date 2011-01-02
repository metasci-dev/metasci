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
    """A lightweight cross-section cache based off of python dictionaries."""

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


##############################
### Partial group collapse ###
##############################

def phi_g(E_g, E_n, phi_n):
    """Calculates the lower resolution flux, phi_g, from the lower resolution group stucture E_g, 
    the higher resolution groups E_n, and the higher resolution flux phi_n.

    Args:
        * E_g (sequence of floats): lower resolution energy group structure [MeV]
          that is of length G+1. Ordered from lowest-to-highest energy.
        * E_n (sequence of floats): higher resolution energy group structure [MeV]
          that is of length N+1. Ordered from lowest-to-highest energy.
        * phi_n (sequence of floats): The high-fidelity flux [n/cm^2/s] to collapse the fission 
          cross-section over.  Length N.  Ordered from lowest-to-highest energy.

    Returns: 
        * phi_g_array (numpy array of floats): The flux collapsed to G energy groups.
    """
    # Some convienence paramters
    G = len(E_g) - 1
    N = len(E_n) - 1
    assert N == len(phi_n)

    index_E = np.arange(N+1)

    # Setup masks
    inner_mask_E = np.array([(E_g[g] <= E_n) & (E_n <= E_g[g+1]) for g in range(G)])
    inner_mask_phi = inner_mask_E[:, :-1]

    # si
    print inner_mask_E.sum(axis=0)

    phi_g = np.array([phi_n[inner_mask_phi[g]].sum() for g in range(G)])
    print phi_g 

    lower_index = np.array([index_E[inner_mask_E[g]][0] for g in range(G)])
    upper_index = np.array([index_E[inner_mask_E[g]][-1] for g in range(G)])
    print lower_index
    print upper_index

#def partial_group_collapse()

################################
### Cross-section generators ###
################################

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
