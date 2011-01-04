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

def partial_energy_matrix(E_g, E_n):
    """Gerenates a matrix of fractional values that may be used to converts a high-resolution 
    flux array with group structure E_n to a low-resolution flux array with group-structure E_g.

    Args:
        * E_n (sequence of floats): higher resolution energy group structure [MeV]
          that is of length N+1. Ordered from lowest-to-highest energy.
        * E_g (sequence of floats): lower resolution energy group structure [MeV]
          that is of length G+1. Ordered from lowest-to-highest energy.

    Returns:
        * pem (array of fractions): This is a GxN sized matrix that when dotted with a 
          high-resolution flux produces a low-resolution flux.
    """
    # Some convienence paramters
    G = len(E_g) - 1
    N = len(E_n) - 1

    index_E_n = np.arange(N+1)

    # Some assertions to ensure that everything is well formed
    assert E_n[0] <= E_g[0]
    assert E_g[-1] <= E_n[-1]

    # Get the interior points for each gth group in n-space
    inner_mask = np.array([(E_g[g] <= E_n) & (E_n <= E_g[g+1]) for g in range(G)])

    # Get the upper and lower nth index for every gth group
    lower_index = np.array([index_E_n[inner_mask[g]][0] for g in range(G)])
    upper_index = np.array([index_E_n[inner_mask[g]][-1] for g in range(G)])

    # Convert the mask to initialize the partial enery matrix
    # Hack off the last index of the mask to make the right size
    pem = np.array(inner_mask[:, :-1], dtype=float)

    # Check for partial contibutions at the edges
    for g in range(G):
        # Lower bound
        lig = lower_index[g]
        if lig != 0:
            pem[g][lig-1] = (E_n[lig] - E_g[g]) / (E_n[lig] - E_n[lig-1])

        # Upper bound
        uig = upper_index[g]
        if uig < N:
            pem[g][uig] = (E_g[g+1] - E_n[uig]) / (E_n[uig+1] - E_n[uig])

    return pem


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
    pem = partial_energy_matrix(E_g, E_n)
    phi_g_array = np.dot(pem, phi_n)
    return phi_g_array


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
