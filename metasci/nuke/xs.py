"""This module provides functionality to automatically extract cross-sections
from the nulear database."""
import re

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

        # Grab fission cross-sections from the file
        m = re.match('sigma_f_n_(\d+)', key)
        if (m is not None) and (key not in self):
            iso = int(m.group(1))
            self[key] = get_sigma_f_n(iso)

        # Calculate the low-res flux as needed
        if (key == 'phi_g') and ('phi_g' not in self):
            self['phi_g'] = phi_g(self['E_g'], self['E_n'], self['phi_n'])

        # Return the value requestion
        return super(XSCache, self).__getitem__(key)


    def __setitem__(self, key, value):
        """Overwrites the dict's key lookup by prodiving some custom setting functionality.
        """
        # Set the E_g
        if (key == 'E_g'):
            value = np.array(value, dtype=float)

            # If the E_gs are the same, don't set anything
            if ('E_g' in self):
                if (len(value) == len(self['E_g'])) and (value == self['E_g']).all():
                    return 

            # Otherwise, preload some stuff.
            self['partial_energy_matrix'] = partial_energy_matrix(value, self['E_n'])

            # And remove any previous paramters dependent on E_g
            dirty_keys = [k for k in self if '_g' in k]
            for dk in dirty_keys:
                del self[dk]


        # Set the E_n
        if (key == 'E_n'):
            value = np.array(value, dtype=float)
            
            # If the E_gs are the same, don't set anything
            if ('E_n' in self):
                if (len(value) == len(self['E_n'])) and (value == self['E_n']).all():
                    return 

            # Otherwise, preload some stuff.
            if 'E_g' in self:
                self['partial_energy_matrix'] = partial_energy_matrix(self['E_g'], value)


        # Set the high resolution flux, phi_n
        if key == 'phi_n':
            value = np.array(value, dtype=float)

            # If the flux is same, don't set or remove anything
            if ('phi_n' in self):
                if (len(value) == len(self['phi_n'])) and (value == self['phi_n']).all():
                    return 

            # And remove any previous paramters dependent on phi_n
            dirty_keys = [k for k in self if ('_g' in k and k != 'E_g')]
            for dk in dirty_keys:
                del self[dk]

        # Set the value normally
        super(XSCache, self).__setitem__(key, value)

xs_cache = XSCache()


def get_sigma_f_n(iso):
    """Grabs an isotope's fission cross-section from the nuc_data library.

    Args:
        * iso (zzaaam): Isotope name, in appropriate form.

    Returns:
        * sigma_f_n (numpy array): This isotope's fission cross-section pulled from the
          database library file.  If not present in the library, a zero-array is returned.
    """
    with tb.openFile(nuc_data, 'r') as f:
        N = f.root.neutron.xs_mg.fission.coldescrs['xs'].shape[0]
        rows = [np.array(row['xs']) for row in f.root.neutron.xs_mg.fission.where('iso_zz == {0}'.format(iso))]

    if len(rows) == 0:
        # Not fissionable, return zero-array
        sigma_f_n = np.zeros(N, dtype=float)
    elif len(rows) == 1:
        # Return cross-section from file
        sigma_f_n = rows[0]
    else:
        iso_LL = isoname.zzaaam_2_LLAAAM(iso)
        err_str = "The database contains multiple entries for the fission cross-section for {0}!".format(iso_LL)
        raise ValueError(err_str)

    return sigma_f_n


def get_sigma_a_n(iso):
    """Grabs an isotope's absorption cross-section from the nuc_data library.

    Args:
        * iso (zzaaam): Isotope name, in appropriate form.

    Returns:
        * sigma_f_n (numpy array): This isotope's absorption cross-section pulled from the
          database library file.  If not present in the library, a zero-array is returned.
    """
    with tb.openFile(nuc_data, 'r') as f:
        N = f.root.neutron.xs_mg.absorption.coldescrs['xs'].shape[0]
        rows = [np.array(row['xs']) for row in 
                f.root.neutron.xs_mg.absorption.where("(from_iso_zz == {0}) & (reaction_type != 'c')".format(iso))]

    if len(rows) == 0:
        # No absportion, return zero-array
        sigma_a_n = np.zeros(N, dtype=float)
    else:
        rows = np.array(rows)
        sigma_a_n = rows.sum(axis=0)

    return sigma_a_n


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


def partial_group_collapse(sigma_n, E_g=None, E_n=None, phi_n=None):
    """Calculates the group cross-sections for an isotope for a new, lower resolution
    group structure using a higher fidelity flux.  Note that g indexes G, n indexes N, and G < N.

    Args:
        * sigma_n (sequence of floats): A high-fidelity cross-section.

    Keyword Args:
        If any of these are None-valued, values from the cache are used.

        * E_g (sequence of floats): New, lower fidelity energy group structure [MeV]
          that is of length G+1. Ordered from lowest-to-highest energy.
        * E_n (sequence of floats): higher resolution energy group structure [MeV]
          that is of length N+1. Ordered from lowest-to-highest energy.
        * phi_n (sequence of floats): The high-fidelity flux [n/cm^2/s] to collapse the fission 
          cross-section over.  Length N.  Ordered from lowest-to-highest energy.

    Returns:
        * sigma_g (numpy array): A numpy array of the collapsed fission cross-section.
    """
    # Load the appropriate values into the cache
    if E_n is not None:
        xs_cache['E_n'] = E_n

    if E_g is not None:
        xs_cache['E_g'] = E_g

    if phi_n is not None:
        xs_cache['phi_n'] = phi_n

    # Get some other stuff from the cache
    pem = xs_cache['partial_energy_matrix']

    # Calulate partial group collapse
    sigma_g_numer = np.dot(pem, sigma_n * xs_cache['phi_n'])
    sigma_g_denom = xs_cache['phi_g']
    sigma_g = sigma_g_numer / sigma_g_denom

    return sigma_g


################################
### Cross-section generators ###
################################

def sigma_f(iso, E_g=None, E_n=None, phi_n=None):
    """Calculates the neutron fission cross-section for an isotope for a new, lower resolution
    group structure using a higher fidelity flux.  Note that g indexes G, n indexes N, and G < N.

    Note: This always pulls the fission cross-section out of the nuc_data library.    

    Args:
        * iso (int or str): An isotope to calculate the fission cross-section for.

    Keyword Args:
        If any of these are None-valued, values from the cache are used.

        * E_g (sequence of floats): New, lower fidelity energy group structure [MeV]
          that is of length G+1. Ordered from lowest-to-highest energy.
        * E_n (sequence of floats): higher resolution energy group structure [MeV]
          that is of length N+1. Ordered from lowest-to-highest energy.
        * phi_n (sequence of floats): The high-fidelity flux [n/cm^2/s] to collapse the fission 
          cross-section over.  Length N.  Ordered from lowest-to-highest energy.

    Returns:
        * sigma_f_g (numpy array): A numpy array of the collapsed fission cross-section.
    """
    # Ensure that the low-fidelity group structure is in the cache
    if E_n is not None:
        xs_cache['E_n'] = E_n

    if E_g is not None:
        xs_cache['E_g'] = E_g

    if phi_n is not None:
        xs_cache['phi_n'] = phi_n

    # Get the fission XS
    iso_zz = isoname.mixed_2_zzaaam(iso)
    sigma_f_n_iso_zz = 'sigma_f_n_{0}'.format(iso_zz)
    sigma_f_g_iso_zz = 'sigma_f_g_{0}'.format(iso_zz)

    # Don't recalculate anything if you don't have to
    if sigma_f_g_iso_zz in xs_cache:
        return xs_cache[sigma_f_g_iso_zz]
    else:
        sigma_f_n = xs_cache[sigma_f_n_iso_zz]

    # Perform the group collapse, knowing that the right data is in the cache
    sigma_f_g = partial_group_collapse(sigma_f_n)

    # Put this value back into the cache, with the appropriate label
    xs_cache[sigma_f_g_iso_zz] = sigma_f_g

    return sigma_f_g
