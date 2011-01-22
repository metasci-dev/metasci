"""This module provides functionality to automatically extract cross-sections
from the nulear database."""
import re
from itertools import product

import numpy as np
import tables as tb
from scipy import integrate
from scipy import constants
from scipy.special import erf

import isoname

from . import nuc_data

# Bolzman's constant in MeV/K
k = constants.physical_constants['Boltzmann constant in eV/K'][0] * (10**-6)

# Neutron mass in amu
m_n = constants.physical_constants['neutron mass in u'][0]

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

        # Grab absorption cross-sections from the file
        m = re.match('sigma_a_n_(\d+)', key)
        if (m is not None) and (key not in self):
            iso = int(m.group(1))
            self[key] = get_sigma_a_n(iso)

        # Grab scattering length from the file
        m = re.match('b_(\d+)', key)
        if (m is not None) and (key not in self):
            iso = int(m.group(1))
            self[key] = get_b(iso)

        # Grab atomic mass from the file
        m = re.match('M_A_(\d+)', key)
        if (m is not None) and (key not in self):
            iso = int(m.group(1))
            self[key] = get_M_A(iso)

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
        * sigma_a_n (numpy array): This isotope's absorption cross-section pulled from the
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

    # Add in the fission cross-section
    sigma_f_n = get_sigma_f_n(iso)
    sigma_a_n += sigma_f_n

    return sigma_a_n


def get_b(iso):
    """Grabs an isotope's neutron scattering length from the nuc_data library.
    Because the coherent and incoherent lenghts are in general complex, the 
    scattering length is computed as follows.

    .. math::
        b = \\sqrt{\\left| b_{\mbox{coh}} \\right|^2 + \\left| b_{\mbox{inc}} \\right|^2}

    Args:
        * iso (zzaaam): Isotope name, in appropriate form.

    Returns:
        * b (float): The scattering length in [cm].
    """
    with tb.openFile(nuc_data, 'r') as f:
        # Try finding the isotope
        rows = [(row['b_coherent'], row['b_incoherent']) for row in 
                f.root.neutron.scattering_lengths.where("(iso_zz == {0})".format(iso))]

        # Try finding the element
        if len(rows) == 0:
            rows = [(row['b_coherent'], row['b_incoherent']) for row in 
                    f.root.neutron.scattering_lengths.where("(iso_zz == {0})".format(10000 * (iso // 10000)))]

        # Guess the next lowest element
        n = 1
        while len(rows) == 0:
            rows = [(row['b_coherent'], row['b_incoherent']) for row in 
                    f.root.neutron.scattering_lengths.where("(iso_zz == {0})".format(10000 * ((iso // 10000) - n)))]

    b_coh, b_inc = rows[0]

    b = np.sqrt(np.abs(b_coh)**2 + np.abs(b_inc)**2)

    return b


def get_M_A(iso):
    """Grabs an isotope's atomic mass from the nuc_data library.

    Args:
        * iso (zzaaam): Isotope name, in appropriate form.

    Returns:
        * M_A (float): The atomic mass [amu].
    """
    with tb.openFile(nuc_data, 'r') as f:
        # Try finding the isotope
        rows = [row['value'] for row in f.root.A.where("(iso_zz == {0})".format(iso))]

        # If the isotope is not found, maybe a stable version is present.
        if len(rows) == 0:
            rows = [row['value'] for row in f.root.A.where("(iso_zz == {0})".format(10 * (iso // 10)))]

    if len(rows) == 0:
        # Take best guess 
        m = (iso // 10) % 1000
        M_A = float(m)
    else:
        M_A = rows[0]

    return M_A

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


def sigma_a(iso, E_g=None, E_n=None, phi_n=None):
    """Calculates the neutron absorption cross-section for an isotope for a new, lower resolution
    group structure using a higher fidelity flux.  Note that g indexes G, n indexes N, and G < N.

    Note: This always pulls the absorption cross-section out of the nuc_data library.    

    Args:
        * iso (int or str): An isotope to calculate the absorption cross-section for.

    Keyword Args:
        If any of these are None-valued, values from the cache are used.

        * E_g (sequence of floats): New, lower fidelity energy group structure [MeV]
          that is of length G+1. Ordered from lowest-to-highest energy.
        * E_n (sequence of floats): higher resolution energy group structure [MeV]
          that is of length N+1. Ordered from lowest-to-highest energy.
        * phi_n (sequence of floats): The high-fidelity flux [n/cm^2/s] to collapse the fission 
          cross-section over.  Length N.  Ordered from lowest-to-highest energy.

    Returns:
        * sigma_a_g (numpy array): A numpy array of the collapsed absorption cross-section.
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
    sigma_a_n_iso_zz = 'sigma_a_n_{0}'.format(iso_zz)
    sigma_a_g_iso_zz = 'sigma_a_g_{0}'.format(iso_zz)

    # Don't recalculate anything if you don't have to
    if sigma_a_g_iso_zz in xs_cache:
        return xs_cache[sigma_a_g_iso_zz]
    else:
        sigma_a_n = xs_cache[sigma_a_n_iso_zz]

    # Perform the group collapse, knowing that the right data is in the cache
    sigma_a_g = partial_group_collapse(sigma_a_n)

    # Put this value back into the cache, with the appropriate label
    xs_cache[sigma_a_g_iso_zz] = sigma_a_g

    return sigma_a_g

#
# Scattering cross section
#

def alpha(E_prime, E, theta, M_A=1.0, T=300.0):
    """Scattering kernel alpha value.

    .. math::

        \\alpha = \\frac{E^\prime + E - 2\sqrt{E^\prime E}\cos\\theta}{\\frac{M_A}{m_n}kT}

    Args:
        * E_prime (float): The exiting energy of the neutron after the 
          scattering event [MeV].
        * E (float): The incident energy of the neutron prior to the 
          scattering event [MeV].
        * theta (float): Scattering angle in [radians].

    Keyword Args:
        * M_A (float): Atomic mass of the target nucleus [amu].
        * T (float): Tempurature of the target material [kelvin].

    Returns:
        * a (float): alpha value.
    """
    a = (E_prime + E - 2 * np.sqrt(E_prime*E) * np.cos(theta)) / (k * T * M_A / m_n)
    return a


def beta(E_prime, E, T=300.0):
    """Scattering kernel beta value.

    .. math::

        \\beta = \\frac{E^\prime - E}{kT}

    Args:
        * E_prime (float): The exiting energy of the neutron after the 
          scattering event [MeV].
        * E (float): The incident energy of the neutron prior to the 
          scattering event [MeV].

    Keyword Args:
        * T (float): Tempurature of the target material [kelvin].

    Returns:
        * b (float): beta value.
    """
    b = (E_prime - E) / (k*T)
    return b


def alpha_given_theta_0(E_prime, E, M_A=1.0, T=300.0):
    """Scattering kernel alpha value at the lower bound of the scattering angle.

    .. math::

        \\alpha_{\\theta=0} = \\frac{E^\prime + E - 2\sqrt{E^\prime E}}{\\frac{M_A}{m_n}kT}

    Args:
        * E_prime (float): The exiting energy of the neutron after the 
          scattering event [MeV].
        * E (float): The incident energy of the neutron prior to the 
          scattering event [MeV].

    Keyword Args:
        * M_A (float): Atomic mass of the target nucleus [amu].
        * T (float): Tempurature of the target material [kelvin].
    """
    a = (E_prime + E - 2 * np.sqrt(E_prime*E)) / (k * T * M_A / m_n)
    return a


def alpha_given_theta_pi(E_prime, E, M_A=1.0, T=300.0):
    """Scattering kernel alpha value at the upper bound of the scattering angle.

    .. math::

        \\alpha_{\\theta=\pi} = \\frac{E^\prime + E + 2\sqrt{E^\prime E}}{\\frac{M_A}{m_n}kT}

    Args:
        * E_prime (float): The exiting energy of the neutron after the 
          scattering event [MeV].
        * E (float): The incident energy of the neutron prior to the 
          scattering event [MeV].

    Keyword Args:
        * M_A (float): Atomic mass of the target nucleus [amu].
        * T (float): Tempurature of the target material [kelvin].
    """
    a = (E_prime + E + 2 * np.sqrt(E_prime*E)) / (k * T * M_A / m_n)
    return a


def one_over_gamma_squared(E):
    """The realitivistic correction factor for the bound scattering length.

    .. math::

        \\frac{1}{\gamma^2} = \left( 1 - \\frac{2E}{931.46 \\cdot m_n} \\right)

    Args:
        * E (float): The incident energy of the neutron prior to the 
          scattering event [MeV].
    """
    v_squared = 2.0 * E / m_n
    rcf = 1.0 - (v_squared / 931.46)
    return rcf


def d2sigma_s_dE_prime_dOmega(E_prime, E, theta, b=1.0, M_A=1.0, T=300.0):
    """Computes the double differential total scattering cross section from the equation.  
    This is likely only valid in the thermal region.

    .. math::

        \\frac{d\sigma_s(E)}{dE^\prime} = \mbox{[1]}

    FIXME: I am untested

    Args:
        * E_prime (float): The exiting energy of the neutron after the 
          scattering event [MeV].
        * E (float): The incident energy of the neutron prior to the 
          scattering event [MeV].
        * theta (float): Scattering angle in [radians].

    Keyword Args:
        * b (float): The bound scattering length of the target nucleus.
        * M_A (float): Atomic mass of the target nucleus [amu].
        * T (float): Tempurature of the target material [kelvin].

    Refs:
        1. Mattes M, Keinert J. Thermal neutron scattering data for the moderator 
           materials H2O, D2O and ZrHx in ENDF-6 format and as ACE library for 
           MCNP (X) codes. IAEA report INDC (NDS)-0470. 2005;(April). Available at: 
           http://200.136.52.101/reports-new/indc-reports/indc-nds/indc-nds-0470.pdf.
    """
    kT = k * T
    rcf = one_over_gamma_squared(E)

    _alpha = alpha(E_prime, E, theta, M_A, T) 
    _beta = beta(E_prime, E, T)

    power_term = (_beta/2.0) + (_alpha/4.0) + (_beta**2 / (4.0 * _alpha))

    return rcf * (b**2 / kT) * np.sqrt((np.pi * E_prime) / (_alpha * E)) * np.exp(-power_term)
    
#
# Failed analytic soltion to the integral of the double differential 
# scattering cross-section over all solid angles.  This is because the 
# S(alpha, beta) term was only valid in the thermal regime but was used 
# even up to fast energies.
#

#def dsigma_s_dE_prime(E_prime, E, b=1.0, M_A=1.0, T=300.0):
#    """Computes the differential total scattering cross section from an analytic
#    solution to the integral of the double-differentional scattering cross section, 
#    integrated over all solid angles.
#
#    .. math::
#
#        \\frac{d\sigma_s(E)}{dE^\prime} = b^2 \left( 1 - \\frac{2E}{931.46 \\cdot m_n} \\right) 
#            \\frac{e^{-\\frac{\\beta + |\\beta|}{2}}}{2E} \\frac{M_A}{m_n} Q
#
#        Q = \left( \mbox{Erf}\left(\\frac{|\\beta| - \\alpha_{\\theta=0}}{2 \sqrt{\\alpha_{\\theta=0}}}\\right) 
#            - \mbox{Erf}\left(\\frac{|\\beta| - \\alpha_{\\theta=\pi}}{2 \sqrt{\\alpha_{\\theta=\pi}}}\\right) \\right) 
#            - e^{-\\frac{|\\beta|}{2}} \left( \mbox{Erf}\left(\\frac{|\\beta| + \\alpha_{\\theta=0}}{2 \sqrt{\\alpha_{\\theta=0}}}\\right) 
#            - \mbox{Erf}\left(\\frac{|\\beta| + \\alpha_{\\theta=\pi}}{2 \sqrt{\\alpha_{\\theta=\pi}}}\\right) \\right)
#
#    Args:
#        * E_prime (float): The exiting energy of the neutron after the 
#          scattering event [MeV].
#        * E (float): The incident energy of the neutron prior to the 
#          scattering event [MeV].
#
#    Keyword Args:
#        * b (float): The bound scattering length of the target nucleus.
#        * M_A (float): Atomic mass of the target nucleus [amu].
#        * T (float): Tempurature of the target material [kelvin].
#    """
#    kT = k * T
#    rcf = one_over_gamma_squared(E)
#
#    alpha_lower = alpha_given_theta_0(E_prime, E, M_A, T)
#    alpha_upper = alpha_given_theta_pi(E_prime, E, M_A, T)
#    _beta = beta(E_prime, E, T)
#    abs_beta = np.abs(_beta)
#
#    Q = erf((abs_beta - alpha_lower) / (2.0 * np.sqrt(alpha_lower))) - \
#        erf((abs_beta - alpha_upper) / (2.0 * np.sqrt(alpha_upper))) + \
#        np.exp(-abs_beta/2.0) * (erf((abs_beta + alpha_lower) / (2.0 * np.sqrt(alpha_lower))) - \
#        erf((abs_beta + alpha_upper) / (2.0 * np.sqrt(alpha_upper))))
#
#    deriv = (rcf * b**2) * (np.exp(-(_beta + abs_beta)/2.0) / (2.0*E)) * (M_A / m_n) * Q * np.exp(-np.sqrt(E))
#
#    return deriv


def E_prime_min(E, M_A=1.0):
    """The minimum possible exiting enegy of a neuron after a scattering collision. 
    This is based on the incident energy and the mass of the target.
    For a proof, use the conservation of energy and momentum.

    .. math::

        \mbox{min}(E^\prime) = \left(\\frac{M_A - m_n}{M_A + m_n}\\right)^2 E

    Args:
        * E (float): The incident energy of the neutron prior to the 
          scattering event [MeV].

    Keyword Args:
        * M_A (float): Atomic mass of the target nucleus [amu].
    """
    alph = ((M_A - m_n)/(M_A + m_n))**2 * E
    min_E = alph * E
    return min_E


def E_prime_max(E, T=300.0):
    """The maximum possible exiting enegy of a neuron after a scattering collision. 
    This is based on the incident energy and the tempurature of the target.
    The neutron can gain no more energy than the kinetic energy of the target.
    In a macroscopic system, this is on average equal to kT.

    .. math::

        \mbox{max}(E^\prime) = E + kT

    Args:
        * E (float): The incident energy of the neutron prior to the 
          scattering event [MeV].

    Keyword Args:
        * T (float): Tempurature of the target material [kelvin].
    """
    max_E = E + (k * T)
    return max_E


def sigma_s_const(b):
    """Computes the constant scattering cross-section based on the 
    scattering length.

    .. math::

        \sigma_s = 4 \pi b^2

    Args:
        * b (float): The bound scattering length of the target nucleus [cm].

    Returns:
        * sig_s (float): the micorscopic scattering cross-section [barns].
    """
    sig_s = 4.0 * np.pi * (b**2) * (10.0**24)

    return sig_s


def sigma_s_E(E, b=1.0, M_A=1.0, T=300.0):
    """Computes the total scattering cross section from an empirical model.

    .. math::

        \\sigma_s(E) = 4 \\pi b^2 \\cdot \\left( 1 - \\frac{2E}{931.46 \\cdot m_n} \\right) \\cdot
                      \\left( 1 + \\frac{m_n}{M_A} \\frac{kT}{E} \\cdot e^{-\\frac{M_A}{m_n}\\frac{E}{kT}} \\right) 
                      \\cdot \\left( 1 - \\mbox{Exp}\\left[-\\sqrt{\\frac{0.1}{E}}\\right] \\right)

    Args:
        * E (float): The incident energy of the neutron prior to the 
          scattering event [MeV].

    Keyword Args:
        * b (float): The bound scattering length of the target nucleus.
        * M_A (float): Atomic mass of the target nucleus [amu].
        * T (float): Tempurature of the target material [kelvin].
    """
    kT_over_AE = k * T / ((M_A / m_n) * E)

    sig_s = sigma_s_const(b)
    rcf = one_over_gamma_squared(E)
    
    sig_s_E = (rcf * sig_s) * (1.0 + kT_over_AE * np.exp(-1.0/kT_over_AE)) * (1.0 - np.exp(-np.sqrt(0.1/E))) 

    return sig_s_E


def P(E, E_prime, M_A=1.0, T=300.0):
    """Calculates the neutron scattering transfer probability from

    .. math::
        P(E \\to E^\\prime) = \\frac{1}{E + kT - \\left( \\frac{M_A - m_n}{M_A + m_n} \\right)^2 E}

    as long as 

    .. math::
        \\left( \\frac{M_A - m_n}{M_A + m_n} \\right)^2 E \\le E^\\prime \\le E + kT

    The probability is zero outside of this range.

    Args:
        * E (float): The incident energy of the neutron prior to the 
          scattering event [MeV].
        * E_prime (float): The exiting energy of the neutron after the 
          scattering event [MeV].

    Keyword Args:
        * M_A (float): Atomic mass of the target nucleus [amu].
        * T (float): Tempurature of the target material [kelvin].

    Returns:
        * p (float): A transfer probablity.
    """
    # Find bounds
    E_prime_lower = E_prime_min(E, M_A)    
    E_prime_upper = E_prime_max(E, T)    

    if E_prime_lower <= E_prime <= E_prime_upper:
        p = 1.0 / (E_prime_upper - E_prime_lower)
    else:
        p = 0.0

    return p
    


def sigma_s_gh(iso, T, E_g=None, E_n=None, phi_n=None):
    """Calculates the neutron scattering cross-section kernel for an isotope for a new, lower resolution
    group structure using a higher fidelity flux.  Note that g, h index G, n indexes N, and G < N.
    g is for the incident energy and h is for the exiting energy.

    .. math::
        \\sigma_{s, g\\to h} = \\frac{\\int_{E_g}^{E_{g+1}} \\int_{E_h}^{E_{h+1}} \\sigma_s(E) P(E \\to E^\\prime) \\phi(E) dE^\\prime dE}
                                {\\int_{E_g}^{E_{g+1}} \\phi(E) dE}

    Note: This always pulls the scattering length out of the nuc_data library.    

    Args:
        * iso (int or str): An isotope to calculate the fission cross-section for.
        * T (float): Tempurature of the target material [kelvin].

    Keyword Args:
        If any of these are None-valued, values from the cache are used.

        * E_g (sequence of floats): New, lower fidelity energy group structure [MeV]
          that is of length G+1. Ordered from lowest-to-highest energy.
        * E_n (sequence of floats): higher resolution energy group structure [MeV]
          that is of length N+1. Ordered from lowest-to-highest energy.
        * phi_n (sequence of floats): The high-fidelity flux [n/cm^2/s] to collapse the fission 
          cross-section over.  Length N.  Ordered from lowest-to-highest energy.

    Returns:
        * sig_s_gh (numpy array): A numpy array of the scattering kernel.
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
    sigma_s_gh_iso_zz_T = 'sigma_s_gh_{0}_{1}'.format(iso_zz, T)

    # Don't recalculate anything if you don't have to
    if sigma_s_gh_iso_zz_T in xs_cache:
        return xs_cache[sigma_s_gh_iso_zz_T]

    # Get some needed data
    G = len(xs_cache['E_g']) - 1
    b = xs_cache['b_{0}'.format(iso_zz)]
    M_A = xs_cache['M_A_{0}'.format(iso_zz)]

    # Initialize the scattering kernel
    sig_s_gh = np.zeros((G, G), dtype=float)


    # Calculate all values of the kernel
    for g, h in product(range(G), range(G)):
        # Numerator inetgration term 
        dnumer = lambda _E_prime, _E: sigma_s_E(_E, b, M_A, T) *  P(_E, _E_prime, M_A, T) * xs_cache['phi_g'][g]

        # Integral
        numer = integrate.dblquad(dnumer, xs_cache['E_g'][g], xs_cache['E_g'][g+1], 
                                          lambda E_h: xs_cache['E_g'][h], lambda E_h: xs_cache['E_g'][h+1])
        numer = numer[0]

        # Denominator term, analytically integrated
        denom = xs_cache['phi_g'][g] * (xs_cache['E_g'][g+1] - xs_cache['E_g'][g])

        # Cross section value
        sig_s_gh[g, h] = numer / denom 

    # Put this value back into the cache, with the appropriate label
    xs_cache[sigma_s_gh_iso_zz_T] = sigma_s_gh

    return sig_s_gh


def sigma_s(iso, T, E_g=None, E_n=None, phi_n=None):
    """Calculates the neutron scattering cross-section for an isotope. 

    .. math::
        \\sigma_{s, g} = \\sum_{h} \\sigma_{s, g\\to h} 

    Args:
        * iso (int or str): An isotope to calculate the fission cross-section for.
        * T (float): Tempurature of the target material [kelvin].

    Keyword Args:
        If any of these are None-valued, values from the cache are used.

        * E_g (sequence of floats): New, lower fidelity energy group structure [MeV]
          that is of length G+1. Ordered from lowest-to-highest energy.
        * E_n (sequence of floats): higher resolution energy group structure [MeV]
          that is of length N+1. Ordered from lowest-to-highest energy.
        * phi_n (sequence of floats): The high-fidelity flux [n/cm^2/s] to collapse the fission 
          cross-section over.  Length N.  Ordered from lowest-to-highest energy.

    Returns:
        * sig_s_g (numpy array): A numpy array of the scattering cross-section.
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
    sigma_s_g_iso_zz_T = 'sigma_s_g_{0}_{1}'.format(iso_zz, T)
    sigma_s_gh_iso_zz_T = 'sigma_s_gh_{0}_{1}'.format(iso_zz, T)

    # Don't recalculate anything if you don't have to
    if sigma_s_g_iso_zz_T in xs_cache:
        return xs_cache[sigma_s_g_iso_zz_T]

    # This calculation requires the scattering kernel
    if sigma_s_gh_iso_zz_T not in xs_cache:
        xs_cache[sigma_s_gh_iso_zz_T] = sigma_s_gh(iso, T, E_g, E_n, phi_n)

    # Sum over all h indeces
    sig_s_g = xs_cache[sigma_s_gh_iso_zz_T].sum(axis=1)

    # Put this value back into the cache, with the appropriate label
    xs_cache[sigma_s_g_iso_zz_T] = sigma_s_g

    return sig_s_g

