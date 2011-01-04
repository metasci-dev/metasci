import os

import numpy as np
import tables as tb

from nose.tools import assert_equal, assert_not_equal
from numpy.testing import assert_array_equal

from metasci.nuke import xs
from metasci.nuke import nuc_data



def test_xs_cache_E_n():
    with tb.openFile(nuc_data, 'r') as f:
        E_n = np.array(f.root.neutron.xs_mg.E_g)

    from_cache = xs.xs_cache['E_n']
    assert_not_equal(id(E_n), id(from_cache))
    assert_equal(id(from_cache), id(xs.xs_cache['E_n']))

    assert_array_equal(E_n, xs.xs_cache['E_n'])


def test_partial_energy_matrix1():
    E_g = np.array([0.0, 10.0])
    E_n = np.array([0.0, 10.0])

    pem = xs.partial_energy_matrix(E_g, E_n)

    expected = np.array([[1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix2():
    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 5.0, 10.0])

    pem = xs.partial_energy_matrix(E_g, E_n)

    expected = np.array([[1.0, 0.0], 
                         [0.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix3():
    E_g = np.array([1.25, 5.0, 7.5])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = xs.partial_energy_matrix(E_g, E_n)

    expected = np.array([[0.5, 1.0, 0.0, 0.0], 
                         [0.0, 0.0, 1.0, 0.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix4():
    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = xs.partial_energy_matrix(E_g, E_n)

    expected = np.array([[1.0, 1.0, 0.0, 0.0], 
                         [0.0, 0.0, 1.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix5():
    E_g = np.array([0.0, 4.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = xs.partial_energy_matrix(E_g, E_n)

    expected = np.array([[1.0, 0.6, 0.0, 0.0], 
                         [0.0, 0.4, 1.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix6():
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = xs.partial_energy_matrix(E_g, E_n)

    expected = np.array([[1.0, 0.6, 0.0, 0.0], 
                         [0.0, 0.4, 1.0, 0.2]])

    assert_array_equal(pem, expected)    


"""def test_phi_g():
    # Set up energies
    G = 10
    E_g = np.logspace(-9, 1, G+1)
    E_n = xs.xs_cache['E_n']

    # setup flux
    N = len(E_n) - 1
    phi_n = np.ones(N)

    # Calc new flux
    phi_g = xs.phi_g(E_g, E_n, phi_n)

    assert False
"""
