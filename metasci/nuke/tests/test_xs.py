import os

import numpy as np
import tables as tb

from nose.tools import assert_equal, assert_not_equal
from numpy.testing import assert_array_equal, assert_array_almost_equal

from metasci.nuke import xs
from metasci.nuke import nuc_data



def test_xs_cache_E_n():
    with tb.openFile(nuc_data, 'r') as f:
        E_n = np.array(f.root.neutron.xs_mg.E_g)

    from_cache = xs.xs_cache['E_n']
    assert_not_equal(id(E_n), id(from_cache))
    assert_equal(id(from_cache), id(xs.xs_cache['E_n']))

    assert_array_equal(E_n, xs.xs_cache['E_n'])


def test_xs_cache_sigma_f_n():
    with tb.openFile(nuc_data, 'r') as f:
        sigma_f_n_U235 = np.array(f.root.neutron.xs_mg.fission[28]['xs'])

    from_cache = xs.xs_cache['sigma_f_n_922350']

    assert_not_equal(id(sigma_f_n_U235), id(from_cache))
    assert_equal(id(from_cache), id(xs.xs_cache['sigma_f_n_922350']))

    assert_array_equal(sigma_f_n_U235, xs.xs_cache['sigma_f_n_922350'])

#
# Test Partial Energy Matrix
#

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


#
# Test Partial Energy Matrix
#

def test_phi_g1():
    E_g = np.array([0.0, 10.0])
    E_n = np.array([0.0, 10.0])

    phi_n = np.ones(1)

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([1.0])

    assert_array_equal(phi_g, expected)    


def test_phi_g2():
    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 5.0, 10.0])

    phi_n = np.ones(2)

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([1.0, 1.0])

    assert_array_equal(phi_g, expected)    


def test_phi_g3():
    E_g = np.array([1.25, 5.0, 7.5])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([1.5, 1.0])

    assert_array_equal(phi_g, expected)    


def test_phi_g4():
    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([2.0, 2.0])

    assert_array_equal(phi_g, expected)    


def test_phi_g5():
    E_g = np.array([0.0, 4.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([1.6, 2.4]) 

    assert_array_equal(phi_g, expected)    


def test_phi_g6():
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([1.6, 1.6])

    # Floating point error here requires 'alomst' equal
    assert_array_almost_equal(phi_g, expected)    


def test_phi_g7():
    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.array([0.0, 2.0, 1.0, 0.5])

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([1.2, 1.9])

    # Floating point error here requires 'alomst' equal
    assert_array_almost_equal(phi_g, expected)    

