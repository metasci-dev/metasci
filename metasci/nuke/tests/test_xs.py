import os

import numpy as np
import tables as tb

from nose.tools import assert_equal, assert_not_equal
from numpy.testing import assert_array_equal, assert_array_almost_equal

from metasci.nuke import xs
from metasci.nuke import nuc_data



def test_xs_cache_E_n():
    xs.xs_cache.clear()

    with tb.openFile(nuc_data, 'r') as f:
        E_n = np.array(f.root.neutron.xs_mg.E_g)

    from_cache = xs.xs_cache['E_n']
    assert_not_equal(id(E_n), id(from_cache))
    assert_equal(id(from_cache), id(xs.xs_cache['E_n']))

    assert_array_equal(E_n, xs.xs_cache['E_n'])


def test_xs_cache_sigma_f_n():
    xs.xs_cache.clear()

    with tb.openFile(nuc_data, 'r') as f:
        sigma_f_n_U235 = np.array(f.root.neutron.xs_mg.fission[28]['xs'])

    from_cache = xs.xs_cache['sigma_f_n_922350']

    assert_not_equal(id(sigma_f_n_U235), id(from_cache))
    assert_equal(id(from_cache), id(xs.xs_cache['sigma_f_n_922350']))

    assert_array_equal(sigma_f_n_U235, xs.xs_cache['sigma_f_n_922350'])


def test_xs_cache_set_E_g():
    xs.xs_cache.clear()

    # Add an energy stucture
    xs.xs_cache['E_g'] = [1.0, 10.0]
    E_g = xs.xs_cache['E_g']

    # Assert that the cache is working
    assert_equal(E_g.shape, (2, ))
    assert_equal(id(E_g), id(xs.xs_cache['E_g']))

    # Assert that the cache has been reloaded
    xs.xs_cache['E_g'] = [1.0, 2.0, 10.0]
    assert_not_equal(id(E_g), id(xs.xs_cache['E_g']))
    assert_equal(len(E_g), 2)
    assert_equal(len(xs.xs_cache['E_g']), 3)

    # Assert that the partial energy matrix is calculated
    assert_equal(len(xs.xs_cache['partial_energy_matrix']), 2)

    # Assert that the reloading is done properly
    xs.xs_cache['has_some_g'] = True
    xs.xs_cache['E_g'] = [1.0, 2.0, 8.0, 10.0]
    assert_equal(len(xs.xs_cache['partial_energy_matrix']), 3)
    assert 'has_some_g' not in xs.xs_cache
    

def test_xs_cache_set_phi_n():
    xs.xs_cache.clear()

    xs.xs_cache['E_n'] = np.array([0.0, 5.0, 10.0])
    xs.xs_cache['E_g'] = np.array([0.0, 5.0, 10.0])
    xs.xs_cache['phi_n'] = [1.0, 10.0]
    assert_array_equal(xs.xs_cache['phi_n'], np.array([1.0, 10.0]))

    # Test that resetting the flux cleans the cache properly
    phi_g = xs.xs_cache['phi_g']
    assert 'phi_g' in xs.xs_cache

    xs.xs_cache['phi_n'] = [1.0, 5.0]

    assert 'E_g' in xs.xs_cache
    assert 'phi_g' not in xs.xs_cache


def test_xs_cache_get_phi_g():
    xs.xs_cache.clear()

    xs.xs_cache['E_n'] = np.array([0.0, 5.0, 10.0])
    xs.xs_cache['E_g'] = np.array([0.0, 5.0, 10.0])

    xs.xs_cache['phi_n'] = [1.0, 1.0]

    phi_g = xs.xs_cache['phi_g']

    expected = np.array([1.0, 1.0])

    assert_array_equal(phi_g, expected)    

#
# Test Partial Energy Matrix
#

def test_partial_energy_matrix1():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 10.0])
    E_n = np.array([0.0, 10.0])

    pem = xs.partial_energy_matrix(E_g, E_n)

    expected = np.array([[1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix2():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 5.0, 10.0])

    pem = xs.partial_energy_matrix(E_g, E_n)

    expected = np.array([[1.0, 0.0], 
                         [0.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix3():
    xs.xs_cache.clear()

    E_g = np.array([1.25, 5.0, 7.5])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = xs.partial_energy_matrix(E_g, E_n)

    expected = np.array([[0.5, 1.0, 0.0, 0.0], 
                         [0.0, 0.0, 1.0, 0.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix4():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = xs.partial_energy_matrix(E_g, E_n)

    expected = np.array([[1.0, 1.0, 0.0, 0.0], 
                         [0.0, 0.0, 1.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix5():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 4.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = xs.partial_energy_matrix(E_g, E_n)

    expected = np.array([[1.0, 0.6, 0.0, 0.0], 
                         [0.0, 0.4, 1.0, 1.0]])

    assert_array_equal(pem, expected)    


def test_partial_energy_matrix6():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    pem = xs.partial_energy_matrix(E_g, E_n)

    expected = np.array([[1.0, 0.6, 0.0, 0.0], 
                         [0.0, 0.4, 1.0, 0.2]])

    assert_array_equal(pem, expected)    


#
# Test cache helper functions.
#
def test_get_sigma_f_n1():
    sigma_f_n = xs.get_sigma_f_n(922350)
    expected = np.array([1.74780000e+03,   1.09570000e+03,   8.54720000e+02,
                         8.21910000e+02,   5.96110000e+02,   6.55820000e+02,
                         4.85430000e+02,   5.24960000e+02,   4.01070000e+02,
                         3.84060000e+02,   8.32680000e+02,   3.68510000e+02,
                         2.66930000e+02,   2.27710000e+02,   1.83750000e+02,
                         1.67020000e+02,   7.96280000e+01,   6.53830000e+01,
                         2.87850000e+01,   1.43510000e+01,   1.87710000e+01,
                         1.92710000e+01,   7.72680000e+01,   4.90740000e+01,
                         5.32240000e+01,   4.62680000e+01,   2.47770000e+01,
                         2.08130000e+01,   2.07720000e+01,   1.36800000e+01,
                         1.30990000e+01,   8.23490000e+00,   6.81700000e+00,
                         8.26300000e+00,   5.23320000e+00,   4.96880000e+00,
                         4.32240000e+00,   3.26220000e+00,   2.71850000e+00,
                         2.31530000e+00,   2.16830000e+00,   1.98670000e+00,
                         1.80300000e+00,   1.61870000e+00,   1.46980000e+00,
                         1.32110000e+00,   1.23810000e+00,   1.18940000e+00,
                         1.15190000e+00,   1.13810000e+00,   1.18470000e+00,
                         1.22020000e+00,   1.25640000e+00,   1.29290000e+00,
                         1.26850000e+00,   1.19820000e+00,   1.12000000e+00,
                         1.06560000e+00,   1.53220000e+00,   2.06170000e+00,
                         2.10070000e+00,   1.96770000e+00,   1.96770000e+00])

    assert_array_equal(sigma_f_n, expected)


def test_get_sigma_f_n2():
    sigma_f_n = xs.get_sigma_f_n(10010)
    expected = np.zeros(63)
    assert_array_equal(sigma_f_n, expected)


#
# Test Partial Energy Matrix
#

def test_phi_g1():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 10.0])
    E_n = np.array([0.0, 10.0])

    phi_n = np.ones(1)

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([1.0])

    assert_array_equal(phi_g, expected)    


def test_phi_g2():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 5.0, 10.0])

    phi_n = np.ones(2)

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([1.0, 1.0])

    assert_array_equal(phi_g, expected)    


def test_phi_g3():
    xs.xs_cache.clear()

    E_g = np.array([1.25, 5.0, 7.5])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([1.5, 1.0])

    assert_array_equal(phi_g, expected)    


def test_phi_g4():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([2.0, 2.0])

    assert_array_equal(phi_g, expected)    


def test_phi_g5():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 4.0, 10.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([1.6, 2.4]) 

    assert_array_equal(phi_g, expected)    


def test_phi_g6():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([1.6, 1.6])

    # Floating point error here requires 'alomst' equal
    assert_array_almost_equal(phi_g, expected)    


def test_phi_g7():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.array([0.0, 2.0, 1.0, 0.5])

    phi_g = xs.phi_g(E_g, E_n, phi_n)

    expected = np.array([1.2, 1.9])

    # Floating point error here requires 'alomst' equal
    assert_array_almost_equal(phi_g, expected)    


#
# Partial Group Collapse Tests
#

def test_partial_group_collapse1():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.ones(4)
    sigma_n = np.ones(4)    

    sigma_g = xs.partial_group_collapse(sigma_n, E_g, E_n, phi_n)

    expected = np.array([1.0, 1.0])

    assert_array_equal(sigma_g, expected)


def test_partial_group_collapse2():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 5.0, 10.0])
    E_n = np.array([0.0, 5.0, 10.0])

    phi_n = np.array([2.0, 1.0])
    sigma_n = np.ones(2)    

    sigma_g = xs.partial_group_collapse(sigma_n, E_g, E_n, phi_n)

    expected = np.array([1.0, 1.0])

    assert_array_equal(sigma_g, expected)


def test_partial_group_collapse3():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.array([2.0, 1.0, 1.0, 1.0])
    sigma_n = np.ones(4)    

    sigma_g = xs.partial_group_collapse(sigma_n, E_g, E_n, phi_n)

    expected = np.array([1.0, 1.0])

    assert_array_equal(sigma_g, expected)


def test_partial_group_collapse4():
    xs.xs_cache.clear()

    E_g = np.array([0.0, 4.0, 8.0])
    E_n = np.array([0.0, 2.5, 5.0, 7.5, 10.0])

    phi_n = np.array([2.0, 1.0, 1.0, 1.0])
    sigma_n = np.array([2.0, 1.0, 1.0, 1.0])    

    sigma_g = xs.partial_group_collapse(sigma_n, E_g, E_n, phi_n)

    expected = np.array([4.6/2.6, 1.0])

    assert_array_equal(sigma_g, expected)


#
# Fission XS Test
#

def test_sigma_f1():
    xs.xs_cache.clear()

    N = len(xs.xs_cache['E_n']) - 1

    E_g = [10.0, 20.0]
    phi_n = np.ones(N)

    sigma_f_g = xs.sigma_f(922350, E_g=E_g, phi_n=phi_n)

    expected = np.array([(2.0617 + 2.1007 + 1.9677)/3.0])

    assert_array_equal(sigma_f_g, expected)


def test_sigma_f2():
    xs.xs_cache.clear()

    N = len(xs.xs_cache['E_n']) - 1

    E_g = [10.0, 20.0]
    phi_n = np.ones(N)
    phi_n[59] = 2.0

    sigma_f_g = xs.sigma_f(922350, E_g=E_g, phi_n=phi_n)

    expected = np.array([(2.0617*2.0 + 2.1007 + 1.9677)/4.0])

    assert_array_equal(sigma_f_g, expected)


def test_sigma_f3():
    xs.xs_cache.clear()

    N = len(xs.xs_cache['E_n']) - 1

    E_g = np.copy(xs.xs_cache['E_n'])
    phi_n = np.ones(N)

    sigma_f_g = xs.sigma_f("H1", E_g=E_g, phi_n=phi_n)

    expected = np.zeros(N)

    assert_array_equal(sigma_f_g, expected)

