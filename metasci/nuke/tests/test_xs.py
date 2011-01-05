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


def test_get_sigma_a_n1():
    # Test example with one entry
    sigma_a_n = xs.get_sigma_a_n(10010)

    expected = np.array([
         9.96360000e-01,   6.07160000e-01,   4.72250000e-01,
         3.99360000e-01,   3.52130000e-01,   3.18550000e-01,
         2.93030000e-01,   2.69290000e-01,   2.46390000e-01,
         2.27390000e-01,   2.11380000e-01,   1.95050000e-01,
         1.76480000e-01,   1.50820000e-01,   1.20850000e-01,
         9.42780000e-02,   7.26880000e-02,   5.65930000e-02,
         4.38660000e-02,   3.42250000e-02,   2.68640000e-02,
         2.08190000e-02,   1.61410000e-02,   1.27460000e-02,
         9.89720000e-03,   7.70440000e-03,   5.88870000e-03,
         4.59450000e-03,   3.57770000e-03,   2.78980000e-03,
         2.18600000e-03,   1.74630000e-03,   1.32340000e-03,
         1.13220000e-03,   1.04010000e-03,   9.54700000e-04,
         7.95210000e-04,   6.03360000e-04,   4.53900000e-04,
         3.60530000e-04,   2.99660000e-04,   2.36340000e-04,
         1.71240000e-04,   1.19680000e-04,   8.19370000e-05,
         5.64050000e-05,   4.48700000e-05,   3.98850000e-05,
         3.71800000e-05,   3.54070000e-05,   3.45980000e-05,
         3.44370000e-05,   3.43310000e-05,   3.43120000e-05,
         3.49470000e-05,   3.58100000e-05,   3.62750000e-05,
         3.60950000e-05,   3.48160000e-05,   2.97130000e-05,
         2.89150000e-05,   2.67690000e-05,   2.60400000e-05])

    assert_array_equal(sigma_a_n, expected)


def test_get_sigma_a_n2():
    # Test example with multiple entries but not that have reaction_type = 'c'
    sigma_a_n = xs.get_sigma_a_n(10020)

    expected = np.array([
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.       ,  0.       ,  0.       ,  0.       ,  0.       ,
        0.0011439,  0.019188 ,  0.047752 ,  0.087012 ,  0.18032  ,
        0.18847  ,  0.20148  ,  0.2053   ])

    expected += np.array([
         1.52240000e-03,   9.30720000e-04,   7.23200000e-04,
         6.13370000e-04,   5.40020000e-04,   4.89210000e-04,
         4.49150000e-04,   4.13780000e-04,   3.78790000e-04,
         3.48940000e-04,   3.24070000e-04,   2.99220000e-04,
         2.70960000e-04,   2.31220000e-04,   1.85480000e-04,
         1.44650000e-04,   1.11350000e-04,   8.69150000e-05,
         6.73760000e-05,   5.26030000e-05,   4.15620000e-05,
         3.20750000e-05,   2.48770000e-05,   1.95970000e-05,
         1.52790000e-05,   1.17130000e-05,   9.09390000e-06,
         7.05040000e-06,   5.53650000e-06,   4.39520000e-06,
         3.43240000e-06,   2.68440000e-06,   2.00390000e-06,
         1.70200000e-06,   1.58790000e-06,   1.46300000e-06,
         1.28160000e-06,   1.09490000e-06,   1.01380000e-06,
         1.04050000e-06,   1.07440000e-06,   1.17800000e-06,
         1.41250000e-06,   1.77080000e-06,   2.30510000e-06,
         2.96130000e-06,   3.51820000e-06,   3.92870000e-06,
         4.43470000e-06,   4.95710000e-06,   5.57910000e-06,
         6.32150000e-06,   7.01920000e-06,   7.69640000e-06,
         8.43000000e-06,   9.14600000e-06,   9.83610000e-06,
         1.03010000e-05,   1.06530000e-05,   9.44650000e-06,
         9.18780000e-06,   8.38510000e-06,   8.10500000e-06])

    assert_array_equal(sigma_a_n, expected)


def test_get_sigma_a_n3():
    # Test example with multiple entries including one that has reaction_type = 'c'
    sigma_a_n = xs.get_sigma_a_n(10020)

    expected = np.array([])

    expected += np.array([
         1.52240000e-03,   9.30720000e-04,   7.23200000e-04,
         6.13370000e-04,   5.40020000e-04,   4.89210000e-04,
         4.49150000e-04,   4.13780000e-04,   3.78790000e-04,
         3.48940000e-04,   3.24070000e-04,   2.99220000e-04,
         2.70960000e-04,   2.31220000e-04,   1.85480000e-04,
         1.44650000e-04,   1.11350000e-04,   8.69150000e-05,
         6.73760000e-05,   5.26030000e-05,   4.15620000e-05,
         3.20750000e-05,   2.48770000e-05,   1.95970000e-05,
         1.52790000e-05,   1.17130000e-05,   9.09390000e-06,
         7.05040000e-06,   5.53650000e-06,   4.39520000e-06,
         3.43240000e-06,   2.68440000e-06,   2.00390000e-06,
         1.70200000e-06,   1.58790000e-06,   1.46300000e-06,
         1.28160000e-06,   1.09490000e-06,   1.01380000e-06,
         1.04050000e-06,   1.07440000e-06,   1.17800000e-06,
         1.41250000e-06,   1.77080000e-06,   2.30510000e-06,
         2.96130000e-06,   3.51820000e-06,   3.92870000e-06,
         4.43470000e-06,   4.95710000e-06,   5.57910000e-06,
         6.32150000e-06,   7.01920000e-06,   7.69640000e-06,
         8.43000000e-06,   9.14600000e-06,   9.83610000e-06,
         1.03010000e-05,   1.06530000e-05,   9.44650000e-06,
         9.18780000e-06,   8.38510000e-06,   8.10500000e-06])

    assert_array_equal(sigma_a_n, expected)



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

