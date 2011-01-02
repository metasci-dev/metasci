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