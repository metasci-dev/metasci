"""Functions to help integration."""

import numpy as np
from scipy import integrate

def dlbtrapz(func, x, y):
    """Performs a double trapazoidal integral for a function or an x and y quadrature.

    Args:
        * func (function): pointer to function to integrate, f(y, x).
        * x (sequence): Sequence of x points to integrate over.
        * y (sequence): Sequence of y points to integrate over.

    Returns:
        * out (float): integration value.
    """

    func_xy = np.array([[func(_y, _x) for _y in y] for _x in x])

    xtrapz = np.array([integrate.trapz(xrow, y) for xrow in func_xy])

    out = integrate.trapz(xtrapz, x)

    return x
