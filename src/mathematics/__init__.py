"""Mathematics Package: functions, objects, & help that cannot be found elsewhere.
"""

def H(x, c=0.0):
    """The Heaviside Step Function. `<http://mathworld.wolfram.com/HeavisideStepFunction.html>`_.

    Args:
       * `x` (numeric): point to evaluate at.

    Keyword Args:
       * `c` (numeric): Heaviside offset, H(x - c).

    Returns:
       0.0, 0.5, or 1.0 (float): Value of Heaviside function.

    .. math::
       :nowrap:

       \\begin{eqnarray*} 
                           & 0            & x < c \\\\ 
       H_c(x) = H(x - c) = & \\frac{1}{2} & x = c \\\\ 
                           & 1            & c < x 
       \end{eqnarray*}
    """

    try:	
        if x < c:
            pass
    except:
        Hlist = []
        for i in range(len(x)):
            Hlist.append( H(x[i], c) )
        return Hlist

    if x < c :
        return 0.0
    elif x == c:
        return 0.5
    elif c < x:
        return 1.0
    else:
        raise ArithmeticError("Heaviside Step Function failed.  x = {0} and c = {1} could not be compared!".format(x,c))

def orderfloor(x, base=10.0):
    """Returns the order of magnitude of base that is just less than x.

    Args:
        * `x` (numeric): point to evaluate floor of.
        * `base` (numeric): Base value for the order of magnitude.
    """

    n = 0

    if x < base**n:
        while x < base**n:
            n = n - 1
    elif base**n < x:
        while base**n < x:
            n = n + 1
        n = n - 1

    return n

def orderceil(x, base=10.0):
    """Returns the order of magnitude of base that is just greater than x.

    Args:
        * `x` (numeric): point to evaluate the ceiling of.
        * `base` (numeric): Base value for the order of magnitude.
    """

    n = 0

    if x < base**n:
        while x < base**n:
            n = n - 1
        n = n + 1
    elif base**n < x:
        while base**n < x:
            n = n + 1

    return n

def min_above(x, p):
    """Finds the entry in a data set that is the lowest value above some threshold point.

    Args:
        * `x` (numeric sequence): The data set to search in.
        * `p` (numeric): The threshold point such that p < x_min.

    Returns:
        * `x_min` (numeric):  The lowest value in x that is still above p.
          for instance min_above([12.0, 0.0, 1,], 0) = 1.  Returns min(x) if no points
          are bove p.
    """

    x_min = min(x)
    if x_min <= p:
        x_min = max(x)
        for n in range(len(x)):
            if (p < x[n] < x_min):
                x_min = x[n]
    return x_min

def max_below(x, p):
    """Finds the entry in a data set that is the highest value below some threshold point.

    Args:
        * `x` (numeric sequence): The data set to search in.
        * `p` (numeric): The threshold point such that x_max < p.

    Returns:
        * `x_max` (numeric):  The highest value in x that is still below p.
          for instance max_below([12.0, 0.0, 1,], 10) = 1.  Returns max(x) if no points
          are bove p.
    """

    x_max = max(x)
    if p <= x_max:
        x_max = min(x)
        for n in range(len(x)):
            if (x_max < x[n] < p):
                x_max = x[n]
    return x_max

import polynomial 
