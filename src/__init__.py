"""MetaSci Basics"""

#__all__ = ["data", "graph"]

import os

import numpy as np

#########################
### Basic Definitions ###
#########################

delims = ".,:;[]{}()<>?!~$%^*\ \t\n"
"""Common delimiters: .,:;[]{}()<>?!~$%^*[space][tab][newline]"""

numbers = "0123456789"
"""Base 10 digits: 0123456789"""

alphaupper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
"""Uppercase letters: ABCDEFGHIJKLMNOPQRSTUVWXYZ"""

alphalower = "abcdefghijklmnopqrstuvwxyz"
"""Lowercase letters: abcdefghijklmnopqrstuvwxyz"""

alphabet = alphaupper + alphalower
"""Upper and lower case letters: ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"""

alphanumeric = alphabet + numbers
"""The alphabet and digits: ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"""

alphadelim = alphabet + delims
"""The alphabet and delimiters: ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz.,:;[]{}()<>?!~$%^*[space][tab][newline]"""

numdelim = numbers + delims
"""Digits and delimiters: 0123456789.,:;[]{}()<>?!~$%^*[space][tab][newline]"""

alphanumericdelim = alphanumeric + delims
"""The alphabet, digits, and delimiters: ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.,:;[]{}()<>?!~$%^*[space][tab][newline]"""


#########################
### General Functions ###
#########################

def OnceYouPopYouJustCantStop(l, s):
    """Pops all elements of a list l that are s and returns new list.

    Args:
        * `l` (list):   List container.
        * `s` (object): Elements that may be in l to be removed.

    Returns:
        * `new_l` (list): New list without any s elements.
    """
    for i in range(l.count(s)):
        l.remove(s)
    return l

def SplitClean(s, ss, cl):
    """Splits as string s by ss and cleans the resulting list by removing the characters in character list cl.

    Args:
        * `s`  (str):  String to be split and cleaned.
        * `ss` (str):  Split string to partition s by.
        * `cl` (list): List of characters to remove from split list, eg [',', '', ' ']

    Returns:
        * `l` (list): List of strings, like built-in function split().
    """
    l = s.split(ss)
    for c in cl:
        l = OnceYouPopYouJustCantStop(l, c)
    return l

def HasSubString(s, sl):
    """Checks if a string contains any of the substrings in a list.

    Args:
        * `s`  (str):  String to check.
        * `sl` (list): List of strings to test if present in s.

    Returns:
        * `contains` (True or False): Membership boolean.
    """

    for sub in sl:
        if 0 <= s.find(sub):
            return True
    return False

def Zeros(N):
    """Returns a list of length N that is filled with zeros.

    Args:
        * `N` (int): Length of new list.

    Returns:
        * `z` (list): List filled with zeros such that len(z) = N.

    .. warning:: This function should be avoided in favor of numpy.zeros().    
    """
    if not type(N) == int:
        N = int(N)

    l = []
    for i in range(N):
        l.append(0.0)
    return l

def make2D(lyst, xsize = 100):
    """Makes a 1D list two dimentional with each row now having length xsize.

    Args:
        * `lyst` (list): List to be reshaped.

    Keyword Args:
        * `xsize` (int): Length of new rows.

    Returns:
        * `newlyst` (list): A list of lists of size xsize.
    """
    if not( len(lyst)%xsize == 0 ):
        print "Sizes all wrong!"
        return

    newlyst = []
    for n in range( len(lyst) ):
        if n%xsize == 0:
            newlyst.append([])
        newlyst[-1].append(lyst[n])
    return newlyst

def ToFloatList(l, multby=1.0):
    """Attempts to convert all elements of a list to floats.  Then multiplies these values
    by multby.

    Args:
        * `l` (list): List of various types.

    Keyword Args:
        * `multby` (float): Value to multiply new floats by.

    Returns:
        * `newl` (list): List of floats.
    """

    newl = []
    if multby == 1.0:
        for el in l:
            newl.append( float(el) )
    else:
        for el in l:
            newl.append( float(el) * multby )

    return newl

def MultiSplit(s, d):
    """Splits a string based on delimiter(s).

    Args:
        * `s` (str): String to be split.
        * `d` (str or list): Delimiters; may either be a string of delimiters or a list of them.

    Returns:
        * `l` (list): List of strings, like built-in function split().
    """

    if isinstance(d, str):
        d = list(d)
    l = [""]

    for n in range(len(s)):
        if s[n] in d:
            l.append("")
        else:
            l[-1] = l[-1] + s[n]
    for n in range(len(l)-1, -1, -1):
        if l[n] == "":
            del l[n]

    return l

def AlphanumericSplit(s):
    """Takes a string and returns a list of strings where all alphabetic and numeric groups 
    are split from one another.  This ignores delimiters.

    Args:
        * `s` (str): String to be split.

    Returns:
        * `l` (list): List of strings, like built-in function split().
    """

    l = [""]
    AlphaOrNum = None
    for n in range(len(s)):
        sn_type = None
        if s[n] in list(delims):
            l[-1] = l[-1] + s[n]
            continue
        elif s[n] in list(numbers):
            sn_type = "Num"
        elif s[n] in list(alphabet):
            sn_type = "Alpha"
        else:
            print "%s in %s is not a number, a letter, or a delimiter."%(s[n], s)
            raise SystemExit


        if AlphaOrNum == None:
            AlphaOrNum = sn_type

        if sn_type == AlphaOrNum:
            l[-1] = l[-1] + s[n]
        else:
            l.append(s[n])
            AlphaOrNum = sn_type

    return l

def ReverseDic(d):
    """Reverses the keys and values of a dictionary.  Values are converted to strings in this process.

    Args:
        * `d` (dict): Original dictionary.

    Returns:
        * `newd` (dict): Dictionary with keys and values swapped from d.
    """
    newd = {}
    for dk in d.keys():
        newd[str(d[dk])] = dk
    return newd

def MatrixSwitchDim(m):
    """Transposes a 2D list, or matrix.
    
    Args:
        * `m` (2D list): Matrix of size AxB.

    Returns:
        * `mstar` (2D list): Matrix of size BxA.
    """

    A = len(m)
    B = len(m[0])

    mstar = []
    for b in range(B):
        mstar.append([])
        for a in range(A):
            mstar[b].append( m[a][b])

    return mstar

def grabFromList(lyst, inORout):
    """Grabs Input or Output data from a list that has the following form::

        lyst = [paramname, 1indata, 1inerr, 1outdata, 1outerr, 2indata, 2inerr, 2outdata, 2outerr, ...]

    This is useful for getting information from data files. 
    Note that it converts data to float before returning!

    Args:
        * `lyst`    (list): List of above form.
        * `inORout` (str):  Grabs output data if equal to 'Out'.  Grabs input data otherwise.

    Returns:
        * `datalyst` (list): List of data types.  See other metasci documentation.
    """

    if isinstance(lyst, str):
        lyst = lyst.split()

    templyst = []
    for n in range( len(lyst) // 4 ):
        if inORout == "Out":
            templyst.append( data.data( lyst[(4*n)+3], lyst[(4*n)+4]) )
        else:
            templyst.append( data.data( lyst[(4*n)+1], lyst[(4*n)+2]) )

    return templyst

##########################
### Physics Phunctions ###
##########################
def slope(x2, y2, x1, y1):
    """Finds the slope of a line as defined by two points.  Useful for interpolation.

    Args:
        * `x2` (numeric): Independent variable value of the second point.
        * `y2` (numeric): Dependent variable value of the second point.
        * `x1` (numeric): Independent variable value of the first point.
        * `y1` (numeric): Dependent variable value of the first point.

    Returns:
        * `slope` (numeric): Slope of the line between points (x1, y1) and (x2, y2)
    """
    return (y2 - y1) / (x2 - x1)

def SolveLine(x, x2, y2, x1, y1):
    """Solves a linear interpolation for y given x between two points (x1, y1) and (x2, y2).

    Args:
        * `x`  (numeric): Independent variable point at which to evaluate the line.
        * `x2` (numeric): Independent variable value of the second point.
        * `y2` (numeric): Dependent variable value of the second point.
        * `x1` (numeric): Independent variable value of the first point.
        * `y1` (numeric): Dependent variable value of the first point.

    Returns:
        * `y` (numeric): Solution of the line y = mx + b for line between points (x1, y1) and (x2, y2)
    """
    return (slope(x2,y2,x1,y1) * (x - x2)) + y2

def time2sec(t, u):
    """Takes a time value with units and returns this time in seconds.

    Args:
        * `t` (float-convertible): Time value.
        * `u` (str): Units of time value.  May be any standad unit expression.  Additionally, time 
          may be expressed in 'MeV'.

    Returns:
        * `tsec` (float): Float value of time in seconds.
    """
    u = u.lower()
    if t == "inf":
        return float(t) 
    elif u in ["fs", "femtosec", "femtosecond", "femtoseconds"]:
        return float(t)*(10.0**(-15))
    elif u in ["ps", "picosec", "picosecond", "picoseconds"]:
        return float(t)*(10.0**(-12))
    elif u in ["ns", "nanosec", "nanosecond", "nanoseconds"]:
        return float(t)*(10.0**(-9))
    elif u in ["us", "microsec", "microsecond", "microseconds"]:
        return float(t)*(10.0**(-6))
    elif u in ["ms", "millisec", "millisecond", "milliseconds"]:
        return float(t)*(10.0**(-3))
    elif u in ["s", "sec", "second", "seconds"]:
        return float(t)
    elif u in ["m", "min", "minute", "minutes"]:
        return float(t)*60.0
    elif u in ["h", "hr", "hour", "hours"]:
        return float(t)*3600.0
    elif u in ["d", "dy", "day", "days"]:
        return float(t)*3600.0*24.0
    elif u in ["y", "yr", "year", "years"]:
        return float(t)*3600.0*24.0*365.0
    elif u ==  "mev":
        return float(t) * (10.0**(-15)) * float(7.6e-08)/6.03
    else:
        print "Time given is not a number!"
        print "Time =", time, "Unit = ", unit
        return "nan"

###############################
### Bin Structure Functions ###
###############################

def LinearUniformBins(a, b, N):
    """Splits the range [a,b] up into N linear-uniform bins defined by N+1 points.
    For example::

        >>> LinearUniformBins(1.0, 2.5, 3)
        [1.0, 1.5, 2.0, 2.5]

    Args: 
        * `a` (float-convertible): Range lower bound.
        * `b` (float-convertible): Range upper bound.
        * `N` (int-convertible):   Number of bins to split the range into.

    Returns:
        * `l` (list): List of floats that define the bin boundaries, len N+1.    
    """

    if not type(a) == float: 
        a = float(a)
    if not type(b) == float: 
        b = float(b)
    if not type(N) == int:
        N = int(N)

    l = [a]
    for i in range(1, N):
        l.append( a + (b-a)*float(i)/float(N) )
    l.append(b)
    return l

def LogUniformBins(a, b, N):
    """Splits the range [a,b] up into N log-uniform bins defined by N+1 points.
    For example::

        >>> LinearUniformBins(0.0001, 0.1, 3)
        [0.0001, 0.001, 0.01, 0.1]

    Args: 
        * `a` (float-convertible): Range lower bound.
        * `b` (float-convertible): Range upper bound.
        * `N` (int-convertible):   Number of bins to split the range into.

    Returns:
        * `l` (list): List of floats that define the bin boundaries, len N+1.    
    """

    if not type(a) == float: 
        a = float(a)
    if not type(b) == float: 
        b = float(b)
    if not type(N) == int:
        N = int(N)

    l = [a]
    for i in range(1, N):
        l.append( (a**(float(N-i)/float(N))) * (b**(float(i)/float(N)))  ) 
    l.append(b)
    return l

def NinesUniformBins(alpha, beta, N):
    """Splits the range [alpha, beta] up into N one-minus-log-uniform bins defined by N+1 points.
    In the vernacular, the range is 'split in the nines'.
    For example::

        >>> LinearUniformBins(0.9, 0.9999, 3)
        [0.9, 0.99, 0.999, 0.9999]

    Args: 
        * `alpha` (float-convertible): Range lower bound.
        * `beta`  (float-convertible): Range upper bound.
        * `N` (int-convertible): Number of bins to split the range into.

    Returns:
        * `l` (list): List of floats that define the bin boundaries, len N+1.    
    """

    if not type(alpha) == float: 
        alpha = float(alpha)
    if not type(beta) == float: 
        beta = float(beta)
    if not type(N) == int: 
        N = int(N)

    l = LogUniformBins(1.0 - alpha, 1.0 - beta, N)
    for i in range(N+1):
        l[i] = 1.0 - l[i]
    return l

def StairStep(xdat, ydat, N):
    """Makes data sets suitable for graphing as a stair step plot.

    Args:
        * `xdat` (sequence): Independent variable data, has length N+1.
        * `ydat` (sequence): Dependent variable data, has length N.
        * `N` (int): Number of stairs, or bins.

    Returns:
        * `xl` (list): Independent variable stair-step data.
        * `yl` (list): Dependent variable stair-step data.
    """

    xl = [xdat[0]]
    yl = [ydat[0]]
    for n in range(N-1):
        xl.append(xdat[n+1])
        yl.append(ydat[n])
        xl.append(xdat[n+1])
        yl.append(ydat[n+1])
    xl.append(xdat[N])
    yl.append(ydat[N-1])
    return xl, yl

########################
### Scale Exchangers ###
########################

def Nines2Log(dat, base=10.0):
    """Converts Nines Data ([0.9, 0.99, ...]) to form suitable for log-based graphing.

    Args:
        * `dat` (list): Nines data.

    Keyword Args:
        * `base` (float): Base of logarithm.

    Returns:
        * `l` (list): Nines data converted to log scale.
    """
    l = []
    for el in dat:
        l.append(base/(1.0 - el))
    return l

def Nines2Lin(dat, base=10.0):
    """Converts Nines Data ([0.9, 0.99, ...]) to form suitable for Linear-based graphing.

    Args:
        * `dat` (list): Nines data.

    Keyword Args:
        * `base` (float): Base of logarithm to convert through. (Note that nines data is inherently log-based.)

    Returns:
        * `lined` (list): Nines data converted to linear scale.
    """
    logged = base/(1.0 - dat)
    lined = np.log10(logged)
    return lined

def LogLabelConvert(xy, plot, base = 10.0):
    """Converts log labels to appropriate nines values.
    xy is a string that contains either or both 'x' and 'y' for each axis that needs to be converted.

    Args:
        * `xy` (str): A string that contains either 'x', 'y', or both for each axis that needs to be converted.
        * `plot` (plot): A matplotlib plot object, ie the 'plt' in 'import matplotlib.pyplot as plt'.

    Keyword Args:
        * `base` (float): Base of logarithm.
    """

    if 'x' in xy:
        xlocs, xlabels = plot.xticks()
        xlabels = []
        for n in xlocs:
            xlabels.append(1.0 - base/n)
        plot.xticks(xlocs, xlabels)
    if 'y' in xy:
        ylocs, ylabels = plot.yticks()
        ylabels = []
        for n in ylocs:
            ylabels.append( 1.0 - base/n )
        plot.yticks(ylocs, ylabels)
    if not ('x' in xy) and not ('y' in xy):
        print "No axes converted because neither \'x\' nor \'y\' appeared in the xy variable."

    return

#######################
### Other Functions ###
#######################

def SafeRemove(p, IsDir = False):
    """Trys to remove file(s), even if non-existent.
    
    Args:
        * `p` (str): Path to file or directorty.

    Keyword Args:
        * `IsDir` (True or False): If True, then recursively removes p, 
          as a directory.
    """
    try:
        if not IsDir:
            os.remove(p)
        else:
            for root, dirs, files in os.walk(p, topdown=False):
                for name in files:
                    os.remove(os.path.join(root, name))
                for name in dirs:
                    os.rmdir(os.path.join(root, name))
    except Exception as e:
        pass
        
    return

#Automatically imported modules
import data
