"""
Defines a new numeric 'data' object that handles associated uncertainty in scientific data automatically.
Additionally a 'vector' object made up of data points is also defined. 
"""

import math
from . import MultiSplit

class SciDataTypeError(Exception):
    """Error for handling when inputs to scientific data type cannot be converted into floats."""

    def __init__(self, dat, sig):
        """
        Args: 
            * `dat` (`float`): data value
            * `sig` (`float`): sigma value
        """
    
        self.dat = dat
        self.sig = sig
        return 

    def __str__(self):
        """Returns string representation of the exception."""
        s = "Input to data model could not be converted into float.\n"
        s = s + "\tData  Value: " + str(self.dat) + "\n"
        s = s + "\tSigma Value: " + str(self.sig)
        return s

class SciDataTypePowerError(Exception):
    """Error for handling when two scientific data points are raised to the power of one another."""

    def __init__(self):
        """Exception takes no arguments."""
        return 

    def __str__(self):
        """Returns string representation of the exception."""
        return  "Scientific Data Types Both Had Associated Error.\n" + \
            "Uncertainty could not be propagated when raising to the powers."

class SciVectorDataLengthError(Exception):
    """Error for when input data and uncertainty are of two different lengths."""

    def __init__(self, d):
        """
        Args:
            `d` (`list`): data list
        """
        self.d = d
        return 

    def __str__(self):
        """Returns string representation of the exception."""
        s =     "Vector could not be initialized.\n"
        s = s + "Data entry is of length: {0}\n".format(len(d))
        s = s + "For instance: {0}\n".format(d[0])
        s = s + "This fails to conform to the proper standard.\n"
        s = s + "dat\t= data list of one of the following types:\n"
        s = s + "\t[data, ...]\n"
        s = s + "\t[float || int || str, ...]\n"
        s = s + "\t[[float || int || str, float || int || str], ...]"
        return s

class SciVectorDataUncertaintyLengthError(Exception):
    """Error for when input data and uncertainty are of two different lengths."""

    def __init__(self, dat, sig):
        """
        Args:
            * `dat` (`list`): data list
            * `sig` (`list`): sigma list
        """
        self.dat = dat
        self.sig = sig
        return

    def __str__(self):
        """Returns string representation of the exception."""
        s =     "Vector could not be initialized.\n"
        s = s + "Data list is of length: {0}\n".format(len(self.dat))
        s = s + "While the uncertainty list is of length:".format(len(self.sig))
        return s

class SciVectorBadLength(Exception):
    """Error for when vector types have different lengths."""

    def __init__(self):
        """Exception takes no arguments."""
        return 

    def __str__(self):
        """Returns string representation of the exception."""
        return "Vectors have different lengths which make math impossible."

class data:
    """
    Data object for scientific computing that handles scientific data and an associated error.
    Standard math functions on this object perform appropriate corresponding calculations on error data.
    All math functions assume that the data is uncorrelated at this point in time.
    (That is to say that all the covariances are zero.)
    """

    dat = None
    """The data value of the `data` object."""
    sig = None
    """The uncertainty associated with the `data` object."""

    def __init__(self, dat=0.0, sig=0.0):
        """
        Args:
            * `dat` (`float`-convertible, data): data value
            * `sig` (`float`-convertible): uncertainty value
        """

        if isinstance(dat, str):
            l = MultiSplit(dat, ",:;[]{}()<>?!~$%^*\ \t\n")
            if len(l) == 1:
                dat = l[0]
            elif 1 < len(l):
                dat = l[0]
                sig = l[1]
        elif isinstance(dat, data):
            self.dat = dat.dat
            self.sig = dat.sig
            return
    
        try:
            self.dat = float( dat )
            self.sig = float( sig )
        except:
            raise SciDataTypeError(dat, sig)
        
    def __str__(self):
        """Returns string representation of the data object."""
        return "{0:G} +/- {1:G}".format(self.dat, self.sig)

    def write(self):
        """
        Similar to __str__, this produces a string of the data.
        However, this string is suitable for writing to an output file.
        The number after the decimal point is the precision (currently 6).
        """
        return "{0:.6E}\t{1:.6E}".format(self.dat, self.sig)

    def __float__(self):
        """Returns the data dat value.  **USE WITH CAUTION!**"""
        return self.dat

    def __getitem__(self, key):
        """
        Enables key-based access.

        Args:
            `key` (`str`): May be either 'dat' or 'sig'.

        Returns:
                    (`float`) of the appropriate attribute.
        """

        if key == 'dat':
            return self.dat
        elif key == 'sig':
            return self.sig
        else:
            return KeyError

    def __lt__(self, other):
        """
        "Less Than" (<) comparison for data object to another numerical type.

        Args:
                    `other` (`data`-convertible): Another numerical type to compare to.

        Returns:
                    (`True` or `False`)
        """

        if not isinstance(other, data):
            other = data(other)

        if self.dat < other.dat:
            return True
        else:
            return False

    def __le__(self, other):
        """
        "Less Than or Equal To" (<=) comparison for data object to another numerical type.

        Args:
                    `other` (`data`-convertible): Another numerical type to compare to.

        Returns:
                    (`True` or `False`)
        """

        if not isinstance(other, data):
            other = data(other)

        if self.dat <= other.dat:
            return True
        else:
            return False

    def __eq__(self, other):
        """
        "Equal To" (==) comparison for data object to another numerical type.

        Args:
                    `other` (`data`-convertible): Another numerical type to compare to.

        Returns:
                    (`True` or `False`)
        """

        if not isinstance(other, data):
            other = data(other)

        if (self.dat == other.dat) and (self.sig == other.sig):
            return True
        else:
            return False

    def __ne__(self, other):
        """
        "Not Equal To" (!=) comparison for data object to another numerical type.

        Args:
                   `other` (`data`-convertible): Another numerical type to compare to.

        Returns:
                    (`True` or `False`)
        """

        if not isinstance(other, data):
            other = data(other)

        if (self.dat != other.dat) and (self.sig != other.sig):
            return True
        else:
            return False
        
    def __gt__(self, other):
        """
        "Greater Than" (>) comparison for data object to another numerical type.

        Args:
                    `other` (`data`-convertible): Another numerical type to compare to.

        Returns:
                    (`True` or `False`)
        """

        if not isinstance(other, data):
            other = data(other)

        if self.dat > other.dat:
            return True
        else:
            return False

    def __ge__(self, other):
        """
        "Greater Than or Equal To" (>=) comparison for data object to another numerical type.

        Args:
                    `other` (`data`-convertible): Another numerical type to compare to.

        Returns:
                    (`True` or `False`)
        """

        if not isinstance(other, data):
            other = data(other)

        if self.dat >= other.dat:
            return True
        else:
            return False

    def __add__(self, other):
        """
        Data & uncertainty addition.

        Args:
            `other` (`data`-convertible): Another numerical type to add with.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = dat_a + dat_b

               \sigma_c = \sqrt{\sigma_a^2 + \sigma_b^2 }
        """

        if isinstance(other, vector):
            return NotImplemented
        if not isinstance(other, data):
            other = data(other)
            
        return data( self.dat + other.dat, math.sqrt(self.sig**2 + other.sig**2) )

    def __sub__(self, other):
        """
        Data & uncertainty subtraction.

        Args:
            `other` (`data`-convertible): Another numerical type to subtract with.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = dat_a - dat_b

               \sigma_c = \sqrt{\sigma_a^2 + \sigma_b^2 }
        """

        if isinstance(other, vector):
            return NotImplemented
        if not isinstance(other, data):
            other = data(other)
            
        return data( self.dat - other.dat, math.sqrt(self.sig**2 + other.sig**2) )

    def __mul__(self, other):
        """
        Data & uncertainty multiplication.

        Args:
            `other` (`data`-convertible): Another numerical type to multiply with.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = dat_a \cdot dat_b

               \sigma_c = \sqrt{ (dat_b \cdot \sigma_a)^2 + (dat_a \cdot \sigma_b)^2 }
        """

        if isinstance(other, vector):
            return NotImplemented
        if not isinstance(other, data):
            other = data(other)
            
        return data( self.dat * other.dat, math.sqrt( (other.dat * self.sig)**2 + (self.dat * other.sig)**2) )

    def __div__(self, other):
        """
        Data & uncertainty division.

        Args:
            `other` (`data`-convertible): Another numerical type to divide with.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = \\frac{dat_a}{dat_b}

               \sigma_c = \sqrt{ \left(\\frac{\sigma_a}{dat_b}\\right)^2 + \left(\\frac{dat_a \cdot \sigma_b}{dat_b^2}\\right)^2 }
        """

        if isinstance(other, vector):
            return NotImplemented
        if not isinstance(other, data):
            other = data(other)
            
        return data( self.dat / other.dat, math.sqrt( (self.sig / other.dat)**2 + (self.dat * other.sig / (other.dat**2) )**2) )

    def __floordiv__(self, other):
        """
        Data & uncertainty floor division.

        Args:
            `other` (`data`-convertible): Another numerical type to divide with.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = \left\lfloor\\frac{dat_a}{dat_b}\\right\\rfloor

               \sigma_c = \sqrt{ \left(\\frac{\sigma_a}{dat_b}\\right)^2 + \left(\\frac{dat_a \cdot \sigma_b}{dat_b^2}\\right)^2 }
        """

        if isinstance(other, vector):
            return NotImplemented
        if not isinstance(other, data):
            other = data(other)
            
        return data( self.dat // other.dat, math.sqrt( (self.sig / other.dat)**2 + (self.dat * other.sig / (other.dat**2) )**2) )

    def __pow__(self, other):
        """
        Data & uncertainty power operator.
        Used only when other is not of data type or has no associated error.

        Args:
            `other` (non-data or other.sig = 0.0): Another (non-data) numerical type to raise `self` to the power of.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = (dat_a)^{dat_b}

               \sigma_c = |dat_b| \cdot (dat_a)^{dat_b-1} \cdot \sigma_a
        """

        if isinstance(other, vector):
            return NotImplemented
        if not isinstance(other, data):
            other = float(other)
        elif other.sig == 0.0:
            other = other.dat
        else:
            raise SciDataTypePowerError
            
        return data( self.dat**other, abs( other ) * self.dat**(other - 1) * self.sig )

    def __radd__(self, other):
        """
        Right data & uncertainty addition is commutative.  Calls `__add__()`.

        Args:
            `other` (`data`-convertible): Another numerical type to add with.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = dat_a + dat_b

               \sigma_c = \sqrt{\sigma_a^2 + \sigma_b^2 }
        """
        return self.__add__(other)

    def __rsub__(self, other):
        """
        Right data & uncertainty subtraction.

        Args:
            `other` (`data`-convertible): Another numerical type to subtract with.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = dat_b - dat_a

               \sigma_c = \sqrt{\sigma_a^2 + \sigma_b^2 }
        """

        if not isinstance(other, data):
            other = data(other)
            
        return data( other.dat - self.dat, math.sqrt(self.sig**2 + other.sig**2) )

    def __rmul__(self, other):
        """
        Right data & uncertainty multiplication is commutative! Calls `__mul__()`.

        Args:
            `other` (`data`-convertible): Another numerical type to multiply with.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = dat_a \cdot dat_b

               \sigma_c = \sqrt{ (dat_b \cdot \sigma_a)^2 + (dat_a \cdot \sigma_b)^2 }
        """
        return self.__mul__(other)

    def __rdiv__(self, other):
        """
        Right data & uncertainty division.

        Args:
            `other` (`data`-convertible): Another numerical type to divide with.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = \\frac{dat_b}{dat_a}

               \sigma_c = \sqrt{ \left(\\frac{\sigma_b}{dat_a}\\right)^2 + \left(\\frac{dat_b \cdot \sigma_a}{dat_a^2}\\right)^2 }
        """

        if not isinstance(other, data):
            other = data(other)
            
        return data( other.dat / self.dat, math.sqrt( (other.sig / self.dat)**2 + (other.dat * self.sig / (self.dat**2) )**2) )

    def __rfloordiv__(self, other):
        """
        Data & uncertainty floor division.

        Args:
            `other` (`data`-convertible): Another numerical type to divide with.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = \left\lfloor\\frac{dat_b}{dat_a}\\right\\rfloor

               \sigma_c = \sqrt{ \left(\\frac{\sigma_b}{dat_a}\\right)^2 + \left(\\frac{dat_b \cdot \sigma_a}{dat_a^2}\\right)^2 }
        """

        if not isinstance(other, data):
            other = data(other)
            
        return data( other.dat // self.dat, math.sqrt( (other.sig / self.dat)**2 + (other.dat * self.sig / (self.dat**2) )**2) )

    def __rpow__(self, other):
        """
        Data & uncertainty power operator.
        Used only when other is not of data type or has no associated error.

        Args:
            `other` (non-data or other.sig = 0.0): Another (non-data) numerical type to raise `self` to the power of.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = (dat_b)^{dat_a}

               \sigma_c = (dat_b)^{dat_a} \cdot \ln(dat_b) \cdot \sigma_a
        """

        if not isinstance(other, data):
            other = float(other)
        elif other.sig == 0.0:
            other = other.dat
        else:
            raise SciDataTypePowerError
            
        return data( other**self.dat, (other**self.dat) * math.log(other) * self.sig )

    def __neg__(self):
        """
        Data & uncertainty negative operator.
        The data value's sign is changed while the uncertainty remains positive.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = - dat_a

               \sigma_c = \sigma_a
        """
        return data( - self.dat, self.sig )

    def __pos__(self):
        """
        Data & uncertainty positive operator.
        The data value's sign is changed while the uncertainty remains positive.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = + dat_a

               \sigma_c = \sigma_a
        """
        return data( + self.dat, self.sig )

    def __abs__(self):
        """
        Data & uncertainty absolute value operator.
        The data value's sign is forced positive while the uncertainty remains positive.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = |dat_a|

               \sigma_c = \sigma_a
        """
        return data( abs( self.dat ), self.sig )

    def log(self, base = math.e):
        """
        Data & uncertainty logarithm operator.

        Args:
            `base` (int or float): base of log.

        Returns:
                    `data_c` (`data`): New data object calculated from, 

            .. math::
               dat_c = \log_b( dat_a )

               \sigma_c = \\frac{\sigma_a}{|\ln(b)| \cdot dat_a}
        """
        return data( math.log( self.dat, base ), self.sig / ( abs(math.log(base)) * self.dat) )

class vector:
    """
    This is a vector class that amounts to a list of `data` types.
    Thus arithmetic on these vectors will automatically carry through all uncertainties.
    """

    vec = None
    """The vector list that contains all of the `data` objects."""

    def __init__(self, dat = [], sigma = [], zeros = 0):
        """
        Args:
            * `dat` (`list`): data list of one of the following forms:
               * [data, ...]
               * [float || int || str, ...]
               * [[float || int || str, float || int || str], ...]

            * `sigma` (`list`): list of uncertainties [float || int || str, ...]. 
                      Only use if dat is of [float || int || str, ...] form as well.
            * `zeros` (`int`): number of empty data points to append to list.

        .. warning:: 
           Vector does not check for consistency in the input parameters. So Be CAREFUL!
        """

        if len(dat) == 0:
            self.vec = []
        elif isinstance(dat[0], data):
            self.vec = dat
        elif not isinstance(dat[0], list):
            self.vec = []
            if len(sigma) == 0:
                for d in dat:
                    self.vec.append( data(d) )
            elif len(dat) == len(sigma):
                for n in range(len(dat)):
                    self.vec.append( data(dat[n], sigma[n]) )
            else:
                raise SciVectorDataUncertaintyLengthError(dat, sigma)
        elif len(dat[0]) == 2:
            self.vec = []
            for d in dat:
                self.vec.append( data(d[0], d[1]) )
        else:
            raise SciVectorDataLengthError(dat)

        for n in range(zeros):
            self.vec.append( data() )
            
    def __str__(self):
        """Returns string representation of vector."""
        s = "["
        for n in range(len(self.vec)):
            if n == 0:
                s = s + str(self.vec[n])
            else:
                s = s + ", " + str(self.vec[n])
        s = s + "]"
        return s

    def __len__(self):
        """Returns the number of elements in vector object."""
        return len(self.vec)

    def __get__(self):
        """Returns the list of elements of the vector object."""
        return self.vec

    def __getitem__(self, key):
        """Returns an element of the vector object."""
        return self.vec[key]

    def __setitem__(self, key, value):
        """Sets an element of the vector object."""
        self.vec[key] = value
        return 

    def __delitem__(self, key):
        """Deletes an element of the vector object."""
        del self.vec[key]
        return

    def __iter__(self):
        """Returns an iterator over the vector object."""
        return self.vec.__iter__()

    def __reversed__(self):
        """Returns a reversed list of elements of the vector object."""
        return reversed( self.vec )

    def __contains__(self, item):
        """Returns a boolean on the membership of `item` to the vector object."""
        if not isinstance(item, data):
            item = data(item)

        for n in self:
            if n == item:
                return True
        return False

    def __add__(self, other):
        """Adds a vector, list, data, or other numerical type to this list."""
        if isinstance(other, vector):
            if not (len(self) == len(other)):
                raise SciVectorBadLength
            l = []
            for n in range(len(self)):
                l.append( self.vec[n] + other.vec[n] )
            return vector(l) 
        elif isinstance(other, list):
            return self.__add__( vector(other) )
        elif isinstance(other, data):
            l = []
            for n in range(len(self)):
                l.append( self.vec[n] + other )
            return vector(l) 
        else:
            return self.__add__( data(other) ) 

    def __sub__(self, other):
        """Subtracts a vector, list, data, or other numerical type from this list."""
        if isinstance(other, vector):
            if not (len(self) == len(other)):
                raise SciVectorBadLength
            l = []
            for n in range(len(self)):
                l.append( self.vec[n] - other.vec[n] )
            return vector(l) 
        elif isinstance(other, list):
            return self.__sub__( vector(other) )
        elif isinstance(other, data):
            l = []
            for n in range(len(self)):
                l.append( self.vec[n] - other )
            return vector(l) 
        else:
            return self.__sub__( data(other) ) 

    def __mul__(self, other):
        """Multiplies a vector, list, data, or other numerical type to this list."""
        if isinstance(other, vector):
            if not (len(self) == len(other)):
                raise SciVectorBadLength
            l = []
            for n in range(len(self)):
                l.append( self.vec[n] * other.vec[n] )
            return vector(l) 
        elif isinstance(other, list):
            return self.__mul__( vector(other) )
        elif isinstance(other, data):
            l = []
            for n in range(len(self)):
                l.append( self.vec[n] * other )
            return vector(l) 
        else:
            return self.__mul__( data(other) ) 

    def __div__(self, other):
        """Divides a vector, list, data, or other numerical type from this list."""
        if isinstance(other, vector):
            if not (len(self) == len(other)):
                raise SciVectorBadLength
            l = []
            for n in range(len(self)):
                l.append( self.vec[n] / other.vec[n] )
            return vector(l) 
        elif isinstance(other, list):
            return self.__div__( vector(other) )
        elif isinstance(other, data):
            l = []
            for n in range(len(self)):
                l.append( self.vec[n] / other )
            return vector(l) 
        else:
            return self.__div__( data(other) ) 

    def __floordiv__(self, other):
        """Floor divides a vector, list, data, or other numerical type from this list."""
        if isinstance(other, vector):
            if not (len(self) == len(other)):
                raise SciVectorBadLength
            l = []
            for n in range(len(self)):
                l.append( self.vec[n] // other.vec[n] )
            return vector(l) 
        elif isinstance(other, list):
            return self.__floordiv__( vector(other) )
        elif isinstance(other, data):
            l = []
            for n in range(len(self)):
                l.append( self.vec[n] // other )
            return vector(l) 
        else:
            return self.__floordiv__( data(other) ) 

    def __pow__(self, other):
        """Raises this vector to the power of another vector, list, data, or other numerical type."""
        if isinstance(other, vector):
            if not (len(self) == len(other)):
                raise SciVectorBadLength
            l = []
            for n in range(len(self)):
                l.append( self.vec[n] ** other.vec[n] )
            return vector(l) 
        elif isinstance(other, list):
            return self.__pow__( vector(other) )
        elif isinstance(other, data):
            l = []
            for n in range(len(self)):
                l.append( self.vec[n] ** other )
            return vector(l) 
        else:
            return self.__pow__( data(other) ) 

    def __radd__(self, other):
        """Right adds a vector, list, data, or other numerical type to this list.  Commutative, so calls `__add__()`."""
        return self.__add__(other)

    def __rsub__(self, other):
        """Right subtracts a vector, list, data, or other numerical type from this list."""
        if isinstance(other, vector):
            if not (len(self) == len(other)):
                raise SciVectorBadLength
            l = []
            for n in range(len(self)):
                l.append( other.vec[n] - self.vec[n] )
            return vector(l) 
        elif isinstance(other, list):
            return self.__rsub__( vector(other) )
        elif isinstance(other, data):
            l = []
            for n in range(len(self)):
                l.append( other - self.vec[n] )
            return vector(l) 
        else:
            return self.__rsub__( data(other) ) 

    def __rmul__(self, other):
        """Right multiplies a vector, list, data, or other numerical type to this list.  Commutative, so calls `__mul__()`."""
        return self.__mul__(other)

    def __rdiv__(self, other):
        """Right divides a vector, list, data, or other numerical type from this list."""
        if isinstance(other, vector):
            if not (len(self) == len(other)):
                raise SciVectorBadLength
            l = []
            for n in range(len(self)):
                l.append( other.vec[n] / self.vec[n] )
            return vector(l) 
        elif isinstance(other, list):
            return self.__rdiv__( vector(other) )
        elif isinstance(other, data):
            l = []
            for n in range(len(self)):
                l.append( other / self.vec[n] )
            return vector(l) 
        else:
            return self.__rdiv__( data(other) ) 

    def __rfloordiv__(self, other):
        """Right floor divides a vector, list, data, or other numerical type from this list."""
        if isinstance(other, vector):
            if not (len(self) == len(other)):
                raise SciVectorBadLength
            l = []
            for n in range(len(self)):
                l.append( other.vec[n] // self.vec[n] )
            return vector(l) 
        elif isinstance(other, list):
            return self.__rfloordiv__( vector(other) )
        elif isinstance(other, data):
            l = []
            for n in range(len(self)):
                l.append( other // self.vec[n] )
            return vector(l) 
        else:
            return self.__rfloordiv__( data(other) ) 

    def __rpow__(self, other):
        """Raises a vector, list, data, or other numerical type to the power of this list."""
        if isinstance(other, vector):
            if not (len(self) == len(other)):
                raise SciVectorBadLength
            l = []
            for n in range(len(self)):
                l.append( other.vec[n] ** self.vec[n] )
            return vector(l) 
        elif isinstance(other, list):
            return self.__rpow__( vector(other) )
        elif isinstance(other, data):
            l = []
            for n in range(len(self)):
                l.append( other ** self.vec[n] )
            return vector(l) 
        else:
            return self.__rpow__( data(other) ) 

    def __neg__(self):
        """Negates all elements of this vector."""
        l = []
        for n in range(len(self)):
            l.append( - self.vec[n] )
        return vector(l) 

    def __pos__(self):
        """Positives all elements of this vector."""
        l = []
        for n in range(len(self)):
            l.append( + self.vec[n] )
        return vector(l) 

    def __abs__(self):
        """Forces positive all elements of this vector."""
        l = []
        for n in range(len(self)):
            l.append( abs( self.vec[n] ) )
        return vector(l) 

    def log(self, base = math.e):
        """
        Computes the logarithm of all elements of this vector.

        Args:
            `base` (`int` or `float`): base of logarithm.

        Returns:
            New `vector` object.
        """
        l = []
        for n in range(len(self)):
            l.append( self.vec[n].log(base)  )
        return vector(l) 

    def dot(self, other):
        """
        Computes the dot product of this vector with another vector or list.

        Args:
            `other` (`vector` or `list`): Another vector to dot with this one.

        Returns:
            (`data`): Result of dot product.
        """
        if isinstance(other, vector):
            if not (len(self) == len(other)):
                raise SciVectorBadLength
            d = data()
            for n in range(len(self)):
                d = d + (self.vec[n] * other.vec[n])
            return d
        else:
            return self.dot( vector(other) )

