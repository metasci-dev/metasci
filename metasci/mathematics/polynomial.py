"""
Classes for handling N-dimensional polynomials of order J.
Polynomials are split into a coefficient class and a functional class.
The first class is only responsible for the heady task of 
coefficient management.  The second class calculates the polynomial results at any 
well-defined point. The function class carries a coefficient instance as an attribute.
"""

from itertools import permutations

class polyNdcoef(object):
    """
    N-dimensional Polynomial Coefficients known to Order J.
    For example with N=3 and J=1:

       f(x, y, z) = p[0,0,0] + p[1,0,0] x + p[0,1,0] y + p[0,0,1] z

    This object is a storage container for the coefficients p[].
    The coefficients require their own container to ease the 
    distinction between coefficient access and polynomial function access.

    A specific coefficient may be accessed via an N-tuple key::

       pc = polyNdcoef(N=2, J=3) #Two dimensional coefficients
       pc[0,0]                   #The zeroth term
       pc[1,0]                   #First order in x, zero order in y
       pc[1,2]                   #First order in x, second order in y

    Additionally, if 1 < N, elements may be accessed by single integer. 
    This integer relates to the position of an N-tuple key in polyNdcoef.indices::

       pc = polyNdcoef(N=2, J=3) #Two dimensional coefficients
       pc.indices = [(0, 0), (0, 1), (1, 0), (2, 0), (0, 2), (1, 1), (1, 2), ...]
       pc[0] = pc[0,0]           #The zeroth term
       pc[2] = pc[1,0]           #First order in x, zero order in y
       pc[6] = pc[1,2]           #First order in x, second order in y

    This notation becomes ambiguous in the 1D case and is not allowed.
    """

    N = None
    """Number of independent variables."""
    J = None
    """Polynomial order."""
    indices = None
    """
    List of integer tuples that are valid polynomial indices.  
    For each order j less than or equal to J, the sum of the tuple values 
    for that index must exactly equal j.  Examples, 

       * For N = 1, J = 3, then indices = [(0), (1), (2), (3), ...]
       * For N = 2, J = 3, then indices = [(0, 0), (0, 1), (1, 0), (2, 0), (0, 2), (1, 1), ...]
       * For N = 3, J = 3, then indices = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), ...]
    """

    def __init__(self, N=1, J=0, default=0.0):
        """
        Keyword Args:
           * `N` (int): The dimensionality or number of independent variables.  E.g. f(x,y) has N = 2.
           * `J` (int): Polynomial coefficient order.
           * `default` (numeric type): The default value to fill the coefficients with.

        """

        self.N = N		#Number of independent variables 
        self.J = J		#Polynomial order
        self.indices = []	#Valid polynomial index list

        #Build coefficients set to default value		
        basis = [0 for n in range(N)]
        for j in range(0, J+1):
            basis = basis + range(0, j+1)
            perms = permutations(basis, N)
            for i in set([p for p in perms if sum(p) == j]):
                self.indices.append(i)
                self.__dict__[i] = default
                
        return

    def __len__(self):
        """Returns the length of the polynomial."""
        return len(self.indices)

    def __getitem__(self, key):
        """Returns the value associated with the key."""
        if (self.N == 1) and (type(key) == int):
            key = tuple([key])
        elif type(key) == int:
            key = self.indices[key]

        if not (len(key) == self.N):
            raise TypeError
        elif (sum(key) < 0.0) and (self.J < sum(key)):
            raise IndexError
        elif not (key in self.__dict__):
            raise KeyError
        else:
            return self.__dict__[key]

    def __setitem__(self, key, value):
        """
        Returns the value associated with the key.
        Only allows instantiated coefficients to be set.
        """
        if (self.N == 1) and (type(key) == int):
            key = tuple([key])
        elif type(key) == int:
            key = self.indices[key]

        if not (len(key) == self.N):
            raise TypeError
        elif (sum(key) < 0.0) and (self.J < sum(key)):
            raise IndexError
        elif not (key in self.__dict__):
            raise KeyError
        else:
            self.__dict__[key] = value

    def __contains__(self, key):
        """Returns whether key is a valid index."""
        if key in self.indices:
            return True
        else:
            return False

    def keys(self, j=None):
        """
        Returns the coefficient indices, which act as keys.

        Keyword Args:
           * `j` (int): Polynomial order such that  0 <= j <= J.

        Returns:
           * A (list) of the polynomial indices.
           * If `j` is defined, returns only the keys of order j.
        """

        if j == None:
            return self.indices
        elif type(j) == int:
            if 0 <= j <= self.J:
                return [k for k in self.indices if sum(k) == j]
            else:
                raise IndexError
        else:
            raise TypeError

    def values(self, j=None):
        """
        Returns the coefficient values associated with the indices.

        Keyword Args:
           * `j` (int): Polynomial order such that  0 <= j <= J.

        Returns:
           * A (list) of the polynomial coefficients.
           * If `j` is defined, returns only the coefficients of order j.
        """

        if j == None:
            return [self.__dict__[k] for k in self.indices]
        elif type(j) == int:
            if 0 <= j <= self.J:
                return [self.__dict__[k] for k in self.indices if sum(k) == j]
            else:
                raise IndexError
        else:
            raise TypeError

    def items(self, j=None):
        """
        Returns the coefficient index-value pairs (index, p[index]). 

        Keyword Args:
           * `j` (int): Polynomial order such that  0 <= j <= J.

        Returns:
           * A (list) of the polynomial indices and coefficients.
           * If `j` is defined, returns only the indices and coefficients of order j.
        """
        if j == None:
            return [(k, self.__dict__[k]) for k in self.indices]
        elif type(j) == int:
            if 0 <= j <= self.J:
                return [(k, self.__dict__[k]) for k in self.indices if sum(k) == j]
            else:
                raise IndexError
        else:
            raise TypeError

    def has_key(self, key):
        """Returns whether key is a valid index."""
        return self.__contains__(key)

    def get(self, key):
        """Tries to get the key, otherwise returns None."""
        try:
            return self.__getitem__[key]
        except:
            return None

    def __str__(self):
        """Returns string representation of polynomial coefficients."""

        s = "{0}D Polynomial of Order J = {1}:\n\n".format(self.N, self.J)
        for j in range(0, self.J+1):
            s = s + "Terms for j = {0}:\n".format(j)
            for i in self.keys(j):
                s = s + "{0} = {1}\n".format(i, self[i])
            s = s + "\n"
        return s

    def LaTeX(self):
        """Returns LaTeX string representation of polynomial coefficients."""

        s =     "\\begin{center}\n"
        s = s + "\\begin{table}[htbp]\n"
        s = s + "\\caption{{{0}D Polynomial of Order J = {1}}}\n".format(self.N, self.J)
        s = s + "\\begin{center}\n"
        s = s + "\\begin{tabular}{|c|c|}\n"
        s = s + "\\hline\n"
                
        for j in range(0, self.J+1):
            s = s + "\multicolumn{{2}}{{|c|}}{{Coefficients of Order j = {0}}}\\\\\n".format(j)
            s = s + "\\hline\n"
            for i in self.keys(j):
                s = s + "{0} & {1}\\\\\n".format(i, self[i])
                s = s + "\\hline\n"

        s = s + "\\end{tabular}\n"
        s = s + "\\end{center}\n"
        s = s + "\\end{table}\n"
        s = s + "\\end{center}\n\n"

        return s

class polyNd(object):
    """
    N-dimensional Polynomial of Order J.
    For example with N=3 and J=1:

       f(x, y, z) = p[0,0,0] + p[1,0,0] x + p[0,1,0] y + p[0,0,1] z

    This is a true representation of the polynomial in that the class 
    is callable, like a function.  Therefore given any (x, y, z), this 
    class will calculate f in the above example.  Naturally, this requires 
    that the polynomial function here be initialized with a polynomial 
    coefficient data set polyNdcoef.
    """

    N = None
    """Number of independent variables."""
    J = None
    """Polynomial order."""
    coef = None
    """Polynomial coefficients for this polynomial object."""
    labels = None
    """String list of labels for printing, length `N`."""

    def __init__(self, N=None, J=None, coef=None, labels=["x"]):
        """
        Keyword Args:
           * `N` (int): The dimensionality or number of independent variables.  E.g. f(x,y) has N = 2.
           * `J` (int): Polynomial order.
           * `coef` (polyNdcoef or None): N-dimensional, Jth order polynomial coefficients.
           * `labels` (list): String list of labels for printing, length N.
        """

        #set N
        if N == None:
            self.N = 1
        else:
            self.N = int(N)

        #set J
        if J == None:
            self.J = 1
        else:
            self.J = int(J)

        #set coefficients
        if coef == None:
            self.coef = polyNdcoef(N=self.N, J=self.J)
        elif type(coef) == polyNdcoef:
            if (N == None) and (J == None):	
                self.N = coef.N
                self.J = coef.J
                self.coef = coef
            elif (N == coef.N):
                self.coef = polyNdcoef(N=self.N, J=self.J)
                for j in range(0, min(J, coef.J)+1):
                    for key in self.coef.keys(j):
                        self.coef[key] = coef[key]
            else:
                raise IndexError("Number of polynomial independent variables N={0} does not equal the number of coefficient independent variables N={1}".format(self.N, coef.N))
        elif type(coef) == float:
            self.coef = polyNdcoef(N=self.N, J=self.J, default=coef)
        elif type(coef) == int:
            self.coef = polyNdcoef(N=self.N, J=self.J, default=float(coef))
        else:
            raise TypeError

        self.labels = labels

        return

    def __str__(self):
        """Returns a string representation of the polynomial function."""
        s = "f(" + self.labels[0]
        for n in range(1, len(self.labels)):
            s = s + ", " + self.labels[n]
        s = s + ") = {0:.6G}".format(self.coef.values(0)[0])

        for j in range(1, self.J+1):
            jtems = self.coef.items(j)
            for jtem in jtems:
                if jtem[1] < 0.0:
                    s = s + "  -  "
                else:
                    s = s + "  +  "
                s = s + "({0:.6G}".format(abs(jtem[1]))
                for n in range(self.N):
                    if 0 < jtem[0][n]:
                        s = s + " {0}^{1}".format(self.labels[n], jtem[0][n])
                s = s + ")"

        return s

    def LaTeX(self):
        """Returns a LaTeX string representation of the polynomial function."""
        s = "f(" + self.labels[0]
        for n in range(1, len(self.labels)):
            s = s + ", " + self.labels[n]
        s = s + ") = {0:.6G}".format(self.coef.values(0)[0])

        for j in range(1, self.J+1):
            jtems = self.coef.items(j)
            for jtem in jtems:
                if jtem[1] < 0.0:
                    s = s + "  -  "
                else:
                    s = s + "  +  "
                s = s + "({0:.6G}".format(abs(jtem[1]))
                for n in range(self.N):
                    if 0 < jtem[0][n]:
                        s = s + " {0}^{{{1}}}".format(self.labels[n], jtem[0][n])
                s = s + ")"

        return s

    def __call__(self, *args, **kwargs):
        """
        Allows the polynomial class to be called as a function (which it totally is)!
        This enables the following behavior::

           p = polyNd(N=2, coef=...)
           result = p(1, 2)

        which evaluates a polynomial object at the x-y point (1,2).

        Args:
           * `args[0:N]` (numeric-type): The first N arguments are treated as the point value 
             at which to evaluate the function.
           * `args[N]` (polyNdcoef or list): If an optional N+1 argument is given, it is 
             treated as a new set of coefficients to overwrite the current ones with before 
             evaluating the function.

        Keyword Args:
           * `coef` (polyNdcoef or list): If an optional N+1 argument is given, it is 
             treated as a new set of coefficients to overwrite the current ones with before 
             evaluating the function.

        Returns:
           * `result` (float): The polynomial function evaluated at the point it was called.
        """
        if len(args) < self.N:
            raise IndexError

        #point to evaluate is the first N arguments called
        p = args[:self.N]

        #If it exists, try assigning the next argument as the coefficient set.
        if self.N < len(args):
            if type(args[self.N]) == polyNdcoef:
                self.coef = args[self.N]
            elif type(args[self.N]) == list:
                if len(self.coef) != len(args[self.N]):
                    raise IndexError
                for i in range(len(args[self.N])):
                    self.coef[i] = args[self.N][i]

        #If it exists, try assigning the keyword argument 'coef' as the coefficient set.
        if 'coef' in kwargs.keys():
            if type(kwargs['coef']) == polyNdcoef:
                self.coef = kwargs['coef']
            elif type(kwargs['coef']) == list:
                if len(self.coef) != len(kwargs['coef']):
                    raise IndexError
                for i in range(len(kwargs['coef'])):
                    self.coef[i] = kwargs['coef'][i]
    
        result = 0.0
        for j in range(0, self.J+1):
            jtems = self.coef.items(j)
            for jtem in jtems:
                term = jtem[1]
                for n in range(self.N):
                    if 0 < jtem[0][n]:
                        term = term * (p[n]**jtem[0][n])
                result = result + term
        return result
