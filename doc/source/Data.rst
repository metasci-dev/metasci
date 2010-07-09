*********************
Scientific Data Model
*********************
Real world data always has errors and uncertainties associated with them.
Computationally managing such errors can become a headache, especially if
anything more than basic arithmetic must be performed on the data set.

To aid in data handling, the `metasci.data` module 
provides a new scientific numeric type that automatically propagates uncertainty
for all arithmetic operations.  The `data` class has two fundamental attributes:

   * `data.dat`: the raw data value
   * `data.sig`: the data value's uncertainty.

Arithmetic operations on other, non-`data` numeric types will always return a `data`-type object.  
For ease of use with data sets, a `vector` class is also provided whose elements are `data`-typed.

.. note::
   Certain private members are displayed explicitly in the documentation.  This is because arithmetic
   on uncertainty is non-obvious.  Therefore `__add__()` is displayed in order to be explicit about how
   the (+) operator works.

=============================================
:mod:`metasci.data` -- Scientific Data Module
=============================================

.. automodule:: metasci.data
   :members:
   :inherited-members:
   :exclude-members: data, vector

   .. autoclass:: data
      :members: 

      .. automethod:: __add__
      .. automethod:: __sub__
      .. automethod:: __mul__
      .. automethod:: __div__
      .. automethod:: __pow__
      .. automethod:: __rpow__

   .. autoclass:: vector
      :members: 

