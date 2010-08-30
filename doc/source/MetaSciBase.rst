***************************
Basic MetaSci Functionality
***************************
MetaSci provides a lot of functionality on its lowest level without the need to delve into sub-packages.
In fact, much of the content of the sub-packages relies on these top-level modules. 
Currently there are three modules in the base metasci package:

   * MetaSci itself (`__init__.py`)
   * A scientific data and data vector type (`data.py`)
   * A set of general plotting helper functions (`graph.py`) [probably broken]

These functions and modules may be accessed via::

   from metasci import *

=====================================
:mod:`metasci` -- MetaSci Base Module
=====================================

.. automodule:: metasci
   :members:
   :inherited-members:

====================================================
:mod:`metasci.colortext` -- MetaSci String Colorizer
====================================================

.. automodule:: metasci.colortext
   :members:
   :inherited-members:
