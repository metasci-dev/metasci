===============
Nuclear Library
===============

The following two modules server to aggregate a variety of nuclear data needs, mostly from different 
websources.  Having this data, an HDF5 database containing some or all of this data may be made for 
easy future use.

The code below is a recipe that uses the two modules to make a database::

    from metasci.nuke import data_grab
    from metasci.nuke import data_make

    # Grab atomic weight information
    data_grab.grab_kaeri_atomic_weights(file_out='atomic_weight.txt')

    # Grab the one-group neutron cross-sections from KAERI
    with open('atomic_weight.txt', 'r') as f:
        nuclist = [line.split()[0] for line in f]
    data_grab.grab_kaeri_neutron_xs(nuclist, dir_out='xs_html/')

    # Make atomic weight data table
    data_make.make_atomic_weight(data_file='atomic_weight.txt')

    # Make one-group data tables
    data_make.make_xs_1g(data_dir='xs_html/')

-------------------------------------------------------------
:mod:`metasci.nuke.data_grab` -- Nuclear Data Grabber Package
-------------------------------------------------------------

.. automodule:: metasci.nuke.data_grab
   :members:


-------------------------------------------------------------------
:mod:`metasci.nuke.data_make` -- Nuclear Database Generator Package
-------------------------------------------------------------------

.. automodule:: metasci.nuke.data_make
   :members:

