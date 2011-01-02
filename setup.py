#!/usr/bin/env python
 
from distutils.core import setup
from distutils.extension import Extension
from setuptools.command.develop import develop


setup(name="metasci",
	version='0.23',
	description='MetaSci',
	author='Anthony Scopatz',
	author_email='scopatz@gmail.com',
	url='http://www.scopatz.com/',
	packages=['metasci', 'metasci.chem', 'metasci.mathematics', 'metasci.nuke', 'metasci.stats'],
	package_dir={'metasci': 'metasci'}, 
	package_data={'metasci': ['nuke/nuc_data.h5']},
	)

