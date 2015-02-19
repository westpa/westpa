===============
WESTPA 1.0 beta
===============


--------
Overview
--------

WESTPA is a package for constructing and running stochastic simulations using the "weighted ensemble" approach 
of Huber and Kim (1996). For further information, please see the following:

Zwier, M.C., Adelman, J.L., Kaus, J.W., Pratt, A.J., Wong, K.F., Rego, N.B., Suarez, E., Lettieri, S.,
Wang, D. W., Grabe, M., Zuckerman, D. M., and Chong, L. T. "WESTPA: An Interoperable, Highly 
Scalable Software Package For Weighted Ensemble Simulation and Analysis," J. Chem. Theory Comput., 11: 800−809 (2015). 

To help us fund development, please cite the article listed above and
consider contributing documentation or code to the WESTPA community.

WESTPA is free software, licensed under the terms of the GNU General Public
License, Version 3. See the file ``COPYING`` for more information.


------------
Requirements
------------

WESTPA is written in Python and requires version 2.7. WESTPA further requires
a large number of scientific software libraries for Python and other
languages. The simplest way to meet these requirements is to download the
Anaconda Python distribution from www.continuum.io (free for all users).

WESTPA currently runs on Unix-like operating systems, including Linux and
Mac OS X. It is developed and tested on x86_64 machines running Linux.


------------
Installation
------------

After obtaining a copy of the code (see
https://chong.chem.pitt.edu/wewiki/Obtaining_the_WESTPA_code for details), run
``setup.sh`` in the ``westpa`` directory. If the version of Python you will
be using to run the code is not first on your $PATH, then set the environment
variable WEST_PYTHON to the Python interpreter you want to use. For example::

    cd westpa
    export WEST_PYTHON=/opt/anaconda/bin/python2.7
    ./setup.sh


---------------
Getting started
---------------

High-level tutorials of how to use the WESTPA software are available from
https://chong.chem.pitt.edu/wewiki/WESTPA_tutorials. Further, all WESTPA
command-line tools (located in ``westpa/bin``) provide detailed help when
given the -h/--help option.


------------
Getting help
------------

Documentation is available from the WESTPA wiki, located at
https://chong.chem.pitt.edu/wewiki.

A mailing list for WESTPA is available, at which one can ask questions (or see
if a question one has was previously addressed). This is the preferred means
for obtaining help and support. See http://groups.google.com/group/westpa-users
to sign up or search archived messages.

-------------------------------------------------------
Copyright, license, and warranty information for WESTPA
-------------------------------------------------------

The WESTPA package is copyright (c) 2013, Matthew C. Zwier and
Lillian T. Chong. (Individual contributions noted in each source file.)

WESTPA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

WESTPA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (see the included file ``COPYING``).  If not,
see <http://www.gnu.org/licenses/>.

Unless otherwise noted, source files included in this distribution and
lacking a more specific attribution are subject to the above copyright,
terms, and conditions.


-------------------------------------------------------
Copyright and license information for included software
-------------------------------------------------------

Distributions of WESTPA include a number of components without modification,
each of which is subject to its own individual terms and conditions. Please
see each package's documentation for the most up-to-date possible information
on authorship and licensing. Such packages include:

  h5py
    See lib/h5py/docs/source/licenses.rst
    
  blessings
    See lib/blessings/LICENSE
    
In addition, the ``wwmgr`` work manager is derived from the
``concurrent.futures`` module (as included in Python 3.2) by Brian Quinlan and
copyright 2011 the Python Software Foundation. See 
http://docs.python.org/3/license.html for more information.
