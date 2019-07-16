===============
WESTPA 1.0 beta
===============


--------
Overview
--------

WESTPA is a package for constructing and running stochastic simulations using the "weighted ensemble" approach 
of Huber and Kim (1996) (see overview_). 

For use of WESTPA please cite the following:

Zwier, M.C., Adelman, J.L., Kaus, J.W., Pratt, A.J., Wong, K.F., Rego, N.B., Suarez, E., Lettieri, S.,
Wang, D. W., Grabe, M., Zuckerman, D. M., and Chong, L. T. "WESTPA: An Interoperable, Highly 
Scalable Software Package For Weighted Ensemble Simulation and Analysis," J. Chem. Theory Comput., 11: 800−809 (2015). 

To help us fund development and improve WESTPA please fill out a one-minute survey_ and consider 
contributing documentation or code to the WESTPA community.

WESTPA is free software, licensed under the terms of the GNU General Public
License, Version 3. See the file ``COPYING`` for more information.

.. _survey: https://docs.google.com/forms/d/e/1FAIpQLSfWaB2aryInU06cXrCyAFmhD_gPibgOfFk-dspLEsXuS9-RGQ/viewform
.. _overview: https://westpa.github.io/westpa/overview.html

--------------------------------
Obtaining and Installing WESTPA
--------------------------------

WESTPA is developed and tested on Unix-like operating systems, including Linux and Mac OS X.

Before installing WESTPA, you will need to first install the Python 2.7 version provided by the latest free `Anaconda Python distribution`_. After installing the Anaconda Python distribution, either add the Python executable to your $PATH or set the environment variable WEST_PYTHON::

    export WEST_PYTHON=/opt/anaconda/bin/python3

We recommend obtaining the latest release of WESTPA by downloading the corresponding tar.gz file from the `releases page`_. After downloading the file, unpack the file and install WESTPA by executing the following::

    tar xvzf westpa-master.tar.gz
    cd westpa
    ./setup.sh

A westpa.sh script is created during installation, and will set the following environment variables::

    WEST_ROOT
    WEST_BIN
    WEST_PYTHON

These environment variables must be set in order to run WESTPA on your computing cluster.

To define environment variables post-installation, simply source the 
``westpa.sh`` script in the ``westpa`` directory from the command line
or your setup scripts.

.. _`releases page`: https://github.com/westpa/westpa/releases
.. _`Anaconda Python distribution`: https://www.continuum.io/downloads 

---------------
Getting started
---------------

A Quickstart guide and tutorials are provided here_. 

.. _here: https://github.com/westpa/westpa/wiki

------------
Getting help
------------

FAQ
###

Responses to frequently asked questions (FAQ) can be found in the following page: 
  
- `Frequently Asked Questions (FAQ)`_
  

A mailing list for WESTPA is available, at which one can ask questions (or see
if a question one has was previously addressed). This is the preferred means
for obtaining help and support. See http://groups.google.com/group/westpa-users
to sign up or search archived messages.

Further, all WESTPA command-line tools (located in ``westpa/bin``) provide detailed help when
given the -h/--help option.

Finally, while WESTPA is a powerful tool that enables expert simulators to access much longer 
timescales than is practical with standard simulations, there can be a steep learning curve to 
figuring out how to effectively run the simulations on your computing resource of choice. 
For serious users who have completed the online tutorials and are ready for production simulations 
of their system, we invite you to contact Lillian Chong (ltchong AT pitt DOT edu) about spending 
a few days with her lab and/or setting up video conferencing sessions to help you get your 
simulations off the ground.

.. _`Frequently Asked Questions (FAQ)`: https://github.com/westpa/westpa/wiki/Frequently-Asked-Questions-%28FAQ%29

-------------------------------------------------------
Copyright, license, and warranty information
-------------------------------------------------------

For WESTPA
###########

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


For included software
######################

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
