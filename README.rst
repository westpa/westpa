===============
WESTPA 1.0
===============


--------
Overview
--------

WESTPA is a package for constructing and running stochastic simulations using the "weighted ensemble" approach 
of Huber and Kim (1996). For use of WESTPA please cite the following:

Zwier, M.C., Adelman, J.L., Kaus, J.W., Pratt, A.J., Wong, K.F., Rego, N.B., Suarez, E., Lettieri, S.,
Wang, D. W., Grabe, M., Zuckerman, D. M., and Chong, L. T. "WESTPA: An Interoperable, Highly 
Scalable Software Package For Weighted Ensemble Simulation and Analysis," J. Chem. Theory Comput., 11: 800−809 (2015). 

See this page_ for an overview of weighted ensemble simulation.

To help us fund development and improve WESTPA please fill out a one-minute survey_ and consider 
contributing documentation or code to the WESTPA community.

WESTPA is free software, licensed under the terms of the GNU General Public
License, Version 3. See the file ``COPYING`` for more information.

.. _survey: https://docs.google.com/forms/d/e/1FAIpQLSfWaB2aryInU06cXrCyAFmhD_gPibgOfFk-dspLEsXuS9-RGQ/viewform
.. _page: https://westpa.github.io/westpa/overview.html

--------------------------------
Obtaining and Installing WESTPA
--------------------------------

First, install the `Anaconda Python distribution`_. Then, make sure you are able to activate conda environments (this is usually taken care of by the Anaconda installer).

WESTPA can then be installed through conda in a dedicated environment with the following.

``conda create -n westpa -c conda-forge westpa``
  
WESTPA will be ready to use after activation with the following command.

``. $(dirname $(dirname `which python2.7`))/$conda_env/westpa-2017.10/westpa.sh``
  
Feel free to install any other conda packages alongside WESTPA in your environment. AmberTools, GROMACS and OpenMM all
provide conda installations of their MD packages. An example command to create an environment containing WESTPA and AmberTools is given below.

``conda create -n westpa -c conda-forge -c ambermd westpa ambertools``
    
.. _`Anaconda Python distribution`: https://www.continuum.io/downloads 

---------------
Getting started
---------------

WESTPA simulation checklist_ 

To define environment variables post-installation, simply source the 
``westpa.sh`` script in the ``westpa`` directory from the command line
or your setup scripts.

High-level tutorials of how to use the WESTPA software can be found here_.
Further, all WESTPA command-line tools (located in ``westpa/bin``) provide detailed help when
given the -h/--help option.

Finally, while WESTPA is a powerful tool that enables expert simulators to access much longer 
timescales than is practical with standard simulations, there can be a steep learning curve to 
figuring out how to effectively run the simulations on your computing resource of choice. 
For serious users who have completed the online tutorials and are ready for production simulations 
of their system, we invite you to contact Lillian Chong (ltchong AT pitt DOT edu) about spending 
a few days with her lab and/or setting up video conferencing sessions to help you get your 
simulations off the ground.

.. _here: https://github.com/westpa/westpa/wiki/WESTPA-Tutorials
.. _checklist: https://github.com/westpa/westpa/wiki/Checklist-for-Setting-Up-a-WESTPA-Simulation

------------
Getting help
------------

WESTPA FAQ_

A mailing list for WESTPA is available, at which one can ask questions (or see
if a question one has was previously addressed). This is the preferred means
for obtaining help and support. See http://groups.google.com/group/westpa-users
to sign up or search archived messages.

.. _FAQ: https://github.com/westpa/westpa/wiki/Frequently-Asked-Questions-%28FAQ%29

-------------------------------------------------------
Copyright, license, and warranty information
-------------------------------------------------------

For WESTPA
###########

The WESTPA package is copyright (c) 2013, Matthew C. Zwier, A. J. Pratt,
Lillian T. Chong, and contributors

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
See ``lib/h5py/docs/source/licenses.rst``
    
blessings
See ``lib/blessings/LICENSE``
    
In addition, the ``wwmgr`` work manager is derived from the
``concurrent.futures`` module (as included in Python 3.2) by Brian Quinlan and
copyright 2011 the Python Software Foundation. See 
http://docs.python.org/3/license.html for more information.
