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

See this page_ and this powerpoint_ for an overview of weighted ensemble simulation.

To help us fund development and improve WESTPA please fill out a one-minute survey_ and consider 
contributing documentation or code to the WESTPA community.

WESTPA is free software, licensed under the terms of the GNU General Public
License, Version 3. See the file ``COPYING`` for more information.

.. _survey: https://docs.google.com/forms/d/e/1FAIpQLSfWaB2aryInU06cXrCyAFmhD_gPibgOfFk-dspLEsXuS9-RGQ/viewform
.. _page: https://westpa.github.io/westpa/overview.html
.. _powerpoint: https://pitt.box.com/s/metui7tsfwh3bcv1xgbbj4g6fe0uokag

------------
Requirements
------------

WESTPA is written in Python and requires version 3. WESTPA further requires
a large number of scientific software libraries for Python and other
languages. The simplest way to meet these requirements is to download the
Anaconda Python distribution from www.continuum.io (free for all users).

WESTPA currently runs on Unix-like operating systems, including Linux and
Mac OS X. It is developed and tested on x86_64 machines running Linux.

--------------------------------
Obtaining and Installing WESTPA
--------------------------------

WESTPA is developed and tested on Unix-like operating systems, including Linux and Mac OS X.

We recommend installing WESTPA through conda. See the quick install instructions on our `wiki`_ for how to do this. To install from source (not recommended), please continue reading.

Before installing WESTPA, you will need to first install the Python 3 version provided by the latest free `Anaconda Python distribution`_. After installing the Anaconda Python distribution, either add the Python executable to your $PATH or set the environment variable WEST_PYTHON::

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

.. _`releases page`: https://github.com/westpa/westpa/releases
.. _`Anaconda Python distribution`: https://www.anaconda.com/products/individual
.. _`wiki`: https://github.com/westpa/westpa/wiki/WESTPA-Quick-Installation

---------------
Getting started
---------------

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

.. _here: https://github.com/westpa/westpa/wiki/Tutorials

------------
Getting help
------------

WESTPA FAQ_

A mailing list for WESTPA is available, at which one can ask questions (or see
if a question one has was previously addressed). This is the preferred means
for obtaining help and support. See http://groups.google.com/group/westpa-users
to sign up or search archived messages.

.. _FAQ: https://westpa.github.io/westpa/users_guide/faq.html

----------
Developers
----------

Search archived messages or post to the westpa-devel Google group: https://groups.google.com/group/westpa-devel.
