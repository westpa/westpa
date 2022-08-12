===============
WESTPA 2.0 
===============

|ghactions| |anaconda| |tutorial1| |tutorial2| 

|docs| |usersgroup| |develgroup| 

.. |ghactions| image:: https://github.com/westpa/westpa/actions/workflows/test.yaml/badge.svg?branch=westpa-2.0-restruct
              :target: https://github.com/westpa/westpa/actions/workflows/test.yaml
              :alt: GitHub Actions

.. |docs| image:: https://readthedocs.org/projects/westpa/badge/?version=latest
         :target: https://westpa.readthedocs.io/en/latest/?badge=latest
         :alt: Documentation Status

.. |tutorial1| image:: https://img.shields.io/badge/WESTPA-tutorials-brightgreen.svg
        :target: https://github.com/westpa/westpa_tutorials 
        :alt: WESTPA Tutorials GitHub

.. |tutorial2| image:: https://img.shields.io/badge/WESTPA2-tutorials-blueviolet.svg
        :target: https://github.com/westpa/westpa2_tutorials
        :alt: WESTPA 2.0 Tutorials GitHub

.. |usersgroup| image:: https://img.shields.io/badge/Google%20Group-Users-lightgrey.svg
        :target: http://groups.google.com/group/westpa-users 
        :alt: Users Google Group

.. |develgroup| image:: https://img.shields.io/badge/Google%20Group-Developers-lightgrey.svg
        :target: https://groups.google.com/g/westpa-devel 
        :alt: Developers Google Group

.. |anaconda| image:: https://anaconda.org/conda-forge/westpa/badges/version.svg
   :alt: Anaconda
   :target: https://anaconda.org/conda-forge/westpa


--------
Overview
--------

WESTPA is a package for constructing and running stochastic simulations using the "weighted ensemble" approach 
of Huber and Kim (1996). For use of WESTPA please cite the following:

Zwier, M.C., Adelman, J.L., Kaus, J.W., Pratt, A.J., Wong, K.F., Rego, N.B., Suarez, E., Lettieri, S.,
Wang, D.W., Grabe, M., Zuckerman, D.M., and Chong, L.T. "WESTPA: An Interoperable, Highly 
Scalable Software Package For Weighted Ensemble Simulation and Analysis," J. Chem. Theory Comput., 11: 800−809 (2015). 

Russo, J. D., Zhang, S., Leung, J.M.G., Bogetti, A.T., Thompson, J.P., DeGrave, A.J., Torrillo, P.A., Pratt, A.J., 
Wong, K.F., Xia, J., Copperman, J., Adelman, J.L., Zwier, M.C., LeBard, D.N., Zuckerman, D.M., Chong, L.T. 
WESTPA 2.0: High-Performance Upgrades for Weighted Ensemble Simulations and Analysis of Longer-Timescale Applications. 
J. Chem. Theory Comput., 18 (2): 638–649 (2022).

See this page_ and this powerpoint_ for an overview of weighted ensemble simulation.

To help us fund development and improve WESTPA please fill out a one-minute survey_ and consider 
contributing documentation or code to the WESTPA community.

WESTPA is free software, licensed under the terms of the MIT License. See the file ``LICENSE`` for more information.

.. _survey: https://docs.google.com/forms/d/e/1FAIpQLSfWaB2aryInU06cXrCyAFmhD_gPibgOfFk-dspLEsXuS9-RGQ/viewform
.. _page: https://westpa.github.io/westpa/overview.html
.. _powerpoint: https://pitt.box.com/s/metui7tsfwh3bcv1xgbbj4g6fe0uokag

------------
Requirements
------------

WESTPA is written in Python and requires version 3.7 or later. WESTPA also requires a number of Python scientific software packages. 
The simplest way to meet these requirements is to download the
Anaconda Python distribution from www.continuum.io (free for all users).

WESTPA currently runs on Unix-like operating systems, including Linux and
Mac OS X. It is developed and tested on x86_64 machines running Linux.

--------------------------------
Obtaining and Installing WESTPA
--------------------------------

WESTPA is developed and tested on Unix-like operating systems, including Linux and Mac OS X.


Regardless of the chosen method of installation, before installing WESTPA, we recommend you to first install the Python 3 version provided by the latest free `Anaconda Python distribution`_. After installing Anaconda, create a new python environment for the WESTPA install with the following::

    conda create -n westpa-2.0 python=3.9
    conda activate westpa-2.0

Then, we recommend installing WESTPA through conda or pip. Execute either of the following::

    conda install -c conda-forge westpa

or::

    python -m pip install westpa
    
See the install instructions on our `wiki`_ for more detailed information. 
    

To install from source (**not recommended**), start by downloading the corresponding tar.gz file from the `releases page`_. After downloading the file, unpack the file and install WESTPA by executing the following::

    tar xvzf westpa-main.tar.gz
    cd westpa
    python -m pip install -e .

.. _`releases page`: https://github.com/westpa/westpa/releases
.. _`Anaconda Python distribution`: https://www.anaconda.com/products/individual
.. _`wiki`: https://github.com/westpa/westpa/wiki/WESTPA-Quick-Installation

---------------
Getting started
---------------

High-level tutorials of how to use the WESTPA software can be found here_.
Further, all WESTPA command-line tools provide detailed help when
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
