Documentation Practices
=======================

Documentation for WESTPA is maintained using `Sphinx <http://sphinx-doc.org/>`_
Docstrings are formatted in the `Numpy style
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_,
which are converted to ReStructuredText using Sphinx' `Napoleon
<http://sphinxcontrib-napoleon.readthedocs.org/en/latest/>`_ plugin, which is
included with Sphinx 1.3.

The documentation may be built by navigating to the ``doc`` folder, and
running::

  make html

to prepare an html version or::

  make latexpdf

To prepare a pdf. The latter requires ``latex`` to be available.
