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

A quick command to update the documentation in gh-pages repo is also available::
  
  make ghpages

This command will run Sphinx html command and change the htmls to fit with the gh-pages
format, it also runs::

  git checkout gh-pages
  git commit -a
  git push

for you. Also note that this will change the current branch you are at to gh-pages
branch. It also leaves behind a doc/_build folder that is no longer useful. Once you run 
ghpages command I suggest going up a folder and removing the unnecessary doc folder 
that is there by::

  cd ../
  rm -r doc

Remeber to make sure you are indeed in gh-pages branch, this branch is not supposed to have
a folder named doc. Sometimes if you are not careful git checkout fails and you might end up
removing the folder you were working with if you are not careful. 
