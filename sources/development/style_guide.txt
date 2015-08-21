Style Guide
===========

Preface
-------

The WESTPA documentation should help the user to understand how WESTPA works
and how to use it. To aid in effective communication, a number of guidelines
appear below.

When writing in the WESTPA documentation, please be:

- Correct
- Clear
- Consistent
- Concise

Articles in this documentation should follow the guidelines on this page.
However, there may be cases when following these guidelines will make an
article confusing: when in doubt, use your best judgment and ask for the
opinions of those around you.

Style and Usage
---------------

Acronyms and abbreviations
~~~~~~~~~~~~~~~~~~~~~~~~~~

- Software documentation often involves extensive use of acronyms and
  abbreviations.

  Acronym: A word formed from the initial letter or letters of each or most of
  the parts of a compound term

  Abbreviation: A shortened form of a written word or name that is used in
  place of the full word or name

- Define non-standard acronyms and abbreviations on their first use by using
  the full-length term, followed by the acronym or abbreviation in parentheses.

  A potential of mean force (PMF) diagram may aid the user in visuallizing the
  energy landscape of the simulation.

- Only use acronyms and abbreviations when they make an idea more clear than
  spelling out the full term. Consider clarity from the point of view of a new
  user who is intelligent but may have little experience with computers.

  Correct: The WESTPA wiki supports HyperText Markup Language (HTML). For
  example, the user may use HTML tags to give text special formatting. However,
  be sure to test that the HTML tag gives the desired effect by previewing
  edits before saving.

  Avoid: The WESTPA wiki supports HyperText Markup Language. For example, the
  user may use HyperText Markup Language tags to give text special formatting.
  However, be sure to test that the HyperText Markup Language tag gives the
  desired effect by previewing edits before saving.

  Avoid: For each iter, make sure to return the pcoord and any auxdata.

- Use all capital letters for abbreviating file types. File extensions should
  be lowercase.

  HDF5, PNG, MP4, GRO, XTC

  west.h5, bound.png, unfolding.mp4, protein.gro, segment.xtc

- Provide pronunciations for acronyms that may be difficult to sound out.
- Do not use periods in acronyms and abbreviations except where it is
  customary:

  Correct: HTML, U.S.

  Avoid: H.T.M.L., US

Capitalization
~~~~~~~~~~~~~~

- Capitalize at the beginning of each sentence.
- Do not capitalize after a semicolon.
- Do not capitalize after a colon, unless multiple sentences follow the colon.
- In this case, capitalize each sentence.
- Preserve the capitalization of computer language elements (commands,
- utilities, variables, modules, classes, and arguments).
- Capitilize generic Python variables according to the
- `PEP 0008 Python Style Guide
  <http://www.python.org/dev/peps/pep-0008/#class-names>`_. For example,
  generic class names should follow the *CapWords* convention, such as
  ``GenericClass``.

Contractions
~~~~~~~~~~~~

- Do not use contractions. Contractions are a shortened version of word
  characterized by the omission of internal letters.

  Avoid: can't, don't, shouldn't

- Possessive nouns are not contractions. Use possessive nouns freely.

Internationalization
~~~~~~~~~~~~~~~~~~~~

- Use short sentences (less than 25 words). Although we do not maintain
  WESTPA documentation in languages other than English, some users may use
  automatic translation programs. These programs function best with short
  sentences.
- Do not use technical terms where a common term would be equally or more
  clear.
- Use multiple simple sentences in place of a single complicated sentence.

Italics
~~~~~~~

-  Use italics (surround the word with * * on each side) to highlight words
   that are not part of a sentence's normal grammer.

   Correct: The word *istates* refers to the initial states that WESTPA uses to
   begin trajectories.

Non-English words
~~~~~~~~~~~~~~~~~

- Avoid Latin words and abbreviations.

  Avoid: etc., et cetera, e.g., i.e.

Specially formatted characters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Never begin a sentence with a specially formatted character. This includes
  abbreviations, variable names, and anything else this guide instructs to use
  with special tags. Sentences may begin with *WESTPA*.

  Correct: The program ``ls`` allows the user to see the contents of a
  directory.

  Avoid: ``ls`` allows the user to see the contents of a directory.

- Use the word *and* rather than an ``&`` ampersand .
- When a special character has a unique meaning to a program, first use the
  character surrounded by \`\` tags and then spell it out.

  Correct: Append an ``&`` ampersand to a command to let it run in the
  background.

  Avoid: Append an "&" to a command... Append an ``&`` to a command... Append
  an ampersand to a command...

- There are many names for the ``#`` hash mark, including hash tag, number
  sign, pound sign, and octothorpe. Refer to this symbol as a "hash mark".

Subject
~~~~~~~

- Refer to the end WESTPA user as *the user* in software documentation.

  Correct: The user should use the ``processes`` work manager to run segments
  in parallel on a single node.

- Refer to the end WESTPA user as *you* in tutorials (you is the implied
  subject of commands). It is also acceptable to use personal pronouns such as
  *we* and *our*. Be consistent within the tutorial.

  Correct: You should have two files in this directory, named ``system.py`` and
  ``west.cfg``.

Tense
~~~~~

- Use *should* to specify proper usage.

  Correct: The user should run ``w_truncate -n <var>iter</var>`` to remove
  iterations after and including iter from the HDF5 file specified in the
  WESTPA configuration file.

- Use *will* to specify expected results and output.

  Correct: WESTPA will create a HDF5 file when the user runs ``w_init``.

Voice
~~~~~

- Use active voice. Passive voice can obscure a sentence and add unnecessary
  words.

  Correct: WESTPA will return an error if the sum of the weights of segments
  does not equal one.

  Avoid: An error will be returned if the sum of the weights of segments does
  not equal one.

Weighted ensemble
~~~~~~~~~~~~~~~~~

- Refer to weighted ensemble in all lowercase, unless at the beginning of a
  sentence. Do not hyphenate.

  Correct: WESTPA is an implementation of the weighted ensemble algorithm.

  Avoid: WESTPA is an implementation of the weighted-ensemble algorithm.

  Avoid: WESTPA is an implementation of the Weighted Ensemble algorithm.

WESTPA
~~~~~~

- Refer to WESTPA in all capitals. Do not use bold, italics, or other special
  formatting except when another guideline from this style guide applies.

  Correct: Install the WESTPA software package.

- The word *WESTPA* may refer to the software package or a entity of running
  software.

  Correct: WESTPA includes a number of analysis utilities.

  Correct: WESTPA will return an error if the user does not supply a
  configuration file.

Computer Language Elements
--------------------------

Classes, modules, and libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Display class names in fixed-width font using the `````` tag.

  Correct: ``WESTPropagator``

  Correct: The ``numpy`` library provides access to various low-level
  mathematical and scientific calculation routines.

- Generic class names should be relevant to the properties of the class; do not
  use *foo* or *bar*

       ``class UserDefinedBinMapper(RectilinearBinMapper)``

Methods and commands
~~~~~~~~~~~~~~~~~~~~

- Refer to a method by its name without parentheses, and without prepending
  the name of its class. Display methods in fixed-width font using the ``````
  tag.

  Correct: the ``arange`` method of the ``numpy`` library

  Avoid: the ``arange()`` method of the ``numpy`` library

  Avoid: the ``numpy.arange`` method

- When referring to the arguments that a method expects, mention the method
  without arguments first, and then use the method's name followed by
  parenthesis and arguments.

  Correct: WESTPA calls the ``assign`` method as assign(coords, mask=None,
  output=None)

-  Never use a method or command as a verb.

   Correct: Run ``cd`` to change the current working directory.

   Avoid: ``cd`` into the main simulation directory.

Programming languages
~~~~~~~~~~~~~~~~~~~~~

- Some programming languages are both a language and a command. When referring
  to the language, capitalize the word and use standard font. When referring
  to the command, preserve capitalization as it would appear in a terminal and
  use the `````` tag.

  Using WESTPA requires some knowledge of Python.

  Run ``python`` to launch an interactive session.

  The Bash shell provides some handy capabilities, such as wildcard matching.

  Use ``bash`` to run ``example.sh``.

Scripts
~~~~~~~

- Use the ``.. code-block::`` directive for short scripts. Options are
  available for some languages, such as ``.. code-block:: bash`` and
  ``.. code-block:: python``.

.. code-block:: bash

  #!/bin/bash
  # This is a generic Bash script. 

  BASHVAR="Hello, world!"
  echo $BASHVAR

.. code-block:: python

  #!/usr/bin/env python
  # This is a generic Python script. 

  def main():
      pythonstr = "Hello, world!"
      print(pythonstr)
      return
  if __name__ == "__main__":
      main()

- Begin a code snippet with a ``#!`` *shebang* (yes, this is the real term),
  followed by the usual path to a program. The line after the shebang should be
  an ellipsis, followed by lines of code. Use ``#!/bin/bash`` for Bash scripts,
  ``#!/bin/sh`` for generic shell scripts, and ``#!/usr/bin/env python`` for
  Python scripts. For Python code snippets that are not a stand-alone script,
  place any import commands between the shebang line and ellipsis.

.. code-block:: python

    #!/usr/bin/env python
    import numpy
    ...
    def some_function(generic_vals):
        return 1 + numpy.mean(generic_vals)

- Follow the `PEP 0008 Python Style Guide
  <http://www.python.org/dev/peps/pep-0008/#class-names>`_ for Python scripts.

  - Indents are four spaces.
  - For comments, use the ``#`` hash mark followed by a single space, and
    then the comment's text.
  - Break lines after 80 characters.

-  For Bash scripts, consider following `Google's Shell Style Guide
   <https://google-styleguide.googlecode.com/svn/trunk/shell.xml>`_

  - Indents are two spaces.
  - Use blank lines to improve readability
  - Use ``; do`` and ``; then`` on the same line as ``while``, ``for``, and
    ``if``.
  -  Break lines after 80 characters.

- For other languages, consider following a logical style guide. At minimum, be
  consistent.

Variables
~~~~~~~~~

- Use the fixed-width `````` tag when referring to a variable.

  the ``ndim`` attribute

- When explicitly referring to an attribute as well as its class, refer to an
  attribute as: the ``attr`` attribute of ``GenericClass``, rather than
  ``GenericClass.attr``
- Use the ``$`` dollar sign before Bash variables.

  WESTPA makes the variable ``$WEST_BSTATE_DATA_REF`` available to new
  trajectories.
