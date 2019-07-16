# Copyright (C) 2013 Joshua L. Adelman
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.

"""westext.stringmethod - Plugin to drive the adaptive evolution of one or more
strings of Voronoi bins

Joshua L. Adelman 2011
"""

from abc import ABCMeta, abstractmethod, abstractproperty


class WESTStringMethod(object):

    ___metaclass__ = ABCMeta

    def __init__(self, centers, **kwargs):
        pass

    @abstractproperty
    def centers(self):
        """ Return the centers of all of the strings
        """
        pass

    @abstractproperty
    def length(self):
        """ Return a list of the lengths of each string
        """
        pass

    @abstractmethod
    def update_string_centers(self, avgcoords, binprob):
        """ Given a set of average coordinates (avgcoords) in each bin
        and the individual probabilities for each bin (binprob), update
        the string centers
        """
        pass

from . import string_method
from .string_method import DefaultStringMethod

from . import string_driver
from .string_driver import StringDriver
