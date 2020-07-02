"""westext.stringmethod - Plugin to drive the adaptive evolution of one or more
strings of Voronoi bins

Joshua L. Adelman 2011
"""

from abc import ABCMeta, abstractmethod, abstractproperty

from . import string_method
from .string_method import DefaultStringMethod

from . import string_driver
from .string_driver import StringDriver


class WESTStringMethod:

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


__all__ = ['string_method', 'DefaultStringMethod', 'string_driver', 'StringDriver']
