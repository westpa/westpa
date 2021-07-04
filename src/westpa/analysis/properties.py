import functools
import numpy as np

from typing import Callable
from westpa.analysis.core import Walker, Iteration


class Property:

    def __set_name__(self, owner, name):
        if hasattr(self, 'owner'):
            raise TypeError("can't assign Property to multiple owners")
        if name != self.name:
            raise AttributeError("can't assign multiple names to Property")

        self.owner = owner

    def __init__(self, fget, *, owner=None, cache=False):
        if fget is None:
            return functools.partial(self.__init__, owner=owner, cache=cache)

        if not isinstance(fget, Callable):
            raise TypeError('fget must be callable')
        if not isinstance(cache, bool):
            raise TypeError('cache must be True or False')

        self.fget = fget
        self.cache = cache

        if owner is not None:
            if hasattr(owner, self.name):
                raise AttributeError(f"class '{owner.__name__}' already has "
                                     f"attribute '{self.name}'")
            setattr(owner, self.name, self)
            self.__set_name__(owner, self.name)

    @property
    def name(self):
        return self.fget.__name__

    @property
    def __doc__(self):
        return self.fget.__doc__

    @property
    def private_name(self):
        return '_' + self.name

    def __get__(self, instance, owner):
        if instance is None:
            return self

        if owner is not self.owner:
            raise TypeError(f'owner must be {self.owner}')

        if hasattr(instance, self.private_name):
            value = getattr(instance, self.private_name)
        else:
            value = self.fget(instance)
            if self.cache:
                setattr(instance, self.private_name, value)

        return value

    def __set__(self, instance, value):
        raise AttributeError("can't set Property attribute")

    def __call__(self, arg):
        if not isinstance(arg, self.owner):
            raise TypeError(f'argument must be an instance of {self.owner}')
        return self.__get__(arg, self.owner)

    def max(self, args):
        return max(map(self, args))

    def min(self, args):
        return min(map(self, args))

    def argmax(self, args):
        return argmax(self, args)

    def argmin(self, args):
        return argmin(self, args)

    def mean(self, args, weights=None):
        return mean(self, args, weights=None)


class WalkerProperty(Property):

    def __init__(self, fget, *, cache=False):
        super().__init__(fget, owner=Walker, cache=cache)

    def weighted_average(self, iterations):
        num_iterations = 0
        for iteration in iterations:
            num_iterations += 1
            values = [self(walker) for walker in iteration]
            if num_iterations == 1:
                value = np.dot(iteration.weights, values)
                continue
            value += np.dot(iteration.weights, values)
        return value / num_iterations


def argmax(func, args):
    argmax = None
    maxval = -np.inf
    for arg, val in zip(args, map(func, args)):
        if val > maxval:
            maxval = val
            argmax = arg
    return argmax


def argmin(func, args):
    return argmax(lambda x: -func(x), args)
