import ast
import math
from dataclasses import dataclass

import numpy as np


VALID_BINARY_OPS = {ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow}
VALID_COMPARISON_OPS = {ast.Lt, ast.LtE, ast.Gt, ast.GtE}
VALID_UNARY_OPS = {ast.UAdd, ast.USub, ast.Not}
VALID_FUNCTIONS = {
    'acos',
    'acosh',
    'asin',
    'asinh',
    'atan',
    'atanh',
    'cos',
    'cosh',
    'degrees',
    'erf',
    'erfc',
    'exp',
    'expm1',
    'gamma',
    'lgamma',
    'log',
    'log10',
    'log1p',
    'log2',
    'pow',
    'radians',
    'sin',
    'sinh',
    'sqrt',
    'tan',
    'tanh',
}


@dataclass
class Sink:
    """A container object representing the subset of progress coordinate space
    designated as the sink."""
    variables: str
    predicate: str

    @classmethod
    def from_string(cls, string):
        """Construct a :class:`Sink` object from its string representation.

        Parameters
        ----------
        string : str
            A string of the form ``'<variables> : <predicate>'``, where
            ``<variables>`` is a comma separated list of variable names, and
            ``<predicate>`` is a boolean expression involving those variables.
            The variables are understood to correspond to consecutive axes of
            progress coordinate space, starting at index 0. (Hence the maximum
            number of variables that may be provided is ``pcoord_ndim``.)
            The predicate may include floating-point constants, comparisons for
            inequality (``>``, ``>=``, ``<``, ``<=``), and algebraic operations
            (``+``, ``-``, ``*``, ``/``, ``**``), as well as power, exponential,
            logarithmic, trigonometric, and hyperbolic functions defined by the
            C standard (e.g., ``exp()``, ``sin()``, ``sqrt()``).

        """
        words = string.split(':')
        if len(words) != 2:
            raise ValueError("invalid syntax: expected '<variables> : <predicate>'")
        variables, predicate = [word.strip() for word in words]
        return cls(variables, predicate)

    def __post_init__(self):
        # Parse the variables (e.g., 'x' or 'x, y').
        expr = ast.parse(self.variables, mode='eval')
        if type(expr.body) is ast.Name:
            variable_names = [expr.body.id]
        elif type(expr.body) is ast.Tuple and all(type(elt) is ast.Name for elt in expr.body.elts):
            variable_names = [elt.id for elt in expr.body.elts]
            if len(variable_names) != len(set(variable_names)):
                raise ValueError('variable names must be unique')
        else:
            raise ValueError('<variables> must be a variable name or a tuple of variable names')

        # Compile the assignment (e.g., 'x, *_ = _x' or 'x, y, *_ = _x').
        for name in variable_names:
            if name in ('_', '_x'):
                raise ValueError(f'variable name {name!r} is reserved')
        self._assignment = compile(', '.join(variable_names) + ', *_ = _x', '<string>', 'exec')

        # Parse and compile the predicate (e.g., 'x > 0' or 'x**2 + y**2 < 1').
        expr = ast.parse(self.predicate, mode='eval')
        validator = PredicateValidator(self.predicate, variable_names)
        try:
            validator.visit(expr.body)
        except (TypeError, ValueError):
            raise
        else:
            # Transform function calls from 'func()' to 'math.func()'.
            transformer = MathFunctionTransformer()
            expr = ast.fix_missing_locations(transformer.visit(expr))
        self._predicate = compile(expr, '<string>', 'eval')

        # Check that the predicate evaluates to a boolean.
        ndim = len(variable_names)
        try:
            result = np.zeros(ndim) in self
        except Exception as e:
            raise RuntimeError(f'an error occurred while evaluating the predicate: {e}')
        if not isinstance(result, (bool, np.bool_)):
            raise TypeError(f'predicate must evaluate to a boolean, not {type(result).__name__}')

    def __contains__(self, _x):
        """Return whether the given progress coordinate value is in the sink.

        Parameters
        ----------
        _x : ndarray, shape=(pcoord_ndim,)
            A point in progress coordinate space.

        Returns
        -------
        bool
            True if `_x` is in the sink, else False.

        """
        exec(self._assignment)
        return eval(self._predicate)

    def __str__(self):
        return f'{self.variables} : {self.predicate}'


@dataclass
class PredicateValidator(ast.NodeVisitor):
    source: str
    variable_names: set[str]

    def get_source_segment(self, node):
        return ast.get_source_segment(self.source, node)

    def visit_Attribute(self, node):
        raise ValueError(f'attribute references are not supported: {self.get_source_segment(node)}')

    def visit_BinOp(self, node):
        if type(node.op) not in VALID_BINARY_OPS:
            raise ValueError(f'invalid binary operation: {self.get_source_segment(node)}')
        self.visit(node.left)
        self.visit(node.right)

    def visit_BoolOp(self, node):
        for value in node.values:
            self.visit(value)

    def visit_Call(self, node):
        if type(node.func) is not ast.Name:
            self.visit(node.func)
        elif node.func.id not in VALID_FUNCTIONS:
            raise ValueError(f'{node.func.id}() is not a recognized function')
        for arg in node.args + node.keywords:
            self.visit(arg)

    def visit_Compare(self, node):
        for op in node.ops:
            if type(op) not in VALID_COMPARISON_OPS:
                raise ValueError(f'invalid comparison: {self.get_source_segment(node)}')
        self.visit(node.left)
        for comparator in node.comparators:
            self.visit(comparator)

    def visit_Constant(self, node):
        if not isinstance(node.value, (float, int)):
            typename = type(node.value).__name__
            raise TypeError(f'constants must be floats or integers, not {typename}: {self.get_source_segment(node)}')

    def visit_Name(self, node):
        if node.id not in self.variable_names:
            raise ValueError(f'{node.id!r} is not a recognized variable name')

    def visit_Subscript(self, node):
        raise ValueError(f'subscripts are not supported: {self.get_source_segment(node)}')

        if type(node.slice) is not ast.Constant:
            raise ValueError('index must be a constant value')
        if not isinstance(node.slice.value, int):
            typename = type(node.slice.value).__name__
            raise TypeError(f'index must be an integer, not {typename}: {self.get_source_segment(node)}')

    def visit_UnaryOp(self, node):
        if type(node.op) not in VALID_UNARY_OPS:
            raise ValueError(f'invalid unary operation: {self.get_source_segment(node)}')
        self.visit(node.operand)

    def generic_visit(self, node):
        raise ValueError(f'invalid input: {self.get_source_segment(node)}')


class MathFunctionTransformer(ast.NodeTransformer):

    def visit_Call(self, node):
        return ast.Call(
            func=ast.Attribute(
                value=ast.Name(id=math.__name__, ctx=ast.Load()),
                attr=node.func.id,
                ctx=ast.Load(),
            ),
            args=node.args,
            keywords=node.keywords,
        )
