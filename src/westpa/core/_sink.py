import ast
import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike


@dataclass
class Sink:
    indicator_function: Callable[[ArrayLike], bool]

    def __contains__(self, x):
        return self.indicator_function(x)

    @classmethod
    def from_string(cls, string):
        """Construct a sink from a string representation.

        Parameters
        ----------
        string : str
            A string of the form ``<variables>: <predicate>``, where
            ``<variables>`` is a variable name or tuple of variable names, and
            ``<predicate>`` is a boolean expression involving those variables.
            The predicate may include arithmetic operations, comparisons, and
            calls to functions provided by the :py:mod:`math` module.

        Returns
        -------
        Sink
            A container object whose membership test evaluates the predicate.

        """
        words = string.split(':')
        if len(words) != 2:
            raise ValueError('invalid syntax: expected <variables>: <predicate>')
        variables, predicate = [word.strip() for word in words]

        # Parse the variables (e.g., 'x' or 'x, y').
        expr = ast.parse(variables, mode='eval')
        if type(expr.body) is ast.Name:
            variable_names = {expr.body.id}
            unpacked = False
        elif type(expr.body) is ast.Tuple and all(type(elt) is ast.Name for elt in expr.body.elts):
            variable_names = {elt.id for elt in expr.body.elts}
            unpacked = True
        else:
            raise ValueError('<variables> must be a variable name or a tuple of variable names')
        if '_x' in variable_names:
            raise ValueError("variable name '_x' is reserved")

        # Parse the predicate (e.g., 'x[0] > 0' or 'x**2 and y**2 < 1').
        expr = ast.parse(predicate, mode='eval')
        validator = PredicateValidator(predicate, variable_names, unpacked)
        try:
            validator.visit(expr.body)
        except (TypeError, ValueError):
            raise
        else:
            # Transform function calls from 'func()' to 'math.func()'.
            transformer = MathFunctionTransformer()
            expr = ast.fix_missing_locations(transformer.visit(expr))

        def indicator_function(_x):
            exec(f'{variables} = _x')
            return eval(compile(expr, '<string>', 'eval'))

        # Infer the minimum dimension of the coordinate space.
        if unpacked:
            ndim = len(variable_names)
        else:
            ndim = 1 + max(node.slice.value for node in ast.walk(expr) if type(node) is ast.Subscript)

        # Check that the predicate evaluates to a boolean.
        try:
            result = indicator_function(np.zeros(ndim))
        except Exception as e:
            raise RuntimeError(f'an error occurred while evaluating the predicate: {e}')
        if not isinstance(result, (bool, np.bool_)):
            raise TypeError(f'predicate must evaluate to a boolean, not {type(result).__name__}')

        if isinstance(result, np.bool_):
            return cls(lambda x: bool(indicator_function(x)))
        else:
            return cls(indicator_function)


@dataclass
class PredicateValidator(ast.NodeVisitor):
    source: str
    variable_names: set[str]
    unpacked: bool

    def get_source_segment(self, node):
        return ast.get_source_segment(self.source, node)

    def visit_Attribute(self, node):
        raise ValueError(f'attribute references are not supported: {self.get_source_segment(node)}')

    def visit_BinOp(self, node):
        if type(node.op) not in (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow):
            raise ValueError(f'invalid binary operation: {self.get_source_segment(node)}')
        self.visit(node.left)
        self.visit(node.right)

    def visit_BoolOp(self, node):
        for value in node.values:
            self.visit(value)

    def visit_Call(self, node):
        if type(node.func) is not ast.Name:
            self.visit(node.func)
        elif node.func.id not in (
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
        ):
            raise ValueError(f'{node.func.id}() is not a recognized function')
        for arg in node.args + node.keywords:
            self.visit(arg)

    def visit_Compare(self, node):
        for op in node.ops:
            if type(op) not in (ast.Lt, ast.LtE, ast.Gt, ast.GtE):
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
        if not self.unpacked:
            raise ValueError(f'array variable {node.id!r} must be subscripted')

    def visit_Subscript(self, node):
        if type(node.value) is not ast.Name:
            raise ValueError(f'subscripted object must be a variable name: {self.get_source_segment(node)}')
        if node.value.id not in self.variable_names:
            raise ValueError(f'{node.value.id!r} is not a recognized variable name')
        if self.unpacked:
            raise ValueError(f'unpacked variable {node.value.id!r} may not be subscripted')

        if type(node.slice) is not ast.Constant:
            raise ValueError('index must be a constant value')
        if not isinstance(node.slice.value, int):
            typename = type(node.slice.value).__name__
            raise TypeError(f'index must be an integer, not {typename}: {self.get_source_segment(node)}')

    def visit_UnaryOp(self, node):
        if type(node.op) not in (ast.UAdd, ast.USub, ast.Not):
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
