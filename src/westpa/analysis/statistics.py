import numpy as np


def time_average(observable, iterations):
    """Compute the time average of an observable.

    Parameters
    ----------
    observable : Callable[[Walker], ArrayLike]
        Function that takes a walker as input and returns a number or
        a fixed-size array of numbers.
    iterations : Sequence[Iteration]
        Sequence of iterations over which to compute the average.

    Returns
    -------
    ArrayLike
        The time average of `observable` over `iterations`.

    """
    for iteration in iterations:
        values = [observable(walker) for walker in iteration]
        try:
            value += np.dot(iteration.weights, values)
        except NameError:
            value = np.dot(iteration.weights, values)
    return value / len(iterations)
