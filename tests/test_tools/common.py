'''Module of common tool test class'''

import sys
import pytest

flaky_on_macos = pytest.mark.flaky(condition=sys.platform.startswith('darwin'), reruns=5, reason='flaky on macos')


# Placeholder class that will set all kwargs as attributes
class MockArgs:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
