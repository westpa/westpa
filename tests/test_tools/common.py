'''Module of common tool test class'''


# Placeholder class that will set all kwargs as attributes
class MockArgs:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
