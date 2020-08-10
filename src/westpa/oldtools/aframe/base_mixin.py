class ArgumentError(RuntimeError):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class AnalysisMixin:
    def __init__(self):
        super().__init__()

    def add_args(self, parser, upcall=True):
        if upcall:
            try:
                upfunc = super().add_args
            except AttributeError:
                pass
            else:
                upfunc(parser)

    def process_args(self, args, upcall=True):
        if upcall:
            try:
                upfunc = super().process_args
            except AttributeError:
                pass
            else:
                upfunc(args)
