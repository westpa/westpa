import westpa
from westpa.core.progress import ProgressIndicator
from westpa.tools.core import WESTToolComponent


class ProgressIndicatorComponent(WESTToolComponent):
    def __init__(self):
        super().__init__()
        self.indicator = None

    def add_args(self, parser):
        pass

    def process_args(self, args):
        self.indicator = ProgressIndicator()
        if westpa.rc.quiet_mode or westpa.rc.verbose_mode or westpa.rc.debug_mode:
            self.indicator.fancy = False
