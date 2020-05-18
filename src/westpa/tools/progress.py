
from westpa.progress import ProgressIndicator
from westtools.core import WESTToolComponent
import westpa

class ProgressIndicatorComponent(WESTToolComponent):
    def __init__(self):
        super(ProgressIndicatorComponent,self).__init__()
        self.indicator = None
        
    def add_args(self, parser):
        pass

    def process_args(self, args):
        self.indicator = ProgressIndicator()
        if westpa.rc.quiet_mode or westpa.rc.verbose_mode or westpa.rc.debug_mode:
            self.indicator.fancy = False
