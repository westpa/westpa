from __future__ import print_function, division; __metaclass__ = type
from westpa.progress import ProgressIndicator
from westtools.tool_classes.core import WESTToolComponent
import westpa

class ProgressIndicatorComponent(WESTToolComponent):
    def __init__(self):
        super(ProgressIndicatorComponent,self).__init__()
        self.indicator = None
        
    def add_args(self, parser):
        pass

    def process_args(self, args):
        self.indicator = ProgressIndicator()
        if westpa.rc.quiet_mode:
            self.indicator.do_fancy = False


