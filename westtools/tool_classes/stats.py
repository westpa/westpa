from core import WESTToolComponent

from westtools.stats.mcbs import get_bssize

class UsesMCBS(WESTToolComponent):
    def __init__(self):
        super(UsesMCBS,self).__init__()
        self.bssize = None
        self.bsconf = None

    def add_args(self, parser):
        group = parser.add_argument_group('bootstrapping options')
        group.add_argument('--confidence', dest='bsconf', metavar='CONFIDENCE', type=float, default=0.95,
                            help='Construct a confidence interval of width CONFIDENCE (default: 0.95=95%%)')
        group.add_argument('--bssize', dest='bssize', metavar='BSSIZE', type=int,
                            help='Use a bootstrap of BSSIZE samples to calculate error (default: chosen from confidence)')
        
    def process_args(self, args):
        self.bsconf = args.bsconf
        self.bssize = args.bssize or get_bssize(1-self.westtools_bsconf)
