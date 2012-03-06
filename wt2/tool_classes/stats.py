from core import WEMDTool

from wt2.stats.mcbs import get_bssize

class UsesMCBS(WEMDTool):
    def __init__(self):
        super(UsesMCBS,self).__init__()
        self.wt2_bssize = None
        self.wt2_bsconf = None

    def add_args(self, parser):
        group = parser.add_argument_group('bootstrapping options')
        group.add_argument('--confidence', dest='wt2_bsconf', metavar='CONFIDENCE', type=float, default=0.95,
                            help='Construct a confidence interval of width CONFIDENCE (default: 0.95=95%%)')
        group.add_argument('--bssize', dest='wt2_bssize', metavar='BSSIZE', type=int,
                            help='Use a bootstrap of BSSIZE samples to calculate error (default: chosen from confidence)')
        
    def process_args(self, args):
        self.wt2_bsconf = args.bsconf
        self.wt2_bssize = args.bssize or get_bssize(1-self.wt2_bsconf)
