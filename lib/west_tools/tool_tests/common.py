'''Module of common tool test class'''

import nose
import tempfile, shutil, os, itertools, re

def get_arg_flag(key):
    return "--{}".format(re.sub('_', '-', key))

def make_args(outfile,**args):
    '''Makes a list of argument strings from a dictionary for command line tools to process'''
    
    arglist = ['{}={}'.format(get_arg_flag(key), value) if type(value) is not bool else '{}'.format(key) for (key, value) in args.items() if value is not None]
    #arglist = []
    #for key, value in args.items():
    #    if value is not None:
    #        if type(value) is not bool:
    #            arglist.append('{}={}'.format(get_arg_flag(key), value))
    #        else:
    #            arglist.append('{}'.format(get_arg_flag(key)))

    arglist.append('-o={}'.format(outfile))

    return arglist

def cycle_args(arg_list):
    '''Generator function to generate different argument combo's - returns as dictionaries'''
    keys = arg_list.keys()
    for values in itertools.product(*arg_list.values()):
        yield {k:v for (k, v) in zip(keys, values)}

class CommonToolTest:

    test_name = None

    @classmethod
    def setUpClass(cls):
        cls.tempdir = tempfile.mkdtemp(prefix='w_tools_test')
        os.chdir(cls.tempdir)
        
    @classmethod
    def tearDownClass(cls):
        os.chdir('{}/lib/west_tools/tests'.format(os.environ['WEST_ROOT']))
        shutil.rmtree(cls.tempdir)

    def mktemp(self, suffix='.h5', prefix = '', dirname = None):
        '''Helper method to make a tempfile (to be written/read by tools)'''

        if dirname is None:
            dirname = self.tempdir

        return tempfile.mktemp(suffix, prefix, dir=dirname)

    def check_args(self, outcome, errmsg, args):
        self.check_args.__func__.description = '{} Args: '.format(self.test_name) + (', '.join(args) or '[defaults]')
        assert outcome, str(errmsg)

    def check_runs_with_args(self, **kwargs):

        try:
            self.w.go() #Run tool with given args

            self.check_output(**kwargs) #Check that output hdf5 file is set up correctly

        except Exception as e:
            return (0, e)
        else:
            return(1, None)

    def check_output(self, **kwargs):
        raise NotImplementedError

    def check_args_processed(self, **kwargs):
        raise NotImplementedError

    def test_args(self, **kwargs):
        raise NotImplementedError
