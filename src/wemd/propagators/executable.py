import os, sys, signal, random, subprocess, time, tempfile
import numpy
import logging
from wemd.states import BasisState, InitialState
log = logging.getLogger(__name__)

# Get a list of user-friendly signal names
SIGNAL_NAMES = {getattr(signal, name): name for name in dir(signal) 
                if name.startswith('SIG') and not name.startswith('SIG_')}

import wemd
from wemd import Segment
from wemd.propagators import WEMDPropagator

def pcoord_loader(fieldname, pcoord_return_filename, destobj, single_point):
    """Read progress coordinate data into the ``pcoord`` field on ``destobj``. 
    An exception will be raised if the data is malformed.  If ``single_point`` is true,
    then only one (N-dimensional) point will be read, otherwise system.pcoord_len points
    will be read.
    """
    
    system = wemd.rc.get_system_driver()
    
    assert fieldname == 'pcoord'
    
    pcoord = numpy.loadtxt(pcoord_return_filename, dtype=system.pcoord_dtype)
    
    if single_point:
        expected_shape = (system.pcoord_ndim,)
        if pcoord.ndim == 0:
            pcoord.shape = (1,)
    else:
        expected_shape = (system.pcoord_len, system.pcoord_ndim)
        if pcoord.ndim == 1:
            pcoord.shape = (len(pcoord),1)
    if pcoord.shape != expected_shape:
        raise ValueError('progress coordinate data has incorrect shape {!r} [expected {!r}]'.format(pcoord.shape,
                                                                                                    expected_shape))
    destobj.pcoord = pcoord

def aux_data_loader(fieldname, data_filename, segment, single_point):
    data = numpy.loadtxt(data_filename)
    segment.data[fieldname] = data
    if data.nbytes == 0:
        raise ValueError('could not read any data for {}'.format(fieldname))
    
    
class ExecutablePropagator(WEMDPropagator):
    ENV_CURRENT_ITER         = 'WEMD_CURRENT_ITER'
    
    # Environment variables set during propagation
    ENV_CURRENT_SEG_ID       = 'WEMD_CURRENT_SEG_ID'
    ENV_CURRENT_SEG_DATA_REF = 'WEMD_CURRENT_SEG_DATA_REF'
    ENV_CURRENT_SEG_INITPOINT= 'WEMD_CURRENT_SEG_INITPOINT_TYPE'
    ENV_PARENT_SEG_ID        = 'WEMD_PARENT_ID'
    ENV_PARENT_DATA_REF      = 'WEMD_PARENT_DATA_REF'
    
    # Environment variables set during propagation and state generation
    ENV_BSTATE_ID            = 'WEMD_BSTATE_ID'
    ENV_BSTATE_DATA_REF      = 'WEMD_BSTATE_DATA_REF'
    ENV_ISTATE_ID            = 'WEMD_ISTATE_ID'
    ENV_ISTATE_DATA_REF      = 'WEMD_ISTATE_DATA_REF'
    
    # Environment variables for progress coordinate calculation
    ENV_STRUCT_DATA_REF      = 'WEMD_STRUCT_DATA_REF'
    
    # Set everywhere a progress coordinate is required
    ENV_PCOORD_RETURN        = 'WEMD_PCOORD_RETURN'
    
    ENV_RAND16               = 'WEMD_RAND16'
    ENV_RAND32               = 'WEMD_RAND32'
    ENV_RAND64               = 'WEMD_RAND64'
    ENV_RAND128              = 'WEMD_RAND128'
    ENV_RANDFLOAT            = 'WEMD_RANDFLOAT'
        
    def __init__(self):
        super(ExecutablePropagator,self).__init__()
            
        # A mapping of environment variables to template strings which will be
        # added to the environment of all children launched.
        self.addtl_child_environ = dict()
        
        # A mapping of executable name ('propagator', 'pre_iteration', 'post_iteration') to 
        # a dictionary of attributes like 'executable', 'stdout', 'stderr', 'environ', etc.
        self.exe_info = {}
        self.exe_info['propagator'] = {}
        self.exe_info['pre_iteration'] = {}
        self.exe_info['post_iteration'] = {}
        self.exe_info['get_pcoord'] = {}
        self.exe_info['gen_istate'] = {}
        
        # A mapping of data set name ('pcoord', 'coord', 'com', etc) to a dictionary of
        # attributes like 'loader', 'dtype', etc
        self.data_info = {}
        self.data_info['pcoord'] = {}
 
        # Process configuration file information
        # Check in advance for keys we know we need
        for key in ('propagator', 'get_pcoord',
                    'segment_data_ref', 'parent_data_ref', 'basis_state_data_ref', 'initial_state_data_ref'):
            wemd.rc.config.require('executable.{}'.format(key))

        #self.exe_info['propagator']['executable'] = wemd.rc.config.get_path('executable.propagator')
        self.segment_ref_template = wemd.rc.config.get('executable.segment_data_ref')
        self.parent_ref_template  = wemd.rc.config.get('executable.parent_data_ref')
        self.initial_state_ref_template = wemd.rc.config.get('executable.initial_state_data_ref')
        self.basis_state_ref_template = wemd.rc.config.get('executable.basis_state_data_ref')
        
        # Load configuration items related to all child processes
        for (key,value) in wemd.rc.config.iteritems():
            if key.startswith('executable.env.'):
                self.addtl_child_environ[key[len('executable.env.'):]] = value
        log.debug('addtl_child_environ: {!r}'.format(self.addtl_child_environ))
        
        # Load configuration items relating to child processes
        for child_type in ('propagator', 'pre_iteration', 'post_iteration', 'get_pcoord', 'gen_istate'):
            child_key_prefix = 'executable.{}'.format(child_type)
            child_info = {key: value for (key,value) in wemd.rc.config.iteritems() if key.startswith(child_key_prefix)}
            if child_key_prefix in wemd.rc.config:
                self.exe_info[child_type]['executable'] = wemd.rc.config.get(child_key_prefix)
                self.exe_info[child_type]['stdin']  = wemd.rc.config.get('{}.stdin'.format(child_key_prefix, os.devnull))
                self.exe_info[child_type]['stdout'] = wemd.rc.config.get('{}.stdout'.format(child_key_prefix), None)
                self.exe_info[child_type]['stderr'] = wemd.rc.config.get('{}.stderr'.format(child_key_prefix), None)
                self.exe_info[child_type]['cwd'] = wemd.rc.config.get('{}.cwd'.format(child_key_prefix), None)
                
                if {key for key in child_info if key.startswith('{}.env.'.format(child_key_prefix))}:
                    # Strip leading 'executable.CHILD_TYPE.env.' from leading edge of entries
                    #log.debug('child_key_prefix: {!r}'.format(child_key_prefix))
                    offset = len('{}.env.'.format(child_key_prefix))
                    #log.debug('offset: {:d}'.format(offset))
                    self.exe_info[child_type]['environ'] = {key[offset:]: value 
                                                            for (key,value) in child_info.iteritems() 
                                                            if key.startswith('{}.env.'.format(child_key_prefix))}                                                            
        log.debug('exe_info: {!r}'.format(self.exe_info))
        
        # Load configuration items relating to data return
        data_info = {key: value for (key, value) in wemd.rc.config.iteritems() if key.startswith('executable.data.')}
        for (key, value) in data_info.iteritems():
            fields = key.split('.')
            try:
                dsname = fields[2]
                spec = fields[3]
            except IndexError:
                raise ValueError('invalid data specifier {!r}'.format(key))
            else:
                if spec not in ('enabled', 'loader'):
                    raise ValueError('invalid dataset option {!r}'.format(spec))
                            
            try:
                self.data_info[dsname][spec] = value
            except KeyError:
                self.data_info[dsname] = {spec: value}
        
        for dsname in self.data_info:
            if 'enabled' not in self.data_info[dsname]:
                self.data_info[dsname]['enabled'] = True
            else:
                self.data_info[dsname]['enabled'] = wemd.rc.config.get_bool('executable.data.{}.enabled'.format(dsname))
                
            if 'loader' not in self.data_info[dsname]:
                self.data_info[dsname]['loader'] = pcoord_loader if dsname == 'pcoord' else aux_data_loader
            else:
                self.data_info[dsname]['loader'] = wemd.rc.config.get_python_callable('executable.data.{}.loader'.format(dsname))
                
        if not self.data_info['pcoord']['enabled']:
            log.warning('configuration file requests disabling pcoord data collection; overriding')
            self.data_info['pcoord']['enabled'] = True
                        
        log.debug('data_info: {!r}'.format(self.data_info))
                
    @staticmethod                        
    def makepath(template, template_args = None,
                  expanduser = True, expandvars = True, abspath = False, realpath = False):
        template_args = template_args or {}
        path = template.format(**template_args)
        if expandvars: path = os.path.expandvars(path)
        if expanduser: path = os.path.expanduser(path)
        if realpath:   path = os.path.realpath(path)
        if abspath:    path = os.path.abspath(path)
        path = os.path.normpath(path)
        return path

    def random_val_env_vars(self):
        '''Return a set of environment variables containing random seeds. These are returned
        as a dictionary, suitable for use in ``os.environ.update()`` or as the ``env`` argument to
        ``subprocess.Popen()``. Every child process executed by ``exec_child()`` gets these.'''
        
        return {self.ENV_RAND16:               str(random.randint(0,2**16)),
                self.ENV_RAND32:               str(random.randint(0,2**32)),
                self.ENV_RAND64:               str(random.randint(0,2**64)),
                self.ENV_RAND128:              str(random.randint(0,2**128)),
                self.ENV_RANDFLOAT:            str(random.random())}
        
    def exec_child(self, executable, environ=None, stdin=None, stdout=None, stderr=None, cwd=None):
        '''Execute a child process with the environment set from the current environment, the
        values of self.addtl_child_environ, the random numbers returned by self.random_val_env_vars, and
        the given ``environ`` (applied in that order). stdin/stdout/stderr are optionally redirected.
        
        This function waits on the child process to finish, then returns
        (rc, rusage), where rc is the child's return code and rusage is the resource usage tuple from os.wait4()'''
        
        all_environ = dict(os.environ)
        all_environ.update(self.addtl_child_environ)
        all_environ.update(self.random_val_env_vars())
        all_environ.update(environ or {})
        
        stdin  = file(stdin, 'rb') if stdin else sys.stdin        
        stdout = file(stdout, 'wb') if stdout else sys.stdout
        if stderr == 'stdout':
            stderr = stdout
        else:
            stderr = file(stderr, 'wb') if stderr else sys.stderr
                
        # close_fds is critical for preventing out-of-file errors
        proc = subprocess.Popen([executable],
                                cwd = cwd,
                                stdin=stdin, stdout=stdout, stderr=stderr if stderr != stdout else subprocess.STDOUT,
                                close_fds=True, env=all_environ)

        # Wait on child and get resource usage
        (_pid, _status, rusage) = os.wait4(proc.pid, 0)
        # Do a subprocess.Popen.wait() to let the Popen instance (and subprocess module) know that
        # we are done with the process, and to get a more friendly return code
        rc = proc.wait()
        return (rc, rusage)
    
    def exec_child_from_child_info(self, child_info, template_args, environ):
        for (key, value) in child_info.get('environ', {}).iteritems():
            environ[key] = self.makepath(value)        
        return self.exec_child(executable = self.makepath(child_info['executable'], template_args),
                               environ = environ,
                               cwd = self.makepath(child_info['cwd'], template_args) if child_info['cwd'] else None,
                               stdin = self.makepath(child_info['stdin'], template_args) if child_info['stdin'] else os.devnull,
                               stdout= self.makepath(child_info['stdout'], template_args) if child_info['stdout'] else None,
                               stderr= self.makepath(child_info['stderr'], template_args) if child_info['stderr'] else None)        
        
    
    # Functions to create template arguments and environment values for child processes
    def update_args_env_basis_state(self, template_args, environ, basis_state):
        new_template_args = {'basis_state': basis_state}
        new_env = {self.ENV_BSTATE_ID: str(basis_state.state_id or -1),
                   self.ENV_BSTATE_DATA_REF: self.makepath(self.basis_state_ref_template, new_template_args)}
        template_args.update(new_template_args)
        environ.update(new_env)
        return template_args, environ
    
    def update_args_env_initial_state(self, template_args, environ, initial_state):
        new_template_args = {'initial_state': initial_state}
        new_env = {self.ENV_ISTATE_ID: str(initial_state.state_id or -1),
                   self.ENV_ISTATE_DATA_REF: self.makepath(self.initial_state_ref_template, new_template_args)}
        
        if initial_state.basis_state is not None:
            basis_state = initial_state.basis_state
        else:
            basis_state = self.basis_states[initial_state.basis_state_id]
          
        self.update_args_env_basis_state(new_template_args, new_env, basis_state)
        
        template_args.update(new_template_args)
        environ.update(new_env)
        return template_args, environ
    
    def update_args_env_iter(self, template_args, environ, n_iter):
        environ[self.ENV_CURRENT_ITER] = template_args['n_iter'] = str(n_iter if n_iter is not None else -1)
        return template_args, n_iter
    
    def update_args_env_segment(self, template_args, environ, segment):
        template_args['segment'] = segment
        
        environ[self.ENV_CURRENT_SEG_INITPOINT] = Segment.initpoint_type_names[segment.initpoint_type]
        
        if segment.initpoint_type == Segment.SEG_INITPOINT_CONTINUES:
            # Could use actual parent object here if the work manager cared to pass that much data
            # to us (we'd need at least the subset of parents for all segments sent in the call to propagate)
            # that may make a good wemd.cfg option for future crazy extensibility, but for now,
            # just populate the bare minimum
            parent = Segment(n_iter=segment.n_iter-1, seg_id=segment.parent_id)
            template_args['parent'] = parent
            
            environ[self.ENV_PARENT_SEG_ID] = str(segment.parent_id)            
            environ[self.ENV_PARENT_DATA_REF] = self.makepath(self.parent_ref_template, template_args)
        elif segment.initpoint_type == Segment.SEG_INITPOINT_NEWTRAJ:
            # This segment is initiated from a basis state; WEMD_PARENT_SEG_ID and WEMD_PARENT_DATA_REF are
            # set to the basis state ID and data ref
            initial_state = self.initial_states[segment.initial_state_id]
            basis_state = self.basis_states[initial_state.basis_state_id]
            
            if self.ENV_BSTATE_ID not in environ:
                self.update_args_env_basis_state(template_args, environ, basis_state)
            if self.ENV_ISTATE_ID not in environ:
                self.update_args_env_initial_state(template_args, environ, initial_state)
            
            assert initial_state.istate_type in (InitialState.ISTATE_TYPE_BASIS, InitialState.ISTATE_TYPE_GENERATED)
            if initial_state.istate_type == InitialState.ISTATE_TYPE_BASIS:
                environ[self.ENV_PARENT_DATA_REF] = environ[self.ENV_BSTATE_DATA_REF]
            else: # initial_state.type == InitialState.ISTATE_TYPE_GENERATED  
                environ[self.ENV_PARENT_DATA_REF] = environ[self.ENV_ISTATE_DATA_REF]
            
        environ[self.ENV_CURRENT_SEG_ID] = str(segment.seg_id or -1)
        environ[self.ENV_CURRENT_SEG_DATA_REF] = self.makepath(self.segment_ref_template, template_args)
        return template_args, environ
    
    def exec_for_segment(self, child_info, segment, addtl_env = None):
        '''Execute a child process with environment and template expansion from the given
        segment.'''
        template_args, environ = {}, {}
        self.update_args_env_iter(template_args, environ, segment.n_iter)
        self.update_args_env_segment(template_args, environ, segment)        
        environ.update(addtl_env or {})
        return self.exec_child_from_child_info(child_info, template_args, environ)
            
    def exec_for_iteration(self, child_info, n_iter, addtl_env = None):
        '''Execute a child process with environment and template expansion from the given
        iteration number.'''
        template_args, environ = {}, {}
        self.update_args_env_iter(template_args, environ, n_iter)
        environ.update(addtl_env or {})
        return self.exec_child_from_child_info(child_info, template_args, environ)

    def exec_for_basis_state(self, child_info, basis_state, addtl_env = None):
        '''Execute a child process with environment and template expansion from the
        given basis state'''
        template_args, environ = {}, {}
        self.update_args_env_basis_state(template_args, environ, basis_state)
        environ.update(addtl_env or {})
        return self.exec_child_from_child_info(child_info, template_args, environ)
        
    def exec_for_initial_state(self, child_info, initial_state,  addtl_env = None):
        '''Execute a child process with environment and template expansion from the given
        initial state.'''
        template_args, environ = {}, {}
        self.update_args_env_initial_state(template_args, environ, initial_state)
        environ.update(addtl_env or {})
        return self.exec_child_from_child_info(child_info, template_args, environ)

    # Specific functions required by the WEMD framework
    def get_pcoord(self, state):
        '''Get the progress coordinate of the given basis or initial state.'''
        
        template_args, environ = {}, {}
        
        if isinstance(state, BasisState):
            execfn = self.exec_for_basis_state
            self.update_args_env_basis_state(template_args, environ, state)
            struct_ref = environ[self.ENV_BSTATE_DATA_REF]
        elif isinstance(state, InitialState):
            execfn = self.exec_for_initial_state
            self.update_args_env_initial_state(template_args, environ, state)
            struct_ref = environ[self.ENV_ISTATE_DATA_REF]
        else:
            raise TypeError('state must be a BasisState or InitialState')
        
        child_info = self.exe_info.get('get_pcoord')
        fd, rfname = tempfile.mkstemp()
        os.close(fd)
        
        addtl_env = {self.ENV_PCOORD_RETURN: rfname,
                     self.ENV_STRUCT_DATA_REF: struct_ref}

        try:
            rc, rusage = execfn(child_info, state, addtl_env)
            if rc != 0:
                log.error('get_pcoord executable {!r} returned {}'.format(child_info['executable'], rc))
                
            loader = self.data_info['pcoord']['loader']
            loader('pcoord', rfname, state, single_point = True)
        finally:
            try:
                os.unlink(rfname)
            except Exception as e:
                log.warning('could not delete progress coordinate return file {!r}: {}'.format(rfname, e))
                
    def gen_istate(self, basis_state, initial_state):
        '''Generate a new initial state from the given basis state.'''
        child_info = self.exe_info.get('gen_istate')
        rc, rusage = self.exec_for_initial_state(child_info, initial_state)
        if rc != 0:
            log.error('gen_istate executable {!r} returned {}'.format(child_info['executable'], rc))
            initial_state.istate_status = InitialState.ISTATE_STATUS_FAILED
            return            
    
        # Determine and load the progress coordinate value for this state
        try:
            self.get_pcoord(initial_state)
        except:
            log.exception('could not get progress coordinate for initial state')
            initial_state.istate_status = InitialState.ISTATE_STATUS_FAILED
        else:
            initial_state.istate_status = InitialState.ISTATE_STATUS_PREPARED
                        
    def prepare_iteration(self, n_iter, segments):
        child_info = self.exe_info.get('pre_iteration')
        if child_info:
            try:
                rc, rusage = self.exec_for_iteration(child_info, n_iter)
            except OSError as e:
                log.warning('could not execute pre-iteration program {!r}: {}'.format(child_info['executable'], e))
            else:
                if rc != 0:
                    log.warning('pre-iteration executable {!r} returned {}'.format(child_info['executable'], rc))
        
    def finalize_iteration(self, n_iter, segments):
        child_info = self.exe_info.get('post_iteration')
        if child_info:
            try:
                rc, rusage = self.exec_for_iteration(child_info, n_iter)
            except OSError as e:
                log.warning('could not execute post-iteration program {!r}: {}'.format(child_info['executable'], e))
            else:
                if rc != 0:
                    log.warning('post-iteration executable {!r} returned {}'.format(child_info['executable'], rc))
        
                
    def propagate(self, segments):
        child_info = self.exe_info['propagator']
        
        for segment in segments:
            starttime = time.time()

            addtl_env = {}
            
            # A mappping of data set name to (filename, delete) pairs, where delete is a bool
            # indicating whether the file should be deleted (i.e. it's a temporary file) or not
            return_files = {}
            for dataset in self.data_info:
                (fd, rfname) = tempfile.mkstemp()
                os.close(fd)
                return_files[dataset] = rfname

                addtl_env['WEMD_{}_RETURN'.format(dataset.upper())] = return_files[dataset]
                                        
            # Spawn propagator and wait for its completion
            rc, rusage = self.exec_for_segment(child_info, segment, addtl_env) 
            
            if rc == 0:
                segment.status = Segment.SEG_STATUS_COMPLETE
            elif rc < 0:
                log.error('child process for segment %d exited on signal %d (%s)' % (segment.seg_id, -rc, SIGNAL_NAMES[-rc]))
                segment.status = Segment.SEG_STATUS_FAILED
                continue
            else:
                log.error('child process for segment %d exited with code %d' % (segment.seg_id, rc))
                segment.status = Segment.SEG_STATUS_FAILED
                continue
            
            # Extract data and store on segment for recording in the master thread/process/node
            for dataset in self.data_info:
                filename = return_files[dataset]
                loader = self.data_info[dataset]['loader']
                try:
                    loader(dataset, filename, segment, single_point=False)
                except Exception as e:
                    log.error('could not read {} from {!r}: {}'.format(dataset, filename, e))
                    segment.status = Segment.SEG_STATUS_FAILED 
                    break
                else:
                    try:
                        os.unlink(filename)
                    except Exception as e:
                        log.warning('could not delete {} file {!r}: {}'.format(dataset, filename, e))
                            
            if segment.status == Segment.SEG_STATUS_FAILED:
                continue
                                        
            # Record timing info
            segment.walltime = time.time() - starttime
            segment.cputime = rusage.ru_utime
        return segments