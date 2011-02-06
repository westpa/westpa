import os, sys, time, contextlib, tempfile
import numpy
import logging
log = logging.getLogger(__name__)

from wemd import Segment
from wemd.propagators import WEMDPropagator
from wemd.util.rtracker import ResourceTracker

@contextlib.contextmanager
def changed_cwd(target_dir):
    if target_dir:
        init_dir = os.getcwd()
        log.debug('chdir(%r)' % target_dir)
        os.chdir(target_dir)
        
    yield # to with block
    
    if target_dir:
        log.debug('chdir(%r)' % init_dir)
        os.chdir(init_dir)
    
class ExecutablePropagator(WEMDPropagator):
    EXTRA_ENVIRONMENT_PREFIX = 'executable.env.'

    ENV_CURRENT_ITER         = 'WEMD_CURRENT_ITER'

    ENV_CURRENT_SEG_ID       = 'WEMD_CURRENT_SEG_ID'
    ENV_CURRENT_SEG_DATA_REF = 'WEMD_CURRENT_SEG_DATA_REF'
    
    ENV_PARENT_SEG_ID        = 'WEMD_PARENT_SEG_ID'
    ENV_PARENT_SEG_DATA_REF  = 'WEMD_PARENT_SEG_DATA_REF'
    
    ENV_PCOORD_RETURN        = 'WEMD_PCOORD_RETURN'
        
    def __init__(self, sim_manager):
        super(ExecutablePropagator,self).__init__(sim_manager)
        
        self.rtracker = ResourceTracker()
        
        self.exename = None
    
        # Common environment variables for all child processes;
        # overridden by those specified per-executable
        self.child_environ = dict()

        # Information about child programs (executables, output redirections,
        # etc)
        self.propagator_info =      {'executable': None,
                                     'environ': dict(),
                                     'cwd': None}
        self.pre_iteration_info =   {'executable': None,
                                     'environ': dict(),
                                     'cwd': None}
        self.post_iteration_info =  {'executable': None,
                                     'environ': dict(),
                                     'cwd': None}
        
        # Process configuration file information
        runtime_config = self.sim_manager.runtime_config
        runtime_config.require('executable.propagator')
        self.segment_dir    = runtime_config.require('executable.segment_dir')
        self.parent_dir     = runtime_config.require('executable.parent_dir')
        self.initial_state_dir = runtime_config.get('executable.initial_state_dir', self.parent_dir)
                
        if 'executable.pcoord_loader' in runtime_config:
            from wemd.util import extloader
            pathinfo = runtime_config.get_pathlist('executable.module_path', default=None)
            self.pcoord_loader = extloader.get_object(runtime_config['executable.pcoord_loader'], 
                                                      pathinfo)
        else:
            self.pcoord_loader = None
        
        
        prefixlen = len(self.EXTRA_ENVIRONMENT_PREFIX)
        for (k,v) in runtime_config.iteritems():
            if k.startswith(self.EXTRA_ENVIRONMENT_PREFIX):
                evname = k[prefixlen:]                
                self.child_environ[evname] = v                
                log.info('including environment variable %s=%r for all child processes' % (evname, v))
        
        for child_type in ('propagator', 'pre_iteration', 'post_iteration'):
            child_info = getattr(self, child_type + '_info')
            child_info['child_type'] = child_type
            executable = child_info['executable'] = runtime_config.get('executable.%s' % child_type, None)            
            if executable:
                log.info('%s executable is %r' % (child_type, executable))

                stdout = child_info['stdout'] = runtime_config.get('executable.%s.stdout' % child_type, None)
                if stdout:
                    log.info('redirecting %s standard output to %r'% (child_type, stdout))
                stderr = child_info['stderr'] = runtime_config.get('executable.%s.stderr' % child_type, None)
                if stderr:
                    if stderr == 'stdout':
                        log.info('merging %s standard error with standard output' % child_type)
                    else:
                        log.info('redirecting %s standard error to %r' % (child_type, stderr))
                        
    def makepath(self, template, template_args = None,
                  expanduser = True, expandvars = True, abspath = False, realpath = False):
        path = template.format(**template_args)
        if expandvars: path = os.path.expandvars(path)
        if expanduser: path = os.path.expanduser(path)
        if realpath:   path = os.path.realpath(path)
        if abspath:    path = os.path.abspath(path)
        path = os.path.normpath(path)
        return path
        
    
    def _popen(self, child_info, addtl_environ = None, template_args = None):
        """Create a subprocess.Popen object for the appropriate child
        process, passing it the appropriate environment and setting up proper
        output redirections
        """
        
        template_args = template_args or dict()
                
        exename = self.makepath(child_info['executable'], template_args)
        child_type = child_info['child_type']
        child_environ = dict(os.environ)
        child_environ.update(self.child_environ)
        child_environ.update(addtl_environ or {})
        child_environ.update(child_info['environ'])
        
        log.debug('preparing to execute %r (%s) in %r' % (exename, child_info['child_type'], 
                                                          os.getcwd()))
        
        stdout = sys.stdout
        stderr = sys.stderr
        if child_info['stdout']:
            stdout_path = self.makepath(child_info['stdout'], template_args)
            log.debug('redirecting stdout to %r' % stdout_path)
            stdout = open(stdout_path, 'wb')
        if child_info['stderr']:
            if child_info['stderr'] == 'stdout':
                stderr = stdout
            else:
                stderr_path = self.makepath(child_info['stderr'], template_args)
                log.debug('redirecting standard error to %r' % stderr_path)
                stderr = open(stderr_path, 'wb')
                    
        log.debug('launching %s executable %r' % (child_type, exename))
        pid = os.fork()
        if pid:
            # in parent process
            #while True:
                #id, rc = os.waitpid(pid, os.WNOHANG)
                #if id == pid:
                #    break
                #time.sleep(1)
            (id,rc) = os.wait()    
            return rc
        else:
            # in child process
            # redirect stdout/stderr
            stdout_fd = stdout.fileno()
            stderr_fd = stderr.fileno()
            os.dup2(stdout_fd,1)
            os.dup2(stderr_fd,2)
                        
            # Execute
            os.execlpe(exename, 'wemd_worker[%s]' % os.path.basename(exename), child_environ)
    
    def _iter_env(self, n_iter):
        addtl_environ = {self.ENV_CURRENT_ITER: str(n_iter)}
        return addtl_environ
    
    def _segment_env(self, segment):
        template_args = self.segment_template_args(segment)
        parent_template = (self.parent_dir if segment.p_parent_id >= 0 
                           else (self.initial_state_dir or self.parent_dir))            
        addtl_environ = {self.ENV_CURRENT_ITER: str(segment.n_iter),
                         self.ENV_CURRENT_SEG_ID: str(segment.seg_id),
                         self.ENV_PARENT_SEG_ID: str(segment.p_parent_id),
                         self.ENV_CURRENT_SEG_DATA_REF: self.makepath(self.segment_dir, template_args),
                         self.ENV_PARENT_SEG_DATA_REF: self.makepath(parent_template, template_args)}
        return addtl_environ
    
    def segment_template_args(self, segment):
        template_args = {'n_iter': segment.n_iter,
                         'segment': segment}
        
        if segment.p_parent_id < 0:
            # (Re)starting from an initial state
            system = self.sim_manager.system
            istate = -segment.p_parent_id - 1
            parent_segment = Segment(seg_id = istate,
                                     n_iter = 0)
            template_args['parent'] = parent_segment
            template_args['initial_region_name'] = system.initial_states[istate][system.INITDIST_NAME]
        else:
            # Continuing from another segment
            parent_segment = Segment(seg_id = segment.p_parent_id,
                                     n_iter = segment.n_iter - 1)
        return template_args
    
    def iter_template_args(self, n_iter):
        return {'n_iter': n_iter}
    
    def _run_pre_post(self, child_info, env_func, template_func, args=(), kwargs={}):
        if child_info['executable']:
            rc = self._popen(child_info, env_func(*args, **kwargs), template_func(*args, **kwargs))
            if rc != 0:
                log.warning('%s executable %r returned %s'
                            % (child_info['child_type'], 
                               child_info['executable'],
                               rc))
            else:
                log.debug('%s executable exited successfully' 
                          % child_info['child_type'])
        
    def pre_iter(self, n_iter):
        self.rtracker.begin('pre_iter')
        with changed_cwd(self.pre_iteration_info['cwd']):
            self._run_pre_post(self.pre_iteration_info, self._iter_env, self.iter_template_args, args=(n_iter,))
        self.rtracker.end('pre_iter')

    def post_iter(self, n_iter):
        self.rtracker.begin('post_iter')
        with changed_cwd(self.post_iteration_info['cwd']):
            self._run_pre_post(self.post_iteration_info, self._iter_env, self.iter_template_args, args=(n_iter,))
        self.rtracker.end('post_iter')
            
    def prepare_iteration(self, n_iter, segments):
        self.pre_iter(n_iter)
        
    def finalize_iteration(self, n_iter, segments):
        self.post_iter(n_iter)
    
    def propagate(self, segments):
        #log.info('propagating %d segment(s)' % len(segments))
        for segment in segments:
            # Record start timing info
            self.rtracker.begin('propagation')

            # Fork the new process
            with changed_cwd(self.propagator_info['cwd']):
                log.debug('propagating segment %r' % segment)
                addtl_env = self._segment_env(segment)
                
                if not self.pcoord_loader:
                    (pc_return_fd, pc_return_filename) = tempfile.mkstemp()
                    log.debug('expecting return information in %r' % pc_return_filename)
                    os.close(pc_return_fd)
                    addtl_env[self.ENV_PCOORD_RETURN] = pc_return_filename
                
                # Spawn propagator and wait for its completion
                rc = self._popen(self.propagator_info, 
                                   addtl_env, 
                                   self.segment_template_args(segment))
                
                if rc == 0:
                    log.debug('child process for segment %d exited successfully'
                              % segment.seg_id)
                    segment.status = Segment.SEG_STATUS_COMPLETE
                else:
                    log.warn('child process for segment %d exited with code %s' 
                             % (segment.seg_id, rc))
                    segment.status = Segment.SEG_STATUS_FAILED
                    return
                
                # Extract progress coordinate
                if self.pcoord_loader:
                    segment.pcoord[...] = self.pcoord_loader(segment, self.sim_manager)
                else:
                    try:
                        segment.pcoord[...] = numpy.loadtxt(pc_return_filename, self.sim_manager.system.pcoord_dtype)
                    except Exception as e:
                        log.error('could not read progress coordinate from %r: %s' % (pc_return_filename, e))
                        segment.status = Segment.SEG_STATUS_FAILED
                    
                    try:
                        os.unlink(pc_return_filename)
                    except OSError as e:
                        log.warning('could not delete progress coordinate return file %r: %s' % (pc_return_filename, e))
                    else:
                        log.debug('deleted %r' % pc_return_filename)
            
            # Record end timing info
            self.rtracker.end('propagation')            
            segment.walltime = self.rtracker.difference['propagation'].walltime
            segment.cputime  = self.rtracker.difference['propagation'].cputime
