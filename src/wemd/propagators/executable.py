import os, sys, time
import logging
log = logging.getLogger(__name__)

from wemd import Segment
from wemd.propagators import WEMDPropagator
from wemd.util.rtracker import ResourceTracker

class ExecutablePropagator(WEMDPropagator):
    EXTRA_ENVIRONMENT_PREFIX = 'executable.env.'

    ENV_CURRENT_ITER         = 'WEMD_CURRENT_ITER'

    ENV_CURRENT_SEG_ID       = 'WEMD_CURRENT_SEG_ID'
    ENV_CURRENT_SEG_DATA_REF = 'WEMD_CURRENT_SEG_DATA_REF'
    
    ENV_PARENT_SEG_ID        = 'WEMD_PARENT_SEG_ID'
    ENV_PARENT_SEG_DATA_REF  = 'WEMD_PARENT_SEG_DATA_REF'
        
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
                                     'environ': dict()}
        self.pre_iteration_info =   {'executable': None,
                                     'environ': dict()}
        self.post_iteration_info =  {'executable': None,
                                     'environ': dict()}
        self.pre_segment_info =     {'executable': None,
                                     'environ': dict()}
        self.post_segment_info =    {'executable': None,
                                     'environ': dict()}
        
        self.sim_manager.runtime_config.require('executable.executable')
    
        self.segref_template = self.sim_manager.runtime_config.get_compiled_template('executable.segref_template', 
                                                                         'traj_segs/${n_iter}/${seg_id}')
        self.istateref_template = self.sim_manager.runtime_config.get_compiled_template('executable.istateref_template', 
                                                                            'initdist/${region_name}')

                                    
        try:
            if self.sim_manager.runtime_config.get_bool('executable.preserve_environment'):
                log.info('including parent environment')
                log.debug('parent environment: %r' % os.environ)
                self.child_environ.update(os.environ)
        except KeyError:
            pass
        
        prefixlen = len(self.EXTRA_ENVIRONMENT_PREFIX)
        for (k,v) in self.sim_manager.runtime_config.iteritems():
            if k.startswith(self.EXTRA_ENVIRONMENT_PREFIX):
                evname = k[prefixlen:]                
                self.child_environ[evname] = v                
                log.info('including environment variable %s=%r for all child processes' % (evname, v))
        
        for child_type in ('executable', 'pre_iteration', 'post_iteration',
                           'pre_segment', 'post_segment'):
            child_info = getattr(self, child_type + '_info')
            child_info['child_type'] = child_type
            executable = child_info['executable'] = self.sim_manager.runtime_config.get('executable.%s' % child_type, None)            
            if executable:
                log.info('%s executable is %r' % (child_type, executable))

                stdout_template = child_info['stdout_template'] \
                                = self.sim_manager.runtime_config.get_compiled_template('backend.executable.%s.stdout_capture' 
                                                                                        % child_type, None)
                if stdout_template:
                    log.info('redirecting %s standard output to %r'
                             % (child_type, stdout_template.template))
                stderr_template = child_info['stderr_template'] \
                                = self.sim_manager.runtime_config.get_compiled_template('backend.executable.%s.stderr_capture' 
                                                                                        % child_type, None)
                if stderr_template:
                    log.info('redirecting %s standard error to %r'
                             % (child_type, stderr_template.template))
                    
                merge_stderr = child_info['merge_stderr'] \
                             = self.sim_manager.runtime_config.get_bool('backend.executable.%s.merge_stderr_to_stdout' 
                                                                        % child_type, False)
                if merge_stderr:
                    log.info('merging %s standard error with standard output'
                             % child_type)
                
                if stderr_template and merge_stderr:
                    log.warning('both standard error redirection and merge specified for %s; standard error will be merged' 
                                % child_type)
                    child_info['stderr_template'] = None

    def make_data_ref(self, segment):
        return self.segref_template.safe_substitute(segment.__dict__)

    def make_istate_data_ref(self, segment):
        # Fetch the name associated with the region
        system = self.sim_manager.system
        assert segment.p_parent_id < 0
        istate = -segment.p_parent_id - 1
        
        subdict = dict(segment.__dict__)
        subdict['region_name'] = system.initial_states[istate][system.INITDIST_NAME]
        return self.istateref_template.safe_substitute(subdict)
    
    def _popen(self, child_info, addtl_environ = None, template_args = None):
        """Create a subprocess.Popen object for the appropriate child
        process, passing it the appropriate environment and setting up proper
        output redirections
        """
        
        template_args = template_args or dict()
        
        exename = child_info['executable']
        child_type = child_info['child_type']
        child_environ = dict(self.child_environ)
        child_environ.update(addtl_environ or {})
        child_environ.update(child_info['environ'])
        
        stdout = None
        stderr = None
        if child_info['stdout_template']:
            stdout = child_info['stdout_template'].safe_substitute(template_args)
            if child_info['merge_stderr']:
                log.debug('redirecting child stdout and stderr to %r' % stdout)
            else:
                log.debug('redirecting child stdout to %r' % stdout)
            stdout = open(stdout, 'wb')
        if child_info['stderr_template']:
            stderr = child_info['stderr_template'].safe_substitute(template_args)
            log.debug('redirecting child stderr to %r' % stderr)
            stderr = open(stderr, 'wb')
        elif child_info['merge_stderr']:
            stderr = sys.stdout
                    
        log.debug('launching %s executable %r' % (child_type, exename))
        pid = os.fork()
        if pid:
            # in parent process
            while True:
                id, rc = os.waitpid(pid, os.WNOHANG)
                if id == pid:
                    break
                
                time.sleep(1)
                
            return rc
        else:
            # in child process
            # redirect stdout/stderr
            stderr_fd = stderr.fileno()
            stdout_fd = stdout.fileno()
            
            os.dup2(stdout_fd, 1)
            os.dup2(stderr_fd, 2)
            
            # Execute
            os.execlpe(exename, child_environ)
    
    def _iter_env(self, n_iter, ):
        addtl_environ = {self.ENV_CURRENT_ITER: str(n_iter)}
        return addtl_environ
    
    def _segment_env(self, segment):

        addtl_environ = {self.ENV_CURRENT_ITER: str(segment.n_iter),
                         self.ENV_CURRENT_SEG_DATA_REF: self.make_data_ref(segment),
                         self.ENV_CURRENT_SEG_ID: str(segment.seg_id),
                         self.ENV_PARENT_SEG_ID: str(segment.p_parent.seg_id),}
        
        if segment.p_parent_id < 0:
            # Restarting from an initial state    
            addtl_environ[self.ENV_PARENT_SEG_DATA_REF] = self.make_istate_data_ref(segment)
            
        return addtl_environ
    
    def _segment_template_args(self, segment):
        return dict(segment.__dict__)
    
    def _iter_template_args(self, n_iter):
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
        self._run_pre_post(self.pre_iteration_info, self._iter_environ, self._iter_template_args, args=(n_iter,))
        self.rtracker.end('pre_iter')

    def post_iter(self, n_iter):
        self.rtracker.begin('post_iter')
        self._run_pre_post(self.post_iteration_info, self._iter_environ, self._iter_template_args, args=(n_iter,))
        self.rtracker.end('post_iter')
    
    def pre_segment(self, segment):
        self.rtracker.begin('pre_segment')
        self._run_pre_post(self.pre_segment_info, self._segment_env, self._segment_template_args, args=(segment,))
        self.rtracker.end('pre_segment')
    
    def post_segment(self, segment):
        self.rtracker.begin('post_segment')
        self._run_pre_post(self.post_segment_info, self._segment_env, self._segment_template_args, args=(segment,))
        self.rtracker.end('post_segment')
        
    def prepare_iteration(self, n_iter, segments):
        self.pre_iter(n_iter)
        
    def finalize_iteration(self, n_iter, segments):
        self.post_iter(n_iter)
    
    def propagate(self, segments):
        #log.info('propagating %d segment(s)' % len(segments))
        for segment in segments:
            # Record start timing info
            self.rtracker.begin('propagation')
            
            self.pre_segment(segment)
            
            # Fork the new process
            log.debug('propagating segment %d' % segment.seg_id)
            addtl_env = self._segment_env(segment)
            

            rc = self._popen(self.propagator_info, 
                               addtl_env, 
                               segment.__dict__)
            
            if rc == 0:
                log.debug('child process for segment %d exited successfully'
                          % segment.seg_id)
                segment.status = Segment.SEG_STATUS_COMPLETE
            else:
                log.warn('child process for segment %d exited with code %s' 
                         % (segment.seg_id, rc))
                segment.status = Segment.SEG_STATUS_FAILED
                return
                
            self.post_segment(segment)

            # Record end timing info
            self.rtracker.end('propagation')            
            segment.walltime = self.rtracker.difference['walltime']
            segment.cputime  = self.rtracker.difference['cputime']
