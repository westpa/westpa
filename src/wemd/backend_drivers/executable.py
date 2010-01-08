import os, sys, subprocess, time, datetime, tempfile, string
from resource import getrusage, RUSAGE_CHILDREN
import logging
log = logging.getLogger(__name__)

import numpy
from wemd.core import Segment
from wemd.backend_drivers import BackendDriver

class ExecutableBackend(BackendDriver):
    EXTRA_ENVIRONMENT_PREFIX = 'backend.executable.env.'

    ENV_CURRENT_ITER         = 'WEMD_CURRENT_ITER'
    ENV_SEG_PCOORD_RETURN    = 'WEMD_SEG_PCOORD_RETURN'

    ENV_CURRENT_SEG_ID       = 'WEMD_CURRENT_SEG_ID'
    ENV_CURRENT_SEG_DATA_REF = 'WEMD_CURRENT_SEG_DATA_REF'
    
    ENV_PARENT_SEG_ID        = 'WEMD_PARENT_SEG_ID'
    ENV_PARENT_SEG_DATA_REF  = 'WEMD_PARENT_SEG_DATA_REF'
    
    def __init__(self, runtime_config):
        super(ExecutableBackend, self).__init__(runtime_config)
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
        
        assert(runtime_config['backend.driver'].lower() == 'executable')
        runtime_config.require('backend.executable.propagator')
        
        drtemplate = self.runtime_config.setdefault('backend.executable.segref_template', 
                                                   'traj_segs/${we_iter}/${seg_id}')
        ctemplate = string.Template(drtemplate)
        try:
            ctemplate.safe_substitute(dict())
        except ValueError, e:
            raise ConfigError('invalid data ref template %r' % drtemplate)
        else:
            self.segref_template = ctemplate
            
        pcoord_file_format = self.runtime_config.get('backend.executable.pcoord_file.format', 'text')
        if pcoord_file_format != 'text':
            raise ConfigError('invalid pcoord file format %r' % pcoord_file_format)

        try:
            if runtime_config\
            .get_bool('backend.executable.preserve_environment'):
                log.info('including parent environment')
                log.debug('parent environment: %r' % os.environ)
                self.child_environ.update(os.environ)
        except KeyError:
            pass
        
        prefixlen = len(self.EXTRA_ENVIRONMENT_PREFIX)
        for (k,v) in runtime_config.iteritems():
            if k.startswith(self.EXTRA_ENVIRONMENT_PREFIX):
                evname = k[prefixlen:]                
                self.child_environ[evname] = v                
                log.debug('including environment variable %s=%r for all child processes' 
                          % (evname, v))
        
        for child_type in ('propagator', 'pre_iteration', 'post_iteration',
                           'pre_segment', 'post_segment'):
            child_info = getattr(self, child_type + '_info')
            child_info['child_type'] = child_type
            executable = child_info['executable'] \
                       = runtime_config.get('backend.executable.%s' 
                                            % child_type, None)            
            if executable:
                log.debug('%s executable is %r' % (child_type, executable))

                stdout_template = child_info['stdout_template'] \
                                = runtime_config.get_compiled_template('backend.executable.%s.stdout_capture' % child_type,
                                                                       None)
                if stdout_template:
                    log.debug('redirecting %s standard output to %r'
                             % (child_type, stdout_template.template))
                stderr_template = child_info['stderr_template'] \
                                = runtime_config.get_compiled_template('backend.executable.%s.stderr_capture' % child_type,
                                                                       None)
                if stderr_template:
                    log.debug('redirecting %s standard error to %r'
                             % (child_type, stderr_template.template))
                    
                merge_stderr = child_info['merge_stderr'] \
                             = runtime_config.get_bool('backend.executable.%s.merge_stderr_to_stdout' % child_type, 
                                                       False)
                if merge_stderr:
                    log.debug('merging %s standard error with standard output'
                             % child_type)
                
                if stderr_template and merge_stderr:
                    log.warning('both standard error redirection and merge specified for %s; standard error will be merged' % child_type)
                    child_info['stderr_template'] = None

    def make_data_ref(self, segment):
        return self.segref_template.safe_substitute(segment.__dict__)

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
            stderr = subprocess.STDOUT
        
        if log.getEffectiveLevel() <= logging.DEBUG:
            log.debug('launching %s executable %r with environment %r' 
                      % (child_type, exename, child_environ))
        else:
            log.debug('launching %s executable %r'
                     % (child_type, exename))
        proc = subprocess.Popen([exename], stdout = stdout, stderr = stderr,
                                env = child_environ)
        return proc
    
    def _iter_environ(self, we_iter):
        addtl_environ = {self.ENV_CURRENT_ITER: str(we_iter.n_iter)}
        return addtl_environ
    
    def _segment_env(self, segment):
        addtl_environ = {self.ENV_CURRENT_ITER: str(segment.n_iter),
                         self.ENV_CURRENT_SEG_DATA_REF: self.make_data_ref(segment),
                         self.ENV_CURRENT_SEG_ID: str(segment.seg_id),}
        if segment.p_parent:
            addtl_environ[self.ENV_PARENT_SEG_ID] = str(segment.p_parent.seg_id)
            addtl_environ[self.ENV_PARENT_SEG_DATA_REF] = self.make_data_ref(segment.p_parent)        
        return addtl_environ
    
    def _run_pre_post(self, child_info, env_func, env_obj):
        if child_info['executable']:
            proc = self._popen(child_info, env_func(env_obj), env_obj.__dict__)
            rc = proc.wait()
            if rc != 0:
                log.warning('%s executable %r returned %s'
                            % (child_info['child_type'], 
                               child_info['executable'],
                               rc))
            else:
                log.debug('%s executable exited successfully' 
                          % child_info['child_type'])
        
    def pre_iter(self, we_iter):
        self._run_pre_post(self.pre_iteration_info, self._iter_environ, we_iter)

    def post_iter(self, we_iter):
        self._run_pre_post(self.post_iteration_info, self._iter_environ, we_iter)
    
    def pre_segment(self, segment):
        self._run_pre_post(self.pre_segment_info, self._segment_env, segment)
    
    def post_segment(self, segment):
        self._run_pre_post(self.post_segment_info, self._segment_env, segment)
    
    def propagate_segments(self, segments):
        log.debug('propagating %d segment(s)' % len(segments))
        for segment in segments:
            self.pre_segment(segment)
            # Create a temporary file for the child process to return 
            # progress coordinate information to us
            (pc_return_fd, pc_return_filename) = tempfile.mkstemp()
            log.debug('expecting return information in %r' % pc_return_filename)
            os.close(pc_return_fd)
            
            # Fork the new process
            addtl_env = self._segment_env(segment)
            addtl_env[self.ENV_SEG_PCOORD_RETURN] = pc_return_filename
            proc = self._popen(self.propagator_info, 
                               addtl_env, 
                               segment.__dict__)
            # Record start timing info
            segment.starttime = datetime.datetime.now()
            init_walltime = time.time()
            init_cputime = getrusage(RUSAGE_CHILDREN).ru_utime
            log.debug('launched at %s' % segment.starttime)
            
            rc = proc.wait()
            
            # Record end timing info
            final_cputime = getrusage(RUSAGE_CHILDREN).ru_utime
            final_walltime = time.time()
            segment.endtime = datetime.datetime.now()
            segment.walltime = final_walltime - init_walltime
            segment.cputime = final_cputime - init_cputime
            log.debug('completed at %s (wallclock %s, cpu %s)' 
                      % (segment.endtime,
                         segment.walltime,
                         segment.cputime))
            
            if rc == 0:
                log.debug('child process exited successfully')
                segment.status = Segment.SEG_STATUS_COMPLETE
            else:
                log.warn('child process exited with code %s' % rc)
                segment.status = Segment.SEG_STATUS_FAILED
                return
                                
            try:
                self.update_pcoord_from_output(segment, pc_return_filename)
            except Exception, e:
                log.error('could not read progress coordinate file: %s' % e)
                segment.status = Segment.SEG_STATUS_FAILED
                
            try:
                os.unlink(pc_return_filename)
            except OSError, e:
                log.warn('could not delete progress coordinate file: %s' % e)
            else:
                log.debug('deleted progress coordinate file %r' % pc_return_filename)
                
            self.post_segment(segment)
            
        
    def update_pcoord_from_output(self, segment, pc_return_filename):
        pcarray =  numpy.loadtxt(pc_return_filename, dtype=numpy.float64)
        
        #FIXME: Need to communicate/store timestep in some elegant way
        #(or really, ANY WAY AT ALL)
        if self.runtime_config.get_bool('backend.executable.pcoord_file.eliminate_time_column', True):
            if len(pcarray) > 1:
                segment.addtl_attrs['t0'] = pcarray[0,0]
                segment.addtl_attrs['dt'] = pcarray[1,0] - pcarray[0,0] 
            segment.pcoord = pcarray[:,1:]
        else:
            segment.pcoord = pcarray
