import os, sys, subprocess, time, datetime, tempfile
from resource import getrusage, RUSAGE_CHILDREN
import logging
log = logging.getLogger(__name__)

import numpy
import simplejson as json
from wemd.core import Segment
from wemd.backend_drivers import BackendDriver

class ExecutableBackend(BackendDriver):
    EXTRA_ENVIRONMENT_PREFIX = 'backend.executable.env.'

    ENV_CURRENT_ITER         = 'WEMD_CURRENT_ITER'
    ENV_SEG_DATA_RETURN      = 'WEMD_SEG_DATA_RETURN'

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
        self.pre_sim_info =         {'executable': None,
                                     'environ': dict()}
        self.post_sim_info =        {'executable': None,
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
                log.info('%s executable is %r' % (child_type, executable))

                stdout_template = child_info['stdout_template'] \
                                = runtime_config.get_compiled_template('backend.executable.%s.stdout_capture' % child_type,
                                                                       None)
                if stdout_template:
                    log.info('redirecting %s standard output to %r'
                             % (child_type, stdout_template.template))
                stderr_template = child_info['stderr_template'] \
                                = runtime_config.get_compiled_template('backend.executable.%s.stderr_capture' % child_type,
                                                                       None)
                if stderr_template:
                    log.info('redirecting %s standard error to %r'
                             % (child_type, stderr_template.template))
                    
                merge_stderr = child_info['merge_stderr'] \
                             = runtime_config.get_bool('backend.executable.%s.merge_stderr_to_stdout' % child_type, 
                                                       False)
                if merge_stderr:
                    log.info('merging %s standard error with standard output'
                             % child_type)
                
                if stderr_template and merge_stderr:
                    log.warning('both standard error redirection and merge specified for %s; standard error will be merged' % child_type)
                    child_info['stderr_template'] = None
    
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
            log.info('redirecting child stdout to %r' % stdout)
            stdout = open(stdout, 'wb')
        if child_info['stderr_template']:
            stderr = child_info['stderr_template'].safe_substitute(template_args)
            log.info('redirecting child stderr to %r' % stderr)
            stderr = open(stderr, 'wb')
        elif child_info['merge_stderr']:
            stderr = subprocess.STDOUT
            log.info('merging child stderr to stdout')
        
        log.debug('launching %s executable %r with environment %r' 
                  % (child_type, exename, child_environ))
        proc = subprocess.Popen([exename], stdout = stdout, stderr = stderr,
                                env = child_environ)
        return proc        

    def pre_sim(self):
        pass
    
    def post_sim(self):
        pass
        
    def pre_iter(self, we_iter):
        pass
    
    def post_iter(self, we_iter):
        pass
    
    def pre_segment(self, segment):
        pass
    
    def post_segment(self, segment):
        pass

    
    def propagate_segments(self, segments):
        log.debug('propagating %d segment(s)' % len(segments))
        for segment in segments:
            # Create a temporary file for the child process to return information
            # to us
            (return_fd, return_filename) = tempfile.mkstemp()
            log.debug('expecting return information in %r' % return_filename)
            os.close(return_fd)
            
            # Set up the child process environment
            addtl_environ = {}
            addtl_environ[self.ENV_CURRENT_ITER] = str(segment.we_iter)
            addtl_environ[self.ENV_SEG_DATA_RETURN] = return_filename
            addtl_environ[self.ENV_CURRENT_SEG_DATA_REF] = segment.data_ref
            addtl_environ[self.ENV_CURRENT_SEG_ID] = str(segment.seg_id)
            if segment.p_parent:
                addtl_environ[self.ENV_PARENT_SEG_ID] = str(segment.p_parent.seg_id)
                addtl_environ[self.ENV_PARENT_SEG_DATA_REF] = segment.p_parent.data_ref or ''
            
            # Fork the new process
            proc = self._popen(self.propagator_info, addtl_environ, 
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
                stream = open(return_filename)
            except OSError, e:
                log.error('could not open output return file: %s' % e)
                segment.status = Segment.SEG_STATUS_FAILED
            else:
                try:            
                    self.update_segment_from_output(segment, stream)
                finally:
                    stream.close()
                log.info('segment %s progress coordinate = %s'
                         % (segment.seg_id, segment.pcoord))
                    
            try:
                os.unlink(return_filename)
            except OSError, e:
                log.warn('could not delete output return file: %s' % e)
            else:
                log.debug('deleted output return file %r' % return_filename)
            
        
    def update_segment_from_output(self, segment, stream):
        try:
            stream.seek(0)
        except AttributeError:
            pass
        
        retdata = json.load(stream)
        try:
            pcoord = retdata['pcoord']
        except KeyError:
            pcoord = retdata
            
        try:
            iter(pcoord)
        except TypeError:
            pcoord = [pcoord]
            
        segment.pcoord = numpy.array(pcoord)
            
        
        
        