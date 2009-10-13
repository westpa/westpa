import os, sys, subprocess, time, datetime, tempfile
from resource import getrusage, RUSAGE_CHILDREN
import logging
log = logging.getLogger(__name__)

import numpy
import simplejson as json

from wemd.core import Segment
from wemd.backend_drivers import BackendDriver

# This class should be instantiated in client processes
class ExecutableBackend(BackendDriver):
    EXTRA_ENVIRONMENT_PREFIX = 'backend.executable.environment.'

    ENV_CURRENT_ITER         = 'WEMD_CURRENT_ITER'
    ENV_SEG_DATA_RETURN      = 'WEMD_SEG_DATA_RETURN'

    ENV_CURRENT_SEG_ID       = 'WEMD_CURRENT_SEG_ID'
    ENV_CURRENT_SEG_DATA_REF = 'WEMD_CURRENT_SEG_DATA_REF'
    
    ENV_PARENT_SEG_ID        = 'WEMD_PARENT_SEG_ID'
    ENV_PARENT_SEG_DATA_REF  = 'WEMD_PARENT_SEG_DATA_REF'
    
    def __init__(self):
        super(ExecutableBackend, self).__init__()
        self.exename = None
        self.child_environ = None
        
    def initialize(self, runtime_config):
        assert(runtime_config['backend.driver'].lower() == 'executable')
        runtime_config.require('backend.executable')
        self.exename = runtime_config['backend.executable']
        log.info('executable for segment propagation: %s' % self.exename)
        self.child_environ = dict()
        
        try:
            if runtime_config.get_bool('backend.executable.preserve_environment'):
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
                log.debug('including environment variable %s=%r' 
                          % (evname, v))
    
    def propagate_segment(self, segment):
        
        # Create a temporary file for the child process to return information
        # to us
        (return_fd, return_filename) = tempfile.mkstemp()
        os.close(return_fd)
        
        # Set up the child process environment
        child_environ = self.child_environ.copy()
        child_environ[self.ENV_CURRENT_ITER] = str(segment.we_iter)
        child_environ[self.ENV_SEG_DATA_RETURN] = return_filename
        child_environ[self.ENV_CURRENT_SEG_DATA_REF] = segment.data_ref
        child_environ[self.ENV_CURRENT_SEG_ID] = str(segment.seg_id)
        if segment.p_parent:
            child_environ[self.ENV_PARENT_SEG_ID] = str(segment.p_parent.seg_id)
            child_environ[self.ENV_PARENT_SEG_DATA_REF] = segment.p_parent.data_ref or ''
        
        # Record start timing info
        log.debug('launching %r with environment %r' % (self.exename, child_environ))
        segment.starttime = datetime.datetime.now()
        log.debug('launched at %s' % segment.starttime)
        init_walltime = time.time()
        init_cputime = getrusage(RUSAGE_CHILDREN).ru_utime
        
        proc = subprocess.Popen([self.exename], env=child_environ)
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
            return
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
            
        
        
        