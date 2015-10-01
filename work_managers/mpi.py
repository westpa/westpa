# TODO: Add license 
#
#
from __future__ import division, print_function; __metaclass__ = type

import logging, sys, threading, time, traceback

from collections import deque
from mpi4py import MPI
from core import WorkManager, WMFuture

log = logging.getLogger( __name__ )

# +------+
# | Task |
# +------+
class Task:
    """Tasks are tuples of (task_id, function, args, keyword args)
    """
    def __init__(self, task_id, fn, args, kwargs):
        self.task_id = task_id
        self.fn = fn
        self.args = args if args is not None else ()
        self.kwargs = kwargs if kwargs is not None else {}
        
    def __repr__(self):
        return '<Task {self.task_id}: {self.fn!r}(*{self.args!r}, **{self.kwargs!r})>'\
               .format(self=self)


# +----------------+
# | MPIWorkManager |
# +----------------+
class MPIWorkManager( WorkManager ):
    """MPIWorkManager factory.
    """

    @classmethod
    def from_environ( cls, wmenv=None ):
        return cls()
        

    def __new__( cls ):
        """Creates a Serial WorkManager if size is 1.  Otherwise creates a 
        single Master and size-1 Slaves.
        """
        log.debug( 'MPIWorkManager.__new__()' )
        assert( MPI.Is_initialized() )
        assert( MPI.Is_thread_main() )
        
        rank = MPI.COMM_WORLD.Get_rank()
        size = MPI.COMM_WORLD.Get_size()
        
        if size == 1:
            return super( MPIWorkManager, cls ).__new__( Serial )
        elif rank == 0:
            return super( MPIWorkManager, cls ).__new__( Master )
        else:
            return super( MPIWorkManager, cls ).__new__( Slave )

    
    def __init__( self ):
        """Initialize info shared by Master and Slave classes.
        """
        log.debug( 'MPIWorkManager.__init__()' )
        
        super( MPIWorkManager, self ).__init__()
        comm = MPI.COMM_WORLD
        self.comm = MPI.COMM_WORLD
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()
        self.name = MPI.Get_processor_name()
        
        # some tags
        self.task_tag     = 110 # tag for server to client msgs
        self.result_tag   = 120 # tag for client to server msgs
        self.shutdown_tag = 130 # tag for server to client to stop
        
        self.masterID = 0
        
        
    def submit( self, fn, args=None, kwargs=None ):
        """Adhere to WorkManager interface.  This method should never be 
        called.
        """
        assert( False )


# +--------+
# | Serial |
# +--------+
# TODO: no need for the code replication here, just use the original serial wm
class Serial( MPIWorkManager ):
    """Replication of the serial work manager.  This is a fallback for MPI runs 
    that request only 1 (size=1) processor.
    """
    
    def __init__( self ):
        super( Serial, self ).__init__()
        log.debug( 'Serial.__init__()' )
        
        
    def submit( self, fn, args=None, kwargs=None ):
        log.debug( 'Serial.__init__()' )
        
        ft = WMFuture()
        try:
            result = fn( *(args if args is not None else ()), 
                        **(kwargs if kwargs is not None else {}) )
        except Exception as e:
            ft._set_exception( e, sys.exc_info()[2] )
        else:
            ft._set_result( result )
        return ft


# +--------+
# | Master |
# +--------+
class Master( MPIWorkManager ):
    """Master of the MPIWorkManage.  Distributes tasks to Slaves as they are 
    received from the sim_manager.  In addition to the main thread, this class 
    spawns two threads, a receiver and a dispatcher.
    """

    def __init__( self ):
        """Initialize different state variables used by Master.
        """
        super( Master, self ).__init__()
        log.debug( 'Master__init__()' )  

        # number of slaves
        self.nslaves = self.size - 1
        
        # list of slave ranks
        self.slaveIDs = range( 1, self.size )
        assert self.nslaves == len( self.slaveIDs )
        
        # deque of idle slaves
        self.dests = deque( self.slaveIDs ) 
        
        # deque of tesks
        self.tasks = deque()
        
        # number of unmatched send/receive pairs
        self.nPending = 0
        
        # thread shutdown sentinel
        self.shutItDown = False
        
        # task_id, future key value pair
        self.pending_futures = dict()
        
        # list of master threads
        self.workers = []
        
        # thread lock
        self.lock = threading.Lock()

        
    def startup( self ):
        """Spawns the dispatcher and receiver threads.
        """
        log.debug( 'Master.startup()' )
        
        if not self.running:
            self.workers.append( threading.Thread( name='dispatcher', 
                                                   target=self._dispatcher) )
            self.workers.append( threading.Thread( name='receiver', 
                                                   target=self._receiver ) )
            
            for t in self.workers:
                t.start()
                log.info( 'Started thread: %s' % t.getName() )
                
            self.running = True
        
        
    def _dispatcher( self ):
        """Continuously dispatches tasks to idle destinations until the 
        shutdown sentinel is set.
        """
        log.debug( 'Master._dispatcher()' )
        assert( MPI.Is_thread_main() == False )
        assert( threading.currentThread().getName() == "dispatcher" )
        
        while not self.shutItDown:
            
            req = []
            # do we have work and somewhere to send it?
            while self.tasks and self.dests:
                
                with self.lock:
                    task = self.tasks.popleft()
                    sendTo = self.dests.popleft()
                    self.nPending += 1
                
                req.append( self.comm.isend( task, 
                                             dest=sendTo, 
                                             tag=self.task_tag ) )
                
            # make sure all sends completed
            MPI.Request.Waitall( requests=req )
            
            # force context switching ( 1ms )
            time.sleep( 0.001 )
        
    
    def _receiver( self ):
        """Continuously receives futures from slaves until the shutdown 
        sentinel is set.
        """
        log.debug( 'Master._receiver()' )
        assert( MPI.Is_thread_main() == False )
        assert( threading.currentThread().getName() == "receiver" )
        
        while not self.shutItDown:
            
            # are we waiting on any results?
            while self.nPending:
            
                stat = MPI.Status()
                ( tid, msg, val ) = self.comm.recv( source=MPI.ANY_SOURCE, 
                                                    tag=self.result_tag, 
                                                    status=stat )
                log.debug( 'Master._receiver received task: %s' % tid )
            
                # update future
                ft = self.pending_futures.pop( tid )
                if msg == 'exception':
                    ft._set_exception( *val )
                else:
                    ft._set_result( val )
            
                with self.lock:
                    self.dests.append( stat.Get_source() )
                    self.nPending -= 1
                    
            # force context switching ( 1ms )
            time.sleep( 0.001 )
    

    def submit( self, fn, args=None, kwargs=None ):
        """Receive task from simulation manager and add it to pending_futures.
        """
        log.debug( 'Master.submit()' )
        
        ft = WMFuture()
        task_id = ft.task_id
        with self.lock:
            self.tasks.append( Task( task_id, fn, args, kwargs ) )
            self.pending_futures[task_id] = ft
        
        return ft


    def shutdown( self ):
        """Send shutdown tag to all slave processes, and set the shutdown 
        sentinel to stop the receiver and dispatcher loops.
        """
        log.debug( 'Master.shutdown()' )
        
        # wait on any unfinished work
        while self.pending_futures:
            pass
        
        # send shutdown msg to all slaves
        req = [ MPI.REQUEST_NULL ]*self.nslaves        
        for rank in self.slaveIDs:
            req[rank-1] = self.comm.isend( MPI.BOTTOM, 
                                           dest=rank, 
                                           tag=self.shutdown_tag )
            
        MPI.Request.Waitall( requests=req )
        
        # stop threads
        self.shutItDown = True
        
        for t in self.workers:
            t.join()
            log.info( 'Stopped thread: %s' % t.getName() )
        
        self.running = False
    

# +-------+
# | Slave |
# +-------+
class Slave( MPIWorkManager ):
    """Client class for executing tasks as distributed by the Master in the
    MPI Work Manager
    """

    def __init__( self ):
        super( Slave, self ).__init__()
        log.debug( 'Slave.__init__() %s' % self.rank )
        
    
    def startup( self ):
        """Clock the slave in for work.
        """
        log.debug( 'Slave.startup() %s' % self.rank )
        if not self.running:
            
            self.clockIn()
            
            self.running = True

    
    def clockIn( self ):
        """Do each task as it comes in.  The completion of a task is
        notice to the master that more work is welcome.
        """
        log.info( 'Slave %s clocking in.' % self.rank )
        
        comm = self.comm
        
        while True:
            
            stat = MPI.Status()
            task = comm.recv( source=self.masterID, 
                             tag=MPI.ANY_TAG,
                             status=stat )
                             
            tag = stat.Get_tag()
            
            if tag == self.task_tag:
                
                log.debug( 'Slave %s received task: %s' % 
                           ( self.rank, task.task_id ) )
                
                # do the work
                try:
                    rv = task.fn( *task.args, **task.kwargs )
                except BaseException as e:
                    ro = ( task.task_id, 
                           'exception', 
                           ( e, traceback.format_exc() ) )
                else:
                    ro = ( task.task_id, 'result', rv )
                
                # send result back to master
                comm.send( ro, dest=self.masterID, tag=self.result_tag )
                
            if tag == self.shutdown_tag:
                log.info( 'Slave %s clocking out.' % self.rank )
                return


    @property
    def is_master( self ):
        """Slave processes need to be marked as not master.  This ensures that
        the proper branching is followed in w_run.py.
        """
        return False
        
