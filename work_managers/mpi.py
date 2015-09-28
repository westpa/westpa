# TODO: License 
#
#
from __future__ import division, print_function; __metaclass__ = type

import logging, sys, re, threading

from collections import deque
from mpi4py import MPI
from core import WorkManager, WMFuture
#from . import WorkManager, WMFuture

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
    """MPI work manager factory.
    """

    @classmethod
    def from_environ( cls, wmenv=None ):
        return cls()
        

    def __new__( cls ):
        print( 'MPIWorkManager.__new__()' )
        print( 'Is_initialized(): %s' % MPI.Is_initialized() )
        print( 'Is_thread_main(): %s' % MPI.Is_thread_main() )
        rank = MPI.COMM_WORLD.Get_rank()
        size = MPI.COMM_WORLD.Get_size()
        
        
        if size == 1:
            return super( MPIWorkManager, cls ).__new__( Serial )
        elif rank == 0:
            return super( MPIWorkManager, cls ).__new__( Master )
        else:
            return super( MPIWorkManager, cls ).__new__( Slave )

    
    def __init__( self ):
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
        pass


# +--------+
# | Serial |
# +--------+
class Serial( MPIWorkManager ):
    """Replication of the serial work manager.  This is a fallback for MPI runs 
    that request only 1 processor.
    """
    
    def __init__( self ):
        super( Serial, self ).__init__()
        print( 'Serial construction' )
        
        
    def submit( self, fn, args=None, kwargs=None ):
        # TODO: add docstring
        """
        """
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
    """Master class of the MPI work manager for distributing tasks as received 
    through the sim manager through submit.
    """

    def __init__( self ):
        super( Master, self ).__init__()
        log.debug( 'Master construction' )
        print( 'Master.__init__()' )

        # TODO: doc        
        self.nslaves   = self.size - 1        
        self.slaveIDs  = range( 1, self.size )
        
        self.dests     = deque( self.slaveIDs ) 
        self.tasks     = deque()
        self.nPending  = 0 # number of unmatched send/recv pairs
        
        self.shutItDown = False
        
        assert self.nslaves == len( self.slaveIDs )
        
        self.pending_futures = dict()
        self.workers = []

        
    def startup( self ):
        # TODO: add docstring
        """
        """
        log.debug( 'Master.startup()' )
        print( 'Master.startup()' )
        
        if not self.running:
            self.workers.append( threading.Thread( name='dispatcher', target=self._dispatcher) )
            self.workers.append( threading.Thread( name='receiver', target=self._receiver ) )
            
            for t in self.workers:
                log.info( 'Starting thread: %s' % t.getName() )
                t.start()
                
            self.running = True
        
        
    def _dispatcher( self ):
        # TODO: add docstring
        """
        """
        log.debug( 'Master._dispatcher()' )
        assert( MPI.Is_thread_main() == False )
        assert( threading.currentThread().getName() == "dispatcher" )
        
        while not self.shutItDown:
            
            req = []
            # do we have work and somewhere to send it?
            while self.tasks and self.dests:
                
                ## TODO: do we need to lock here
                task = self.tasks.popleft()
                sendTo = self.dests.popleft()
                self.nPending += 1
                print( "_dispatcher.nPending: %s" % self.nPending )
                ## unlock
                
                req.append( self.comm.isend( task, dest=sendTo, tag=self.task_tag ) )
                
            # make sure all sends completed
            MPI.Request.Waitall( requests=req )
        
    
    def _receiver( self ):
        # TODO: add docstring
        """
        """
        log.debug( 'Master._receiver()' )
        assert( MPI.Is_thread_main() == False )
        assert( threading.currentThread().getName() == "receiver" )
        
        while not self.shutItDown:
            
            # are we waiting on any results?
            while self.nPending:
            
                stat = MPI.Status()
                ( tid, msg, val ) = self.comm.recv( source=MPI.ANY_SOURCE, 
                                        tag=self.result_tag, status=stat )
                print( "Master._receiver() received result" )
            
                ft = self.pending_futures.pop( tid )
                if msg == 'exception':
                    ft._set_exception( val )
                else:
                    ft._set_result( val )
            
                ## TODO: do we need a lock here
                self.dests.append( stat.Get_source() )
                self.nPending -= 1
                print( "_receiver.nPending: %s" % self.nPending )
                ## unlock
        
    

    def submit( self, fn, args=None, kwargs=None ):
        """Receive task from simulation manager.
        """
        log.debug( 'Master.submit()' )
        print( 'Master.submit()' )
        
        ft = WMFuture()
        task_id = ft.task_id
        self.tasks.append( Task( task_id, fn, args, kwargs ) )
        self.pending_futures[task_id] = ft
        
        return ft


    def shutdown( self ):
        """Send shutdown tag to all slave processes.
        """
        log.debug( 'Master.shutdown()' )
        print( 'Master.shutdown()' )
        
        # wait on any unfinished work
        while self.pending_futures:
            pass
        
        # send shutdown msg to all slaves
        req = [ MPI.REQUEST_NULL ]*self.nslaves        
        for rank in self.slaveIDs:
            req[rank-1] = self.comm.isend( MPI.BOTTOM, dest=rank, tag=self.shutdown_tag )
            
        MPI.Request.Waitall( requests=req )
        
        # stop threads
        self.shutItDown = True
        
        for t in self.workers:
            print( 'Stopping thread: %s' % t.getName() )
            log.info( 'Stopping thread: %s' % t.getName() )
            t.join()
        
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
        log.debug( 'Slave construction %s' % self.rank )
        print("I am the slave. %s" % self.rank )
        
        # get ready to work
        self.clockIn()

    
    def clockIn( self ):
        """Do each task as it comes in.  The completion of a task is
        notice to the master that more work is welcome.
        """
        log.debug( 'Slave.clockIn()' )
        
        log.info( 'Slave %s clocking in.' % self.rank )
        comm = self.comm
        
        while True:
            
            print( "Slave waiting for work..." )
            stat = MPI.Status()
            task = comm.recv( source=self.masterID, 
                             tag=MPI.ANY_TAG,
                             status=stat )
            print ( "Slave got work" )

            log.info( "Slave.clockIn() received task." )
                             
            tag = stat.Get_tag()
            print( "Slave msg tag: %s" % tag )
            
            if tag == self.task_tag:
                
                # do the work
                try:
                    rv = task.fn( *task.args, **task.kwargs )
                    log.debug( "Slave.clockIn() does this function evaluate?" )
                except Exception:
                    # TODO: better return value?
                    ro = ( task.task_id, 'exception', None )
                else:
                    ro = ( task.task_id, 'result', rv )
                
                # send result back to master
                comm.isend( ro, dest=self.masterID, tag=self.result_tag )
                
            if tag == self.shutdown_tag:
                print( 'Slave %s clocking out' % self.rank )
                log.info( 'Slave %s clocking out' % self.rank )
                return


    @property
    def is_master( self ):
        """Slave processes need to be marked as not master.  This ensures that
        the proper branching is followed in w_run.py.
        """
        return False



"""
m = MPIWorkManager()
    
with m:
    
    print( '%s' % MPI.Is_thread_main() )
    # submit a series of tasks as in the code
    for task_id in range( 5 ):
        t = Task( task_id, max, None, None )
        m.submit( t.fn, t.args, t.kwargs )
"""