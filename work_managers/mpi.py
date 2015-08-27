import logging, sys

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
    """MPI work manager factory.
    """


    def __new__( cls ):
        rank = MPI.COMM_WORLD.Get_rank()
        size = MPI.COMM_WORLD.Get_size()
        if size == 1:
            return super( MPIWorkManager, cls ).__new__( Serial )
        elif rank == 0:
            return super( MPIWorkManager, cls ).__new__( Master )
        else:
            return super( MPIWorkManager, cls ).__new__( Slave )

    
    def __init__( self ):
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
        print( 'Master construction')

        # doc        
        self.nslaves   = self.size - 1        
        self.slaveIDs  = range( 1, self.size )
        self.init_dests = self.slaveIDs[:]
        assert self.nslaves == len( self.slaveIDs )
        

    def _dispatch( self, task ):
        """The first self.nslaves tasks are dispatched to the individual
        slaves.  Subsequently, slaves get more work when they are finished 
        with a given task.
        """
        log.debug( 'Master._dispatch() %s' % self.rank )
        
        # the first self.nclient jobs
        if self.init_dests:
            # send work
            self.comm.isend( task, dest=self.init_dests.pop(), tag=self.task_tag )
        else:
            # who is done?
            stat = MPI.Status()
            ( tid, msg, val ) = self.comm.recv( source=MPI.ANY_SOURCE, tag=self.result_tag, status=stat )
            
            # send the new task to the now idle client
            dest = stat.Get_source()
            self.comm.isend( task, dest=dest, tag=self.task_tag )
            
            # TODO: update the future of the received result
            pass
            
        # TODO: check the status of the non-blocking sends
        
    

    def submit( self, fn, args=None, kwargs=None ):
        """Receive task from simulation manager.
        """
        log.debug( 'Master.submit()' )
        print( 'Master.submit()' )
        
        ft = WMFuture()
        task_id = ft.task_id
        task = Task( task_id, fn, args, kwargs )
        
        self._dispatch( task )


    def shutdown( self ):
        """Send shutdown tag to all slave processes.
        """
        log.debug( 'Master.shutdown()' )
        print( 'Master.shutdown()' )
        
        req = [ MPI.REQUEST_NULL ]*self.nslaves       
        
        for rank in self.slaveIDs:
            req[rank-1] = self.comm.isend( MPI.BOTTOM, dest=rank, tag=self.shutdown_tag )
            
        MPI.Request.Waitall( requests=req )
    

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

        # TODO: remove - temporary for testing        
        import time        
        
        log.debug( 'Slave.clockIn() %s' % self.rank )
        log.info( 'Slave %s clocking in' % self.rank )
        comm = self.comm
        
        while True:
            
            stat = MPI.Status()
            task = comm.recv( source=self.masterID, 
                             tag=MPI.ANY_TAG,
                             status=stat )
                             
            tag = stat.Get_tag()
            
            if tag == self.task_tag:
                
                # do the work
                try:
                    rv = task.fn( *task.args, **task.kwargs )
                except Exception:
                    # TODO: better return value?
                    ro = ( task.task_id, 'exception', None )
                else:
                    ro = ( task.task_id, 'result', rv )
                
                # TODO: remove - temporary for testing
                # simulate some work done
                time.sleep( 5 )
                print( "I am client %s.  I just did some work." % self.rank )
                
                # send result back to master
                comm.isend( ro, dest=self.masterID, tag=self.result_tag )
                
            if tag == self.shutdown_tag:
                log.info( 'Slave %s clocking out' % self.rank )
                return
            









#
m = MPIWorkManager()
    
with m:
    # submit a series of tasks as in the code
    for task_id in range( 10 ):
        t = Task( task_id, max, None, None )
        m.submit( t.fn, t.args, t.kwargs )


