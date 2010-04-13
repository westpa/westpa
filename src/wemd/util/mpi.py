import logging, os
log = logging.getLogger(__name__)

def init_mpi():
    from wemd.util.extlogger import ExtendedLogger
    try:
        from mpi4py import rc
    except ImportError:
        return
    
    rc.initialize = False
    rc.finalize = False
    rc.threaded = False
    
    from mpi4py import MPI
    if not MPI.Is_initialized():
        log.info('initializing MPI environment')
        MPI.Init()
        
    ExtendedLogger.clsextra['nodename'] = getnodename()
    ExtendedLogger.clsextra['proc_rank'] = getrank()
    log.info('MPI environment initialized; rank %d is %s (process %d)'
             % (getrank(), getnodename(), os.getpid()))

def finalize_mpi():
    try:
        from mpi4py import MPI
    except ImportError:
        pass
    else:
        if not MPI.Is_finalized():
            log.info('finalizing MPI environment')
            MPI.Finalize()

def is_mpi_active():
    try:
        from mpi4py import MPI
    except ImportError:
        return False
    
    if MPI.COMM_WORLD.size == 1:
        return False
    else:
        return True

def is_rank_0():
    try:
        from mpi4py import MPI
    except ImportError:
        return None
    else:
        return bool(MPI.COMM_WORLD.rank == 0)
    
def getnodename():
    try:
        from mpi4py import MPI
    except ImportError:
        from socket import gethostname
        return gethostname()
    else:
        return MPI.Get_processor_name()

def getrank():
    try:
        from mpi4py import MPI
    except ImportError:
        return 0
    else:
        return MPI.COMM_WORLD.rank
