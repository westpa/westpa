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
