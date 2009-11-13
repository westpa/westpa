import logging
import mpi
from logging import Logger, LogRecord


class ExtendedLogger(logging.getLoggerClass()):
    clsextra = {}

    def __init__(self, name):
        Logger.__init__(self,name)
        self.extra = dict(self.clsextra)
        
    def makeRecord(self, name, level, fn, lno, msg, args, exc_info, func=None,
                   extra=None):
        allextra = dict(self.extra)
        allextra.update(extra or {})
        rec = LogRecord(name, level, fn, lno, msg, args, exc_info, func)
        rec.__dict__.update(allextra)
        return rec


ExtendedLogger.clsextra['nodename'] = mpi.getnodename()
ExtendedLogger.clsextra['proc_rank'] = mpi.getrank()
logging.setLoggerClass(ExtendedLogger)
