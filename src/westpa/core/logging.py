"""Logging utilities"""

import logging


class ConsecutiveDuplicateFilter(logging.Filter):
    def filter(self, record):
        current_log = (record.module, record.levelno, record.getMessage())
        if current_log != getattr(self, "last_log", None):
            self.last_log = current_log
            return True
        return False
