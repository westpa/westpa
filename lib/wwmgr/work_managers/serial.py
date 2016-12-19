# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
#
# This file is part of WESTPA.
#
# WESTPA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# WESTPA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with WESTPA.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division; __metaclass__ = type

import logging, sys

log = logging.getLogger(__name__)

from . import WorkManager, WMFuture

class SerialWorkManager(WorkManager):
    @classmethod
    def from_environ(cls, wmenv=None):
        return cls()
    
    def __init__(self):
        log.debug('initializing serial work manager')
        super(SerialWorkManager,self).__init__()
        self.n_workers = 1
        
    def submit(self, fn, args=None, kwargs=None):
        ft = WMFuture()
        try:
            result = fn(*(args if args is not None else ()), **(kwargs if kwargs is not None else {}))
        except Exception as e:
            ft._set_exception(e, sys.exc_info()[2])
        else:
            ft._set_result(result)
        return ft
    