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

from __future__ import print_function, division; __metaclass__ = type
import threading, time, blessings #@UnresolvedImport
from collections import deque
import numpy
from scipy.stats import linregress

def nop():
    pass

class ProgressIndicator:
    def __init__(self, stream=None, interval = 1):
        self.terminal = blessings.Terminal(stream=stream)
        self.interval = interval # how often the status is updated
        
        self._operation = None # string describing operation
        self._extent = None # how far we have to go, total
        self._progress = 0 # how far we've gone
        
        self._endloop = False
        self._startup_event = threading.Event()
        self._event = threading.Event()
        
        self._last_update = None
        self._last_operation = None
        self._operation_start = None
        self._progress_history = deque(maxlen=100)
        
        try:
            self.flush_output = self.terminal.stream.flush
        except AttributeError:
            self.flush_output = nop
            
        self.fancy = self.terminal.is_a_tty
        
    def _completion_est(self):
        if self._extent is None:
            return 'unknown'

        history = numpy.array(self._progress_history)
        if not len(history):
            return 'unknown'
        history[:,1] -= self._extent
        (_slope, intercept, _r, _p, _stderr) = linregress(history[:,1], history[:,0])
        if not numpy.isfinite(intercept):
            return 'unknown'
        else:
            total_seconds_remaining = intercept - time.time()
            
            if total_seconds_remaining < 60:
                return 'less than 1 minute'
            
            # 1 minute or more remains
            minutes = total_seconds_remaining / 60
            if minutes < 60:
                minutes = int(round(minutes))
                return 'about {:d} {}'.format(minutes, 'minutes' if minutes > 1 else 'minute')
            
            # 1 hour or more remains
            hours, minutes = divmod(minutes, 60)
            days, hours = divmod(hours, 24)
            days = int(days)
            hours = int(round(hours))
            minutes = int(round(minutes))
            if days > 0:
                return 'about {:d} {}, {:d} {}'.format(days, 'days' if days > 1 else 'day', 
                                                       hours, 'hours' if hours > 1 else 'hour')
            else:
                return 'about {:d} {}, {:d} {}'.format(hours, 'hours' if hours > 1 else 'hour',
                                                       minutes, 'minutes' if minutes > 1 else 'minute') 
    def draw_fancy(self):
        operation_text = 'Operation:      '
        progress_text  = 'Progress:       '
        remaining_text = 'Time remaining: '

        term = self.terminal
        stream = self.terminal.stream
        
        stream.write('{t.clear_eol}{t.bold}{}{t.normal}{}\n'.format(operation_text,self._operation or '', t=term))
            
        if self._extent:
            width = term.width 
            pct_done = self.progress / self._extent
            pct_part = '{:<4.0%} '.format(pct_done)
            
            barwidth = width - len(progress_text) - 2 - len(pct_part)
            neqs = int(round(pct_done*barwidth))
            eqs = '=' * neqs
            spaces = ' ' * (barwidth - neqs)
            
            stream.write('{t.clear_eol}{t.bold}{progress_text}{t.normal}{pct_part}{t.bold}[{t.normal}{eqs}{spaces}{t.bold}]{t.normal}\n'
                         .format(t=term,progress_text=progress_text,pct_part=pct_part,eqs=eqs,spaces=spaces))

            completion_est = self._completion_est()
            stream.write('{t.clear_eol}{t.bold}{}{t.normal}{}\n'.format(remaining_text,completion_est,t=term))
            stream.write('{t.move_up}'.format(t=term)*2)
        stream.write('{t.move_up}'.format(t=term))
    
        self._last_update = time.time()

        
    def draw_simple(self):
        if self._operation != self._last_operation:
            self.terminal.stream.write((self._operation or '(unknown)')+'...\n')
            self._last_operation = self._operation
            self._last_update = time.time()
            
        
    def draw(self):
        if not self._operation: return
        if self.fancy:
            self.draw_fancy()
        else:
            self.draw_simple()
            
    def clear(self):
        if self.fancy:
            if self._extent:
                nlines = 3
            else:
                nlines = 1

            for _i in xrange(nlines):
                self.terminal.stream.write('{t.clear_eol}{t.move_down}'.format(t=self.terminal))
            for _i in xrange(nlines):
                self.terminal.stream.write('{t.move_up}'.format(t=self.terminal))

    @property
    def operation(self):
        return self._operation
        
    @operation.setter
    def operation(self, op):
        if self._operation is not None: self.clear()
        self._operation = op
        self._operation_start = time.time()
        self._progress_history.clear()
        self._event.set()
        
    @property
    def extent(self):
        return self._extent

    @extent.setter
    def extent(self, ext):
        self._extent = ext
        self._progress_history.clear()
        self._event.set()
        
    @property
    def progress(self):
        return self._progress
    
    @progress.setter
    def progress(self, p):
        self._progress = p
        self._progress_history.append((time.time(), p))
    
    def new_operation(self, operation, extent=None, progress=0):
        self.operation = operation
        self.progress = progress
        self._extent = extent
        self._event.set()

    def _reporter_loop(self):
        self._startup_event.set()
        while not self._endloop:
            self._event.wait(self.interval)
            self._event.clear()
            if self._last_operation != self._operation or time.time() - self._last_update >= self.interval:
                self.draw()
                self.flush_output()
            

    def start(self):
        if self.fancy:
            self.terminal.stream.write(self.terminal.hide_cursor)
        self._endloop = False
        t = threading.Thread(target=self._reporter_loop)
        t.daemon = True
        t.start()
        self._startup_event.wait()
        self.flush_output()
        
    def stop(self):
        self._endloop = True
        self._event.set()
        
        self.clear()
        if self.fancy:
            self.terminal.stream.write(self.terminal.normal_cursor)
        self.flush_output()
        
        
    def __enter__(self):
        self.start()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop()
        return False

if __name__ == '__main__':
    with ProgressIndicator() as pi:
        pi.operation='Test 1'
        pi.extent = 10
        for i in xrange(10):
            pi.progress = i+1
            time.sleep(2)
            
        time.sleep(0.2)


    
