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

'''Miscellaneous routines to help with input and output of WEST-related data in text format'''

import time, sys, getpass, socket
        
class NumericTextOutputFormatter(object):
    comment_string = '# '
    emit_header = True

    def __init__(self, output_file, mode='wt', emit_header=None):
        if hasattr(output_file, 'write'):
            self._file = output_file
        else:
            self._file = open(output_file, mode)
            
        self._header_written = False
        self._header_lines = []

        self.write_header('created by program: {}'.format(sys.argv[0] or 'unknown program'))
        self.write_header('created by user:    {}'.format(getpass.getuser()))
        self.write_header('created on host:    {}'.format(socket.gethostname()))
        self.write_header('created at:         {}'.format(time.strftime('%c', time.localtime(time.time()))))
        
        if emit_header is not None:
            self.emit_header = emit_header
        else:
            self.emit_header = self.__class__.emit_header
        
    def __getattr__(self, attr):
        return getattr(self._file,attr)

    def close(self):
        if not self._header_written:
            self._write_header()
        self._file.close()
    
    def write(self, str):
        if not self._header_written:
            self._write_header()
        self._file.write(str)
    
    def writelines(self, sequence):
        if not self._header_written:
            self._write_header()
        self._file.writelines(sequence)
        
    def write_comment(self, line):
        '''Writes a line beginning with the comment string'''
        self._file.write('{}{}\n'.format(self.comment_string, line))
        
    def write_header(self, line):
        '''Appends a line to those written when the file header is written. The
        appropriate comment string will be prepended, so ``line`` should not include
        a comment character.'''
        if self._header_written:
            raise EnvironmentError('cannot append header lines after header has been written')
        else:
            self._header_lines.append(line)        
        
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._file.close()
        return False
        
    def _write_header(self):
        if self.emit_header and not self._header_written:
            for line in self._header_lines:
                self.write_comment(line)
        self._header_written = True
