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



import logging

log = logging.getLogger(__name__)

import west
from oldtools.aframe import AnalysisMixin

class CommonOutputMixin(AnalysisMixin):
    def __init__(self):
        super(CommonOutputMixin,self).__init__()
        
        include_args = self.include_args.setdefault('CommonOutputMixin',{})
        include_args.setdefault('suppress_headers', True)
        include_args.setdefault('print_bin_labels', True)
        
        self.output_suppress_headers = False
        self.output_print_bin_labels = False
        
    def add_common_output_args(self, parser_or_group):
        if self.include_args['CommonOutputMixin']['suppress_headers']:
            parser_or_group.add_argument('--noheaders', dest='suppress_headers', action='store_true',
                                help='Do not include headers in text output files (default: include headers)')
        if self.include_args['CommonOutputMixin']['print_bin_labels']:
            parser_or_group.add_argument('--binlabels', dest='print_bin_labels', action='store_true',
                                help='Print bin labels in output files, if available (default: do not print bin labels)')
                
        
    def process_common_output_args(self, args):
        if self.include_args['CommonOutputMixin']['suppress_headers']:
            self.output_suppress_headers = bool(args.suppress_headers)
        if self.include_args['CommonOutputMixin']['print_bin_labels']:
            self.output_print_bin_labels = bool(args.print_bin_labels)
        