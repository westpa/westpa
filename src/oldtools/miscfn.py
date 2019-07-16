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


'''Miscellaneous support functions for WEST and WEST tools'''

import re

def parse_int_list(list_string):
    '''Parse a simple list consisting of integers or ranges of integers separated by commas. Ranges are specified
    as min:max, and include the maximum value (unlike Python's ``range``).  Duplicate values are ignored.
    Returns the result as a sorted list.  Raises ValueError if the list cannot be parsed.'''
    
    try:
        entries = set()
        fields = re.split(r'\s*[;,]\s*', list_string)
        for field in fields:
            if ':' in field:
                lb, ub = list(map(int,re.split(r'\s*:\s*', field)))
                entries.update(list(range(lb,ub+1)))
            else:
                entries.add(int(field))
    except (ValueError,TypeError):
        raise ValueError('invalid integer range string {!r}'.format(list_string))
    else:
        return sorted(entries)
    
    
