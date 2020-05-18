

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
    
    
