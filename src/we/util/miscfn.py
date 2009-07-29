import re

reIsSeconds = re.compile(r'^\d+(:?\.\d*)?$')

def parse_elapsed_time(timestr):
    timestr = timestr.strip()
    
    if reIsSeconds.match(timestr):
        return float(timestr)
    else:
        raise ValueError('invalid time string %r' % timestr)
    
