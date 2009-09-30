import re

reIsSeconds = re.compile(r'^\d+(:?\.\d*)?$')
reISODate = re.compile(r'(\d+)-(\d+)-(\d+)[ T](\d+):(\d+):(\d+)\.?(d+)?')

def parse_elapsed_time(timestr):
    timestr = timestr.strip()
    
    if reIsSeconds.match(timestr):
        return float(timestr)
    else:
        raise ValueError('invalid time string %r' % timestr)
    
def datetime_from_iso(dtstr):
    from datetime import datetime
    m = reISODate.match(dtstr)
    ints = [int(field) for field in m.groups() if field is not None]
    return datetime(*ints)