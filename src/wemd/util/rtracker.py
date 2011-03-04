from __future__ import division; __metaclass__ = type

import time, resource, operator
from itertools import izip, starmap
from collections import namedtuple
from contextlib import contextmanager

ResourceUsage = namedtuple('ResourceUsage', ['walltime', 'cputime', 'usertime', 'systime',
                                             'c_usertime', 'c_systime'] )

class ResourceTracker:
    '''A simple class to automate the tracking of consumable resources (particularly CPU and wallclock time)
    for a number of different tasks'''
    
    def __init__(self):
        self.difference = dict()
        self.initial = dict()
        self.final = dict()
    
    @contextmanager
    def tracking(self, label):
        self.begin(label)
        yield
        self.end(label)
    
    def reset(self, label=None):
        if label is None:
            self.initial = dict()
            self.final = dict()
            self.difference = dict()
        else:
            for collection in (self.initial, self.final, self.difference):
                try:
                    del collection[label]
                except KeyError:
                    pass
    
    def begin(self, label):
        walltime, cputime = time.time(), time.clock()
        ru_self = resource.getrusage(resource.RUSAGE_SELF)
        ru_children = resource.getrusage(resource.RUSAGE_CHILDREN)
        
        self.initial[label] = ResourceUsage(walltime, cputime, ru_self.ru_utime, ru_self.ru_stime,
                                             ru_children.ru_utime, ru_children.ru_stime)
    
    def end(self, label):
        walltime, cputime = time.time(), time.clock()
        ru_self = resource.getrusage(resource.RUSAGE_SELF)
        ru_children = resource.getrusage(resource.RUSAGE_CHILDREN)
        
        final = self.final[label] = ResourceUsage(walltime, cputime, ru_self.ru_utime, ru_self.ru_stime,
                                                  ru_children.ru_utime, ru_children.ru_stime)
        self.difference[label] = ResourceUsage(*starmap(operator.sub, izip(final, self.initial[label])))
        
    def dump_differences(self, filename_or_stream):
        if filename_or_stream is None:
            # return a string
            from cStringIO import StringIO
            stream = StringIO()
        else:        
            try:
                filename_or_stream.write
            except AttributeError:
                stream = open(filename_or_stream, 'wt')
            else:
                stream = filename_or_stream
        
        max_label_len = max(map(len, self.difference.viewkeys()))
        field_width=14
        stream.write((' '.join(['#{:<{max_label_width}}'] + 6*['{:>{field_width}}'])+'\n'
                      ).format('Label', 'walltime', 'cputime', 'usertime', 'systime','c_usertime', 'c_systime',
                               max_label_width = max_label_len,
                               field_width = field_width))
        for (label, usage) in self.difference.viewitems():
            stream.write((' '.join([' {label:<{max_label_width}}', 
                                   '{diff.walltime:{field_width}.6f}',
                                   '{diff.cputime:{field_width}.6f}',
                                   '{diff.usertime:{field_width}.6f}',
                                   '{diff.systime:{field_width}.6f}',
                                   '{diff.c_usertime:{field_width}.6f}',
                                   '{diff.c_systime:{field_width}.6f}'])+'\n').format(label=label, diff=usage,
                                                                   max_label_width=max_label_len,
                                                                   field_width=field_width))
        try:
            stream.flush()
        except AttributeError:
            pass
        
        if filename_or_stream is None:
            return stream.getvalue()
    