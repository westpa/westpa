from __future__ import division
import os, sys, re, math
from itertools import izip
from optparse import OptionParser
import cPickle as pickle

class ArraySaver(object):
    text_cell_format = '%-21.17g'
    pickle_protocol = pickle.HIGHEST_PROTOCOL
    default_file_formats = ('text', 'pickle')
    format_exts = {'text': '.txt',
                   'pickle': '.dat'}
    format_fileformats = {'text': 't',
                          'pickle': 'b'}
    
    def __init__(self, file_formats = None,
                 ignore_errors = False):
        self.file_formats = file_formats or self.default_file_formats
        self.ignore_errors = ignore_errors
        
    def save_array(self, data, filename_template, template_params, 
                   file_formats = None, ignore_errors = False):

        file_formats = file_formats or self.file_formats
        tparams = dict(template_params.iteritems())
        
        for fmt in self.file_formats:
            tparams['file_format'] = fmt
            saver = getattr(self, 'save_fmt_%s' % fmt)
            try:
                tparams['format_ext'] = self.format_exts[fmt]
            except IndexError:
                pass
            filename = filename_template % tparams
            
            try:
                stream = open(filename, 'w' + self.format_fileformats[fmt])
            except OSError, e:
                if ignore_errors:
                    continue
                else:
                    raise
            else:
                saver(data, stream)
                stream.close()
        
    def save_fmt_text(self, data, stream):
        if len(data.shape) > 2:
            raise TypeError('cannot save array of rank > 2 as text')
        elif len(data.shape) == 0:
            stream.write((self.text_cell_format % data) + '\n')
            return
        elif len(data.shape) == 1:
            for row in data:
                stream.write((self.text_cell_format % row) + '\n')
            return
        else:
            cols = data.shape[1]
            rowfmt = ' '.join([self.text_cell_format] * cols) + '\n'
            for row in data:
                stream.write(rowfmt % tuple(row))
            return
        
    def save_fmt_pickle(self, data, stream):
        pickle.dump(data, stream, self.pickle_protocol)
        
        

def write_array_as_text(data, stream, format='-15.8g'):
    if len(data.shape) > 2:
        raise TypeError('cannot save array of rank > 2 as text')
    rows = data.shape[0]
    cols = data.shape[1]
    
    rowfmt = ' '.join([format] * cols) + '\n'
    
    for row in data:
        stream.write(rowfmt % row)

parser = OptionParser(usage = '%prog PC_ARRAY TRANS_LB TRANS_UB')
(opts, args) = parser.parse_args()

if len(args) != 3:
    parser.print_help(sys.stderr)
    sys.exit(2)
    
import numpy, wemd

dat = numpy.load(args[0])
if len(dat.shape) < 2:
    Q = dat
    t = numpy.fromiter(xrange(0,dat.shape[0]), numpy.float_)
    sys.stderr.write('Warning: no time values available\n')
else:
    t = dat[:,0]
    Q = dat[:,1]
assert len(Q) == len(t)

Q_trans_lb = float(args[1])
Q_trans_ub = float(args[2])

te_finder = wemd.analysis.transitions.PC1DCrossingFinder(Q_trans_lb, Q_trans_ub)
te_finder.assign_regions(Q)
te_finder.find_crossings()
te_finder.find_transitions()

asaver = ArraySaver()

for (region1, region2) in sorted(te_finder.transition_indices):
    filename_params = {'region1': region1,
                       'region2': region2}
    ti_array = te_finder.transition_indices[region1, region2]
    sys.stdout.write('\n%s->%s Transitions\n' % (region1, region2))
    trans_times = t[ti_array[:,1]] - t[ti_array[:,0]]
    sys.stdout.write('Number of transition events: %d\n' % len(trans_times))
    sys.stdout.write('Minimum event duration:      %g\n' % min(trans_times))
    sys.stdout.write('Maximum event duration:      %g\n' % max(trans_times))
    tt_mean = trans_times.mean()
    se_tt_mean = trans_times.std() / (len(trans_times) ** 0.5)
    sys.stdout.write('Mean event duration:         %g\n' % tt_mean)
    sys.stdout.write('Std. err. of mean duration:  %g\n' % se_tt_mean)
    asaver.save_array(trans_times, '%(region1)s%(region2)s_trans%(format_ext)s',
                      filename_params)
    quantiles = [0.01*x for x in xrange(1,100)]
    tt_qs = numpy.empty((len(quantiles), 2), numpy.float_)
    tt_qs[:,0] = quantiles
    tt_qs[:,1] = wemd.analysis.probability.quantiles(trans_times, quantiles)
    asaver.save_array(tt_qs,
                      '%(region1)s%(region2)s_trans_quantiles%(format_ext)s',
                      filename_params)
    
    try:
        pi_array = te_finder.passage_indices[region1, region2]
    except KeyError:
        pass
    else:
        passage_times = t[pi_array[:,1]] - t[pi_array[:,0]]
        mfpt = passage_times.mean()
        se_mfpt = passage_times.std() / (len(passage_times) ** 0.5)
        k = 1.0/mfpt
        se_k = 1.0/mfpt**2 * se_mfpt
        sys.stdout.write('Minimum first passage time:  %g\n' % min(passage_times))
        sys.stdout.write('Maximum first passage time:  %g\n' % max(passage_times))
        sys.stdout.write('Mean first passage time:     %g\n' % mfpt)
        sys.stdout.write('StdErr of MFPT:              %g\n' % se_mfpt)
        sys.stdout.write('Reaction rate:               %g\n' % k)
        sys.stdout.write('StdErr of reation rate:      %g\n' % se_k)
        asaver.save_array(passage_times, 
                          '%(region1)s%(region2)s_passage%(format_ext)s',
                          filename_params)
        pt_qs = numpy.empty((len(quantiles), 2), numpy.float_)
        pt_qs[:,0] = quantiles
        pt_qs[:,1] = wemd.analysis.probability.quantiles(passage_times, quantiles)
        asaver.save_array(pt_qs, 
                          '%(region1)s%(region2)s_passage_quantiles%(format_ext)s',
                          filename_params)
        
