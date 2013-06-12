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
import sys
import re
from westtools import WESTTool, WESTDataReader
import numpy, h5py, operator, time
import westpa
from westpa import h5io

from west import Segment
from west.states import InitialState
from west.data_manager import (weight_dtype, n_iter_dtype, seg_id_dtype, utime_dtype, vstr_dtype, 
                               istate_type_dtype, istate_status_dtype)

class Trace:
    '''A class representing a trace of a certain trajectory segment back to its origin.'''
    
    def __init__(self, summary, endpoint_type, basis_state, initial_state, data_manager = None):
        self.summary = summary
        self.endpoint_type = endpoint_type
        self.basis_state = basis_state
        self.initial_state = initial_state
        self.data_manager = data_manager or westpa.rc.get_data_manager()
        
        # A mapping from aux file names to open h5py.File objects, to minimize time
        
        self._auxfiles = {}
            
    def __len__(self):
        try:
            return len(self.summary)
        except TypeError:
            return 0
        
    def __getitem__(self, sl):
        return self.summary[sl]
    
    def __iter__(self):
        return iter(self.summary)
            
    @classmethod
    def from_data_manager(cls, n_iter, seg_id, data_manager = None):
        '''Construct and return a trajectory trace whose last segment is identified
        by ``seg_id`` in the iteration number ``n_iter``.'''
        
        data_manager = data_manager or westpa.rc.get_data_manager()
        
        # These values are used later on
        endpoint_type = None
        pcoord_dtype = None
        pcoord_pt_shape = None
        
        seginfo = []
        parent_id = seg_id
        
        while n_iter > 0 and parent_id >= 0:
            seg_id = parent_id
            iter_group = data_manager.get_iter_group(n_iter)
            pcoord_ds = iter_group['pcoord']
            seg_index = iter_group['seg_index']
            n_segs = pcoord_ds.shape[0]
            pcoord_len = pcoord_ds.shape[1]
            
            assert seg_id < n_segs
            
            indexrow = seg_index[seg_id]
            final_pcoord = pcoord_ds[seg_id, pcoord_len-1]
            weight = indexrow['weight']
            cputime = indexrow['cputime']
            walltime = indexrow['walltime']
            
            try:
                parent_id = long(indexrow['parent_id'])
            except IndexError:
                # old HDF5 version
                parent_id = long(iter_group['parents'][indexrow['parents_offset']])
                
            if endpoint_type is None:
                endpoint_type = indexrow['endpoint_type']
                pcoord_pt_shape = pcoord_ds.shape[2:]
                pcoord_dtype = pcoord_ds.dtype
                
            seginfo.append((n_iter, seg_id, weight, walltime, cputime, final_pcoord))
            
            del iter_group, pcoord_ds, seg_index
            n_iter -= 1
            
        # loop terminates with parent_id set to the identifier of the initial state, 
        # seg_id set to the identifier of the first segment in the trajectory, and
        # n_iter set to one less than the iteration of the first segment
        first_iter = n_iter + 1
        first_seg_id = seg_id
        first_parent_id = parent_id

        # Initial segment (for fetching initial state)
        first_segment = Segment(n_iter=first_iter, seg_id=first_seg_id, parent_id=first_parent_id)
        
        seginfo.reverse()
        
        summary_dtype = numpy.dtype([('n_iter', n_iter_dtype),
                                     ('seg_id', seg_id_dtype),
                                     ('weight', weight_dtype),
                                     ('walltime', utime_dtype),
                                     ('cputime', utime_dtype),
                                     ('final_pcoord', pcoord_dtype, pcoord_pt_shape),
                                     ])
        
        summary = numpy.array(seginfo, dtype=summary_dtype)
        
        try:
            initial_state = data_manager.get_segment_initial_states([first_segment], first_iter)[0]
        except KeyError:
            # old HDF5 version
            assert parent_id < 0
            istate_pcoord = data_manager.get_iter_group(first_iter)['pcoord'][first_seg_id,0]
            istate_id = -(first_parent_id+1)
            basis_state = None
            initial_state = InitialState(istate_id, None, iter_created=0, pcoord=istate_pcoord)
            
        else:
            basis_state = data_manager.get_basis_states(first_iter)[initial_state.basis_state_id]
            
        return cls(summary, endpoint_type, basis_state, initial_state, data_manager)
            
    def get_segment_data_slice(self, datafile, dsname, n_iter, seg_id, slice_=None, index_data=None,
                               iter_prec=None):
        '''Return the data from the dataset named ``dsname`` within the given ``datafile`` (an open
        h5py.File object) for the given iteration and segment. By default, it is assumed that the
        dataset is stored in the iteration group for iteration ``n_iter``, but if ``index_data``
        is provided, it must be an iterable (preferably a simple array) of (n_iter,seg_id) pairs,
        and the index in the ``index_data`` iterable of the matching n_iter/seg_id pair is used as
        the index of the data to retrieve.
        
        If an optional ``slice_`` is provided, then the given slicing tuple is appended to that
        used to retrieve the segment-specific data (i.e. it can be used to pluck a subset of the
        data that would otherwise be returned).
        '''

        if slice_ is None:
            slice_ = numpy.s_[...] 
                    
        if index_data is not None:
            dataset = datafile[dsname]

            for i, (i_n_iter,i_seg_id) in enumerate(index_data):
                if (i_n_iter,i_seg_id) == (n_iter,seg_id):
                    break
            else:
                raise KeyError((n_iter,seg_id))
            
            itpl = (i,) + slice_
            return dataset[itpl]
        else:
            if not iter_prec:
                iter_prec = datafile.attrs.get('west_iter_prec', self.data_manager.default_iter_prec)
            igname_tail = 'iter_{:0{iter_prec:d}d}'.format(int(n_iter),iter_prec=int(iter_prec))
            try:
                iter_group = datafile['/iterations/' + igname_tail]
            except KeyError:
                iter_group = datafile[igname_tail]
            
            dataset = iter_group[dsname]
            itpl = (seg_id,) + slice_

            return dataset[itpl]

    def trace_timepoint_dataset(self, dsname, slice_=None, auxfile=None,index_ds=None):
        '''Return a trace along this trajectory over a dataset which is layed out as [seg_id][timepoint][...].
        Overlapping values at segment boundaries are accounted for.  Returns (data_trace, weight), where 
        data_trace is a time series of the dataset along this trajectory, and weight is the corresponding
        trajectory weight at each time point.
        
        If ``auxfile`` is given, then load the dataset from the given HDF5 file, which must be 
        layed out the same way as the main HDF5 file (e.g. iterations arranged as
        iterations/iter_*).
        
        If index_ds is given, instead of reading data per-iteration from iter_* groups, then the
        given index_ds is used as an index of n_iter,seg_id pairs into ``dsname``. In this case,
        the target data set need not exist on a per-iteration basis inside iter_* groups.
        
        If ``slice_`` is given, then *further* slice the data returned from the HDF5 dataset. This can
        minimize I/O if it is known (and specified) that only a subset of the data along the
        trajectory is needed.
        '''
        
        # Figure out where to look for the dataset
        if isinstance(auxfile, basestring):
            datafile = h5py.File(auxfile, 'r')
            close_datafile = True
        elif auxfile is not None:
            datafile = auxfile
            close_datafile = False
        else:
            datafile = self.data_manager.we_h5file
            close_datafile = False
            
        iter_prec = self.data_manager.iter_prec
        get_data_slice = self.get_segment_data_slice
            
        # Load the index if we use it
        if index_ds is not None:
            if isinstance(index_ds,basestring):
                index_ds = datafile[index_ds]
            index_data = index_ds[...]
        else:
            index_data = None
            
        # Be sure to retrieve the time series
        if not slice_:
            first_sl = numpy.index_exp[:, ...]
            other_sl = numpy.index_exp[1:,...]
        else:
            first_sl = numpy.index_exp[:] + slice_
            other_sl = numpy.index_exp[1:] + slice_ 
        
        # Retrieve the first segment's data
        first_n_iter, first_seg_id = self.summary[0]['n_iter'], self.summary[0]['seg_id']
        first_iter_data = get_data_slice(datafile, dsname, first_n_iter, first_seg_id, first_sl, index_data, iter_prec)

        n_segs = len(self)
        n_points_per_seg = len(first_iter_data)
        
        length = n_points_per_seg + (n_segs-1)*(n_points_per_seg-1)
        tracedata = numpy.empty((length,) + first_iter_data.shape[1:], dtype=first_iter_data.dtype)
        traceweight = numpy.empty((length,), weight_dtype)
        
        # Store first segment data
        tracedata[0:n_points_per_seg] = first_iter_data
        traceweight[0:n_points_per_seg] = self.summary[0]['weight']
        del first_iter_data
        
        # Store remainder of data
        
        for iseg, summary_item in enumerate(self.summary[1:]):
            n_iter = summary_item['n_iter']
            seg_id = summary_item['seg_id']
            weight = summary_item['weight']
            
            offset = n_points_per_seg + iseg*(n_points_per_seg-1)
            length = n_points_per_seg - 1
            seg_data = get_data_slice(datafile, dsname, n_iter, seg_id, other_sl, index_data, iter_prec)
            
            tracedata[offset:offset+length] = seg_data 
            traceweight[offset:offset+length] = weight
            del seg_data
        
        if close_datafile:
            datafile.close()
        
        return tracedata, traceweight

    """
    # This is disabled until there is a real use for it; the following code is 
    # outdated
    def trace_perseg_dataset(self, dsname):
        '''Return a trace along this trajectory over a dataset which is layed out as [seg_id][...].
        Returns (data_trace, weight), where  data_trace is a time series of the dataset along this
        trajectory, and weight is the corresponding trajectory weight at each time point.'''

        first_n_iter, first_seg_id = self.summary[0]['n_iter'], self.summary[0]['seg_id']
        first_iter_group = self.data_manager.get_iter_group(first_n_iter)
        first_iter_ds = first_iter_group[dsname]
        n_segs = len(self)
        tracedata = numpy.empty((n_segs,) + first_iter_ds.shape[1:], dtype=first_iter_ds.dtype)
        traceweight = numpy.empty((n_segs,), weight_dtype)
        tracedata[0] = first_iter_ds[first_seg_id]
        traceweight[0] = self.summary[0]['weight']
        for isegm1, summary_item in enumerate(self.summary[1:]):
            iseg = isegm1 + 1
            n_iter = summary_item['n_iter']
            seg_id = summary_item['seg_id']            
            iter_group = self.data_manager.get_iter_group(n_iter)
            seg_data = iter_group[dsname][seg_id]
            tracedata[iseg] = seg_data
            traceweight[iseg] = summary_item['weight']
            del seg_data
            
        return tracedata, traceweight
    """
        
class WTraceTool(WESTTool):
    prog='w_trace'
    description = '''\
Trace individual WEST trajectories and emit (or calculate) quantities along the
trajectory.

Trajectories are specified as N_ITER:SEG_ID pairs. Each segment is traced back
to its initial point, and then various quantities (notably n_iter and seg_id)
are printed in order from initial point up until the given segment in the given
iteration.

Output is stored in several files, all named according to the pattern given by
the -o/--output-pattern parameter. The default output pattern is "traj_%d_%d",
where the printf-style format codes are replaced by the iteration number and
segment ID of the terminal segment of the trajectory being traced.

Individual datasets can be selected for writing using the -d/--dataset option
(which may be specified more than once). The simplest form is ``-d dsname``,
which causes data from dataset ``dsname`` along the trace to be stored to
HDF5.  The dataset is assumed to be stored on a per-iteration basis, with
the first dimension corresponding to seg_id and the second dimension
corresponding to time within the segment.  Further options are specified
as comma-separated key=value pairs after the data set name, as in

    -d dsname,alias=newname,index=idsname,file=otherfile.h5,slice=[100,...]
    
The following options for datasets are supported:

    alias=newname
        When writing this data to HDF5 or text files, use ``newname``
        instead of ``dsname`` to identify the dataset. This is mostly of
        use in conjunction with the ``slice`` option in order, e.g., to
        retrieve two different slices of a dataset and store then with
        different names for future use.

    index=idsname
        The dataset is not stored on a per-iteration basis for all
        segments, but instead is stored as a single dataset whose
        first dimension indexes n_iter/seg_id pairs. The index to
        these n_iter/seg_id pairs is ``idsname``.
    
    file=otherfile.h5
        Instead of reading data from the main WEST HDF5 file (usually
        ``west.h5``), read data from ``otherfile.h5``.
        
    slice=[100,...]
        Retrieve only the given slice from the dataset. This can be
        used to pick a subset of interest to minimize I/O.
        
-------------------------------------------------------------------------------
'''


    pcoord_formats = {'u8': '%20d',
                      'i8': '%20d',
                      'u4': '%10d',
                      'i4': '%11d',
                      'u2': '%5d',
                      'i2': '%6d',
                      'f4': '%14.7g',
                      'f8': '%023.15g'}
    
    def __init__(self):
        super(WTraceTool,self).__init__()
        
        self.data_reader = WESTDataReader()
        #self.h5storage = HDF5Storage()
        self.output_file = None
        self.output_pattern = None
        self.endpoints = None
        self.datasets = []

        
    # Interface for command-line tools
    def add_args(self, parser):
        self.data_reader.add_args(parser)
        #self.h5storage.add_args(parser)
        parser.add_argument('-d', '--dataset', dest='datasets',
                            #this breaks argparse (see http://bugs.python.org/issue11874) 
                            #metavar='DSNAME[,alias=ALIAS][,index=INDEX][,file=FILE][,slice=SLICE]',
                            metavar='DSNAME',
                            action='append',
                            help='''Include the dataset named DSNAME in trace output. An extended form like
                            DSNAME[,alias=ALIAS][,index=INDEX][,file=FILE][,slice=SLICE] will
                            obtain the dataset from the given FILE instead of the main WEST HDF5 file,
                            slice it by SLICE, call it ALIAS in output, and/or access per-segment data by a n_iter,seg_id
                            INDEX instead of a seg_id indexed dataset in the group for n_iter.''')
        parser.add_argument('endpoints',  metavar='N_ITER:SEG_ID', nargs='+',
                            help='''Trace trajectory ending (or at least alive at) N_ITER:SEG_ID.''')
        
        #tgroup = parser.add_argument_group('trace options')
        ogroup = parser.add_argument_group('output options')
        ogroup.add_argument('--output-pattern', default='traj_%d_%d',
                            help='''Write per-trajectory data to output files/HDF5 groups whose names begin with OUTPUT_PATTERN,
                                 which must contain two printf-style format flags which will be replaced with the iteration number
                                 and segment ID of the terminal segment of the trajectory being traced.
                                 (Default: %(default)s.)''')
        ogroup.add_argument('-o', '--output', default='trajs.h5',
                            help='Store intermediate data and analysis results to OUTPUT (default: %(default)s).')
        
    
    def process_args(self, args):
        self.data_reader.process_args(args)
        #self.h5storage.process_args(args)
        self.endpoints = [map(long,endpoint.split(':')) for endpoint in args.endpoints]
        self.output_pattern = args.output_pattern
        
        for dsstr in args.datasets or []:
            self.datasets.append(self.parse_dataset_string(dsstr))        
        
        #self.h5storage.open_analysis_h5file()
        self.output_file = h5py.File(args.output)
        
    def parse_dataset_string(self, dsstr):
        dsinfo = {}

        r = re.compile(r',(?=[^\]]*(?:\[|$))')
        fields = r.split(dsstr)

        dsinfo['dsname'] = fields[0]

        for field in (field.strip() for field in fields[1:]):
            k,v = field.split('=')
            k = k.lower()
            if k in ('alias', 'file', 'index'):
                dsinfo[k] = v
            elif k == 'slice':
                try:
                    dsinfo['slice'] = eval('numpy.index_exp' + v)
                except SyntaxError:
                    raise SyntaxError('invalid index expression {!r}'.format(v))
            else:
                raise ValueError('invalid dataset option {!r}'.format(k))
            
        return dsinfo
    
    def go(self):
        self.data_reader.open('r')
        
        #Create a new 'trajectories' group if this is the first trace
        try:
            trajs_group = h5io.create_hdf5_group(self.output_file, 'trajectories', replace=False, creating_program=self.prog)
        except ValueError:
            trajs_group = self.output_file['trajectories']
        
        for n_iter, seg_id in self.endpoints:
            trajname = self.output_pattern % (n_iter,seg_id)
            trajgroup = trajs_group.create_group(trajname)

            trace = Trace.from_data_manager(n_iter,seg_id, self.data_reader.data_manager)
            
            with open(trajname + '_trace.txt', 'wt') as trace_output:
                self.emit_trace_text(trace, trace_output)
                
            self.emit_trace_h5(trace, trajgroup)
            
            aux_h5files = {}
            for dsinfo in self.datasets:
                dsname = dsinfo['dsname']
                filename = dsinfo.get('file')
                if filename:
                    try:
                        aux_h5file = aux_h5files[filename]
                    except KeyError:
                        aux_h5file = aux_h5files[filename] = h5py.File(filename, 'r')
                else:
                    aux_h5file = None
                    
                slice_ = dsinfo.get('slice')
                alias = dsinfo.get('alias', dsname)
                index = dsinfo.get('index')
                
                data, weights = trace.trace_timepoint_dataset(dsname, auxfile=aux_h5file, slice_=slice_,index_ds=index)
                
                # Save data to HDF5
                try:
                    del trajgroup[alias]
                except KeyError:
                    pass
                trajgroup[alias] = data
                
                # All weight vectors will be the same length, so only store in HDF5 once
                if not ('weights' in trajgroup and trajgroup['weights'].shape == weights.shape):
                    try:
                        del trajgroup['weights']
                    except KeyError:
                        pass    
                    trajgroup['weights'] = weights
                            
    def emit_trace_h5(self, trace, output_group):
        for dsname in ('basis_state', 'initial_state', 'segments'):
            try:
                del output_group[dsname]
            except KeyError:
                pass
        
        if trace.basis_state:
            output_group['basis_state'] = trace.basis_state.as_numpy_record()
        output_group['initial_state'] = trace.initial_state.as_numpy_record()
        output_group['segments'] = trace.summary
            
    def emit_trace_text(self, trace, output_file):
        '''Dump summary information about each segment in the given trace to the given output_file,
        which must be opened for writing in text mode.  Output columns are separated by at least
        one space.'''
        
        if not trace:
            return
                
        pcoord_ndim = trace[0]['final_pcoord'].shape[0]
        lastseg = trace[-1]
        len_n_iter = max(6, len(str(lastseg['n_iter'])))
        len_seg_id = max(6, max(len(str(seg_id)) for seg_id in trace['seg_id']))
        seg_pattern = '    '.join(['{n_iter:{len_n_iter}d}',
                                   '{seg_id:{len_seg_id}d}',
                                   '{weight:22.17e}',
                                   '{walltime:10.6g}',
                                   '{cputime:10.6g}',
                                   '{pcoord_str:s}'
                                   ]) + '\n'
                                  
        
        output_file.write('''\
# Trace of trajectory ending in n_iter:seg_id {n_iter:d}:{seg_id:d} (endpoint type {endpoint_type_text:s})   
# column  0: iteration (0 => initial state)
# column  1: seg_id (or initial state ID)
# column  2: weight
# column  3: wallclock time (s)
# column  4: CPU time (s)
'''.format(n_iter = long(lastseg['n_iter']), 
           seg_id = long(lastseg['seg_id']), 
           endpoint_type_text = Segment.endpoint_type_names[trace.endpoint_type]))
        
        
        if pcoord_ndim == 1:
            output_file.write('''\
# column  5: final progress coordinate value            
''')
        else:
            fpcbegin = 5
            fpcend = fpcbegin + pcoord_ndim - 1
            output_file.write('''\
# columns {fpcbegin:d} -- {fpcend:d}: final progress coordinate value
'''.format(fpcbegin=fpcbegin,fpcend=fpcend))
        
        
        pcoord_formats = self.pcoord_formats
        
        # Output row for initial state
        initial_state = trace.initial_state
        pcoord_str = '    '.join(pcoord_formats.get(pcfield.dtype.str[1:], '%s') % pcfield 
                                 for pcfield in initial_state.pcoord)
        output_file.write(seg_pattern.format(n_iter=0, seg_id=initial_state.state_id,
                                             weight=0.0, walltime=0, cputime=0, pcoord_str=pcoord_str,
                                             len_n_iter=len_n_iter,len_seg_id=len_seg_id))
        
        # Output rows for segments
        for segment in trace:            
            pcoord_str = '    '.join(pcoord_formats.get(pcfield.dtype.str[1:], '%s') % pcfield 
                                     for pcfield in segment['final_pcoord'])
            output_file.write(seg_pattern.format(n_iter = long(segment['n_iter']), 
                                                 seg_id = long(segment['seg_id']),
                                                 weight = float(segment['weight']),
                                                 walltime = float(segment['walltime']),
                                                 cputime = float(segment['cputime']), 
                                                 pcoord_str=pcoord_str,
                                                 len_n_iter=len_n_iter,
                                                 len_seg_id=len_seg_id))
        
    
        

if __name__ == '__main__':
    WTraceTool().main()
    
