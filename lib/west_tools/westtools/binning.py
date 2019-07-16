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

from .core import WESTToolComponent
import sys,logging, pickle, math
from itertools import count
import numpy
import westpa
import westpa.binning
from westpa.binning import RectilinearBinMapper
from westpa.extloader import get_object
from pickle import PickleError

log = logging.getLogger(__name__)

from west.data_manager import weight_dtype 
EPS = numpy.finfo(weight_dtype).eps


def mapper_from_expr(expr):
    namespace = {'numpy': numpy,
                 'inf': float('inf')}
    try:
        mapper = RectilinearBinMapper(eval(expr,namespace))
    except TypeError as e:
        if 'has no len' in str(e):
            raise ValueError('invalid bin boundary specification (a list of lists is required)')
        else:
            raise
    else:
        log.debug('loaded {!r} from expression {!r}'.format(mapper,expr))
        return mapper

def mapper_from_system():
    system = westpa.rc.get_system_driver()
    log.debug('loaded {!r} from {!r}'.format(system.bin_mapper,system))
    return system.bin_mapper

def mapper_from_function(funcspec):
    '''Return a mapper constructed by calling a function in a named module.
    ``funcspec`` should be formatted as ``[PATH]:MODULE.FUNC``. This function 
    loads MODULE, optionally adding PATH to the search path, then returns MODULE.FUNC()'''
    if ':' in funcspec:
        (pathpart, funcpart) = funcspec.rsplit(':')
        pathinfo = ['.'] + pathpart.split(':')
    else:
        funcpart = funcspec
        pathinfo = ['.']
        
    fn = get_object(funcpart,['.'] + pathinfo)
    mapper = fn()
    log.debug('loaded {!r} from {!r}'.format(mapper,fn))
    return mapper

def mapper_from_hdf5(topol_group, hashval):
    '''Retrieve the mapper identified by ``hashval`` from the given bin topology group
    ``topol_group``. Returns ``(mapper, pickle, hashval)``'''
    
    try:
        index_ds = topol_group['index']
        pickle_ds = topol_group['pickles']
    except KeyError:
        log.debug('binhash {} not found'.format(hashval), exc_info=True)
        raise KeyError('binhash {} not found'.format(hashval))
    
    n_entries = len(index_ds)
    if n_entries == 0:
        raise KeyError('hash {} not found'.format(hashval))
    
    chunksize=256
    for istart in range(0,n_entries,chunksize):
        chunk = index_ds[istart:min(istart+chunksize,n_entries)]
        for i in range(len(chunk)):
            if chunk[i]['hash'] == hashval:
                pkldat = bytes(pickle_ds[istart+i,0:chunk[i]['pickle_len']].data)
                mapper = pickle.loads(pkldat) 
                log.debug('loaded {!r} from {!r}'.format(mapper, topol_group))
                log.debug('hash value {!r}'.format(hashval))
                return mapper, pkldat, hashval
    
    raise KeyError('hash {} not found'.format(hashval))

def mapper_from_yaml(yamlfilename):
    import yaml
    ydict = yaml.load(open(yamlfilename, 'rt'))
    ybins = ydict['bins']
    from westpa._rc import bins_from_yaml_dict
    #return mapper_from_dict(ybins)
    return bins_from_yaml_dict(ybins)

# We want this function to live on...
def mapper_from_dict(ybins):
    typename = ybins.pop('type')
    kwargs = ybins
    
    try:
        mapper_type = getattr(sys.modules['westpa.binning'], typename)
    except AttributeError:
        raise KeyError('unknown bin mapper type {!r} in YAML file {!r}'.format(typename, yamlfilename))
    
    if typename == 'RectilinearBinMapper':
        boundary_lists = kwargs.pop('boundaries')
        for ilist, boundaries in enumerate(boundary_lists):
            boundary_lists[ilist] = list(map((lambda x: 
                                           float('inf') 
                                           if (x if isinstance(x, str) else '').lower() == 'inf' 
                                           else x), boundaries))
        return mapper_type(boundary_lists)
    else:
        try:
            return mapper_type(**kwargs)
        except Exception:
            log.exception('exception instantiating mapper')
            raise

def write_bin_info(mapper, assignments, weights, n_target_states, outfile=sys.stdout, detailed=False):
    '''Write information about binning to ``outfile``, given a mapper (``mapper``) and the weights
    (``weights``) and bin assignments (``assignments``) of a set of segments, along with a target state
    count (``n_target_states``). If ``detailed`` is true, then per-bin information is written as well as
    summary information about all bins.'''
     
    norm = weights.sum()
    enorm = norm - 1.0
    enormeps = abs(enorm / EPS)
    bincounts = numpy.bincount(assignments, minlength=len(assignments))
    n_occupied = numpy.count_nonzero(bincounts)
    binweights = numpy.bincount(assignments,weights,minlength=len(assignments))
    nonzero_counts = (bincounts > 0)
    n_active = mapper.nbins - n_target_states
    weights_by_bin = [[] for _i in range(mapper.nbins)]
            
    min_bin_weight = binweights[nonzero_counts].min()
    max_bin_weight = binweights.max()
    min_seg_weight = weights.min()
    max_seg_weight = weights.max()
    
    ndec = int(math.ceil(-math.log10(1/n_active)))
    
    outfile.write('{:d} segments\n'.format(len(weights)))
    outfile.write('{:d} bins total, {:d} targets, {:d} ({:.{ndec}%}) occupied\n'
            .format(mapper.nbins, n_target_states, n_occupied, n_occupied/n_active,ndec=ndec))
    outfile.write('Minimum probability by bin:     {:23.17e}\n'.format(min_bin_weight))
    outfile.write('Maximum probability by bin:     {:23.17e}\n'.format(max_bin_weight))
    outfile.write('Dynamic range (by bin):         {:g} kT\n'.format(-math.log(min_bin_weight/max_bin_weight)))
    outfile.write('Minimum probability by segment: {:23.17e}\n'.format(min_seg_weight))
    outfile.write('Maximum probability by segment: {:23.17e}\n'.format(max_seg_weight))
    outfile.write('Dynamic range (by segment):     {:g} kT\n'.format(-math.log(min_seg_weight/max_seg_weight)))
    outfile.write('Norm = {:g}, error in norm = {:g} ({:.0f} epsilon)\n'
            .format(norm, enorm, enormeps))
    
    if not detailed:
        return
    
    outfile.write('\n')
    for iseg, weight in enumerate(weights):
        weights_by_bin[assignments[iseg]].append(weight)

    mw = max(6,max(len(str(mapper.nbins-1)), len(str(bincounts.max()))))
    fmt = '    '.join(['{ibin:{mw}d}','{nwalkers:{mw}d}','{weight:23.17e}',
                       '{min_weight:23.17e}','{max_weight:23.17e}','{weight_ratio:12.5f}','{label}\n'])
    outfile.write('{:>{mw}s}    {:>{mw}s}    {:23s}    {:23s}    {:23s}    {:12s}    {}\n'
            .format('Index', 'Count', 'Total weight', 'Min seg weight', 'Max seg weight', 'Weight ratio', 'Label',
                    mw=mw))
                
    for ibin, label, nwalkers, total_weight, seg_weights \
    in zip(count(), mapper.labels, bincounts, binweights, weights_by_bin):
        if nwalkers > 0:
            min_seg_weight = min(seg_weights)
            max_seg_weight = max(seg_weights)
            weight_ratio = max_seg_weight/min_seg_weight
        else:
            min_seg_weight = 0
            max_seg_weight = 0
            weight_ratio = 0
        outfile.write(fmt.format(ibin=ibin, label=label, nwalkers=nwalkers, weight=total_weight, min_weight=min_seg_weight,
                           max_weight=max_seg_weight,weight_ratio=weight_ratio,mw=mw))
    
    
    
def write_bin_labels(mapper, dest, 
                     header='# bin labels:\n',
                     fmt='# bin {index:{max_iwidth}d} -- {label!s}\n'):
    
    '''Print labels for all bins in ``mapper`` to the file-like object``dest``.
    
    If provided, ``header`` is printed prior to any labels. A number of expansions
    are available in ``header``:
      * ``mapper`` -- the mapper itself (from which most of the following can be obtained)
      * ``classname``  -- the class name of the mapper
      * ``nbins`` -- number of bins in the mapper 
      
    The ``fmt`` string specifies how bin labels are to be printed. A number of
    expansions are available in ``fmt``:
      * ``index`` -- the zero-based index of the bin
      * ``label`` -- the label of the bin
      * ``max_iwidth`` -- the maximum width (in characters) of the bin index, for pretty alignment
    '''
    
    if header:
        dest.write(header.format(mapper=mapper,classname=mapper.__class__.__name__, nbins=mapper.bins))
         
    max_iwidth = len(str(mapper.nbins-1))
    for (ibin,label) in enumerate(mapper.labels):
        dest.write(fmt.format(index=ibin, label=label, max_iwidth=max_iwidth))


class BinMappingComponent(WESTToolComponent):
    '''Component for obtaining a bin mapper from one of several places based on
    command-line arguments. Such locations include an HDF5 file that contains
    pickled mappers (including the primary WEST HDF5 file), the system object,
    an external function, or (in the common case of rectilinear bins) a
    list of lists of bin boundaries.
    
    Some configuration is necessary prior to calling process_args() if loading a
    mapper from HDF5. Specifically, either set_we_h5file_info() or 
    set_other_h5file_info() must be called to describe where to find the
    appropriate mapper. In the case of set_we_h5file_info(), the mapper used for
    WE at the end of a given iteration will be loaded. In the case of
    set_other_h5file_info(), an arbitrary group and hash value are specified;
    the mapper corresponding to that hash in the given group will be returned.
    
    In the absence of arguments, the mapper contained in an existing HDF5 file
    is preferred; if that is not available, the mapper from the system driver
    is used.
        
    This component adds the following arguments to argument parsers:
    
      --bins-from-system
          Obtain bins from the system driver
        
      --bins-from-expr=EXPR
          Construct rectilinear bins by parsing EXPR and calling
          RectilinearBinMapper() with the result. EXPR must therefore be a list of
          lists.
    
      --bins-from-function=[PATH:]MODULE.FUNC
          Call an external function FUNC in module MODULE (optionally adding PATH
          to the search path when loading MODULE) which, when called, returns a
          fully-constructed bin mapper.
      
      --bins-from-file
          Load bin definitions from a YAML configuration file.
          
      --bins-from-h5file
          Load bins from the file being considered; this is intended to mean the
          master WEST HDF5 file or results of other binning calculations, as
          appropriate.
    '''
    
    def __init__(self):
        
        # The final mapper 
        self.mapper = None
        self.mapper_hash = None
        self.mapper_pickle = None
        self.mapper_source = None
        self.mapper_source_desc = None
        self.bin_target_counts = None
        self.target_counts_required = False
        
        # Cues for where to load a mapper from, for the external file case
        self.mapper_source_group = None
        self.mapper_source_n_iter = None
        self.mapper_source_hash = None
        
        self._parse_target_count_args = False
        
        
    def add_args(self, parser, description='binning options', suppress=[]):
        suppressed_options = set(suppress)
        
        group = parser.add_argument_group(description)
        egroup = group.add_mutually_exclusive_group()
        if '--bins-from-system' not in suppressed_options:
            egroup.add_argument('--bins-from-system', action='store_true',
                                help='''Bins are constructed by the system driver specified in the WEST configuration file
                                (default where stored bin definitions not available).''')
        if '--bins-from-expr' not in suppressed_options:
            egroup.add_argument('--bins-from-expr', '--binbounds', dest='bins_from_expr',
                                help='''Construct bins on a rectilinear grid according to the given BINEXPR. This must
                                be a list of lists of bin boundaries (one list of bin boundaries for each dimension
                                of the progress coordinate), formatted as a Python expression. E.g. "[[0,1,2,4,inf],[-inf,0,inf]]".
                                The numpy module and the special symbol "inf" (for floating-point infinity) are available
                                for use within BINEXPR.''')
        if '--bins-from-function' not in suppressed_options:
            egroup.add_argument('--bins-from-function', '--binfunc', dest='bins_from_function',
                                help='''Supply an external function which, when called, returns a properly constructed
                                bin mapper which will then be used for bin assignments. This should be formatted as
                                "[PATH:]MODULE.FUNC", where the function FUNC in module MODULE will be used; the optional
                                PATH will be prepended to the module search path when loading MODULE.''')
        if '--bins-from-file' not in suppressed_options:
            egroup.add_argument('--bins-from-file', '--binfile', dest='bins_from_file', metavar='BINFILE',
                                help='''Load bin specification from the YAML file BINFILE. This currently
                                takes the form {'bins': {'type': 'RectilinearBinMapper', 'boundaries':
                                [[boundset1], [boundset2], ... ]}}; only rectilinear bin bounds are supported.''')
            
        if '--bins-from-h5file' not in suppressed_options:
            egroup.add_argument('--bins-from-h5file', action='store_true',
                                help='''Load bin specification from the data file being examined
                                (default where stored bin definitions available).''')
            
    def add_target_count_args(self, parser, description='bin target count options'):
        '''Add options to the given parser corresponding to target counts.'''

        self._parse_target_count_args = True        
        group = parser.add_argument_group(description)
        egroup = group.add_mutually_exclusive_group()
        
        egroup.add_argument('--target-counts',
                            help='''Use TARGET_COUNTS instead of stored or system driver target counts.
                            TARGET_COUNTS is a comma-separated list of integers. As a special case, a single
                            integer is acceptable, in which case the same target count is used for all bins.''')
        egroup.add_argument('--target-counts-from', metavar='FILENAME',
                            help='''Read target counts from the text file FILENAME instead of using stored or system
                            driver target counts. FILENAME must contain a list of integers, separated by arbitrary
                            whitespace (including newlines).''')
        
    def process_args(self, args):
        
        # User may have suppressed any of these arguments
        bins_from_system = getattr(args, 'bins_from_system', None)
        bins_from_expr = getattr(args,'bins_from_expr', None)
        bins_from_function = getattr(args, 'bins_from_function', None)
        bins_from_file = getattr(args, 'bins_from_file', None)
        bins_from_h5file = getattr(args, 'bins_from_h5file', None) 
        
        if not any([bins_from_system, bins_from_expr, bins_from_function, bins_from_h5file]):
            log.debug('no argument provided')
            if self.mapper_source_group is None:
                log.info('using bins from system file')
                bins_from_system = True
            else:
                log.info('using bins from HDF5 file')
                bins_from_h5file = True
        
        if bins_from_h5file:
            # it is an error to call process_args() prior to setting source group and hash,
            # preferably using the set_*_h5file_info() functions.
            assert self.mapper_source_group is not None
            assert self.mapper_source_hash is not None
            self.mapper, self.mapper_pickle, self.mapper_hash = mapper_from_hdf5(self.mapper_source_group, self.mapper_source_hash)
            self.mapper_source = 'file'
            self.mapper_source_desc = 'HDF5 file {} (group "{}")'.format(self.mapper_source_group.file.filename,
                                                                         self.mapper_source_group.name)
        elif bins_from_file:
            self.mapper = mapper_from_yaml(args.bins_from_file)
            self.mapper_source = args.bins_from_file
            self.mapper_source_desc = 'YAML file {!r}'.format(args.bins_from_file)
        elif bins_from_system:
            self.mapper = mapper_from_system()
            self.mapper_source = 'system'
            self.mapper_source_desc = 'system driver'
        elif bins_from_expr:
            self.mapper = mapper_from_expr(args.bins_from_expr)
            self.mapper_source = 'expr'
            self.mapper_source_desc = 'boundary expression'
        elif bins_from_function:
            self.mapper = mapper_from_function(args.bins_from_function)
            self.mapper_source = 'func'
            self.mapper_source_desc = 'external function'

        if self.mapper and not self.mapper_hash:
            try:
                self.mapper_pickle, self.mapper_hash = self.mapper.pickle_and_hash()
            except PickleError as e:
                log.debug('could not pickle bin mapper: {}'.format(e))
                self.mapper_pickle = self.mapper_hash = None
                
        log.info('loaded mapper {!r} from {}'.format(self.mapper, self.mapper_source_desc))
        
        if self._parse_target_count_args and self.target_counts_required:
            import re
            if args.target_counts is not None:
                self.bin_target_counts = numpy.array(list(map(int,re.split(r'\s*,\s*', args.target_counts))))
            elif args.target_counts_from is not None:
                self.bin_target_counts = numpy.array(list(map(int,re.split(r'\s*,\s*', open(args.target_counts_from,'rt').read()))))
            else:
                # if target counts required, they will have already been loaded from a master
                if self.bin_target_counts is None and self.target_counts_required:
                    raise EnvironmentError('target counts are required but none have been provided')
                
            if len(self.bin_target_counts) == 1:
                flat_target_count = self.bin_target_counts
                self.bin_target_counts = numpy.empty((self.mapper.nbins,), numpy.int)
                self.bin_target_counts[:] = flat_target_count
            log.debug('bin target counts = {!r}'.format(self.bin_target_counts))
                
    
    def set_we_h5file_info(self, n_iter=None, data_manager=None, required=False):
        '''Set up to load a bin mapper from the master WEST HDF5 file. The mapper is actually loaded
        from the file when self.load_bin_mapper() is called, if and only if command line arguments
        direct this. If ``required`` is true, then a mapper must be available at iteration ``n_iter``,
        or else an exception will be raised.'''
                
        data_manager = data_manager or westpa.rc.get_data_manager()
        n_iter = n_iter or data_manager.current_iteration
        iter_group = data_manager.get_iter_group(n_iter)
        
        try:
            self.mapper_source_group = data_manager.we_h5file['/bin_topologies']
            self.mapper_source_group['index']   # raises KeyError if missing
            self.mapper_source_group['pickles'] # raises KeyError if missing
            self.mapper_source_hash = iter_group.attrs['binhash']
        except TypeError:
            if required:
                raise EnvironmentError('HDF5 file not open')
        except KeyError:
            if required:
                raise EnvironmentError('bin topologies not available')
        
        try:
            self.bin_target_counts = iter_group['bin_target_counts'][...]
        except KeyError:
            pass
        
                    
    def set_other_h5file_info(self, topology_group, hashval):
        '''Set up to load a bin mapper from (any) open HDF5 file, where bin topologies are
        stored in ``topology_group`` (an h5py Group object) and the desired mapper has hash
        value ``hashval``. The mapper itself is loaded when self.load_bin_mapper() is called.'''
        
        self.mapper_source_group = topology_group
        self.mapper_source_hash = hashval
     
        
