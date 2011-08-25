'''
Tools for dealing with bins separate independently of a system file.
'''

from math import ceil, log10
import re

re_transform_inf = re.compile(r'inf', re.IGNORECASE)
inf = float('inf')

def region_set_from_code(code, namespace = None):
    '''Get a RegionSet from a file or string containing Python source. 
    The RegionSet may be constructed in any convenient manner, but there 
    must be a top-level variable called `region_set`` containing it by the
    end of the code's execution.
    
    If provided, ``namespace`` is the namespace in which the code is executed.
    The symbols ``wemd``, ``numpy``, ``RegionSet``, ``RectilinearRegionSet``,
    and ``PiecewiseRegionSet`` will be added to this namespace prior to
    execution. 
    '''
    import numpy
    import wemd
    from wemd.pcoords import RegionSet, PiecewiseRegionSet, RectilinearRegionSet

    namespace = dict(namespace) if namespace is not None else dict()
    namespace.update({'wemd': wemd,
                      'numpy': numpy,
                      'RegionSet': RegionSet,
                      'PiecewiseRegionSet': PiecewiseRegionSet,
                      'RectilinearRegionSet': RectilinearRegionSet,
                      })
    
    try:
        # This will raise SyntaxError as appropriate
        exec code in namespace
        region_set = namespace['region_set']
    except KeyError:
        raise KeyError('code must assign to a variable named "region_set"')
    else:
        return region_set

def rectilinear_region_set_from_expr(expr):
    '''Construct a RectilinearRegionSet from a string representing a list
    of lists of boundaries, one list for each dimension.  The syntax is that of
    Python.  The string 'inf' will be converted to ``float('inf')`` prior to 
    evaluation.  For convenience, this routine ensures that the list is
    surrounded by double brackets ('[[a, b, c, ...]]'), adding one or two pairs
    of brackets if necessary. 
    '''
    from wemd.pcoords import RectilinearRegionSet
    
    transformed_expr = '{}\n'.format(re_transform_inf.sub("float('inf')", expr))
    if not transformed_expr.startswith('['):
        transformed_expr = '[{}]'.format(transformed_expr)
    if not transformed_expr.startswith('[['):
        transformed_expr = '[{}]'.format(transformed_expr)
    binbounds = eval(transformed_expr)
    return RectilinearRegionSet(binbounds)

def add_region_set_options(parser):
    '''Add options relating to obtaining region sets to the given argparse
    parser.  The resulting arguments can be passed to get_region_set_from_args.'''
    
    bingroup = parser.add_argument_group('binning options')
    egroup = bingroup.add_mutually_exclusive_group()
    egroup.add_argument('--binfile', dest='binfile',
                        help='''Execute the Python script BINFILE to construct bins for analysis; must contain an assignment to
                        a variable called "region_set" (default: load from system file).''')
    egroup.add_argument('--binexpr', dest='binexpr',
                         help='''Use the Python code BINEXPR to construct bins for analysis; must contain an assignment to 
                         a variable called "region_set" (default: load from system file).''')
    egroup.add_argument('--binbounds', dest='binbounds',
                        help='''Construct rectilinear bins from BINBOUNDS, which will be parsed as a list of lists
                        of bin boundaries.  The text 'inf' will be replaced by infinity.
                        (Default: load from system file.)''')

def get_region_set_from_args(args, status_stream = None, print_labels=False, print_labels_dest=None):
    '''Get a region set from parsed command-line arguments [as added to an argparse parser by
    add_region_set_options()].  If status_stream is provided, then information about the bins
    is printed thereto.  If print_labels is True, then labels of all bins will be printed to
    the status stream (or print_labels_dest, if provided).'''
    
    print_labels_dest = print_labels_dest or status_stream
    
    args.using_external_bins = True
    if args.binfile:
        if status_stream:
            status_stream.write('loading bin boundaries from {}\n'.format(args.binfile))        
        region_set = region_set_from_code(open(args.binfile, 'rt'))
    elif args.binexpr:
        if status_stream:
            status_stream.write('constructing bin boundaries from the following code:\n{}\n'.format(args.binexpr))        
        region_set = region_set_from_code(args.binexpr)
    elif args.binbounds:
        if status_stream:
            status_stream.write('constructing bins with the following boundaries:\n  {}\n'.format(args.binbounds))        
        region_set = rectilinear_region_set_from_expr(args.binbounds)
    else:
        # Load system file
        if status_stream:
            status_stream.write('loading bin boundaries from WEMD system\n')
            
        import wemd
        
        runtime_config = wemd.rc.read_config(args.run_config_file)
        runtime_config.update_from_object(args)
        sim_manager = wemd.rc.load_sim_manager(runtime_config)
        sim_manager.load_system_driver()
        region_set = sim_manager.system.region_set
        
        args.using_external_bins = False
    
    # No point in printing if we're not connected to an output stream
    if print_labels and print_labels_dest:
        bins = region_set.get_all_bins()        
        maxwidth = int(ceil(log10(len(bins))))
        print_labels_dest.write('bin labels:\n')
        for (ibin, bin) in enumerate(bins):
            print_labels_dest.write('  bin {0:<{maxwidth}d}: {1!s}\n'.format(ibin, bin.label, maxwidth=maxwidth))
    
    return region_set
        
def print_labels(region_set, dest, header='# bin labels:\n', format='#  bin {bin_index:{max_iwidth}d} -- {label!s}\n'):
    '''Print labels for all bins in the given ``region_set`` to ``dest``.  If provided, ``header`` is printed
    before any labels.   The ``format`` string specifies how bin labels are to be printed.  Valid entries are:
      * ``bin_index`` -- the zero-based index of the bin
      * ``label`` -- the label, as obtained by ``bin.label``
      * ``max_iwidth`` -- the maximum width (in characters) of the bin index, for pretty alignment
    '''
    
    dest.write(header)
    bins = region_set.get_all_bins()
    max_iwidth = len(str(len(bins)-1))
    for (ibin, bin) in enumerate(bins):
        dest.write(format.format(bin_index=ibin, label=bin.label, max_iwidth=max_iwidth))
    
