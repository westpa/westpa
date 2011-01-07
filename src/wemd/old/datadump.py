import time
import wemd
from wemd import Segment
import numpy
import cPickle

DO_PRINT = False

runtime_config = wemd.rc.read_config()
data_manager = wemd.data_manager.make_data_manager(runtime_config)

min_iter = 100
max_iter = 250

itime = time.time()
# This returns a list of named-field tuple containing rows from the database; accessing
# fields from the rows is shown below
segments = data_manager.get_segments(Segment.n_iter.between(min_iter,max_iter),
                                     result_format='rows')

# This returns a list of (parent_id, child_id) pairs; merges are not included
conn_pairs = data_manager.get_connectivity(Segment.n_iter.between(min_iter,max_iter))
ftime = time.time()

print "data read in %g seconds" % (ftime-itime)

# This indexes the above result by parent_id for pretty printing below
conn_dict = dict()
for (parent_id, seg_id) in conn_pairs:
    try:
        conn_dict[parent_id].append(seg_id)
    except KeyError:
        conn_dict[parent_id] = [seg_id]

# This demonstrates how to access fields from the rows
if DO_PRINT:
    for segment in segments:
        print "segment %d" % segment.seg_id
        print "  p_parent_id: %s" % segment.p_parent_id
        print "  n_iter:      %d" % segment.n_iter
        print "  weight:      %g" % segment.weight
        print "  pcoord:      %r" % segment.pcoord
        print "  children:    %r" % conn_dict.get(segment.seg_id, [])

# This demonstrates how to reduce things to numpy arrays if you need to
all_segids = numpy.fromiter((seg.seg_id for seg in segments), numpy.uint64) 
all_weights = numpy.fromiter((seg.weight for seg in segments), numpy.float64)

# This demonstrates how to save selected fields from the result list
# We use simple tuples so that SQLAlchemy isn't necessary to unpickle
pickle_data = [(seg.seg_id, seg.n_iter, seg.weight, seg.pcoord) for seg in segments]
cPickle.dump(pickle_data, open('datadump.pickle', 'w'))



