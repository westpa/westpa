import numpy as np 

def get_parents(walker_tuple, h5_file):
    it, wlk = walker_tuple
    parent = h5_file["iterations/iter_{:08d}".format(it)]["seg_index"]["parent_id"][wlk]
    return it-1, parent

def trace_walker(walker_tuple, h5_file):
    # Unroll the tuple into iteration/walker 
    it, wlk = walker_tuple
    # Initialize our path
    path = [(it,wlk)]
    # And trace it
    while it > 1: 
        it, wlk = get_parents((it, wlk), h5_file)
        path.append((it,wlk))
    return np.array(sorted(path, key=lambda x: x[0]))

def get_pcoords(path, h5_file):
    # Initialize a list for the pcoords
    pcoords = []
    # Loop over the path and get the pcoords for each walker
    for it, wlk in path:
        # Here we are taking every 10 time points, feel free to adjust to see what that does
        pcoords.append(h5_file['iterations/iter_{:08d}'.format(it)]['pcoord'][wlk][::10,:])
    return np.array(pcoords)

def adjust_plot(hist, midpoints, binbounds):
    import matplotlib.pyplot as pyplot
    import h5py
    # First adjust axis labels
    pyplot.xlabel("First dimension")
    pyplot.ylabel("Second dimension")

    # Now we can play with some fun stuff
    # First we'll pull which iteration we are on using the 
    # file created by our movie loop
    with open("cur_iter.txt", 'r') as f:
        cur_iter = int(f.readline())

    # Now we want to pull our tranjectory and  plot it 
    # until the current iteration. If you have long 
    # trajectories or a lot of data it's a good idea
    # to calculate this beforehand and just open a pickled 
    # file or something instead of re-calculating it every 
    # iteration like it's done here.
    with h5py.File("west.h5", "r") as w:
        # We will trace a specific trajectory from the latest 
        # iteration to the beginning and then use only the part 
        # that we currently need
        last_iter = w.attrs['west_current_iteration'] - 1
        # Say we pull the 0th walker and trace that 
        path = trace_walker((last_iter, 279), w)
        # And pull pcoords for the path calculated
        pcoords = get_pcoords(path, w)
        cur_pcoords = pcoords[:cur_iter,:,:]
        shp = cur_pcoords.shape
        cur_pcoords = cur_pcoords.reshape(shp[0]*shp[1], shp[2])
        # We got the full path, now plot up to the part we need
        pyplot.plot(cur_pcoords[:,0], cur_pcoords[:,1], c="black", lw=3)
        pyplot.plot(cur_pcoords[:,0], cur_pcoords[:,1], c="white", lw=2)
        # Say we pull the 300th walker and trace that 
        path = trace_walker((last_iter, 435), w)
        # And pull pcoords for the path calculated
        pcoords = get_pcoords(path, w)
        cur_pcoords = pcoords[:cur_iter,:,:]
        shp = cur_pcoords.shape
        cur_pcoords = cur_pcoords.reshape(shp[0]*shp[1], shp[2])
        print("Current pcoords shape")
        print(cur_pcoords.shape)
        # We got the full path, now plot up to the part we need
        #import IPython
        #IPython.embed()
        pyplot.plot(cur_pcoords[:,0], cur_pcoords[:,1], c="black", lw=3)
        pyplot.plot(cur_pcoords[:,0], cur_pcoords[:,1], c="cyan", lw=2)
