# Copyright (C) 2017 Matthew C. Zwier and Lillian T. Chong
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

#import warnings
#warnings.filterwarnings('ignore', category=DeprecationWarning)
#warnings.filterwarnings('ignore', category=RuntimeWarning)
#warnings.filterwarnings('ignore', category=FutureWarning)
import numpy as np
import h5py
# We don't care if we're reimporting blessings, under this context.
#warnings.filterwarnings('ignore')

class Plotter(object):
    '''
    This is a semi-generic plotting interface that has a built in curses based terminal plotter.
    It's fairly specific to what we're using it for here, but we could (and maybe should) build it out into
    a little library that we can use via the command line to plot things.  Might be useful for looking at data later.
    That would also cut the size of this tool down by a good bit.
    '''
    #def __init__(self, kinavg, kinrw, iteration, bin_labels, state_labels, state_pops, bin_pops, interface='matplotlib'):
    def __init__(self, h5file, h5key, iteration=-1, interface='matplotlib'):
        # Need to sort through and fix all this, but hey.
        self.iteration = iteration
        # These two are important for... reasons.
        try:
            self.bin_labels = list(bin_labels[...])
            self.state_labels = list(state_labels[...]) + ['unknown']
        except:
            try:
                self.state_labels = list(h5file['state_labels'][...]) + ['unknown']
            except:
                self.state_labels = None
        # unless we totally fail out.
        self.interface = interface
        # What we should ACTUALLY do is just... yeah, just have it sub in what we need.
        # We'll need to throw in the state labels or whatever, but.
        self.h5file = h5file
        self.h5key = h5key
        # We should determine the number of dimensions of our dataset...
        # This has time data, so an i to j is a 3 dim, and an i is 2.
        try:
            self.dim = len(h5file[h5key].shape)
        except:
            self.dim = 1
        try:
            # does the ci exist?
            a = h5file[h5key]['expected']
        except:
            self.dim = 1

    def plot(self, i=0, j=1, tau=1, iteration=None, dim=0, interface=None):
        if iteration == None:
            iteration = self.iteration
            self.__generic_ci__(self.h5file, iteration, i, j, tau=tau, h5key=self.h5key, dim=dim, interface=interface)

    def __generic_ci__(self, h5file, iteration, i, j, tau, h5key='rate_evolution', dim=0, interface=None):
        # This function just calls the appropriate plot function for our available
        # interface.
        if (interface == None and self.interface == 'text') or interface == 'text':
            if self.dim > 1:
                self.__terminal_ci__(h5file, iteration, i, j, tau, h5key)
            else:
                self.__terminal_expected__(h5file, iteration, i, j, tau, h5key, dim)
        else:
            try:
                import matplotlib
                matplotlib.use('TkAgg')
                from matplotlib import pyplot as plt
                if self.dim == 3:
                    plt.plot(h5file[h5key]['expected'][:iteration, i, j] / tau, color='black')
                    plt.plot(h5file[h5key]['ci_ubound'][:iteration, i, j] / tau, color='grey')
                    plt.plot(h5file[h5key]['ci_lbound'][:iteration, i, j] / tau, color='grey')
                else:
                    plt.plot(h5file[h5key]['expected'][:iteration, i] / tau, color='black')
                    plt.plot(h5file[h5key]['ci_ubound'][:iteration, i] / tau, color='grey')
                    plt.plot(h5file[h5key]['ci_lbound'][:iteration, i] / tau, color='grey')
                plt.show()
            except:
                print('Unable to import plotting interface.  An X server ($DISPLAY) is required.')
                if self.dim > 1:
                    self.__terminal_ci__(h5file, iteration, i, j, tau)
                else:
                    self.__terminal_expected__(h5file, iteration, i, j, tau, h5key, dim)
                return 1

    def __generic_histo__(self, vector, labels):
        # This function just calls the appropriate plot function for our available
        # interface.  Same thing as generic_ci, but for a histogram.
        if self.interface == 'text':
            self.__terminal_histo__(vector, labels)
        else:
            try:
                import matplotlib
                matplotlib.use('TkAgg')
                from matplotlib import pyplot as plt
                plt.bar(list(range(0, np.array(vector).shape[0])), vector, linewidth=0, align='center', color='gold', tick_label=labels)
                plt.show()
            except:
                print('Unable to import plotting interface.  An X server ($DISPLAY) is required.')
                self.__terminal_histo__(h5file, vector, labels)
                return 1

    def __terminal_histo__(self, vector, labels, fullscreen_mode=True):
        from blessings import Terminal

        self.t = Terminal()
        h = int(self.t.height / 4) * 3
        w = self.t.width
        cols = np.array(vector).shape[0]
        # Let's print this business!

        colwidth = w / cols
        with self.t.fullscreen():
            for y in range(0, h):
                for x in range(0, cols):
                    if x == 0:
                        with self.t.location(0, y):
                            print(self.t.red('{0:.4f}|'.format(float(h-y)/float(h))))
                    with self.t.location((x*colwidth)+8+len(labels[x])/2, y):
                        if vector[x] >= (float(h-y)/float(h)):
                            #print(float(h-y)/float(h))
                            print(self.t.on_blue(' '))
            for x in range(0, cols):
                if x == 0:
                    with self.t.location(x, h):
                        print('States| ')
                with self.t.location((x*colwidth)+8, h):
                    print(self.t.blue(labels[x]))

            if fullscreen_mode:
                input("Press enter to continue.")

    def __terminal_ci__(self, h5file, iteration, si, sj, tau, h5key):
        from blessings import Terminal

        self.t = Terminal()
        h = int(self.t.height / 4 * 3.75)
        # We'll figure out how to subsample the timepoints...
        w = self.t.width
        if self.dim == 3:
            in_tup = (iteration-1, si, sj)
        else:
            in_tup = (iteration-1, si)
        yupper = (h5file[h5key]['ci_ubound'][in_tup] / tau) * 2
        ylower = (h5file[h5key]['ci_lbound'][in_tup] / tau) / 2
        # Here are points pertaining to height.
        scale = np.array([0.0] + [ylower+i*(yupper-ylower)/np.float(h) for i in range(0, h)])[::-1]
        if iteration > w:
            block_size = iteration / w
        else:
            block_size = 1

        with self.t.fullscreen():
            try:
                for x in range(0, w-12):
                    iter = x * block_size
                    if self.dim == 3:
                        in_tup = (iter-1, si, sj)
                    else:
                        in_tup = (iter-1, si)
                    yupper = (h5file[h5key]['ci_ubound'][in_tup] / tau)
                    ylower = (h5file[h5key]['ci_lbound'][in_tup] / tau)
                    ci = np.digitize([yupper, ylower], scale)
                    if x == 0:
                        for y in range(0, h+1):
                            with self.t.location(0, y):
                                print(self.t.bold(self.t.red('{0:.7f}|'.format(scale[y]))))
                    for y in range(ci[0], ci[1]):
                        #with self.t.location(x+12, y):
                        print(self.t.move(y, x+12) + self.t.on_blue(' '))
                            #print(self.t.on_blue(' '))
                    #with self.t.location(x+12, np.digitize(h5file['rate_evolution']['expected'][iter-1, si, sj]/tau, scale)):
                    #        print(self.t.on_blue('-'))
                    print(self.t.move(np.digitize(h5file[h5key]['expected'][in_tup]/tau, scale), x+12) + self.t.on_blue('-'))

                for x in range(0, w-12, w/10):
                    if x == 0:
                        with self.t.location(x, h+1):
                            print('Iteration| ')
                    with self.t.location(x+12, h+1):
                        iter = x * block_size
                        print(self.t.blue(str(iter)))
            except:
                pass

            with self.t.location(0, h+2):
                # We need to improve this.
                #if h5key == 'rate_evolution':
                #    print("k_ij from {} to {} from iter 1 to {}".format(self.state_labels[si], self.state_labels[sj], self.iteration))
                #elif h5key == 'conditional_flux_evolution':
                #    print("i->j flux from {} to {} from iter 1 to {}".format(self.state_labels[si], self.state_labels[sj], self.iteration))
                if self.dim == 3:
                    print("{} from {} to {} from iter 1 to {}".format(h5key, self.state_labels[si], self.state_labels[sj], self.iteration))
                else:
                    print("{} of state {} from iter 1 to {}".format(h5key, self.state_labels[si], self.iteration))
            with self.t.location(0, h+3):
                input("Press enter to continue.")

    def __terminal_expected__(self, h5file, iteration, si, sj, tau, h5key, dim):
        from blessings import Terminal

        self.t = Terminal()
        h = int(self.t.height / 4 * 3.75)
        # We'll figure out how to subsample the timepoints...
        w = self.t.width
        if self.dim == 3:
            in_tup = (iteration-1, si, sj)
        else:
            in_tup = (iteration-1, si)
        in_tup = (iteration-1, dim)
        try:
            yupper = (np.max(h5file) / tau) * 2
        except:
            in_tup = (iteration-1)
            yupper = (np.max(h5file) / tau) * 2
        ylower = (np.min(h5file) / tau) * 2
        # Here are points pertaining to height.
        if yupper > 0:
            yupper = (np.max(h5file) / tau) * 1.2
            ylower = (np.min(h5file) / tau) / 2
            scale = np.array([0.0] + [ylower+i*(yupper-ylower)/np.float(h) for i in range(0, h)])[::-1]
        else:
            yupper = (np.max(h5file) / tau) / 2
            ylower = (np.min(h5file) / tau) * 1.2
            scale = np.array([ylower+i*(yupper-ylower)/np.float(h) for i in range(0, h)] + [0.0])[::-1]
        if iteration > w:
            block_size = iteration / w
        else:
            block_size = 1

        with self.t.fullscreen():
            try:
                for x in range(0, w-12):
                    iter = x * block_size
                    if self.dim == 3:
                        in_tup = (iter-1, si, sj)
                    else:
                        in_tup = (iter-1, si)
                    in_tup = (iter-1, dim)
                    try:
                        yupper = (h5file[in_tup] / tau)
                    except:
                        in_tup = (iter-1)
                        yupper = (h5file[in_tup] / tau)
                    ylower = (h5file[in_tup] / tau)
                    ci = np.digitize([yupper, ylower], scale)
                    if x == 0:
                        for y in range(0, h+1):
                            with self.t.location(0, y):
                                print(self.t.bold(self.t.red('{0:.7f}|'.format(scale[y]))))
                    for y in range(ci[0], ci[1]):
                        print(self.t.move(y, x+12) + self.t.on_blue(' '))
                            #print(self.t.on_blue(' '))
                    print(self.t.move(np.digitize(h5file[in_tup]/tau, scale), x+12) + self.t.on_blue('-'))

                for x in range(0, w-12, w/10):
                    if x == 0:
                        with self.t.location(x, h+1):
                            print('Iteration| ')
                    with self.t.location(x+12, h+1):
                        iter = x * block_size
                        print(self.t.blue(str(iter)))
            except:
                pass

            with self.t.location(0, h+2):
                # We need to improve this.
                print("{} from iter 1 to {}".format(h5key, self.iteration))
            with self.t.location(0, h+3):
                input("Press enter to continue.")
