import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=FutureWarning)
import numpy as np
import h5py
# We don't care if we're reimporting blessings, under this context.
warnings.filterwarnings('ignore')

class Plotter():
    '''
    This is a semi-generic plotting interface that has a built in curses based terminal plotter.
    It's fairly specific to what we're using it for here, but we could (and maybe should) build it out into
    a little library that we can use via the command line to plot things.  Might be useful for looking at data later.
    That would also cut the size of this tool down by a good bit.
    '''
    def __init__(self, kinavg, kinrw, iteration, bin_labels, state_labels, state_pops, bin_pops, interface='matplotlib'):
        # Need to sort through and fix all this, but hey.
        self.kinavg_file = kinavg
        self.kinrw_file = kinrw
        self.stateprobs_file = kinavg
        self.iteration = iteration
        self.bin_labels = list(bin_labels[...])
        self.state_labels = list(state_labels[...]) + ['unknown']
        self.state_populations = state_pops
        self.bin_populations = bin_pops
        self.interface = interface

    def __generic_ci__(self, h5file, iteration, i, j, tau, h5key='rate_evolution'):
        if self.interface == 'text':
            self.__terminal_ci__(h5file, iteration, i, j, tau, h5key)
        else:
            try:
                import matplotlib
                matplotlib.use('TkAgg')
                from matplotlib import pyplot as plt
                plt.plot(h5file[h5key]['expected'][:iteration, i, j] / tau, color='black')
                plt.plot(h5file[h5key]['ci_ubound'][:iteration, i, j] / tau, color='grey')
                plt.plot(h5file[h5key]['ci_lbound'][:iteration, i, j] / tau, color='grey')
                plt.show()
            except:
                print('Unable to import plotting interface.  An X server ($DISPLAY) is required.')
                self.__terminal_ci__(h5file, iteration, i, j, tau)
                return 1

    def __generic_histo__(self, vector, labels):
        if self.interface == 'text':
            self.__terminal_histo__(vector, labels)
        else:
            try:
                import matplotlib
                matplotlib.use('TkAgg')
                from matplotlib import pyplot as plt
                plt.bar(range(0, np.array(vector).shape[0]), vector, linewidth=0, align='center', color='gold', tick_label=labels)
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
                raw_input("Press enter to continue.")

    def __terminal_ci__(self, h5file, iteration, si, sj, tau, h5key='rate_evolution'):
        from blessings import Terminal

        self.t = Terminal()
        h = int(self.t.height / 4) * 3
        # We'll figure out how to subsample the timepoints...
        w = self.t.width
        yupper = (h5file[h5key]['ci_ubound'][iteration-1, si, sj] / tau) * 2
        ylower = (h5file[h5key]['ci_lbound'][iteration-1, si, sj] / tau) / 2
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
                    yupper = (h5file[h5key]['ci_ubound'][iter-1, si, sj] / tau)
                    ylower = (h5file[h5key]['ci_lbound'][iter-1, si, sj] / tau)
                    ci = np.digitize([yupper, ylower], scale)
                    if x == 0:
                        for y in range(0, h+1):
                            with self.t.location(0, y):
                                print(self.t.bold(self.t.red('{0:.7f}|'.format(scale[y]))))
                    for y in range(ci[0], ci[1]):
                        #with self.t.location(x+12, y):
                        print self.t.move(y, x+12) + self.t.on_blue(' ')
                            #print(self.t.on_blue(' '))
                    #with self.t.location(x+12, np.digitize(h5file['rate_evolution']['expected'][iter-1, si, sj]/tau, scale)):
                    #        print(self.t.on_blue('-'))
                    print self.t.move(np.digitize(h5file[h5key]['expected'][iter-1, si, sj]/tau, scale), x+12) + self.t.on_blue('-')

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
                if h5key == 'rate_evolution':
                    print("k_ij from {} to {} from iter 1 to {}".format(self.state_labels[si], self.state_labels[sj], self.iteration))
                else:
                    print("{} state population from iter 1 to {}".format(self.state_labels[si], self.iteration))
            with self.t.location(0, h+3):
                raw_input("Press enter to continue.")
            


    def kinavg(self, i=0, j=1, tau=1):
        self.__generic_ci__(self.kinavg_file, self.iteration, i, j, tau)

    def kinrw(self, i=0, j=1, tau=1):
        self.__generic_ci__(self.kinrw_file, self.iteration, i, j, tau)

    def rwstateprobs(self, i=0, j=1, tau=1):
        self.__generic_ci__(self.kinrw_file, self.iteration, i, None, 1, h5key='state_pop_evolution')

    def stateprobs(self, i=0, j=1, tau=1):
        self.__generic_ci__(self.stateprobs_file, self.iteration, i, None, 1, h5key='state_pop_evolution')

    def states(self):
        self.__generic_histo__(self.state_populations, self.state_labels)

    def bins(self):
        self.__generic_histo__(self.bin_populations, self.bin_labels)
