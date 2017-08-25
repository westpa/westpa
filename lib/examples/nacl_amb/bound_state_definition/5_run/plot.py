#!/usr/bin/env python
import numpy
import scipy.stats
import matplotlib
import matplotlib.pyplot as pyplot

class DistancePDF(object):
    '''
    Plot an estimate of the probability density function describing the
    distribution of distances between Na+ and Cl- over the course of the
    simulation.  Perform the estimate using kernel density estimation,
    with the bandwidth set 0.1 Angstroms.
    '''
    def __init__(self):
        matplotlib.rcParams['font.size'] = 7
        self.distance_arr = numpy.loadtxt('./distance.dat', skiprows=1, 
                                          usecols=(1,))
        n_datapoints = self.distance_arr.shape[0]

        # Use the second half of the simulation
        self.distance_arr = self.distance_arr[n_datapoints/2:]
        self.kde = scipy.stats.gaussian_kde(self.distance_arr,
                                            bw_method=0.1)
        xs = numpy.arange(0, 20.1, 0.1, dtype=float)
        ys = self.kde.evaluate(xs)
        self.fig = pyplot.gcf()
        self.ax = pyplot.gca()
        self.fig.set_size_inches(8.8/2.54,6/2.54)
        self.ax.plot(xs, ys, color='black')

        bound_state_idx = numpy.argmax(ys)
        bound_state = xs[bound_state_idx]
        self.ax.plot((bound_state, bound_state), 
                     (0, 5),
                     color='black',
                     linestyle='dashed')
        self.ax.text(bound_state, ys[bound_state_idx], 
                     u' Bound state (separation = {:.01f} \u00c5)'.format(bound_state),
                     ha='left',
                     va='bottom')

        self.ax.set_xlim(0,20)
        self.ax.set_ylim(0,4)
        self.ax.set_xlabel(u'Separation (\u00c5)')
        self.ax.set_ylabel(u'Probability density')
        self.ax.set_xticks(numpy.arange(0,25,5))
        self.ax.set_xticks(numpy.arange(0,21,1), minor=True)
        self.ax.tick_params(direction='out')
        self.fig.subplots_adjust(left=0.2, bottom=0.2)
        pyplot.savefig('probability_distribution.pdf')

if __name__ == "__main__":
    DistancePDF()

        
