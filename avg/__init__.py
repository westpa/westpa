'''
Set of tools for use with averaging together multiple WESTPA simulations.
'''

import logging
log = logging.getLogger(__name__)


def _average_2D_data_slice(slice, dataset):
    # A function which takes in a slice, and a dataset, and adds on 
    # to it.
