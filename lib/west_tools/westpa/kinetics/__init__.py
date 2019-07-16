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

'''
Kinetics analysis library
'''

import logging
log = logging.getLogger(__name__)


from .rate_averaging import RateAverager

from . import _kinetics
from ._kinetics import (calculate_labeled_fluxes, labeled_flux_to_rate,
                        calculate_labeled_fluxes_alllags,
                        nested_to_flat_matrix, nested_to_flat_vector,
                        flat_to_nested_matrix, flat_to_nested_vector, find_macrostate_transitions,
                        sequence_macro_flux_to_rate)
from .events import WKinetics

