import os
import time
from dataclasses import dataclass, field
from datetime import datetime

import numpy as np

from westpa.core.h5io import WESTPAH5File
from westpa.core.run_status import RUN_STATE_COMPLETE
from westpa.core.segment import Segment


@dataclass
class SegmentStatusCounts:
    total: int = 0
    unset: int = 0
    prepared: int = 0
    complete: int = 0
    failed: int = 0
    other: int = 0


