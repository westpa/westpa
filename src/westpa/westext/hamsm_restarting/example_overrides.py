import logging

from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()])
log = logging.getLogger("msm_we")
log.setLevel(logging.DEBUG)

import numpy as np


def processCoordinates(self, coords):
    log.debug("Processing coordinates")

    if self.dimReduceMethod == "none":
        nC = np.shape(coords)
        nC = nC[0]
        data = coords.reshape(nC, 3 * self.nAtoms)
        return data

    if self.dimReduceMethod == "pca" or self.dimReduceMethod == "vamp":

        ### NaCl dimensionality reduction
        log.warning("Hardcoded selection: Doing dim reduction for Na, Cl. This is only for testing!")
        indNA = self.reference_structure.topology.select("element Na")
        indCL = self.reference_structure.topology.select("element Cl")

        diff = np.subtract(coords[:, indNA], coords[:, indCL])

        dist = np.array(np.sqrt(np.mean(np.power(diff, 2), axis=-1)))

        return dist


log.info("Loading user-override functions for modelWE")
