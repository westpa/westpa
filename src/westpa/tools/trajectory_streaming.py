import logging
import MDAnalysis as mda
from imdclient.IMD import IMDReader
import subprocess
import time
import re
import socket
import os
import sys

# TODO: Add options for LAMMPS and NAMD
ACCEPTABLE_MD_ENGINES = ["gromacs"]
IMD_FLAGS = {"gromacs": {"-imdwait": None, "-imdport": "0"}}
IMD_PORT_OUTPUT = {"gromacs": r"IMD connection on port (\d+)"}

log = logging.getLogger("TrajectoryStreamer")
log.setLevel(logging.INFO)
if not log.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)


class TrajectoryStreamer:
    """
    A class for streaming trajectory data with user-defined simulation and analysis functions.
    Utilizes the IMDv3 protocol for communication between the simulation engine and the analysis tools.
    See https://imdclient.readthedocs.io/en/latest/protocol_v3.html for more details.

    The user provides:
    - A simulation function that generates trajectory data
    - Topology file for the simulation
    """

    def __init__(self, md_engine: str, topology: str, simulation_string: str, host: str = "localhost"):
        """
        Initialize the trajectory streamer.

        Parameters
        ----------
            md_engine : str
                Name of the molecular dynamics engine (e.g., 'gromacs')
            topology : str
                Path to the topology file for the simulation. Any format supported by MDAnalysis.
            simulation_string : str
                String representing the simulation command to be executed.
            host : str, optional
                Hostname or IP address of the machine running the simulation. Currently only localhost is supported. (default=localhost)
        """
        self.md_engine = md_engine.lower()
        if self.md_engine not in ACCEPTABLE_MD_ENGINES:
            raise ValueError(f"Unsupported MD engine: {self.md_engine}. " f"Supported engines: {', '.join(ACCEPTABLE_MD_ENGINES)}")
        self.topology = topology
        # Check for the existence of the topology files
        if not os.path.exists(self.topology):
            raise FileNotFoundError(f"Topology file not found: {self.topology}")
        self.set_simulation_function(simulation_string)
        self.host = host

    def set_simulation_function(self, sim_func: str):
        """
        Sets the simulation function.
        """
        self.simulation_function = sim_func.split()
        # Check and fix IMD flags
        md_engine_flags = IMD_FLAGS[self.md_engine]
        for flag in md_engine_flags:
            if md_engine_flags[flag] is None:
                if flag not in self.simulation_function:
                    self.simulation_function.append(flag)
                    log.warning(f"Adding {flag} to simulation function")
            else:
                if flag not in self.simulation_function:
                    self.simulation_function.append(f"{flag}")
                    self.simulation_function.append(md_engine_flags[flag])
                    log.warning(f"Adding {flag} with value {md_engine_flags[flag]} to simulation function")
                elif self.simulation_function[self.simulation_function.index(flag) + 1] != md_engine_flags[flag]:
                    self.simulation_function[self.simulation_function.index(flag) + 1] = md_engine_flags[flag]
                    log.warning(f"Updating {flag} with value {md_engine_flags[flag]} in simulation function")

    def start_sim_and_get_universe(self, stream_timeout: float = 5.0) -> mda.Universe:
        """
        Start the simulation and return the MDAnalysis universe.

        Parameters
        ----------
            stream_timeout : float, optional
                Timeout for the IMD connection in seconds.
                Important if the time between messages from the engine is long. (default=5.0)
        """
        if self.simulation_function is None:
            raise ValueError("No simulation function has been set")

        log.info(f"Launching simulation with {self.md_engine} engine")

        proc = subprocess.Popen(self.simulation_function, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)

        assigned_port = None
        retcode = proc.poll()
        if retcode is not None and retcode != 0:
            raise RuntimeError(f"Simulation returned with error code {retcode}. Check the simulation  for details.")

        start_time = time.time()
        for line in proc.stdout:
            log.info(line.strip())
            m = re.search(IMD_PORT_OUTPUT[self.md_engine], line)
            if m:
                assigned_port = int(m.group(1))
                break
            if time.time() - start_time > 60:  # 1 minute timeout
                raise RuntimeError("IMD port assignment was not printed within 1 minute. Make sure an IMD simulation is being run.")
        else:
            raise RuntimeError(
                f"{self.md_engine.upper()} output did not contain expected '{IMD_PORT_OUTPUT[self.md_engine]}' pattern. Check the simulation  for details."
            )
        log.info(f"Assigned IMD port: {assigned_port}")

        port = assigned_port
        port_open = False
        timeout = 0.2
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        while not port_open:
            #     sock.settimeout(timeout)
            try:
                sock.connect((self.host, port))
            except ConnectionRefusedError:
                time.sleep(timeout)
            else:
                log.info(f"Port {port} on {self.host} is now open!")
                port_open = True

        u = mda.Universe(self.topology, f"imd://{self.host}:{port}", timeout=stream_timeout)
        self.sim_process = proc
        return u

    def end_sim(self):
        """
        Return the remaining output from the simulation process and the end the simulation.
        """
        log.info("Dumping remaining output from the simulation...")
        for line in iter(self.sim_process.stdout.readline, ''):
            log.info(line.strip())
        if self.sim_process.poll() is None:
            log.warning(
                "Simulation is still running. Something likely went wrong. Make sure that the whole simulation was analysed. Sending a termination signal."
            )
            self.sim_process.terminate()

            # Wait a bit for the process to terminate gracefully
            try:
                self.sim_process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                log.warning("Process did not terminate in time. Forcibly killing it.")
                self.sim_process.kill()
                self.sim_process.wait()  # Ensure the process is fully gone

        log.info("Segment Complete")

    def find_port(interface='localhost'):
        """Generic function to find an open port on the local machine."""
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.bind((interface, 0))  # Bind to an ephemeral port
        port = sock.getsockname()[1]
        sock.close()
        return port
