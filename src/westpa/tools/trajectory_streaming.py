import logging
import MDAnalysis as mda
from imdclient.IMD import IMDReader
import subprocess
import time
import re
from typing import Optional
import socket
 
log = logging.getLogger(__name__)

# TODO: Add options for LAMMPS and NAMD
ACCEPTABLE_MD_ENGINES = ["gromacs"]
IMD_FLAGS = {"gromacs": {"-imdwait": None, "-imdport": "0"}}
IMD_PORT_OUTPUT = {"gromacs": r"IMD connection on port (\d+)"}



class TrajectoryStreamer:
    """
    A class for streaming trajectory data with user-defined simulation and analysis functions.

    The user provides:
    - A simulation function that generates trajectory data
    - Topology file for the simulation
    """

    def __init__(self, md_engine: str, topology: str, simulation_string: str):
        """
        Initialize the trajectory streamer.

        Args:
            md_engine: Name of the molecular dynamics engine (e.g., 'gromacs')
            topology: Path to the topology file for the simulation. Any format supported by MDAnalysis.
            simulation_string: String representing the simulation command to be executed.
        """
        self.md_engine = md_engine.lower()
        if self.md_engine not in ACCEPTABLE_MD_ENGINES:
            raise ValueError(f"Unsupported MD engine: {self.md_engine}. " f"Supported engines: {', '.join(ACCEPTABLE_MD_ENGINES)}")
        self.topology = topology
        self.set_simulation_function(simulation_string)

    def set_simulation_function(self, sim_func: str):
        """
        Set the simulation function.

        Args:
            sim_func: String representing the simulation function to be used.
        """
        self.simulation_function = sim_func.split()
        # Check and fix IMD flags
        md_engine_flags = IMD_FLAGS[self.md_engine]
        for flag in md_engine_flags:
            if md_engine_flags[flag] is None:
                if flag not in self.simulation_function:
                    self.simulation_function.append(flag)
                    print(f"Adding {flag} to simulation function")
            else:
                if flag not in self.simulation_function:
                    self.simulation_function.append(f"{flag}")
                    self.simulation_function.append(md_engine_flags[flag])
                    print(f"Adding {flag} with value {md_engine_flags[flag]} to simulation function")
                elif self.simulation_function[self.simulation_function.index(flag) + 1] != md_engine_flags[flag]:
                    self.simulation_function[self.simulation_function.index(flag) + 1] = md_engine_flags[flag]
                    print(f"Updating {flag} with value {md_engine_flags[flag]} in simulation function")

    def start_sim_and_get_universe(self):
        """
        Start the simulation and return the MDAnalysis universe.
        """
        if self.simulation_function is None:
            raise ValueError("No simulation function has been set")


        print(f"Launching simulation with {self.md_engine} engine")

        proc = subprocess.Popen(self.simulation_function, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)

        assigned_port = None
        retcode = proc.poll()
        if retcode is not None and retcode != 0:
            raise RuntimeError(f"Simulation returned with error code {retcode}. Check the simulation log for details.")

        start_time = time.time()
        for line in proc.stdout:
            print(line, end="")
            m = re.search(IMD_PORT_OUTPUT[self.md_engine], line)
            if m:
                assigned_port = int(m.group(1))
                break
            if time.time() - start_time > 60:  # 1 minute timeout
                raise RuntimeError(
                    "IMD port assignment was not printed within 1 minute. Make sure an IMD simulation is being run."
                )
        else:
            raise RuntimeError(
                f"{self.md_engine.upper()} output did not contain expected '{IMD_PORT_OUTPUT[self.md_engine]}' pattern. Check the simulation log for details."
            )
        print(f"Assigned IMD port: {assigned_port}")

        port = assigned_port
        port_open = False
        timeout = 0.2
        host = "localhost"
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        while not port_open:
            #     sock.settimeout(timeout)
            try:
                sock.connect((host, port))
            except ConnectionRefusedError:
                time.sleep(timeout)
            else:
                print(f"Port {port} on {host} is now open!")
                port_open = True

        u = mda.Universe(self.topology, f"imd://{host}:{port}")
        return u
    
    def find_port():
        """Generic function to find an open port on the local machine."""
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.bind(('', 0))  # Bind to an ephemeral port
        port = sock.getsockname()[1]
        sock.close()
        return port