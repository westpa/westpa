from __future__ import division 
__metaclass__ = type

import os, sys
import socket
import threading
import SocketServer
import time
import multiprocessing
import cPickle as pickle
from socket import timeout

from default import DefaultWEMaster
from wemd.sim_managers import WESimManagerBase, WESimClient

from wemd.core.we_sim import WESimIter
from wemd.core.particles import Particle, ParticleCollection
from wemd.core.segments import Segment
from wemd.core.errors import PropagationIncompleteError
from wemd.rc import EX_ERROR

import logging
log = logging.getLogger(__name__)

#For the watchdog (server side)
class ThreadedTCPRequestHandler(SocketServer.BaseRequestHandler):

    def handle(self):

        data = self.request.recv(1024)

        try:
            cid, cmsg = pickle.loads(data)
        except pickle.UnpicklingError:
            log.error(('Error in watchdog unpickling %r\n' % data))
            raise
 
        assert(cmsg == 'client here')

        self.request.send('server here')

        if cid >= self.server.nclients:
            log.error(('Watchdog cid out of range %d\n' % cid))
            raise ValueError

        #update last time client was heard from
        self.server.client_list[cid] = time.time()

class ThreadedTCPServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer):

    def __init__(self, addr, handlerclass, nclients):

        for superclass in self.__class__.__bases__:
            initfunc = getattr(superclass, '__init__', lambda x, y, z: None)
            initfunc(self, addr, handlerclass)

        self.nclients = nclients
        cur_time = time.time()
        self.client_list = [cur_time for i in range(0, nclients)]

def watchdog_check(server, flag, timeout):
    while True:
        curtime = time.time()

        for cid in range(0,len(server.client_list)):

            ctime = server.client_list[cid]

            if curtime - ctime > timeout:
                flag.set()
                log.debug('watchdog flag set '+repr(flag.is_set())+' \n')       
                return

        time.sleep(timeout)

class TCPWEMaster(DefaultWEMaster):
    def __init__(self):
        super(TCPWEMaster,self).__init__()
        self.hostname = None
        self.dport = None #data
        self.sport = None #sync
        self.wport = None #watchdog
        self.nclients = None
        self.sync_sock = None
        self.sock = None
        self.timeout = None #watchdog timeout in s
        #self.cores[cid] = number of cores for client 'cid'
        self.cores = None
        self.watchdog_server = None
        self.watchdog_server_thread = None
        self.watchdog_check_thread = None
        self.watchdog_flag = threading.Event() #set to true when shutdown should occur

    def set_hostname(self, hostname):
        self.hostname = hostname
        
    def runtime_init(self, runtime_config, load_sim_config = True):
        super(TCPWEMaster, self).runtime_init(runtime_config, load_sim_config)
        self.timeout = runtime_config.get_int('server.timeout')
        self.dport = runtime_config.get_int('server.dport')
        self.sport = runtime_config.get_int('server.sport')
        self.wport = runtime_config.get_int('server.wport')
        self.nclients = runtime_config.get_int('server.nclients')  

    def init_server(self):
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        #So server can be restarted quickly
        self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.sock.bind((self.hostname, self.dport))

        #queue up to nclients
        self.sock.listen(self.nclients)
        self.sock.settimeout(self.timeout)

    def init_watchdog_server(self):
        self.watchdog_server = ThreadedTCPServer((self.hostname, self.wport), ThreadedTCPRequestHandler, self.nclients)
        self.watchdog_server_thread = threading.Thread(target=self.watchdog_server.serve_forever)
        self.watchdog_server_thread.daemon = True
        self.watchdog_server_thread.start()

        self.watchdog_check_thread = threading.Thread(target=watchdog_check, args=(self.watchdog_server, self.watchdog_flag, self.timeout))
        self.watchdog_check_thread.daemon = True
        self.watchdog_check_thread.start()

    def init_sync_server(self):
        self.sync_sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        #So server can be restarted quickly
        self.sync_sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.sync_sock.bind((self.hostname, self.sport))

        #queue up to nclients
        self.sync_sock.listen(self.nclients)
        self.sync_sock.settimeout(self.timeout)
    
    def run(self, hostname = None):
        if self.backend_driver is None:
            self.load_backend_driver()
            
        if self.we_driver is None:
            self.load_we_driver() #runtime_init

        self.init_server()
        
        if self.send_cid() == False:
            self.shutdown(EX_ERROR)
            sys.exit(EX_ERROR)

        self.init_watchdog_server()
        
        max_wallclock = self.max_wallclock   
        if( max_wallclock is not None):     
            we_cur_wallclock = time.time() - self.start_wallclock
            loop_start_time = loop_end_time = None
                    
        while self.continue_simulation():
            if( max_wallclock is not None):
                if( loop_end_time is not None):
                    loop_duration = loop_end_time - loop_start_time
                    we_cur_wallclock += loop_duration
                    if( we_cur_wallclock + loop_duration * 2.0 > max_wallclock ):
                        log.info('Shutdown so walltime does not exceed max wallclock:%r'%(max_wallclock))                        
                        self.shutdown(0)
                        sys.exit(0)

                loop_start_time = time.time()                        

            if self.watchdog_flag.is_set():
                self.shutdown(EX_ERROR)
                sys.exit(EX_ERROR)
                                        
            self.prepare_iteration()
                        
            self.backend_driver.pre_iter(self.we_iter)
            if self.propagate_particles() == False:
                self.shutdown(EX_ERROR)
                sys.exit(EX_ERROR)
            
            self.run_we()
            self.backend_driver.post_iter(self.we_iter)
            self.finalize_iteration()
            
            if self.sync() == False:
                self.shutdown(EX_ERROR)
                sys.exit(EX_ERROR)
            
            if( max_wallclock is not None):
                loop_end_time = time.time()

    def send_cid(self):
        '''assigns an id to each client, and gets the number of cores available
           returns True on success, False on failure'''
        
        
        msgs = []
        for i in range(0, self.nclients):
            msg = (i, 'cid')
            msg_pickle = pickle.dumps(msg, pickle.HIGHEST_PROTOCOL)
            msgs.append(msg_pickle)

        data = self.send_to_clients(msgs, timeout = True)
        
        if data is None:
            return False

        msgs = []
        for i in range(0, self.nclients):
            msg = (i, 'cores')
            msg_pickle = pickle.dumps(msg, pickle.HIGHEST_PROTOCOL)
            msgs.append(msg_pickle)

        self.cores = self.send_to_clients(msgs, timeout = True)
        
        if self.cores is None:
            return False

        #this can happen if the # of clients in run.cfg is wrong 
        if 'cid' in self.cores:
            return False

        return True

    def propagate_particles(self):
        current_iteration = self.we_iter.n_iter
        log.info('WE iteration %d (of %d requested)'
                 % (current_iteration, self.max_iterations))
        n_inc = self.data_manager.num_incomplete_segments(self.we_iter)
        log.info('%d segments remaining in WE iteration %d'
                 % (n_inc, current_iteration))
        #segments = self.data_manager.get_prepared_segments(self.we_iter)
        segments = self.data_manager.get_segments((Segment.n_iter == self.we_iter.n_iter)
                                                  &(Segment.status == Segment.SEG_STATUS_PREPARED),
                                                  load_p_parent = True)

        if self.send_segments(segments) == False:
            return False
        
        segments = self.get_segments()
        
        if segments is None:
            return False
        
        #segments[cid][task#] -> segments[i]
        segments = [j for k in segments for j in k if j is not None]
        self.data_manager.update_segments(self.we_iter, segments)
        
        return True

    def send_segments(self, segments):

        cmds = []
        
        #for now, just assume a homogeneous resource
        i = 0
        tasks = [[] for i in xrange(0, self.nclients)]
        
        for segment in map(None, *(iter(segments),) * self.cores[0]):

            if type(segment) is tuple:
                tasks[i].append(segment)                
            else:    
                tasks[i].append((segment,))
                        
            i += 1
            i %= self.nclients
            
        for i in xrange(0, self.nclients):
            #flatten list of tuples into a list
            task = [j for k in tasks[i] for j in k if j is not None]
            command = (i, task)
            work_pickle = pickle.dumps(command, pickle.HIGHEST_PROTOCOL)
            cmds.append(work_pickle)
            
        ret = self.send_to_clients(cmds)

        if ret is None:
            return False

        return True
                        
    def get_segments(self):

        msgs = []
        for i in range(0, self.nclients):
            msg = (i, 'prepare for next commnand')
            msg_pickle = pickle.dumps(msg, pickle.HIGHEST_PROTOCOL)
            msgs.append(msg_pickle)

        segs = self.send_to_clients(msgs, block = False)
        if segs is None:
            return None

        for cid in range(0, len(segs)):

            if segs[cid] is None:
                log.error(('client: %d returned %r\n' % (cid, segs[cid])))
                return None

            log.debug(('client: %d returned %r\n' % (cid, segs[cid])))
             
        return segs

    def sync(self, shutdown = False):
        '''sends sync command to clients on sport
           returns False on timeout, True otherwise'''

        self.init_sync_server()

        sock = self.sync_sock

        clients_alive = []
        while len(clients_alive) < self.nclients:

            if self.watchdog_flag.is_set():
                return False

            try:
                newSocket, address = sock.accept()

                data = newSocket.recv(1024)

                try:
                    cid, cmsg = pickle.loads(data)
                except (pickle.UnpicklingError, EOFError):
                    newSocket.close()
                    self.shutdown_sync_server()
                    return False
                    

                if shutdown == False:
                    newSocket.sendall('sync')
                else:
                    newSocket.sendall('shutdown')

                newSocket.close()

                #only count sync clients
                if cmsg == 'sync' and cid not in clients_alive:
                    clients_alive.append(cid)

            except (socket.timeout, socket.error):
                pass

        self.shutdown_sync_server()
        return True

    def send_to_clients(self, messages, block = True, timeout = False):
        '''Send message to clients
           Returns received data[client_id] or None on timeout
           messages[client_id] is message for client client_id
           set block to False to send message as clients connect
           rather then after all clients have connected'''

        sock = self.sock
        client_sockets = []
        client_sockets_cid = []

        #block until all clients have connected
        client_data = [0 for x in range(0, self.nclients)]

        nclients_alive = 0
        while nclients_alive < self.nclients:

            if self.watchdog_flag.is_set():

                if block == True:
                    for i in range(0, nclients_alive):
                        client_sockets[i].close()

                return None

            try:
                newSocket, address = sock.accept()
                
                log.debug('connected to:'+repr(address)+'\n')
                client_sockets.append(newSocket)

                buf = []
                data_len = None
                #[32 bytes pickled msg length (padded)][pickled msg]
                while True:

                    if self.watchdog_flag.is_set():
                        if block == True:
                            for i in range(0, nclients_alive):
                                client_sockets[i].close()
                        return None

                    data = client_sockets[-1].recv(1024)
                    
                    if data_len is None:
                        data_len = pickle.loads(data)
                        log.debug('recv %d bytes from client' % data_len)
                    
                    buf.append(data)    
                    tmp_buf = ''.join(buf)
                      
                    if len(tmp_buf) - 32 >= data_len:
                        break

                cpickle = ''.join(buf)
                data = cpickle[32:] #pickled msg

                cid, cmsg = pickle.loads(data)

                if cid is None:
                    cid = nclients_alive

                client_data[cid] = cmsg

                client_sockets_cid.append(cid)

                if block == False:
                    client_sockets[-1].sendall(messages[cid])
                    client_sockets[-1].close()

                nclients_alive += 1

            except (pickle.UnpicklingError, EOFError):

                if block == True:
                    for i in range(0, nclients_alive):
                        client_sockets[i].close()

                log.error('Error unpickling client data')

                return None

            except socket.timeout:
                if timeout == True:
                    return None
                #otherwise just ignore error

        if block == True:
            for i in range(0, len(client_sockets)):

                cid = client_sockets_cid[i]
                #send message that corresponds to client
                if cid is not None:
                    client_sockets[i].sendall(messages[cid])
                else:
                    client_sockets[i].sendall(messages[i])

                client_sockets[i].close()

        return client_data


    def shutdown_sync_server(self):
        self.sync_sock.close()

    def shutdown(self, exit_code=0):
        log.debug(('server: Shutdown %d\n' % exit_code))

        if self.sock is not None:
            self.sock.close()

        if self.watchdog_server is not None:
            self.watchdog_server.shutdown()
            #wait until both threads have exited
            self.watchdog_server_thread.join()
            self.watchdog_check_thread.join()
       
    
#client     
def watchdog_check_client(hostname, wport, flag, cid, timeout):
    last_connect_time = time.time()
    #use a shorter timeout so we connect to the server more than needed
    effective_timeout = timeout/4
    while True:
        try:
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.settimeout(effective_timeout)
            sock.connect((hostname, wport))

            msg = (cid, 'client here')
            msg_pickle = pickle.dumps(msg, pickle.HIGHEST_PROTOCOL)   

            sock.sendall(msg_pickle)
            data = sock.recv(1024)

            assert(data == 'server here')

            last_connect_time = time.time()

            sock.close()
        except (socket.timeout, socket.error):
            pass       

        if time.time() - last_connect_time > timeout:
            flag.set()
            log.debug('watchdog flag set '+repr(flag.is_set())+' \n')   
            return

        #if shutdown is called somewhere else
        #exit this thread
        if flag.is_set():
            return

        time.sleep(effective_timeout)

def propagate_segment_thread(backend_driver, segment):
    backend_driver.propagate_segments(segment)
            
class TCPWEWorker(WESimClient):
    def __init__(self):
        WESimClient.__init__(self)
        self.hostname = None
        self.dport = None
        self.sport = None
        self.wport = None
        self.work_pickle = None
        self.segments = None
        self.shutdown_flag = False
        self.cid = None
        self.watchdog_flag = threading.Event()
        self.watchdog_check_thread = None
        self.timeout = None
        # The driver that actually propagates segments
        self.backend_driver = None

    def set_hostname(self, hostname):
        self.hostname = hostname
             
    def runtime_init(self, runtime_config, load_sim_config = True):
        super(TCPWEWorker, self).runtime_init(runtime_config, load_sim_config)
        self.timeout = runtime_config.get_int('server.timeout')
        self.dport = runtime_config.get_int('server.dport')
        self.sport = runtime_config.get_int('server.sport')
        self.wport = runtime_config.get_int('server.wport') 

    def init_watchdog(self):
        self.watchdog_check_thread = threading.Thread(target=watchdog_check_client, args=(self.hostname, self.wport, self.watchdog_flag, self.cid, self.timeout))
        self.watchdog_check_thread.daemon = True
        self.watchdog_check_thread.start()

    def run(self):
        if self.backend_driver is None:
            self.load_backend_driver()
            
        log.debug('entering receive loop')

        if self.get_cid() == False:
            self.shutdown(EX_ERROR)
            sys.exit(EX_ERROR)     

        self.init_watchdog()
                
        while True:
            if self.get_segments() == False:
                self.shutdown(EX_ERROR)
                sys.exit(EX_ERROR)

            if self.propogate_segments() == False:
                self.shutdown(EX_ERROR)
                sys.exit(EX_ERROR)
            
            if self.send_segments() == False:
                self.shutdown(EX_ERROR)
                sys.exit(EX_ERROR)

            if self.sync_clients() == False:
                self.shutdown(0)
                sys.exit(0)

            if self.shutdown_flag == True:
                self.shutdown(EX_ERROR)
                sys.exit(EX_ERROR)
                        
    def get_cid(self):
        '''gets the client id and sends the # of cores available'''
        try:
            #get cid

            data = self.send_to_server((None, 'cid'))
            
            if data is None:
                log.error('problem connecting to server-> exit\n')                
                return False
                        
            self.cid, tmp = pickle.loads(data)

            assert(tmp == 'cid')
   
        except (pickle.UnpicklingError, EOFError):
            log.error('problem connecting to server-> exit\n')
            return False

        try:
            #send number of cores

            data = self.send_to_server((self.cid, multiprocessing.cpu_count()))
    
            if data is None:
                log.error('problem connecting to server-> exit\n')                
                return False
    
            cid, tmp = pickle.loads(data)

            assert(tmp == 'cores')
            assert(self.cid == cid)

        except (pickle.UnpicklingError, EOFError):
            log.error('problem connecting to server-> exit\n')
            return False

        return True              

    def get_segments(self):

        self.work_pickle = None
        
        work_pickle_str = self.send_to_server((self.cid, 'ready for work'))

        if work_pickle_str is None:
            return False

        self.work_pickle = work_pickle_str

        return True    

    def propogate_segments(self):

        self.rc = None

        if self.work_pickle is None:
            log.error('Work pickle is not defined')
            self.rc = None
            return False

        try:
            cid, segs = pickle.loads(self.work_pickle)
            assert(self.cid == cid)
        except (pickle.UnpicklingError, EOFError):
            log.error('Error unpickling work')
            log.error(repr(self.work_pickle))
            self.segments = None
            return False
         
        nsegs = len(segs)
        
        #if there aren't any segs for this client
        #return [] (== segs)
        if nsegs == 0:
            self.segments = segs
            return True
        
        ncpus = multiprocessing.cpu_count()
        returncodes = []
        
        for cur_seg in xrange(0, nsegs, ncpus):

            seg_threads = []

            if cur_seg + ncpus > nsegs:
                seg_list = segs[cur_seg:]
            else:
                seg_list = segs[cur_seg:cur_seg + ncpus]
                
            for seg in seg_list:
                p = self.do_propagate_segment([seg])
                seg_threads.append(p)

            nrunning_tasks = len(seg_list)                
            complete_threads = []

            shutdown_flag = False
            while True:

                #check to see if threads are finished
                for i in [x for x in range(0, nrunning_tasks) if x not in complete_threads]:
                    if seg_threads[i].is_alive() == False:
                        complete_threads.append(i)

                if len(complete_threads) == nrunning_tasks:
                    break

                if self.watchdog_flag.is_set():
                    shutdown_flag = True
 
                    #wait on all incomplete threads
                    for i in [x for x in range(0, nrunning_tasks) if x not in complete_threads]:
                        if seg_threads[i].is_alive() == True:
                            log.debug('Join '+repr(i)+'\n')
                            seg_threads[i].join()
                    break

                time.sleep(5)
                
            if shutdown_flag == True:
                break
            
        self.segments = segs

        if shutdown_flag == True:
            return False
                    
        for seg in segs:
            if seg.status == Segment.SEG_STATUS_FAILED:
                return False
            
        return True

    def do_propagate_segment(self, segment):
        thread = threading.Thread(target=propagate_segment_thread, args=(self.backend_driver, segment,))
        thread.start()        
        return thread

    def send_segments(self):

        try:
            data_pickle = self.send_to_server((self.cid, self.segments))          
            
            if data_pickle is None:
                return False

            cid, data = pickle.loads(data_pickle)
            
            if data != 'prepare for next commnand':
                return False

        except (pickle.UnpicklingError, EOFError):
            log.error('Error unpickling rc data')
            log.error(repr(data_pickle))
            return False

        return True

    def sync_clients(self):

        shutdown_flag = False
        while True:

            if self.watchdog_flag.is_set():
                shutdown_flag = True
                break

            try:
                sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                #do not set too low, otherwise it will never connect with a large # of clients
                sock.settimeout(self.timeout)
                sock.connect((self.hostname, self.sport))

                msg = (self.cid, 'sync')
                msg_pickle = pickle.dumps(msg, pickle.HIGHEST_PROTOCOL)

                sock.sendall(msg_pickle)

                data = sock.recv(1024)

                if data == 'shutdown':
                    shutdown_flag = True
                    break
                #sync means all clients have checked in and are ready to continue
                elif data == 'sync':
                    break
                else:                
                    shutdown_flag = True
                    break                  

            except (socket.timeout, socket.error):
                pass
            
            finally:
                sock.close()

            time.sleep(self.timeout)

        if shutdown_flag == True:
            return False

        return True
    
    def send_to_server(self, msg):
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        msg_pickle = pickle.dumps(msg, pickle.HIGHEST_PROTOCOL)
        len_pickle = pickle.dumps(len(msg_pickle), pickle.HIGHEST_PROTOCOL)
        
        #pad len_pickle
        if len(len_pickle) < 32:
            len_pickle = len_pickle + ' ' * ( 32 - len(len_pickle) )
        else:
            return None
        
        complete_pickle = len_pickle + msg_pickle
        
        try:
            sock.connect((self.hostname, self.dport))
            
            sock.sendall(complete_pickle)
        except socket.error, e:
            log.error('problem connecting to server %r-> exit\n' % e)
            sock.close()
            return None

        buf = []

        while True:

            if self.watchdog_flag.is_set():
                sock.close()
                return None

            try:
                data = sock.recv(1024)
                #stop when server closes connection
                if not data:
                    break
  
                buf.append(data)
            except socket.error, e:
                log.error('problem connecting to server %r-> exit\n' % e)
                sock.close()
                return None

        buf_str = ''.join(buf)
        sock.close()  
              
        return buf_str                   

            
    def shutdown(self, exit_code=0):
        self.watchdog_flag.set()   
