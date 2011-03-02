from __future__ import division 
__metaclass__ = type

import os, stat, sys, socket, threading, time, multiprocessing
import argparse
import cPickle as pickle

from wemd.work_managers import WEMDWorkManager

from wemd.types import Segment

import logging
log = logging.getLogger(__name__)

class TCPWorkManager(WEMDWorkManager):
    def __init__(self, sim_manager):
        super(TCPWorkManager,self).__init__(sim_manager)
        self.worker = None
        
    def read_key_file(self, keyfilename):
        '''Read communications "secret key" from keyfile, which must be between 16 and 512 bytes
        in length (inclusive) and must not be world-readable.'''
        kfstats = os.stat(keyfilename) 
        if (kfstats.st_mode & stat.S_IROTH):
            raise RuntimeError('key file (%r) must not be world-readable' % keyfilename)
        elif not 16 <= kfstats.st_size <= 512:
            raise RuntimeError('key file (%r) must be between 16 and 512 bytes')
        keyfile = file(os.path.expanduser(os.path.expandvars(keyfilename)), 'rb')
        key = keyfile.read()
        return key
        
    def parse_aux_args(self, aux_args):
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument('--help-workmanager', '--help-tcpip', dest='do_help', action='store_true',
                            help='show this help message and exit')
        parser.add_argument('--debug-tcp', dest='do_tcp_debug', action='store_true',
                            help='enable TCP driver-specific debug messages')
        parser.add_argument('-s', '--servername', dest='servername', 
                            help='Hostname of server to connect to (clients) or the hostname to serve on (server) '
                                +'(default: %(default)s)',
                            default=socket.getfqdn())
        parser.add_argument('-k', '--keyfile', dest='keyfile',
                            help='Read communication authentication key from KEYFILE '
                                +'(default: %(default)s)',
                            default='~/.wemd_key')
        parser.add_argument('-d', '--sport', '--dport', dest='dport', type=int,
                            help='Server port (default: %(default)d)',
                            default=5231)
        subparsers = parser.add_subparsers()
        
        client_parser = subparsers.add_parser('client')        
        client_parser.add_argument('-c', '--cport', dest='cport', type=int,
                                   help='Client-side port. Zero chooses an unused port. (Default: 0)',
                                   default=0)
        client_parser.set_defaults(mode='client')
        server_parser = subparsers.add_parser('server')
        server_parser.add_argument('-n', '--nclients', dest='nclients', type=int,
                                   help='Number of clients (hosts) to expect', required=True)
        server_parser.set_defaults(mode='server')

        if '--help-workmanager' in aux_args:
            parser.print_help()
            sys.exit(0)
        args, extra_args = parser.parse_known_args(aux_args)
    
        secret_key = self.read_key_file(args.keyfile)
        
        if args.mode == 'client':
            self.worker = TCPWorkerClient(self.sim_manager, args.servername, secret_key, args.dport, args.cport, 
                                          args.do_tcp_debug)
        else:
            self.worker = TCPWorkerServer(self.sim_manager, args.servername, secret_key, args.dport, args.nclients, 
                                          args.do_tcp_debug)
            
        return extra_args

    def prepare_workers(self):
        return self.worker.prepare_workers()
    
    def prepare_iteration(self, n_iter, segments):
        if self.is_server():
            super(TCPWorkManager,self).prepare_iteration(n_iter, segments)
            self.worker.prepare_iteration(n_iter, segments)    

    def propagate(self, segments):
        if self.is_server():
            self.worker.propagate(segments)
    
    def finalize_iteration(self, n_iter, segments):
        if self.is_server():
            super(TCPWorkManager,self).finalize_iteration(n_iter, segments)
            self.worker.finalize_iteration(n_iter, segments)
    
    def shutdown(self, exit_code=0):
        return self.worker.shutdown(exit_code)
    
    def is_server(self):
        return self.worker.is_server()
    
    def run(self):
        if self.is_server() == False:
            self.worker.run()

class TCPWorkerBase():
    def __init__(self, sim_manager, hostname, secret_key, dport, enable_debug=False):
        log.debug('initializing tcpip work manager')
        self.set_server_hostname(hostname)
        self.secret_key = None
        self.set_secret_key(secret_key)
        
        self.runtime_config = sim_manager.runtime_config
        self.sim_manager = sim_manager
        
        self.timeout = self.runtime_config.get_int('server.timeout')
        self.dport = dport 
        self.sock = None
        
        self.enable_debug = enable_debug

    def run(self):
        pass
                                            
    def set_server_hostname(self, hostname):
        self.hostname = socket.gethostbyname(hostname)
    
    def get_local_hostname(self):
        return socket.gethostbyname(socket.gethostname())
    
    def set_secret_key(self, secret_key):
        self.secret_key = secret_key

    def init_dserver(self, port, nqueue = 1, use_timeout = True):
        self.debug('init_dserver')
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        #So server can be restarted quickly
        self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        if use_timeout:
            self.sock.settimeout(self.timeout)
        self.sock.bind((self.get_local_hostname(), port))
        self.sock.listen(nqueue)
        return self.sock.getsockname()[1] #return the port used
                
    def shutdown_dserver(self):
        if self.sock is not None:
            self.sock.close()        
            self.sock = None
        
    def send_message(self, message, sock):
        '''adds secret key and length to message
           returns [128 bytes secret key][32 bytes pickled msg length (padded)][pickled msg]
           or None on error '''
        #self.debug('send_message' % repr(message))
        
        msg_pickle = pickle.dumps(message, pickle.HIGHEST_PROTOCOL)
        len_pickle = pickle.dumps(len(msg_pickle), pickle.HIGHEST_PROTOCOL)
        
        #pad len_pickle
        if len(len_pickle) < 32:
            len_pickle = len_pickle + ' ' * ( 32 - len(len_pickle) )
        else:
            return False #message too long --> unlikely unless message is very long
        
        assert(self.secret_key is not None)
        assert(len(self.secret_key) >= 64)
        
        data = self.secret_key + len_pickle + msg_pickle
   
        sock.sendall(data)
        
        return True
    
    def send_message_thread(self, message, sock, error):
        '''calls send_message, but sets error on error'''
        if self.send_message(message, sock) == False:
            error.set()
               
    def recv_message(self, sock):
        '''verifies key and recvs on sock until message complete
           returns message or None on error
        '''
        #self.debug('recv_message')
        
        #[128 bytes secret key][32 bytes pickled msg length (padded)][pickled msg]
        key_check = True
        data_len = None
        buf = []
        key_len = len(self.secret_key)
        
        while True:            
            data = sock.recv(1024) 
            
            buf.append(data)    
            tmp_buf = ''.join(buf)
            
            if key_check == False:
                #need to receive at least secret_key_len+32 bytes
                if len(tmp_buf) < (key_len+32):
                    continue
                            
                if tmp_buf[0:key_len] != self.secret_key[0:key_len]:
                    log.error('Invalid Key: %r %r %r' %(repr(sock),repr(data)))
                    return None
                
                key_check = True

            if data_len is None:
                self.debug('pickle.loads data_len start')
                data_len = pickle.loads(data[key_len:key_len+32])
                self.debug('pickle.loads data_len end')
                self.debug('recv %d bytes from client' % data_len)

            if (len(tmp_buf) - (key_len+32)) >= data_len:
                break                
        
        pickled_data = ''.join(buf)

        try:
            self.debug('pickle.loads data start')
            data = pickle.loads(pickled_data[key_len+32:])
            self.debug('pickle.loads data end')
        except (pickle.UnpicklingError, EOFError):
            raise ValueError('problem decoding data pickle')  
                
        return data

    def debug(self, string):
        if self.enable_debug is True and log.isEnabledFor(logging.DEBUG):
            log.debug('DEBUG: Base: ' + string + ' ' + repr(time.time()))
                                                       
class TCPWorkerServer(TCPWorkerBase):
    def __init__(self, sim_manager, hostname, secret_key, dport, nclients, enable_debug=False):
        super(TCPWorkerServer,self).__init__(sim_manager, hostname, secret_key, dport)
        self.nclients = None
        #client info
        #client_info[cid]['key']
        #key = cores -> # cores
        #key = ip -> ip addr
        #key = port -> client port
        self.client_info = None
        self.stop_ping_event = None
        self.ping_thread = None
        self.nclients = nclients
        self.enable_debug=enable_debug
        
    def is_server(self):
        return True    
        
    def prepare_workers(self):
        #only for initial connection (to get client info)
        self.init_dserver(self.dport, nqueue = self.nclients)
        
        if self.send_cid() == False:
            self.shutdown(os.EX_NOHOST)
            sys.exit(os.EX_NOHOST)
            
        self.shutdown_dserver()

    def prepare_iteration(self, n_iter, segments):
        pass

    def debug(self, string):
        if self.enable_debug is True and log.isEnabledFor(logging.DEBUG):
            log.debug('DEBUG: Server: ' + string + ' ' + repr(time.time()))

    def send_cid(self):
        '''assigns an id to each client, and gets the number of cores available
           returns True on success, False on failure'''
        self.debug('send_cid')

        msgs = []
        for i in range(0, self.nclients):
            msgs.append((i,'cid'))

        self.client_info = self.recv_from_clients(msgs)

        if self.client_info is None:
            self.debug('self.client_info is None')
            return False
        
        #cmsg -> cores
        for info in self.client_info:
            info.update({'cores':info['cmsg']})
            info.pop('cmsg')
        
        self.debug('send_cid client_info:'+ repr(self.client_info))
        return True

    def recv_from_clients(self, messages):
        '''Send message to clients -> They must connect to server
        Returns (ip, port, client_id, client 'message') as dictionary
        '''       
        self.debug('recv_from_clients') 
        sock = self.sock
        client_data = [0 for x in range(0, self.nclients)]
        
        nclients_alive = 0
        while nclients_alive < self.nclients:
            new_socket, address = sock.accept()
            self.debug('connected to %r' % (repr(address)))
            data = self.recv_message(new_socket)
                  
            if data is None:
                new_socket.close()
                return None
            
            cid, cmsg, port = data
            self.debug('cid, cmsg, port: %r, %r, %r' %(cid, cmsg, port))

            #assign cid on first come basis
            if cid is None:
                cid = nclients_alive

            client_data[cid] = {'ip':address[0], 'port':port, 'cid':cid, 'cmsg':cmsg}
            self.debug('recv_from_clients: %r' % client_data[cid])

            if self.send_message(messages[cid], new_socket) == False:
                new_socket.close()                
                return None

            new_socket.close()

            nclients_alive += 1
                                
        return client_data

    
    def send_to_clients(self, messages):
        '''Like recv_from_clients but sends then recvs()
           Also initiates connection
           Assumes that self.client_info is known'''
        self.debug('send_to_clients %r ' % messages)

        info = self.client_info
        
        cdata = [0 for x in range(0, self.nclients)]
        
        sock_threads = []
        sock_error = threading.Event()
        socks = []
        for cid in xrange(0, self.nclients):
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)            
            sock.connect((info[cid]['ip'], info[cid]['port']))

            socks.append(sock)
            
            sock_threads.append(threading.Thread(target=self.send_message_thread, args=(messages[cid], sock, sock_error)))
            sock_threads[-1].start()
        
        self.debug('send_to_clients sending data to clients (threads)')

        complete_threads = []
        while True: 
            for cid in xrange(0, self.nclients):
                if cid not in complete_threads:
                    if sock_threads[cid].is_alive() == False:
                        complete_threads.append(cid)
                  
            if len(complete_threads) == len(sock_threads):
                break
        self.debug('send_to_clients finished sending data')
        
        if sock_error.is_set():
            for sock in socks:
                sock.close()
            log.error('TCPServer: socket error in send message')
            return None
        
        for sock in socks:
            data = self.recv_message(sock)
            cid = data[0]
            cdata[cid] = data
            
        for sock in socks:
            self.debug('send_to_clients closing socket %r' % cid)
            sock.close()
        
        self.debug('send_to_clients recv %r' % repr(cdata))
        return cdata
    
    def send_command_to_clients(self, command):
        '''Like send_to_clients but sends the same message to each one
        '''
        self.debug('send_command_to_clients %r' % command)
        msgs = []
        for i in range(0, self.nclients):
            msgs.append((i, command))
        
        return self.send_to_clients(msgs)
                
    def get_clients_status(self):
        self.debug('get_clients_status')
        return self.send_command_to_clients('status')
    
    def shutdown_clients(self):
        self.debug('shutdown clients')
        return self.send_command_to_clients('shutdown')
    
    def send_to_client_by_id(self, message, cid):
        #self.debug('send_to_client_by_id %r %r' % (message, cid))

        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        info = self.client_info
        
        sock.connect((info[cid]['ip'], info[cid]['port']))
        
        if self.send_message(message, sock) == False:
            sock.close()
            return None
        
        data = self.recv_message(sock)
        
        sock.close()
        
        #self.debug('send_to_client_by_id:%r' % repr(data))
        return data
    
    def propagate(self, segments):
        self.debug('propagate')
        self.stop_ping_thread()
        
        status = self.get_clients_status()
        state = self.reduce_client_status(status)
            
        if state != 'ready':
            log.error('Clients not ready for next iter, current state %r' % state)
            raise ValueError('Clients not ready')
        
        if self.send_segments(segments) == False:
            raise ValueError('Segments could not be sent')
        
        
        segs_by_id = {segment.seg_id: segment for segment in segments}
        new_segments = self.get_segments()
            
        if new_segments is None:
            log.error('No segments returned')
            raise ValueError('No segments returned')    

        self.start_ping_thread()
                
        #segments[cid][task#] -> segments[i]
        propagated_segments = [j for k in new_segments for j in k if j is not None]
        
        # Double check that we got back what we sent out
        assert (set(segment.seg_id for segment in propagated_segments if segment is not None) 
                == set(segment.seg_id for segment in segments))
        
        # Update in place
        for propagated_seg in propagated_segments:
            orig_segment = segs_by_id[propagated_seg.seg_id]
            orig_segment.status = propagated_seg.status
            orig_segment.walltime = propagated_seg.walltime
            orig_segment.cputime  = propagated_seg.cputime
            orig_segment.pcoord[...] = propagated_seg.pcoord[...]
        

    def finalize_iteration(self, n_iter, segments):
        pass
    
    def send_segments(self, segments):
        self.debug('send_segments')

        cmds = []
        
        #for now, just assume a homogeneous resource
        i = 0
        tasks = [[] for i in xrange(0, self.nclients)]
        
        for segment in map(None, *(iter(segments),) * self.client_info[0]['cores']):

            if type(segment) is tuple:
                tasks[i].append(segment)                
            else:    
                tasks[i].append((segment,))
                        
            i += 1
            i %= self.nclients
            
        for i in xrange(0, self.nclients):
            #flatten list of tuples into a list
            task = [j for k in tasks[i] for j in k if j is not None]
            cmds.append((i, 'propagate_segs', task))
            
        ret = self.send_to_clients(cmds)

        if ret is None:
            return False

        return True
    
    def reduce_client_status(self, status):
        '''reduces list of status to one status
           All Clients Ready -> ready
           Any Clients Busy -> busy
           All Clients Have Data -> data
           Otherwise (some ready, some data) -> busy
        '''
        self.debug('reduce_client_status: %r' % repr(status))

        if status is None:
            return 'error'
        
        if None in status:
            return 'error'
        
        data = False
        ready = False
        reduced_status = None
        for cid in xrange(0, self.nclients):
            cstatus = status[cid][1]
                
            if cstatus == 'busy':
                return 'busy'
            elif cstatus == 'data':
                if ready == False:
                    data = True
                else:
                    return 'busy'
            elif cstatus == 'ready':
                if data == False:
                    ready = True
                else:
                    return 'busy'
            else:
                raise ValueError('Unknown status %r' % cstatus)
            
        if data == True:
            return 'data'
        elif ready == True:
            return 'ready'
        else:
            raise ValueError('Inconsistent status')
    
    def start_ping_thread(self, cid = None):
        self.debug('PING start_ping_thread %r %r' % (self.ping_thread,self.stop_ping_event))
        
        if self.ping_thread is None:
            self.debug('start_ping_thread')
            self.stop_ping_event = threading.Event()
            self.ping_thread = threading.Thread(target=self.ping_clients, args=(self.stop_ping_event, cid))
            self.ping_thread.start()
    
    def stop_ping_thread(self):
        self.debug('PING stop_ping_thread %r %r' % (self.ping_thread,self.stop_ping_event))
        if self.ping_thread is not None:
            self.debug('stop_ping_thread')
            self.stop_ping_event.set()

            while True: #wait until ping_thread has stopped
                if self.ping_thread.is_alive() == False:
                    break   
                 
                time.sleep(0.01)
                
            self.ping_thread = None
            self.stop_ping_event = None                        
            
    def ping_clients(self, stop, cid = None):
        '''ping all clients except cid if set'''
        while True:
            if stop.is_set():
                self.debug('ending ping clients')
                break
            self.debug('pinging clients')

            for i in xrange(0, self.nclients):
                if cid is not None:
                    if i == cid:
                        self.debug('skipping ping %r'%cid)
                        continue
                self.debug('ping %r' % i)

                self.send_to_client_by_id((i,'status'), i)
            
            time.sleep(0.01)
    
    def get_segments(self):
        self.debug('get_segments')

        getsegs_msgs = []
        for i in xrange(0, self.nclients):
            getsegs_msgs.append((i, 'get_segs'))

        segs = [] #NOTE: these are not sorted by cid
        cids = [] #recv'd segs from these clients

        while len(segs) < self.nclients:
            status = self.get_clients_status()
    
            if self.reduce_client_status(status) == 'ready':
                self.debug('all clients ready, exiting loop')
                assert(len(segs) == self.nclients)
                break
            
            #occurs if client has died (connection refused)
            if status is None:
                log.error('Error connecting to client')
                return None
            
            for cid in xrange(0, self.nclients):
                self.debug('get_segments: status[%r] %r' % (cid,repr(status[cid])))
                if status[cid][1] == 'data':
                    self.debug('get_segments: sending %r to %r' % (getsegs_msgs[cid], cid))

                    self.start_ping_thread(cid)     

                    #this will cause the client to switch to ready
                    data_tuple = self.send_to_client_by_id(getsegs_msgs[cid], cid)
                    
                    self.stop_ping_thread()
                    
                    if data_tuple is None:
                        return None
                    
                    ccid, seg_data, seg = data_tuple
                    
                    assert(cid == ccid)
                    assert(seg_data == 'segs')
                    
                    if seg is None:
                        log.error('Error retrieving seg from client %r' % cid)
                        return None
                    
                    segs.append(seg)
                    cids.append(cid)
                    
            time.sleep(0.01)

        self.debug('segs retrieved %r %r %r'%(segs,cids,len(segs)))
        return segs

    def shutdown(self, exit_code=0):
        self.debug(('server: Shutdown %d\n' % exit_code))
        self.stop_ping_thread()
        self.shutdown_dserver()
        self.shutdown_clients()

#client
def propagate_segment_thread(propagator, segment):
    propagator.propagate(segment)


def propagate_particles(propagator, segments, shutdown_flag, error_flag):
    
    nsegs = len(segments)
   
    #if there aren't any segs for this client
    if nsegs == 0:
        return
   
    ncpus = multiprocessing.cpu_count()
   
    for cur_seg in xrange(0, nsegs, ncpus):

        seg_threads = []

        if cur_seg + ncpus > nsegs:
            seg_list = segments[cur_seg:]
        else:
            seg_list = segments[cur_seg:cur_seg + ncpus]
           
        for seg in seg_list:
            p = threading.Thread(target=propagate_segment_thread, args=(propagator, [seg],))
            p.start()        
            seg_threads.append(p)

        nrunning_tasks = len(seg_list)
                       
        complete_threads = []

        while True:
            #check to see if threads are finished
            for i in [x for x in range(0, nrunning_tasks) if x not in complete_threads]:
                if seg_threads[i].is_alive() == False:
                    complete_threads.append(i)

            if len(complete_threads) == nrunning_tasks:
                break

            if shutdown_flag.is_set():

                #wait on all incomplete threads
                for i in [x for x in range(0, nrunning_tasks) if x not in complete_threads]:
                    if seg_threads[i].is_alive() == True:
                        log.info('Join '+repr(i)+'\n')
                        seg_threads[i].join()
                break

            time.sleep(0.01)
            
        if shutdown_flag.is_set():
            return
       
        for seg in segments:
            if seg.status == Segment.SEG_STATUS_FAILED:
                error_flag.set()
                shutdown_flag.set()
                raise ValueError('Segment(s) Failed')
                
class TCPWorkerClient(TCPWorkerBase):
    def __init__(self, sim_manager, hostname, secret_key, dport, cport, enable_debug=False):
        super(TCPWorkerClient, self).__init__(sim_manager, hostname, secret_key, dport)
        self.work_pickle = None
        self.segments = None
        self.cid = None
        self.shutdown_flag = threading.Event()
        self.error_flag = threading.Event()
        self.state = 'busy'
        self.pp_thread = None
        self.cport = cport
        self.enable_debug=enable_debug
        
        try:
            self.max_attempt = self.runtime_config.get_int('server.max_attempt')
        except KeyError:
            self.max_attempt = 100                       
        
    def prepare_workers(self):
        if self.cport is None:
            self.cport = 0
        
        log.info('TCPClient runtime init')

        self.cport = self.init_dserver(self.cport, use_timeout = False) #ie if cport == 0, have the os choose which port to use
        self.get_cid()        
        
        #propagator isn't loaded until after init
        self.propagator = self.sim_manager.propagator        
   
    def debug(self, string):
        if self.enable_debug is True and log.isEnabledFor(logging.DEBUG):
            log.debug('DEBUG: Client: ' + string + ' ' + repr(time.time()))
      
    def get_command(self):
        self.debug('get_command')

        sock = self.sock
           
        new_socket, address = sock.accept()
        self.debug('get_command connected to %r' % repr(address))
        
        dec_data = self.recv_message(new_socket)

        if dec_data is None:
            new_socket.close()
            raise ValueError('Invalid command from server None')

        if len(dec_data) == 2:
            cid, cmsg = dec_data
            data = None
        elif len(dec_data) == 3:
            cid, cmsg, data = dec_data
        else:
            raise ValueError('Invalid command from server %r' % (repr(dec_data)))

        assert(cid == self.cid)
        
        return cmsg, data, new_socket

    def run(self):
        self.debug('run loop')

        self.state = 'ready'
        segs = []
        while True:
            if self.error_flag.is_set():
                log.error('segments have not completed - client exiting')
                raise ValueError('Incomplete Segments')
                        
            if self.shutdown_flag.is_set():
                return
            
            cmd, data, data_sock = self.get_command()
            self.debug('cmd %r' % (cmd))

            if cmd == 'propagate_segs':
                self.debug('propagate_segs cmd')
                self.state = 'busy'
                
                segs = data
                
                if segs:
                    self.propagate_particles(segs) #nonblocking
                else:
                    self.state = 'data'
                self.debug('send running segs')    
                self.send_message((self.cid, 'running segs'), data_sock)
                self.debug('done send running segs')

            elif cmd == 'status':
                self.debug('sending status')
                self.send_message((self.cid, self.state), data_sock)
            elif cmd == 'get_segs':
                if self.state == 'data':                    
                    self.send_message((self.cid, 'segs', segs), data_sock)
                    segs = []
                    self.state = 'ready'
                else:
                    self.send_message((self.cid, None), data_sock)
            
            elif cmd == 'shutdown' or cmd is None:
                self.send_message((self.cid, 'shutdown'), data_sock)
                data_sock.close()
                self.shutdown()
                return
            
            data_sock.close() #close connection after replying to server
    
            if self.propagate_particles_complete() == True and self.state == 'busy': #ie segs finished                                  
                self.state = 'data'

                
    def propagate_particles(self, segments):
        self.debug('propagate_particles')
        self.pp_thread = threading.Thread(target=propagate_particles, args=(self.propagator, segments, self.shutdown_flag, self.error_flag))
        self.pp_thread.start()    
    
    def propagate_particles_complete(self):
        if self.pp_thread is None:
            return False
        
        if self.pp_thread.is_alive():
            return False
        else:
            return True
                                    
    def is_server(self):
        return False
        
    def get_cid(self):
        '''gets the client id and sends the # of cores available'''
        self.debug('get_cid')
        tries = 0
        while tries < self.max_attempt:
            try:
                sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                sock.connect((self.hostname, self.dport))
                break
            except (socket.error,socket.timeout):
                tries += 1
                
        self.debug('connected to: %r' % repr((self.hostname, self.dport)))
        
        if self.send_message((None, multiprocessing.cpu_count(), self.cport), sock) == False:
            sock.close()
            raise ValueError('could not retrieve cid')

        self.debug('get_cid recv_message')
        cid, msg = self.recv_message(sock)
        sock.close()

        self.cid = cid
        assert(msg == 'cid')
        self.debug('hello from client %r' % self.cid)

    def shutdown(self, exit_code=0):
        self.debug('shutdown')
        self.shutdown_flag.set()   
        self.shutdown_dserver()
    