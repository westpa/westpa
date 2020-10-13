import logging

from westpa.tools import WESTParallelTool, WESTDataReader, IterRangeSelection, ProgressIndicatorComponent
import westpa
from westpa.core.extloader import get_object

log = logging.getLogger('w_crawl')


class WESTPACrawler:
    '''Base class for general crawling execution. This class
    only exists on the master.'''

    def initialize(self, iter_start, iter_stop):
        '''Initialize this crawling process.'''
        pass

    def finalize(self):
        '''Finalize this crawling process.'''
        pass

    def process_iter_result(self, n_iter, result):
        '''Process the result of a per-iteration task.'''
        pass


def _remote_task(n_iter, taskfn):
    data_manager = westpa.rc.get_data_manager()  # gaahhh...globals
    data_manager.open_backing(mode='r')
    return n_iter, taskfn(n_iter, data_manager.get_iter_group(n_iter))


class WCrawl(WESTParallelTool):
    prog = 'w_crawl'
    description = '''\
Crawl a weighted ensemble dataset, executing a function for each iteration.
This can be used for postprocessing of trajectories, cleanup of datasets,
or anything else that can be expressed as "do X for iteration N, then do
something with the result". Tasks are parallelized by iteration, and
no guarantees are made about evaluation order.


-----------------------------------------------------------------------------
Command-line options
-----------------------------------------------------------------------------

'''

    def __init__(self):
        super().__init__()

        # These are used throughout
        self.progress = ProgressIndicatorComponent()
        self.data_reader = WESTDataReader()
        self.iter_range = IterRangeSelection(self.data_reader)

        self.crawler = None
        self.task_callable = None

    def add_args(self, parser):
        self.data_reader.add_args(parser)
        self.iter_range.add_args(parser)

        tgroup = parser.add_argument_group('task options')
        tgroup.add_argument(
            '-c',
            '--crawler-instance',
            help='''Use CRAWLER_INSTANCE (specified as module.instance) as an instance of
                            WESTPACrawler to coordinate the calculation. Required only if initialization,
                            finalization, or task result processing is required.''',
        )
        tgroup.add_argument(
            'task_callable',
            help='''Run TASK_CALLABLE (specified as module.function) on each iteration.
                            Required.''',
        )
        self.progress.add_args(parser)

    def process_args(self, args):
        self.progress.process_args(args)
        self.data_reader.process_args(args)
        with self.data_reader:
            self.iter_range.process_args(args)

        self.task_callable = get_object(args.task_callable, path=['.'])
        if args.crawler_instance is not None:
            self.crawler = get_object(args.crawler_instance, path=['.'])
        else:
            self.crawler = WESTPACrawler()

    def go(self):
        iter_start = self.iter_range.iter_start
        iter_stop = self.iter_range.iter_stop
        iter_count = iter_stop - iter_start
        self.data_reader.open('r')
        pi = self.progress.indicator
        with pi:
            pi.operation = 'Initializing'
            self.crawler.initialize(iter_start, iter_stop)

            try:
                pi.new_operation('Dispatching tasks & processing results', iter_count)
                task_gen = ((_remote_task, (n_iter, self.task_callable), {}) for n_iter in range(iter_start, iter_stop))
                for future in self.work_manager.submit_as_completed(task_gen, self.max_queue_len):
                    n_iter, result = future.get_result(discard=True)
                    if self.crawler is not None:
                        self.crawler.process_iter_result(n_iter, result)
                    pi.progress += 1
            finally:
                pi.new_operation('Finalizing')
                self.crawler.finalize()


def entry_point():
    WCrawl().main()


if __name__ == '__main__':
    entry_point()
