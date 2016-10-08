# Copyright (C) 2013 Matthew C. Zwier and Lillian T. Chong
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

from __future__ import division, print_function

import logging, argparse, traceback, os
log = logging.getLogger('errors.py')
import westpa
import sys
sys.tracebacklimit=5

class WESTErrorReporting:
    """
    General class for handling more user-friendly error messages.
    Programs should import this and create an instance of it; then,
    instead of raising specific errors through the Python interface, they
    should just call "report error" and format the arguments accordingly.
    This way, the users see far more specific information about the actual
    error.

    """

    def __init__(self, cp=''):
        self.pstatus = westpa.rc.pstatus
        self.bash_variables = None
        self.script = None
        log.debug('initializing error handling')
        self.config = westpa.rc.config
        self.system = westpa.rc.get_system_driver()
        try:
            self.report_all_errors = self.config['west']['report_all_errors']
        except:
            self.report_all_errors = False
        self.reported_errors = {}
        # Calling program
        self.cp = cp

        # westpa.rc.pstatus('interrupted; shutting down')

        wiki = "https://chong.chem.pitt.edu/wewiki/WESTPA_Error_Handling"

        executable = os.path.expandvars(self.config['west']['executable']['propagator']['executable'])
        rcfile = self.config['args']['rcfile']
        logfile = self.config['west']['executable']['propagator']['stdout']
        pcoord_len = self.system.pcoord_len
        pcoord_ndim = self.system.pcoord_ndim
        self.llinebreak = "----------------------------------------------------------------"
        self.linebreak = "-------------------------------------------"
        self.format_kwargs = { 'executable': executable, 'rcfile': rcfile, 'pcoord_ndim': pcoord_ndim, 'pcoord_len': pcoord_len, 'logfile': logfile,
                               'wiki': wiki, 'linebreak': self.linebreak, 'cp': cp, 'llinebreak': self.llinebreak }

        self.SEG_ERROR            = """
        {llinebreak}{linebreak}
        ERROR # {id} ON Iteration: {segment.n_iter}, Segment: {segment.seg_id}"""

        self.ITER_ERROR = """
        {llinebreak}{linebreak}
        ERROR # {id} ON Iteration: {iteration}"""

        self.RUNSEG_GENERAL_ERROR = { 'msg': """
        A general error has been caught from the {executable} propagator.
        You should check the indicated log file for more specific errors,
        or see below.

        FILES TO CHECK

        {logfile}
        {executable}


        LAST 10 LINES OF STDERR
        {linebreak}
        {err}
        {linebreak}
        """,
        'id': 0 }

        self.RUNSEG_SHAPE_ERROR = { 'msg': """
        The shape of your progress coordinate return value is {shape},
        which is different from what is specified in your {rcfile}: ({pcoord_len}, {pcoord_ndim}).  

        FILES TO CHECK

        {executable}
        {rcfile}
        Your dynamics engine configuration file.

        They should all agree on number of values/timepoints/progress coordinate values
        that you are returning.

        See {logfile}
        """,
        'id': 1 }

        self.RUNSEG_TMP_ERROR = { 'msg': """
        Could not read the auxdata {dataset} return value from {filename} for segment {segment.seg_id} in iteration {segment.n_iter}.

        POSSIBLE REASONS

        {executable} is not returning anything into the {dataset} return.
            - This could be the result of a failed run, or
            - {executable} is not returning data into the {dataset} return
              (typically, cat/paste into the WEST_DATASET_RETURN variable)
        {filename} is not writable.
            - The space {filename} exists on could be full.  Try cleaning it.

        FILES TO CHECK

        {logfile}
        {executable}
        {rcfile} - did you want to return an auxiliary dataset?

        Specific exception:

        {linebreak}
        {e}
        {linebreak}
        """,
        'id': 2 }

        self.RUNSEG_AUX_ERROR = { 'msg': """
        Your auxiliary data return is empty.  This typically
        means that your {executable} propagator has failed.  Please
        check the indicated log file for more specific
        errors:

        FILES TO CHECK
        {linebreak}
        {logfile}

        LAST 10 LINES OF STDERR
        {linebreak}
        {err}
        {linebreak}

        Also, has anyone ever seen me?
        """,
        'id': 3 }

        self.RUNSEG_PROP_ERROR = { 'msg': """
        Propagation has failed for {failed_segments} segments:
        {linebreak}
        {failed_ids}
        {linebreak}

        Check the corresponding log files for each ID.
        """,
        'id': 4 }

        self.WRUN_INTERRUPTED = { 'msg': """
        INTERRUPTION

        An interruption has been sent to {cp}.
        This has either been done manually (such as the break command or the killing of a queue script),
        or by the local sysadmin.
        """,
        'id': 5 }

        self.RUNSEG_EMPTY_VARIABLES = { 'msg': """
        NOTICE

        Empty variables exist in your {executable}.  This could be a problem.

        {linebreak}
        {empties}
        {linebreak}

        """,
        'id': 99 }

        self.REPORT_ONCE = """
        NOTICE

        The configuration has been set such that each error type is caught only once; all other
        segments which report the same error will have their output suppressed.  This can be disabled.
        """

        self.SEE_WIKI = """
        Check the wiki for more information
        {wiki}

        {llinebreak}{linebreak}
        """ 


    def report_segment_error(self, error, segment, **kwargs):
        # We'll want to pass in the segment object, actually.  But we can't call that from here...
        # ... but, this should still work, for the moment.
        self.format_kwargs.update(kwargs)
        # Pull in the ID.
        self.format_kwargs.update(error)
        self.format_kwargs.update({'segment': segment})
        # Testing for istate/bstates.
        try:
            test = segment.n_iter
        except:
            segment.n_iter = 0
            segment.seg_id = 0
        try:
            # These ones need modifying, as the segment.niter doesn't work with format, directly.
            self.format_kwargs['logfile'] = os.path.expandvars(self.format_kwargs['logfile'].format(segment=segment))
        except:
            pass

        # How can we enable it such that we report one 'type' of error only once?
        # Often, we repeat many errors and it's a pain.  Sometimes, this is useful information,
        # but most of the time it's just indicative of a general problem.
        # In the typical python fashion, we ask forgiveness, not permission.
        if self.report_all_errors == False:
            try:
                if self.reported_errors[self.REPORT_ONCE] == False:
                    self.pstatus(self.REPORT_ONCE.format(**self.format_kwargs))
                    self.reported_errors[self.REPORT_ONCE] = True
            except:
                self.pstatus(self.REPORT_ONCE.format(**self.format_kwargs))
                self.reported_errors[self.REPORT_ONCE] = True

        try:
            if self.reported_errors[error['msg']] == False:
                self.pstatus(self.SEG_ERROR.format(**self.format_kwargs))
                self.pstatus(error['msg'].format(**self.format_kwargs))
                self.pstatus(self.SEE_WIKI.format(**self.format_kwargs))
                if self.report_all_errors == False:
                    self.reported_errors[error['msg']] = True
        except:
            self.pstatus(self.SEG_ERROR.format(**self.format_kwargs))
            self.pstatus(error['msg'].format(**self.format_kwargs))
            self.pstatus(self.SEE_WIKI.format(**self.format_kwargs))
            if self.report_all_errors == False:
                self.reported_errors[error['msg']] = True

    def report_error(self, error, **kwargs):
        self.format_kwargs.update(kwargs)
        # Pull in the ID.
        self.format_kwargs.update(error)
        self.pstatus(self.ITER_ERROR.format(**self.format_kwargs))
        self.pstatus(error['msg'].format(**self.format_kwargs))
        self.pstatus(self.SEE_WIKI.format(**self.format_kwargs))

    def scan_bash_variables(self, script):
        if self.bash_variables == None or self.script != script:
            # Let's not run this more than once, shall we?
            self.bash_variables = []
            import re
            with open(script, 'r') as runseg:
                for line in runseg:
                    # Find the variables!
                    var_group = re.findall('\$\w+', line)
                    if len(var_group) > 0:
                        self.bash_variables += var_group
            self.bash_variables = [self.remove_variable_symbol(s) for s in self.bash_variables]

    def scan_bash_empty_variables(self, out):
        import re
        empties = []
        statedv = []
        # Let's place all these in the empty variables...
        # ... but also sort the others into the 'good' variable
        # section.
        for line in out.splitlines():
            empty = re.findall('\w+=$', line)
            # This should match everything BUT a carriage return, newline, or space.
            # Empties just sort of works, I think.
            filled = re.findall('^\w+=(?!\&$| $|\n$).*', line)
            empties += empty
            statedv += filled
        return [self.remove_equals_symbol(s) for s in empties], statedv

    def remove_variable_symbol(self, s):
        return s[1:]

    def remove_equals_symbol(self, s):
        return s[:-1]

    def does_not_exist_in_list(self, l1, l2):
        # We want to see if any string in l1 is a part of string
        # in l2.
        # We'll return all values that are in the second list.
        rl = []
        for s1 in l1:
            exists = False
            for s2 in l2:
                if s1 in s2:
                    exists = True
            if exists == False:
                rl.append(s1)
        return rl

    def report_general_error_once(self, error, **kwargs):
        # This is a function that respects the 'run only once' setting,
        # but doesn't require extensive iteration.  It's useful for printing a
        # warning a during simulation.

        self.format_kwargs.update(kwargs)
        # Pull in the ID.
        self.format_kwargs.update(error)
        # How can we enable it such that we report one 'type' of error only once?
        # Often, we repeat many errors and it's a pain.  Sometimes, this is useful information,
        # but most of the time it's just indicative of a general problem.
        # In the typical python fashion, we ask forgiveness, not permission.
        if self.report_all_errors == False:
            try:
                if self.reported_errors[self.REPORT_ONCE] == False:
                    self.pstatus(self.REPORT_ONCE.format(**self.format_kwargs))
                    self.reported_errors[self.REPORT_ONCE] = True
            except:
                self.pstatus(self.REPORT_ONCE.format(**self.format_kwargs))
                self.reported_errors[self.REPORT_ONCE] = True

        try:
            if self.reported_errors[error['msg']] == False:
                self.pstatus(error['msg'].format(**self.format_kwargs))
                if self.report_all_errors == False:
                    self.reported_errors[error['msg']] = True
        except:
            self.pstatus(error['msg'].format(**self.format_kwargs))
            if self.report_all_errors == False:
                self.reported_errors[error['msg']] = True
    def raise_exception(self):
        raise Exception('Error reported from {}'.format(self.cp))
