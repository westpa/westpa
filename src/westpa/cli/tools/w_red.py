from h5py import File as H5File
import numpy as np
from westpa import rc
from westpa.tools import WESTParallelTool


class DurationCorrector(object):
    @staticmethod
    def from_kinetics_file(directh5, istate, fstate, dtau, n_iters=None):
        iter_slice = slice(n_iters)

        if isinstance(directh5, H5File):
            dataset = directh5['durations'][iter_slice]
        else:
            with H5File(directh5, 'r') as directh5:
                dataset = directh5['durations'][iter_slice]

        torf = np.logical_and(dataset['istate'] == istate, dataset['fstate'] == fstate)
        torf = np.logical_and(torf, dataset['weight'] > 0)

        durations = dataset['duration']
        weights = dataset['weight']

        weights[~torf] = 0.0  # mask off irrelevant flux

        return DurationCorrector(durations, weights, dtau)

    def __init__(self, durations, weights, dtau, maxduration=None):
        self.weights = np.array(weights)
        self.durations = np.array(durations)
        self.dtau = dtau
        self._f_tilde = None
        self._f_int1 = None

        if maxduration is None:
            self.maxduration = self.durations.shape[0]
        else:
            self.maxduration = maxduration

        if dtau is None:
            all_durations = []
            all_durations.extend(durations)
            all_durations.extend(np.arange(maxduration))
            uniq_durations = np.unique(all_durations)  # unique sorts automatically
            self.dtau = np.min(np.diff(uniq_durations))

        self._build_map()

    @property
    def event_duration_histogram(self):
        return self._f_tilde

    @property
    def cumulative_event_duration_histogram(self):
        return self._f_int1

    def _build_map(self):
        weights = self.weights
        durations = self.durations
        maxduration = self.maxduration
        dtau = self.dtau

        taugrid = np.arange(0, maxduration, dtau, dtype=float)
        f_map = np.zeros(weights.shape, dtype=int) - 1
        for i, tau in enumerate(taugrid):
            matches = np.logical_and(durations >= tau, durations < tau + dtau)
            f_map[matches] = i

        self.taugrid = taugrid
        self.f_map = f_map

    def correction(self, iters, freqs=None):
        r"""
        Return the correction factor

        __                              __  -1
        |  t=theta  tau=t                |
        |    |\      |\                  |
        |    |       | ~                 |
        |    |       | f(tau) dtau dt    |      * maxduration
        |   \|      \|                   |
        |   t=0    tau=0                 |
        |_                              _|

        where
            ~`                        ^
            f(tau) is proportional to f(tau)/(theta-tau), and is normalized to
                            ^
        integrate to 1, and f(tau) is sum of the weights of walkers with
        duration time tau.

        ---------
        Arguments
        ---------
        maxduration: the maximum duration time that could have been observed in
            the simulation, which is usually equal to the length of the
            simulation. This should be in units of tau.
        """

        if iters is None:
            iters = np.arange(len(self.weights))

        if freqs is None:
            freqs = np.ones(len(iters), dtype=float)

        maxduration = np.max(iters) + 1

        f_map = self.f_map[iters]
        weights = self.weights[iters]
        taugrid = self.taugrid  # [self.taugrid < maxduration]

        weights *= freqs[:, None]

        dtau = self.dtau

        f_tilde = np.zeros(len(taugrid), dtype=float)
        for i, tau in enumerate(taugrid):
            if tau < maxduration:
                f_tilde[i] = weights[f_map == i].sum() / (maxduration - tau + 1)

        if f_tilde.sum() != 0:
            f_tilde /= f_tilde.sum() * dtau

        self._f_tilde = f_tilde
        # now integrate f_tilde twice
        # integral1[t/dtau] gives the integral of f_tilde(tau) dtau from 0 to t
        self._f_int1 = integral1 = np.zeros(f_tilde.shape)

        for i, tau in enumerate(taugrid):
            if i > 0 and tau < maxduration:
                integral1[i] = np.trapz(f_tilde[: i + 1], taugrid[: i + 1])

        integral2 = np.trapz(integral1, taugrid)

        if integral2 == 0:
            return 0.0
        return maxduration / integral2


def get_raw_rates(directh5, istate, fstate, n_iters=None):
    rate_evol = directh5['rate_evolution'][slice(n_iters), istate, fstate]
    avg = rate_evol['expected']

    return avg


def calc_avg_rate(directh5_path, istate, fstate, **kwargs):
    """
    Return the raw or RED-corrected rate constant with the confidence interval.

    ---------
    Arguments
    ---------
    dt: timestep (ps)
    nstiter: duration of each iteration (number of steps)
    ntpr: report inteval (number of steps)

    """

    n_iters = kwargs.pop("n_iters", None)

    ntpr = kwargs.pop("report_interval", 20)
    nstiter = kwargs.pop("n_steps_iter", 1000)
    callback = kwargs.pop("callback", None)

    red = kwargs.pop("red", False)

    if len(kwargs) > 0:
        raise ValueError("unparsed kwargs")

    dtau = float(ntpr) / nstiter
    dc = None

    with H5File(directh5_path, 'r') as directh5:
        if n_iters is None:
            n_iters = directh5['rate_evolution'].shape[0]

        rate_evol = directh5['rate_evolution'][n_iters - 1, istate, fstate]
        rate = rate_evol['expected']

        if red:
            dc = DurationCorrector.from_kinetics_file(directh5, istate, fstate, dtau, n_iters)

    if callback is not None:
        kw = {"correction": dc}
        callback(**kw)

    iters = np.arange(n_iters)

    correction = dc.correction(iters) if dc else 1.0

    rate *= correction

    return rate


def calc_rates(directh5_path, istate, fstate, **kwargs):
    """
    Return the raw and RED-corrected rate constants vs. iterations.
    This code is faster than calling calc_rate() iteratively

    ---------
    Arguments
    ---------
    dt: timestep (ps)
    nstiter: duration of each iteration (number of steps)
    ntpr: report inteval (number of steps)

    """

    n_iters = kwargs.pop("n_iters", None)

    ntpr = kwargs.pop("report_interval", 20)
    nstiter = kwargs.pop("n_steps_iter", 1000)
    callback = kwargs.pop("callback", None)

    red = kwargs.pop("red", False)

    if len(kwargs) > 0:
        raise ValueError("unparsed kwargs")

    dtau = float(ntpr) / nstiter
    dc = None

    with H5File(directh5_path, 'r') as directh5:
        rate_evol, cilb, ciub = get_raw_rates(directh5, istate, fstate, n_iters)
        if n_iters is None:
            n_iters = len(rate_evol)
        if red:
            dc = DurationCorrector.from_kinetics_file(directh5, istate, fstate, dtau, n_iters)

    if callback is not None:
        kw = {"correction": dc}
        callback(**kw)

    raw_rates = np.zeros(n_iters)

    rates = np.zeros(n_iters)

    for i in range(n_iters):
        i_iter = i + 1
        print("\riter %d/%d (%3.0f%%)" % (i_iter, n_iters, i_iter * 100.0 / n_iters), end="")

        r = rate_evol[i]

        iters = np.arange(i_iter)

        correction = dc.correction(iters) if dc else 1.0

        raw_rates[i] = r
        rates[i] = raw_rates[i] * correction

    print("\n")

    return rates


class RateCalculator:
    def __init__(self, directh5, istate, fstate, assignh5=None, **kwargs):
        n_iters = kwargs.pop("n_iters", None)
        ntpr = kwargs.pop("report_interval", 20)
        nstiter = kwargs.pop("n_steps_iter", 1000)

        if len(kwargs) > 0:
            for k in kwargs:
                print(k)
            raise ValueError("unparsed kwargs")

        dtau = float(ntpr) / nstiter

        with H5File(directh5, 'r') as f:
            state_labels = {}
            for i, raw_label in enumerate(f['state_labels']):
                label = raw_label.decode() if isinstance(raw_label, bytes) else raw_label
                state_labels[label] = i
            if istate not in state_labels:
                raise ValueError(f"istate not found: {istate}, available options are {list(state_labels.keys())}")
            if fstate not in state_labels:
                raise ValueError(f"istate not found: {fstate}, available options are {list(state_labels.keys())}")
            istate = state_labels[istate]
            fstate = state_labels[fstate]
            cond_fluxes = f['conditional_fluxes'][slice(n_iters), istate, fstate]

        if assignh5 is not None:
            with H5File(assignh5, 'r') as f:
                pops = f['labeled_populations'][slice(n_iters)]
                pops = pops.sum(axis=2)
        else:
            pops = None

        self._dc = None
        self._pops = pops
        self._cond_fluxes = cond_fluxes
        self._dtau = dtau
        self._directh5 = directh5
        self._assignh5 = assignh5
        self._istate = istate
        self._fstate = fstate

    @property
    def conditional_fluxes(self):
        return self._cond_fluxes

    @property
    def populations(self):
        return self._pops

    @property
    def tau(self):
        return self._tau

    @property
    def dtau(self):
        return self._dtau

    @property
    def istate(self):
        return self._istate

    @property
    def fstate(self):
        return self._fstate

    @property
    def n_iters(self):
        return len(self.conditional_fluxes)

    def _get_corrector(self):
        if self._dc is None:
            with H5File(self._directh5, 'r') as f:
                self._dc = DurationCorrector.from_kinetics_file(f, self.istate, self.fstate, self.dtau, self.n_iters)

        return self._dc

    def calc_rate(self, i_iter=None, red=False, **kwargs):

        if i_iter is None:
            i_iter = self.n_iters

        dc = self._get_corrector() if red else None
        found = False
        with H5File(self._directh5, 'r') as f:
            for i in range(f['rate_evolution'].shape[0]):
                rate_evol = f['rate_evolution'][i, self.istate, self.fstate]
                start = rate_evol['iter_start']
                stop = rate_evol['iter_stop']

                if i_iter >= start and i_iter < stop:
                    rate = rate_evol['expected']
                    found = True
                    break

            if not found:
                self.log.error("Can't find rate evolution data for iteration %d!" % i_iter)

        if dc:
            iters = np.arange(i_iter)
            correction = dc.correction(iters)
            rate *= correction

        return rate

    def calc_rates(self, n_iters=None, **kwargs):
        if n_iters is None:
            n_iters = self.n_iters

        rates = np.zeros(n_iters)

        for i in range(n_iters):
            i_iter = i + 1
            print("\riter %d/%d (%3.0f%%)" % (i_iter, n_iters, i_iter * 100.0 / n_iters), end="")

            r = self.calc_rate(i_iter, **kwargs)

            rates[i] = r

        print("\n")

        return rates


class WRed(WESTParallelTool):
    prog = 'w_red'
    description = '''\
Apply the RED scheme to estimate steady-state WE fluxes from
shorter trajectories.

-----------------------------------------------------------------------------
Source data
-----------------------------------------------------------------------------

Source data is provided as a w_ipa "scheme" which is typically defined
in the west.cfg file. For instance, if a user wishes to estimate RED
fluxes for a scheme named "DEFAULT" that argument would be provided
to w_red and WRed would estimate RED fluxes based off of the data
contained in the assign.h5 and direct.h5 files in ANALYSIS/DEFAULT.

'''

    def __init__(self):
        super().__init__()

    def go(self):
        try:
            rc.config['west']['analysis']['red']
        except Exception:
            raise ValueError('No RED parameters are specified in west.cfg.')
        try:
            rc.config['west']['analysis']['red']['scheme']
        except Exception:
            raise ValueError('No scheme specified for RED calculation in west.cfg.')
        try:
            rc.config['west']['analysis']['red']['istate_label']
        except Exception:
            raise ValueError('No intial state label specified for RED calculation in west.cfg.')
        try:
            rc.config['west']['analysis']['red']['fstate_label']
        except Exception:
            raise ValueError('No final state label specified for RED calculation in west.cfg.')
        try:
            rc.config['west']['analysis']['red']['nstiter']
        except Exception:
            raise ValueError('Time step not specified in west.cfg.')
        try:
            rc.config['west']['analysis']['red']['nstrep']
        except Exception:
            raise ValueError('Time step not specified in west.cfg.')

        if rc.config['west']['analysis']['kinetics']['evolution'] == "cumulative":
            pass
        else:
            print("Only RED estimates with cumulative averaging are supported at this time.")
            exit()

        config = rc.config
        adir = config.get(['west', 'analysis', 'directory'])
        name = config.get(['west', 'analysis', 'red', 'scheme'])
        istate = config.get(['west', 'analysis', 'red', 'istate_label'])
        fstate = config.get(['west', 'analysis', 'red', 'fstate_label'])
        n_steps_per_iter = config.get(['west', 'analysis', 'red', 'nstiter'])
        n_steps_per_report = config.get(['west', 'analysis', 'red', 'nstrep'])

        directh5path = '%s/%s/direct.h5' % (adir, name)
        assignh5path = '%s/%s/assign.h5' % (adir, name)

        print('\nConfig successfully read from west.cfg:')
        print('--------------------------------------')
        print('scheme: %s' % name)
        print('states: %s -> %s' % (istate, fstate))
        print('nstiter: %s' % n_steps_per_iter)
        print('nstrep: %s' % n_steps_per_report)
        print('--------------------------------------')
        print('\nEstimating RED fluxes...')
        rater = RateCalculator(
            directh5path,
            istate,
            fstate,
            n_steps_iter=n_steps_per_iter,
            report_interval=n_steps_per_report,
            assignh5=assignh5path,
        )

        rates = rater.calc_rates(red=True, callback=None)

        with H5File(directh5path, "r+") as dest_file:
            try:
                dest_file.create_dataset('red_flux_evolution', data=rates)
                print('saved RED fluxes to red_flux_evolution in ANALYSIS/%s/direct.h5' % name)
            except Exception:
                warning = input('Dataset already exists! Overwrite? (y/n)')
                if warning == "y":
                    dest_file['red_flux_evolution'][...] = rates
                    print('saved RED fluxes to red_flux_evolution in ANALYSIS/%s/direct.h5' % name)
                elif warning == "n":
                    np.save('ANALYSIS/%s/red.npy' % name, rates)
                    print('saved RED fluxes to red_flux_evolution.npy in ANALYSIS/%s' % name)
                else:
                    print('red rates not saved. Exiting...')
                    exit


def entry_point():
    WRed().main()


if __name__ == '__main__':
    entry_point()
