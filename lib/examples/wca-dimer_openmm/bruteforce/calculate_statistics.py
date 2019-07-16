import numpy as np

a_bound = 4.2
b_bound = 7.0

def summary_stats(time_in_a, time_in_b, n_a2b, n_b2a):
    total_time = time_in_a + time_in_b

    frac_a = (1.0 * time_in_a) / total_time
    frac_b = (1.0 * time_in_b) / total_time

    flux_a2b = (1.0 * n_a2b) / total_time
    flux_b2a = (1.0 * n_b2a) / total_time

    try:
        mfpt_a2b = frac_a / flux_a2b
    except:
        mfpt_a2b = np.nan

    try:
        mfpt_b2a = frac_b / flux_b2a
    except:
        mfpt_b2a = np.nan

    return mfpt_a2b, mfpt_b2a


def calc_stats(d):
    niters = d.shape[0]
    
    time_in_a = 0
    time_in_b = 0
    n_a2b = 0
    n_b2a = 0

    tt_a2b = []
    tt_b2a = []

    if d[0] < a_bound:
        time_in_a += 1
        curr_state = 0
    else:
        time_in_b += 1
        curr_state = 1

    for k in range(1, niters):
        if d[k] < a_bound:
            state = 0
        elif d[k] > b_bound:
            state = 1
        else:
            state = curr_state

        if state == 0:
            time_in_a += 1
        else:
            time_in_b += 1

        if state != curr_state:
            if state == 0:
                n_b2a += 1
                tt_b2a.append(k)
            else:
                n_a2b += 1
                tt_a2b.append(k)

        curr_state = state

    return n_a2b, n_b2a, time_in_a, time_in_b, tt_a2b, tt_b2a


if __name__ == '__main__':
    d = np.load('data/md-solvent-langevin-distance.npy')
    n_a2b, n_b2a, time_in_a, time_in_b, tt_a2b, tt_b2a = calc_stats(d)

    #print tt_a2b
    #print '--------'
    #print tt_b2a

    print('ntranstions a->b: ', n_a2b)
    print('ntranstions b->a: ', n_b2a)

    print('time in a: ', time_in_a)
    print('time in b: ', time_in_b)

    total_time = time_in_a + time_in_b

    frac_a = (1.0 * time_in_a) / total_time
    frac_b = (1.0 * time_in_b) / total_time
    print('Frac a: ', frac_a)
    print('Frac b: ', frac_b)

    flux_a2b = (1.0 * n_a2b) / total_time
    flux_b2a = (1.0 * n_b2a) / total_time

    mfpt_a2b = frac_a / flux_a2b
    mfpt_b2a = frac_b / flux_b2a

    print('MFPT a->b: ', mfpt_a2b)
    print('MFPT b->a: ', mfpt_b2a)
