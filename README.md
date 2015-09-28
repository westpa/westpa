# wwmgr
Parallel work manager for WESTPA

## MPIWorkManager

### Test Suite
```
patrick@sami:~/workspace/westpa/lib/wwmgr$ nosetests test_work_managers.test_mpi
........
----------------------------------------------------------------------
Ran 8 tests in 2.752s

OK
```

### Validation
All validation to date was done against the odld case.

#### Serial run:
```
patrick@sami:~/workspace/westpa/lib/examples/odld$ ./run.sh --work-manager=serial

...

Mon Sep 28 17:05:32 2015
Iteration 1000 (1000 requested)
Beginning iteration 1000
940 segments remain in iteration 1000 (940 total)
94 of 100 (94.000000%) active bins are populated
per-bin minimum non-zero probability:       1.98167e-05
per-bin maximum probability:                0.0281443
per-bin probability dynamic range (kT):     7.25857
per-segment minimum non-zero probability:   1.37656e-06
per-segment maximum non-zero probability:   0.00396506
per-segment probability dynamic range (kT): 7.96569
norm = 1, error in norm = 0 (0*epsilon)
Iteration wallclock: 0:00:00.187738, cputime: 0:00:00


Mon Sep 28 17:05:33 2015
WEST run complete.
```

#### MPI run ( single proc )
```
patrick@sami:~/workspace/westpa/lib/examples/odld$ mpirun -np 1 ./run.sh --work-manager=mpi

...

Mon Sep 28 17:09:52 2015
Iteration 1000 (1000 requested)
Beginning iteration 1000
940 segments remain in iteration 1000 (940 total)
94 of 100 (94.000000%) active bins are populated
per-bin minimum non-zero probability:       2.69966e-05
per-bin maximum probability:                0.0381878
per-bin probability dynamic range (kT):     7.25456
per-segment minimum non-zero probability:   2.44167e-06
per-segment maximum non-zero probability:   0.00479355
per-segment probability dynamic range (kT): 7.58234
norm = 1, error in norm = 0 (0*epsilon)
Iteration wallclock: 0:00:00.185674, cputime: 0:00:00


Mon Sep 28 17:09:52 2015
WEST run complete.
```

#### MPI run ( 2 procs )
Currently seq faults after some time

#### MPI run ( 8 procs )
Currently seq faults after some time
