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

#### Serial run (baseline):
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
```
patrick@sami:~/workspace/westpa/lib/examples/odld$ mpirun -np 2 ./run.sh --work-manager=mpi

...

Wed Sep 30 16:58:58 2015
Iteration 1000 (1000 requested)
Beginning iteration 1000
940 segments remain in iteration 1000 (940 total)
94 of 100 (94.000000%) active bins are populated
per-bin minimum non-zero probability:       2.34428e-05
per-bin maximum probability:                0.0314202
per-bin probability dynamic range (kT):     7.20064
per-segment minimum non-zero probability:   5.10842e-07
per-segment maximum non-zero probability:   0.00433192
per-segment probability dynamic range (kT): 9.04546
norm = 1, error in norm = 0 (0*epsilon)
Iteration wallclock: 0:00:00.281643, cputime: 0:00:00


Wed Sep 30 16:58:59 2015
WEST run complete.
```

#### MPI run ( 8 procs )
```
patrick@sami:~/workspace/westpa/lib/examples/odld$ mpirun -np 8 ./run.sh --work-manager=mpi

...

Wed Sep 30 17:05:41 2015
Iteration 1000 (1000 requested)
Beginning iteration 1000
940 segments remain in iteration 1000 (940 total)
94 of 100 (94.000000%) active bins are populated
per-bin minimum non-zero probability:       4.47964e-05
per-bin maximum probability:                0.0315613
per-bin probability dynamic range (kT):     6.55756
per-segment minimum non-zero probability:   2.3645e-06
per-segment maximum non-zero probability:   0.00456024
per-segment probability dynamic range (kT): 7.56457
norm = 1, error in norm = 0 (0*epsilon)
Iteration wallclock: 0:00:00.330630, cputime: 0:00:00


Wed Sep 30 17:05:41 2015
WEST run complete.
```

#### MPI run ( 40 procs - over two physically different nodes )
```
php8@login0b:~/westpa/lib/examples/odld$ srun --nodes=2 --tasks-per-node=20 ./run.sh --work-manager=mpi

...

Wed Sep 30 17:21:36 2015
Iteration 1000 (1000 requested)
Beginning iteration 1000
940 segments remain in iteration 1000 (940 total)
94 of 100 (94.000000%) active bins are populated
per-bin minimum non-zero probability:       1.48948e-05
per-bin maximum probability:                0.0272635
per-bin probability dynamic range (kT):     7.5123
per-segment minimum non-zero probability:   1.23687e-06
per-segment maximum non-zero probability:   0.00401684
per-segment probability dynamic range (kT): 8.08567
norm = 1, error in norm = 0 (0*epsilon)
Iteration wallclock: 0:00:00.425205, cputime: 0:00:00


Wed Sep 30 17:21:36 2015
WEST run complete.
```
