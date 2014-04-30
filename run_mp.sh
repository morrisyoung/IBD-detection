#!/bin/bash
#$ -cwd -S /bin/bash -N IBDReport_mp -j y
#$ -pe orte 2 -l mem=4G,time=4::
#$ -M sy2515@columbia.edu -m bes
LD_LIBRARY_PATH=/opt/OFED/current/lib64:$LD_LIBRARY_PATH
/opt/OFED/current/mpi/gcc/openmpi-1.4.2/bin/mpirun -n 2 IBDReport_mp -f test_10000N_250C_50M.trees_20 -l 50000000 -m 50000 -e 0.01
