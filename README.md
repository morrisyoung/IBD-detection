IBD-detection
=============

Introduction
========
This is the IBD detection algorithm and implementations. One version for multi-threads, one version for multi-processes (OpenMPI). Matlab script (and its result with cutoff=50,000) is also attached for your sanity checks. Also, if you are in Columbia, you can use the run\_mp.sh to submit your OpenMPI jobs in C2B2 clusters.

Compile
========
* multi-threads:

gcc -o IBDReport\_mt IBDReport\_mt.c -lpthread -lm

* multi-processes:

1. install the OpenMPI into your computer;

2. compile:

PC: mpicc -o IBDReport\_mp IBDReport\_mp.c -lm

C2B2 cluster: /opt/OFED/current/mpi/gcc/openmpi-1.4.2/bin/mpicc -o IBDReport\_mp IBDReport\_mp.c

Run
========
* multi-threads:

./IBDReport -f TREE\_FILE -m CUTOFF -d DISCRETIZATION -l LENGTH\_OF\_CHROMOSOME -t NUMBER\_OF\_THREADS

-f: file name, in present working directory

-m: cutoff value; defining the effective IBD segment length

-d: discretization value; processing one tree at least every d trees

-l: length of chromosome; we need enter this manually because we can't get the total length of chromosome in Nexus tree

-t: number of working threads for IBD detection; as there is a main thread and a IBD report thread, so you'd better use [your computer's core number minus 2] as the number of working threads

notes:

1. this version is not fully optimized for memory usage (access), so the running speed may be a little slow

2. this version doesn't support Newick format trees and reading trees from STDIN presently, which are all supported in the next multi-processes version

3. the result is saved locally as "result\_yourtreefilename"


* multi-processes:

mpirun -n 16 IBDReport\_mp -f test\_1000G\_50I\_100Million.trees -F 1 -t 1 -e 0.01 -m 50000 -d 0 -l 100000000

or

"other program generating stdout" | mpirun -n 1 IBDReport\_mp -F 2 -t 2 -e 0.000001 -m 0 -d 0 -l 100000000

-f: file name, in present working directory

-F: format of trees; 1 -> Nexus; 2-> Newick

-t: input type; 1 -> tree file; 2 -> stdin

-m: cutoff value

-e: epsilon value; within which (+-epsilon) we still regard the tMRCA as the same

-d: discretization value; for speeding up

-l: length of chromosomes; needed for Nexus trees

notes:

1. this version is fully optimized both in memory usage and in parallelism, which should be faster; but as we use temporary tree files and temporary result files, so extra consumption (running time) is obvious

2. for the Newick trees, as presently I only try the trees from MaCS's stdout, so actually "-F 2" and " -t 2" are bound together; so if you want try Newick trees, please use MaCS and read them from its stdin; later on the "-F" will be fixed as a separate parameter, and also for "-t"

3. as I just finished reading newick trees from MaCS's stdout, no more sanity checks for it have been performed by now, although it should be right logically; so reading nexus trees from outside files is more recommanded presently

4. as the MPI program is invoked by the "mpirun" excutable, if you need n working processing, please invoke n+1 processes; but this is not for reading trees from stdin, in which you should invoke only 1 process, because there is only serially processing mode for trees read from stdin

5. the result is saved locally as "result\_yourtreefilename"

Contact
========
If you have any problems during using this program, please contact me: sdmorrisys@gmail.com
