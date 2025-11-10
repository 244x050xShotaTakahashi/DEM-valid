#!/bin/bash
#============ LSF Options ============
#QSUB -ug gr10171
#QSUB -q gr10171a
#QSUB -rn
#QSUB -W 48:00
#QSUB -A p=1:t=1:c=1:m=1920M

#============ Shel Script ============

cd $LS_SUBCWD

aprun -n $LSB_PROCS -d $LSB_CPUS -N $LSB_PPN ./dmotion9 160 160 800 2.0e-8 1.0 15 60 80 80 25 30 10000 100 100000 1.956951 7.372429E+05 0 10 1.0e-3