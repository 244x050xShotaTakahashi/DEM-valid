#!/bin/bash
#SBATCH -p gr10451a
#SBATCH --rsc p=32:t=1:c=1
#SBATCH -t 120:00:00
#SBATCH -o stdout.%J.log
#SBATCH -e stderr.%J.log

# set -x

module load intel intelmpi
module load hdf5/1.12.2_intel-2023.2-impi
module load fftw

export EMSES_DEBUG=no

date

rm *_0000.h5
srun ./mpiemses3D plasma.inp

mypython generate_xdmf3.py nd1p00_0000.h5 nd2p00_0000.h5 nd3p00_0000.h5 phisp00_0000.h5
mypython generate_xdmf3.py rho00_0000.h5 rhobk00_0000.h5
mypython generate_xdmf3.py ex00_0000.h5 ey00_0000.h5 ez00_0000.h5
mypython generate_xdmf3.py j1x00_0000.h5 j1y00_0000.h5 j1z00_0000.h5
mypython generate_xdmf3.py j2x00_0000.h5 j2y00_0000.h5 j2z00_0000.h5
mypython generate_xdmf3.py j3x00_0000.h5 j3y00_0000.h5 j3z00_0000.h5

date

# Postprocessing(visualization code, etc.)
mypython plot.py ./
