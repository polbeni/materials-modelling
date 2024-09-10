#!/bin/bash
#SBATCH --account=
#SBATCH --qos=
#SBATCH --job-name=
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --ntasks=112

module load fftw/3.3.10
module load hdf5/1.14.1-2
module load quantumespresso/7.2

time srun   pw.x -input scf.in > scf.out
time srun   ph.x -input ph.in > ph.out
