#!/bin/bash
#SBATCH --account=
#SBATCH --qos=
#SBATCH --job-name=
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --ntasks=112

module load fftw/3.3.10
module load hdf5/1.14.1-2
module load quantumespresso/6.5

time srun   pw.x -input scf.in > scf.out
time srun   ph.x -input ph.in > ph.out
time srun   q2r.x -input q2r.in > q2r.out
time srun   matdyn.x -input matdyn.in > matdyn.out
