#!/bin/bash

#SBATCH --job-name=SF6
#SBATCH --output="%j.o"
#SBATCH --error="%j.e"
#SBATCH --account="thonhauserGrp"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --partition="small"
#SBATCH --mem=120Gb

ulimit -s unlimited

module load openmpi/3.1.0-intel-2018  

#mpirun /deac/thonhauserGrp/chakrad/quantum-espresso-6.2.1/qe-6.2.1-openmpi3.0-intel2018/bin/pw.x < SF6.in > SF6.out
#mpirun /deac/thonhauserGrp/chakrad/quantum-espresso-6.2.1/qe-6.2.1-openmpi3.0-intel2018/bin/pp.x < SF6-pp.in > SF6-pp.out
python C6_calc.py SF6-vdw-df2.cube > SF6-C6.txt

echo "DONE"

