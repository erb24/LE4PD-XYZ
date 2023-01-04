#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -A uoo104
#SBATCH -p compute
#SBATCH --job-name="eq"

BD=`cat BD.txt`
protname=`cat ${BD}/em/protname.txt`
GRO="${BD}/em/after_em_${protname}.gro"

module load openmpi/4.0.4
module load gromacs/2020.4

# npt equilibration
gmx_mpi grompp -v -f npt.mdp -c ${GRO} -o npt.tpr -p ${BD}/em/sys.top -r ${GRO} -maxwarn 5
mpirun -n $SLURM_NTASKS gmx_mpi mdrun -v -s npt.tpr -deffnm npt -ntomp 1

# production run
cd ${BD}
mkdir -v prod
cp -v ${BD}/prod.sh ${BD}/pro.mdp ${BD}/em/*.itp ${BD}/em/protname.txt ${BD}/BD.txt prod/
cd prod
rm -rfv runno
sbatch prod.sh

