#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -A uoo104
#SBATCH -p compute
#SBATCH --job-name="md"

protname=`cat protname.txt`
BD=`cat BD.txt`
module load openmpi/4.0.4
module load gromacs/2020.4

echo $SLURM_NTASKS
# production run
MAX_TIME=24
cycles=10
NAME="pro" #${protname}
PDB="npt.gro"
############################################################################
# Run
############################################################################
if [ -e runno ] ; then
   #########################################################################
   # Restart runs
   #########################################################################
   nn=`tail -n 1 runno | awk '{print $1}'`
   mpirun -n $SLURM_NTASKS gmx_mpi mdrun -maxh ${MAX_TIME} -s ${NAME}.tpr -cpi ${NAME}.cpt -deffnm ${NAME} -ntomp 1 -append
   #########################################################################
else
   #########################################################################
   # First run
   #########################################################################
   nn=1
   gmx_mpi grompp -v -f pro.mdp -c ${BD}/eq/${PDB} -r ${BD}/eq/${PDB} -p ${BD}/em/sys.top -o ${NAME}.tpr
   mpirun -n $SLURM_NTASKS gmx_mpi mdrun -maxh ${MAX_TIME} -s ${NAME}.tpr -deffnm ${NAME} -ntomp 1
   #########################################################################
fi
############################################################################


############################################################################
# Check number of cycles
############################################################################
mm=$((nn+1))
echo ${mm} > runno
#cheking number of cycles
if [ ${nn} -ge ${cycles} ]; then
  exit
fi
############################################################################

############################################################################
# Resubmitting again
############################################################################
sbatch < run.sh

