#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --job-name="md"

PDB=$1
FF=$2
H2O=$3
# assume the desired NAMEOFPROT is the prefix of the PDB file
protname=`echo ${PDB} | sed "s#\.# #g" | awk '{ print $1}'`
echo $protname
echo $protname > protname.txt
module load openmpi/4.0.4
module load gromacs/2020.4

# prepare box
gmx_mpi pdb2gmx -f $PDB -o out.gro -p sys.top -ignh < inp
gmx_mpi editconf -f out.gro -bt cubic -d 0.9 -o boxed.gro
# solvate
gmx_mpi solvate -cp boxed.gro -p sys.top -o solvated.pdb
# add ions
touch ions.mdp
gmx_mpi grompp -c solvated.pdb -r solvated.pdb -f ions.mdp -p sys.top -o ions.tpr -maxwarn 5
echo "SOL" | gmx_mpi genion -s ions.tpr -p sys.top -o ${protname}.gro -neutral -conc 0.1
# energy minimization
gmx_mpi grompp -v -f em.mdp -c ${protname}.gro -o em_${protname}.tpr -p sys.top
#mpirun -n $SLURM_NTASKS
gmx_mpi mdrun -v -s em_${protname}.tpr -o em_${protname}.trr -c after_em_${protname}.gro -g em_${protname}.log -ntomp 1

