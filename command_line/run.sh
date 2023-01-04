#!/bin/bash -l
#SBATCH --job-name="velacc"
#SBATCH --output="./velacc.log"
#SBATCH -A uoo104
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --export=ALL
#SBATCH -t 2-0:00:00

# This file is a script for running GROMACS simulations of a protein in solution, starting from
# just the PDB structure in vacuum. Other initial conditions (e.g. non-canonical 
# amino acids, ligands, missing resiudes, etc.) will almost certainly require
# adjustment to both this script and the force field in use.

BD=$PWD
echo ${BD} > BD.txt
cd $BD

PDB="1UBQ.pdb"

# load modules required for running MPI-enabled GROMACS simulations;
# the modules below are what are used don Expanse as of December 2022 
# and may be to be adjusted depedning on 1) the hardware used and 2)
# the date of use.

module load openmpi/4.0.4
module load gromacs/2020.4

# run the energy minimization
mkdir -v em
cp -v em.mdp em.sh $PDB em/
# Select force field by assigning an integer to the constant $FF.
# Here is the default listing from GROMACS. N.B. for a force field 
# directory in the run foder, enter "1" for $FF. The default GROMACS listing is as follows,
# for GROMACS 2020.4 on Expanse as of December 2022:
# 1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)
# 2: AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)
# 3: AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29, 461-469, 1996)
# 4: AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21, 1049-1074, 2000)
# 5: AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006)
# 6: AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)
# 7: AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)
# 8: CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)
# 9: GROMOS96 43a1 force field
#10: GROMOS96 43a2 force field (improved alkane dihedrals)
#11: GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)
#12: GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)
#13: GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)
#14: GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856, DOI: 10.1007/s00249-011-0700-9)
#15: OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)

FF=1

# Select a water model. The options for GROMACS 2020.4 on Expanse as of December 2022 are as follows;
# selecting 7 will give a vacuum simulation:
#Select the Water Model:
# 1: TIP3P     TIP 3-point, recommended
# 2: TIP4P     TIP 4-point
# 3: TIP4P-Ew  TIP 4-point optimized with Ewald
# 4: TIP5P     TIP 5-point (see http://redmine.gromacs.org/issues/1348 for issues)
# 5: SPC       simple point charge
# 6: SPC/E     extended simple point charge
# 7: None

H2O=1

# run the energy minimization:
cd em/
echo "${FF}" > inp
echo "${H2O}" >> inp
sh em.sh ${PDB} ${FF} ${H2O}
rm -rfv inp
cd $BD

# run equilibration
mkdir -v eq
cp -v eq.sh npt.mdp nvt.mdp ${BD}/BD.txt eq/
cd eq/
sbatch eq.sh 

# production run
#cd ${BD}
#mkdir -v prod
#cp -v ${BD}/prod.sh ${BD}/pro.mdp ${BD}/em/*.itp ${BD}/em/protname.txt ${BD}/BD.txt prod/
#cd prod
#sbatch prod.sh

