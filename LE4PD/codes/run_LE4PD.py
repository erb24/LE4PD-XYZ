# coding: utf-8

import numpy as np
import sys
import subprocess
from LE4PD import *

try:
	PROTNAME = sys.argv[1]
	PDB = sys.argv[2]
	TOP = sys.argv[3]

	#subprocess.call('echo "3" | gmx_mpi trjconv -f '+str(PROTNAME)+'.xtc -s '+str(TOP)+' -o '+str(G96),shell=True)
	#subprocess.call('echo "3" | gmx_mpi trjconv -f '+str(PROTNAME)+'.xtc -s '+str(TOP)+' -o '+str(PDB),shell=True)
	#subprocess.call('echo "1" | gmx_mpi trjconv -f '+str(PROTNAME)+'.xtc -s '+str(TOP)+' -o '+str(PROTNAME)+'_first.pdb -dump 0',shell=True)
	#subprocess.call("sh unformat_traj.sh "+str(PDB),shell=True)

	TOP=str(PROTNAME)+"_first.pdb"
	#gen_protinfo(PROTNAME,PDB,TOP)
	#make_unformatted_traj(G96)
	#convert_traj("unformatted_traj")
	#Umatrix('unformatted_traj.npy')
	#fric_calc(PROTNAME,TOP)
	#LUI_calc(float(np.loadtxt('fratio')),float(np.loadtxt('avblsq')),float(np.loadtxt('sigma.dat')),np.loadtxt('fric'),np.load('Rinv.npy'),T=298)
	mode_mad('unformatted_traj.npy',np.load('Qmatrix.npy'),np.load('QINVmatrix.npy'))
	tau_convert(np.loadtxt('lambda_eig'),float(np.loadtxt('sigma.dat')),np.loadtxt('barriers.dat'))
except IndexError:
	print("Sript requires three command line inputs: PROTNAME, PDB file name, and topology name.")
