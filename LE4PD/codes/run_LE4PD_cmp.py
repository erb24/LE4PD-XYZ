# coding: utf-8

import numpy as np
import sys
import subprocess
from LE4PD_cmp import *

try:
	PROTNAME = sys.argv[1]
	PDB = sys.argv[2]
	TOP = sys.argv[3]
	T = sys.argv[4]

	#subprocess.call('echo "3" | gmx_mpi trjconv -f '+str(PROTNAME)+'.xtc -s '+str(TOP)+' -o '+str(G96),shell=True)
	#subprocess.call('echo "3" | gmx_mpi trjconv -f '+str(PROTNAME)+'.xtc -s '+str(TOP)+' -o '+str(PDB),shell=True)
	#subprocess.call('echo "1" | gmx_mpi trjconv -f '+str(PROTNAME)+'.xtc -s '+str(TOP)+' -o '+str(PROTNAME)+'_first.pdb -dump 0',shell=True)
	#subprocess.call("sh unformat_traj.sh "+str(PDB),shell=True)

	TOP=str(PROTNAME)+"_first.pdb"
	G96=str(PROTNAME)+".g96"
	gen_protinfo(PROTNAME,G96,TOP)
	get_chain_info(PROTNAME,TOP)

	#make_unformatted_traj(G96)
	convert_traj(G96)
	Umatrix('unformatted_traj.npy')
	subprocess.call('echo "'+str(T)+'" > temp', shell=True)
	fric_calc(PROTNAME,TOP)
	subprocess.call('echo "0.306" > visc.txt', shell=True)
	LUI_calc(float(np.loadtxt('fratio')),float(np.loadtxt('avblsq')),float(np.loadtxt('sigma.dat')),np.loadtxt('fric'),np.load('Rinv.npy'),T=310)
	mode_mad('unformatted_traj.npy',np.load('Qmatrix.npy'),np.load('QINVmatrix.npy'))
	tau_convert(np.loadtxt('lambda_eig'),float(np.loadtxt('sigma.dat')),np.loadtxt('barriers.dat'))
except IndexError:
	print("Sript requires three command line inputs: PROTNAME, PDB file name, topology name, and temperature (in Kelvin).")
