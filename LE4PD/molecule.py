import warnings
import numpy as np
import subprocess
import os 

import mdtraj as md
from LE4PD.util import prepare
from LE4PD.codes import LE4PD_XYZ as LE4PD

'''Module to manage topology of varying molecule types. Currently only proteins
are implemented. Future versions will contain nucleic acids such as RNA and DNA.
'''

class protein(object):
	'''

	'''
	def __init__(self, method='simulation', comp = False, isotropic = True, T = 298):
		self._method = method
		self.comp = comp 
		self.isotropic = isotropic
		self.temp = T 


	def load(self, traj, top = None, skip_atoms = None, skip_residues = None, chunk_size = 1000):
		if self._method == 'ensemble':
			self._topfile = traj
			self._skip_atoms = skip_atoms
			self._skip_residues = skip_residues
			self._mdtraj = prepare.fetch(traj, skip_atoms, skip_residues)
			prepare.topology(self)

		elif self._method == 'simulation':
			#It isn't worth it right now to allow the trajectory conversion here. Probably better to give
			#the user a .sh file to do the conversion themselves. The modularized Python yahoos can get over it.

			self._trajfile = traj
			self._topfile = top
			self.protname = traj.split('.')[0]
			N, NFRS, NATOMS = LE4PD.gen_protinfo(self.protname, self._trajfile, self._topfile)
			self.nres = N
			self.nframes = NFRS
			self.natoms = NATOMS


	def prepare_traj(self):
		self.unformatted_traj = LE4PD.convert_traj(self._trajfile)
		self.traj = LE4PD.format_traj(self.unformatted_traj, self.nres, self.nframes)

	def calc_umatrix(self):
		self.Cmatrix, self.Q, self.mu_eig = LE4PD.calc_Cmatrix(self.traj, self.protname, self.nres, self.nframes, self.natoms)

	def Rij(self, path = './'):
		self.Rinv = LE4PD.Rij(self.traj, path)

	def fric_calc(self, intv = 2.71828, viscosity = 1e-3, fd20 = 0.0, path_to_resarea = './'):
		self.fratio, self.sigma, self.fric, self.avfr = LE4PD.fric_calc(self._topfile, self.protname, self.nres, self.nframes, 
																		self.natoms, self.avblsq, self.temp, intv = intv, viscosity = viscosity, fd20 = fd20,
																		path_to_resarea = path_to_resarea)

	def hydrodynamics(self, path = './'):
		self.Hmatrix, self.AIHI, self.QHA, self.QHAINV, self.lambda_eig = LE4PD.eigendecomp(self.Cmatrix, self.Rinv, 
																				self.avfr, self.fratio, self.fric, path = path)
