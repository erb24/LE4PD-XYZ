# coding: utf-8

def gen_protinfo(PROTNAME,G96,TOP):
	import subprocess
	import numpy as np

	protname = str(PROTNAME)
	N = int(subprocess.check_output('grep -c "CA" '+str(TOP),shell=True))
	NFRS = int(subprocess.check_output('grep -c "TIMESTEP" '+str(G96),shell=True))
	NATOMS = int(subprocess.check_output('grep -c "ATOM" '+str(TOP),shell=True))

	return N, NFRS, NATOMS
	#array = np.array([protname,N,NFRS,NATOMS],dtype=str)
	#np.savetxt("protname.txt",array.T,fmt="%s")

def convert_traj(G96):
	import numpy as np
	import subprocess
	import platform

	#Check which system the code is run on. If Linux, I can use 'sed'. For Darwin (macOS) I need 'gsed'.
	if platform.system() == 'Linux':
		status = subprocess.call("sed '/BOX/, +1 d' " + str(G96) + " | sed '/TITLE/, +1 d' | awk 'NF==3' > tmp",shell=True)
	elif platform.system() == 'Darwin':
		status = subprocess.call("gsed '/BOX/, +1 d' " + str(G96) + " | gsed '/TITLE/, +1 d' | awk 'NF==3' > tmp",shell=True)
	else:
		raise OSError("System platform not recognized.")
	if status == 0:
		traj = np.loadtxt('tmp')

		#np.save('unformatted_traj.npy',traj)
		subprocess.call('rm -rfv tmp', shell=True)
		return traj
	else:
		raise OSError('''Something has gone incorrectly and the unformatted trajectory was not generated.
		Please check where the .g96 file is located and make sure the correct PATH is specified 
		in the call to this function.''')

def format_traj(traj, N, NFRS):
	import numpy as np

	ftraj = np.zeros((3*N, NFRS))
	for numba, k in enumerate(range(0, N*NFRS, N)):
		ftraj[::3, numba] = traj[k:k+N, 0]
		ftraj[1::3, numba] = traj[k:k+N, 1]
		ftraj[2::3, numba] = traj[k:k+N, 2]

	return ftraj

def Umatrix(traj, protname, N, nfrs, natoms):
	import numpy as np

	#Old definitions of parameters when reading in from file
	#txt = np.genfromtxt('protname.txt',dtype=str)
	#protname = txt[0]
	#N = int(txt[1])
	#nfrs = int(txt[2])
	#natoms = int(txt[3])

	print(protname,N,nfrs,natoms)

	rx = np.zeros((N,nfrs))
	ry = np.zeros((N,nfrs))
	rz = np.zeros((N,nfrs))
	Rinv = np.zeros((N,N))

	#If the trajectory is unformatted
	#traj = np.load(str(traj))
	#for numba,k in enumerate(range(0,N*nfrs,N)):
	#	rx[:,numba] = traj[k:k+N,0]
	#	ry[:,numba] = traj[k:k+N,1]
	#	rz[:,numba] = traj[k:k+N,2]

	#Already saved the formatted trajectory
	for numba, i in enumerate(range(0, 3*N, 3)):
		rx[numba, :] = traj[i, :]
		ry[numba, :] = traj[i + 1, :]
		rz[numba, :] = traj[i + 2, :]

	#Define bond vectors
	lx = rx[:-1,:] - rx[1:,:] 
	ly = ry[:-1,:] - ry[1:,:]
	lz = rz[:-1,:] - rz[1:,:]

	#Calculate average inverse distances for the
	#hydrodynamic interaction matrix
	for i in range(N):
			for j in range(i,N):
				if i == j:
					Rinv[i,i] = np.nan
				else:
					Rinv[i,j] = (1/np.sqrt((rx[i,:] - rx[j,:])**2 + (ry[i,:] - ry[j,:])**2 + (rz[i,:] - rz[j,:])**2)).sum()
				Rinv[j,i] = Rinv[i,j]
	Rinv = Rinv/nfrs
	lavm = np.zeros(N-1)
	lavmsq = np.zeros(N-1)
	avdot = np.zeros(N-1)
	avblsq = 0
	avbl = 0

	avgx = lx.mean(1)
	avgy = ly.mean(1)
	avgz = lz.mean(1)

	for i in range(N-1):
		for k in range(nfrs):
			dummy = lx[i,k]**2 + ly[i,k]**2 + lz[i,k]**2
			lavm[i] += np.sqrt(dummy)
			lavmsq[i] += dummy
			avblsq += dummy 
		avbl += lavm[i]
		avdot[i] = avgx[i]**2 + avgy[i]**2 + avgz[i]**2
	lavm = lavm/nfrs
	lavmsq = lavmsq/nfrs
	avblsq = (avblsq/(N-1))/nfrs
	avbl = (avbl/(N-1))/nfrs

	print(avbl,avblsq)

	Umat = np.zeros((N-1,N-1))
	for i in range(N-1):
		for j in range(N-1):
			Umat[i,j] = (np.dot(lx[i,:],lx[j,:]) + np.dot(ly[i,:],ly[j,:]) + np.dot(lz[i,:],lz[j,:]))/(lavm[i]*lavm[j]*nfrs)

	#np.save('Umatrix.npy',Umat)
	#np.savetxt('Umatrix',np.insert(np.ravel(Umat),N-1,0))
	#np.save('Rinv.npy',Rinv)
	#np.savetxt('Rij',np.ravel(Rinv))
	#np.savetxt('length',lavm)
	#np.savetxt('lengthsq',lavmsq)
	#np.savetxt('avldot.dat',avdot)
	#np.savetxt('avbl',np.array([avbl]))
	#np.savetxt('avblsq',np.array([avblsq]))

	return Umat, Rinv, avbl, avblsq

def fric_calc(TOP, protname, N, nfrs, natoms, avblsq, T, intv = 2.71828, viscosity = 1e-3, fd20 = 0.0, path_to_resarea = './'):
	import numpy as np
	import sys
	import os
	import subprocess

	#TOP = str(sys.argv[1])
	pi = np.pi

	print(protname, N, nfrs, natoms)

	#Calculate the Miller radius per bead
	mradlist = []
	with open(TOP) as f:
		for line in f:
			if line[0:4] != 'ATOM':
				#print(line)
				pass
			elif line[0:4] == 'ATOM' and line.split()[2] == "CA":
				dummy = line.split()

				#Really horrendous and vestigial; probably smoother
				#to make a dictionary with the residue names plus their
				#assoicated Miller radii
				if dummy[3] == "ALA": mradlist.append((113.0/(4*pi))**.5)
				elif dummy[3] == "ARG" : mradlist.append((241.0/(4*pi))**.5)
				elif dummy[3] == "ASN" : mradlist.append((158.0/(4*pi))**.5)
				elif dummy[3] == "ASP" : mradlist.append((151.0/(4*pi))**.5)
				elif dummy[3] == "CYS" : mradlist.append((140.0/(4*pi))**.5)
				elif dummy[3] == "GLN" : mradlist.append((189.0/(4*pi))**.5)
				elif dummy[3] == "GLU" : mradlist.append((113.0/(4*pi))**.5)
				elif dummy[3] == "GLY" : mradlist.append((85.0/(4*pi))**.5)
				elif dummy[3] == "HIS" : mradlist.append((194.0/(4*pi))**.5)
				elif dummy[3] == "ILE" : mradlist.append((182.0/(4*pi))**.5)
				elif dummy[3] == "LEU" : mradlist.append((180.0/(4*pi))**.5)
				elif dummy[3] == "LYS" : mradlist.append((211.0/(4*pi))**.5)
				elif dummy[3] == "MET" : mradlist.append((204.0/(4*pi))**.5)
				elif dummy[3] == "PHE" : mradlist.append((218.0/(4*pi))**.5)
				elif dummy[3] == "PRO" : mradlist.append((143.0/(4*pi))**.5)
				elif dummy[3] == "SER" : mradlist.append((122.0/(4*pi))**.5)
				elif dummy[3] == "THR" : mradlist.append((146.0/(4*pi))**.5)
				elif dummy[3] == "TRP" : mradlist.append((259.0/(4*pi))**.5)
				elif dummy[3] == "TYR" : mradlist.append((229.0/(4*pi))**.5)
				elif dummy[3] == "VAL" : mradlist.append((160.0/(4*pi))**.5)

	#mrad_array = np.array(mradlist,dtype=str)
	#np.savetxt('mrad.dat',np.array(mradlist).T,fmt='%s')


	#Calculate the average solvent-exposed surface area per bead
	if os.path.exists(path_to_resarea + "resarea.xvg"):
		pass
	else:
		raise FileNotFoundError('''I can't find the resarea.xvg file containing the solvent-exposed surface area of each residue.
								Please either run the process.sh file, if you have not already done so, and move the resarea.xvg 
								file into the current working directory.''')
	resarea = []
	with open(path_to_resarea + 'resarea.xvg') as f:
		for line in f:
			if line[0] == '#' or line[0] == '@':
				pass
			else:
				resarea.append(float(line.split()[1]))

	rad = []
	for area in resarea:
		rad.append(((area/(4*np.pi))**0.5)*10)

	#np.savetxt('avresrad',np.array(rad),fmt="%f")
	fratio = (np.array(rad).sum()/N)/10
	print('fratio: ', fratio)

	#np.savetxt('fratio',np.array([fratio]))

	#Calculate the friction coefficients

	kB = 1.38066E-23
	print('Temperature (K): ',T)
	print('Internal viscosity factor: ',intv)

	#Use NIST formula for viscosity -- good NEAR room temperature and physiological.
	#Won't work higher than, say, 360 K.

	if viscosity == 0:
		print('No viscosity given. Using the NIST formula, which is only valid for physiological conditions,\n')
		print('i.e. between about 273 and 310 K.')
		viscosity = (.2131590-1.96290E-3*T+(.00246411*T)**2+(-.0018462*T)**3)

	print("Viscosity (Pa s): ", viscosity)
	print("fd20", fd20)

	rv = np.array(mradlist)
	rw = np.array(rad)
	rp = np.zeros(N)
	friw = np.zeros(N)
	fri = np.zeros(N)
	friwt = 0
	frit = 0
	for i in range(N):
		if rw[i] < rv[i]: 
			rp[i] = (rv[i]**2 - rw[i]**2)**0.5
		else:
			rp[i] = 0

		friw[i] = 6.0*pi*(rw[i]/10)*viscosity
		fri[i] = 6.0*pi*(rp[i]/10)*(intv*viscosity) + 6.0*pi*(rw[i]/10)*viscosity
		friwt += friw[i]
		frit += fri[i]

	avfr = frit/float(N)
	avfrw = friwt/float(N)
	#np.savetxt('avfr',np.array([avfr*1.0E-9]))

	#avblsq = float(np.loadtxt('avblsq'))
	sigma = (3*kB*T*1E15)/(avblsq*avfr)

	#with open('sigma','w') as f:
	#	f.write('sigma, 1/ps\n')
	#	f.write(str(sigma)+'\n')

	#with open('sigma.dat','w') as f:
	#	f.write(str(sigma)+'\n')

	fric = np.zeros((N+1,2))

	fric[0,0] = avfrw
	fric[0,1] = avfr
	for i in range(N):
		fric[i+1,:] = np.column_stack([friw[i],fri[i]])

	#np.savetxt('fric',fric)

	return fratio, sigma, fric, avfr

def LUI_calc(protname, N, nfrs, U, fratio, avblsq, sigma, fric, Rinv, T):
	import numpy as np
	#Get basic information from the protname.txt file
	#txt = np.genfromtxt('protname.txt',dtype=str)
	#protname = txt[0]
	#N = int(txt[1])
	#nfrs = int(txt[2])
	#natoms = int(txt[3])

	print(protname, N, nfrs)

	avfr = fric[0,1]

	M = np.zeros((N,N))

	for i in range(N):
		for j in range(N):
			if i == 0: 
				M[i,j] = 1/N
			else:
				if i == j + 1:
					M[i,j] = -1
				elif i == j:
					M[i,j] = 1

	a = M[1:,:]

	H = np.zeros((N,N))

	for i in range(N):
		for j in range(i,N):
			if i == j:
				H[i,i] = avfr/fric[i+1,1]
			else:
				H[i,j] = fratio*Rinv[i,j]
				H[j,i] = H[i,j]

	L = np.matmul(a,np.matmul(H,a.T))

	LINV = np.linalg.inv(L)
	L = np.matmul(a,np.matmul(H,a.T))
	LINV = np.linalg.inv(L)
	UINV = U #np.load("Umatrix.npy")
	UILI = np.matmul(UINV,LINV)

	#Eigendecomposition of UILI to find eigenvalues and eigenvectors 
	eigval,Q = np.linalg.eig(UILI)
	eigval = 1/eigval
	QINV = np.linalg.inv(Q)

	perm = np.argsort(np.abs(eigval))

	Q_sorted = np.copy(Q)[:,perm]
	QINV_sorted = np.linalg.inv(Q_sorted)

	eigval_sorted = abs(eigval)[perm]

	mu = 1/(np.diag(np.matmul(QINV_sorted,np.matmul(UINV,QINV_sorted.T))))

	#np.save("UILImatrix",UILI)
	#np.savetxt('UILImatrix',np.ravel(UILI))
	#np.save('Lmatrix',L)
	#np.save("Hmatrix",H)
	#np.save("Qmatrix.npy",Q_sorted)
	#np.savetxt("Qmatrix",np.ravel(Q_sorted))
	#np.save("QINVmatrix.npy",QINV_sorted)
	#np.savetxt("QINVmatrix",np.ravel(QINV_sorted))
	#np.save("lambda_eig.npy",eigval_sorted)
	#np.savetxt("lambda_eig",eigval_sorted)
	#np.save("mu_eig.npy",mu)
	#np.savetxt("mu_eig",mu)

	return UILI, H, Q_sorted, QINV_sorted, eigval_sorted, mu

def mode_mad(traj, protname, N, nfrs, Q, QINV, T, nmodes = 10):
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	import subprocess
	import physt

	pi = np.pi

	print(protname, N, nfrs)

	rx = np.zeros((N,nfrs))
	ry = np.zeros((N,nfrs))
	rz = np.zeros((N,nfrs))

	#Already saved the formatted trajectory
	for numba, i in enumerate(range(0, 3*N, 3)):
		rx[numba, :] = traj[i, :]
		ry[numba, :] = traj[i + 1, :]
		rz[numba, :] = traj[i + 2, :]

	#Define bond vectors
	lx = rx[:-1,:] - rx[1:,:] 
	ly = ry[:-1,:] - ry[1:,:]
	lz = rz[:-1,:] - rz[1:,:]

	xix = np.matmul(QINV,lx)
	xiy = np.matmul(QINV,ly)
	xiz = np.matmul(QINV,lz)
	xim = np.sqrt(xix**2 + xiy**2 + xiz**2)
	xi = xix + xiy +xiz

	theta = np.arccos(xiz/xim)
	phi = np.arctan(xiy/xix)

	if nmodes > N - 1:
		raise ValueError('''Specified more modes than are generated by the analysis.\n
			The isotropic LE4PD generates N - 1 modes, with N the number of residues in the protein.''')
	else:

		for a in range(N - 1):
			for k in range(nfrs):
				if xix[a,k] <= 0.0:
					phi[a,k] += pi
				if phi[a,k] <= 0.0:
					phi[a,k] += 2*pi

		#theta = np.rad2deg(theta)
		#phi = np.rad2deg(phi)

		#Make histogram
		kT = 0.00198*T
		fmadlist = []
		for a in range(N - 1):
			x=xim[a,:]*np.sin(theta[a,:])*np.cos(phi[a,:])
			y=xim[a,:]*np.sin(theta[a,:])*np.sin(phi[a,:])
			z=xim[a,:]*np.cos(theta[a,:])
			h=physt.special.spherical_histogram(np.column_stack([x,y,z]),theta_bins=50,phi_bins=50,radial_bins=1)
			his=(h.densities[0]/h.densities[0].sum()).T
			fes = -kT*np.log(his)
			#fes -= fes.min()
			femax = -kT*np.log(1/nfrs)
			for j in range(fes.shape[0]):
				for k in range(fes.shape[1]):
					if np.isinf(fes[j,k]) == True:
						fes[j,k] = femax
						
			#xx,yy=np.meshgrid(np.linspace(0,180,his.shape[1]),np.linspace(0,360,his.shape[0]))
			#im=plt.contourf(xx,yy,fes,levels=np.linspace(fes.min(),fes.max(),25),cmap='gnuplot')
			#cbar=plt.colorbar(im)
			#cbar.set_label(r'Free Energy $(k_BT)$')
			#plt.contour(xx,yy,fes,levels=np.linspace(fes.min(),fes.max(),25),colors='k')
			#plt.xlabel(r'$\theta$ (deg)')
			#plt.ylabel(r'$\phi$ (deg)')
			#plt.savefig('fes'+str(a+1)+'.eps',dpi=300)
			#plt.show()
			#plt.close()
			dummy = list(np.ravel(fes))
			#Cut-off for outliers 
			cut = -kT*np.log(2/nfrs)
			for num,i in enumerate(dummy):
				if i >= cut: dummy.remove(i)
			fmad = np.median(abs(np.array(dummy) - fes.min()))

			fmadlist.append(fmad)

		#Make a directory for this mode-dependent analysis
		if os.path.exists('mode_analysis'):
			pass
		else:
			os.mkdir('mode_analysis/')
		path = 'mode_analysis/'

		for a in range(nmodes):
			np.save(path + 'fes'+str(a+1)+'.npy', fes)
			np.save(path + 'theta_phi_'+str(a+1)+'.npy', np.column_stack([np.rad2deg(theta[a,:]),np.rad2deg(phi[a,:]), xim[a,:]]))
			np.save(path + "xi_"+str(a+1)+'.npy', xi[a,:])
		np.savetxt(path + 'barriers_kcal.dat', np.array(fmadlist))
		np.savetxt(path + 'barriers.dat', np.array(fmadlist)/kT)

		return np.array(fmadlist)/kT, xi, np.column_stack([np.rad2deg(theta[a,:]),np.rad2deg(phi[a,:]), xim[a,:]])

def tau_convert(eigvals, sigma, bar, T):
	import numpy as np
	import os
	
	if os.path.exists('mode_analysis'):
		pass
	else:
		os.mkdir('mode_analysis/')
	path = 'mode_analysis/'

	tau = (eigvals*sigma)**-1
	tau_scaled = tau*np.exp(bar)
	np.savetxt(path + 'tau.dat', np.column_stack([np.arange(1,len(tau)+1),tau]))
	np.savetxt(path + 'tau_scaled.dat', np.column_stack([np.arange(1,len(tau)+1),tau_scaled]))

	return tau, tau_scaled

def LML(Qmatrix, avbl, mu):
	import numpy as np
	import os
	
	LML = np.zeros_like(Qmatrix)
	for a in range(Qmatrix.shape[0]):
		for i in range(Qmatrix.shape[1]):
			LML[i,a] = (Qmatrix[i,a]**2*avbl**2)/mu[a]
	LML = np.sqrt(LML)

	if os.path.exists('mode_analysis'):
		pass
	else:
		os.mkdir('mode_analysis/')
	path = 'mode_analysis/'

	np.save(path + "LML.npy", LML)

	return LML

#Need to run m1_f2py.sh first to generate the m1int 'module' file
#def calc_m1(N):
#	import numpy as np
#	import m1int

#	m1int.p2m1(N)
