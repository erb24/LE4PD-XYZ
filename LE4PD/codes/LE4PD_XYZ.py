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
	subprocess.call("sed '/BOX/, +1 d' "+str(G96)+" | sed '/TITLE/, +1 d' | awk 'NF==3' > tmp",shell=True)
	traj = np.loadtxt('tmp')

	np.save('unformatted_traj.npy',traj)
	subprocess.call('rm -rfv tmp',shell=True)

def format_traj(traj, N, NFRS):
        import numpy as np

        ftraj = np.zeros((3*N, NFRS))
        for numba, k in enumerate(range(0, N*NFRS, N)):
                ftraj[::3, numba] = traj[k:k+N, 0]
                ftraj[1::3, numba] = traj[k:k+N, 1]
                ftraj[2::3, numba] = traj[k:k+N, 2]

        return ftraj

def calc_Cmatrix(traj, path = './'):
	#traj = np.load(path + 'unformatted_traj.npy')
	avg_coor = traj.mean(1)
	dtraj = np.copy(traj)
	for i in range(traj.shape[0]):
		dtraj[i,:] = traj[i,:] - avg_coor[i]


	#Likely need a better way to perform the calculation of the covariance matrix
	#because this 'brute force' approach requires a lot of memory (but can be run
	#using a large-shared node on Comet).
	covar = np.matmul(dtraj, dtraj.T) / dtraj.shape[1]
	eigvals, R = np.linalg.eigh(covar)
	R_sorted = R[:,np.argsort(1/abs(eigvals))]
	#np.save(path + 'com_covariance_matrix.npy', covar)
	#np.savetxt(path + 'com_covariance_matrix.dat', covar.ravel())
	#np.save(path + 'com_Rmatrix.npy', R_sorted)
	#np.save(path + 'com_RTmatrix.npy', R_sorted.T)
	#np.savetxt(path + 'com_covar_eig', eigvals[np.argsort(1/abs(eigvals))])

	return covar, R_sorted, eigvals[np.argsort(1/abs(eigvals))]

def Rij(traj, path = './'):
	import numpy as np
	N = traj.shape[0] // 3
	nfrs = traj.shape[1]

	print(N, nfrs)

	rx = np.zeros((N,nfrs))
	ry = np.zeros((N,nfrs))
	rz = np.zeros((N,nfrs))
	Rinv = np.zeros((N,N))
	for i,n in enumerate(range(0, N*3, 3)):
		rx[i, :] = traj[n, :]
		ry[i, :] = traj[n + 1, :]
		rz[i, :] = traj[n + 2, :]

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

	np.save(path + 'Rinv.npy', Rinv)
	#np.savetxt(path + 'Rij', np.ravel(Rinv))
	return Rinv

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


def eigendecomp(covar, rij, avfr, avresrad, fric, path = './'):
	import numpy as np
	import matplotlib.pyplot as plt
	import sys
	import os

	N = covar.shape[0] // 3

	#covar = np.reshape(np.loadtxt(path + 'com_covariance_matrix.dat'),(3*N,3*N))
	cov,R = np.linalg.eigh(covar)

	#Sort
	eiglam=np.copy(abs(cov))
	order=np.zeros(3*N, dtype=int)
	eiglamO=np.zeros(3*N)
	for i in range(3*N):
		eiglamO[i]=max(eiglam)
		mineig=np.argmax(eiglam)
		order[i]=mineig
		eiglam[order[i]]=1E-20
		
	#Order the eigenvectors
	QMO=np.zeros((3*N,3*N))
	eigvec=np.copy(R)
	for i in range(3*N):
			QMO[:,i]=eigvec[:,order[i]]
	Q=QMO

	cov = eiglamO

	R = Q

	'''np.savetxt(path + 'com_Rmatrix',np.ravel(R))
	np.savetxt(path + 'com_RTmatrix',np.ravel(R.T))
	np.savetxt(path + 'com_covar_eig',cov)'''

	#Use Kronecker product to get H

	H = np.zeros((N,N))

	#avfr = float(np.loadtxt(path + 'avfr'))
	#avresrad = np.loadtxt(path + 'avresrad').mean()/10
	#avresrad = 0.22
	#fric = np.loadtxt(path + 'fric')[1:,1]*1e-9

	for i in range(H.shape[0]):
		for j in range(i,H.shape[1]):
			if i == j:
				H[i,j] = (avfr)/fric[i]
			else:
				H[i,j] = avresrad*rij[i,j]
				
			H[j,i] = H[i,j]
			
	H3N = np.kron(H,np.eye(3))

	vals,vecs = np.linalg.eig(np.matmul(covar,np.linalg.inv(H3N)))

	#Sort
	eiglam=np.copy(abs(vals))
	order=np.zeros(3*N,dtype=int)
	eiglamO=np.zeros(3*N)
	for i in range(3*N):
		eiglamO[i]=max(eiglam)
		mineig=np.argmax(eiglam)
		order[i]=mineig
		eiglam[order[i]]=1E-20
		
	#Order the eigenvectors
	QMO=np.zeros((3*N,3*N))
	eigvec=np.copy(vecs)
	for i in range(3*N):
			QMO[:,i]=eigvec[:,order[i]]
	Q=QMO
	try:
		QINV = np.linalg.inv(Q)
	except np.linalg.LinAlgError:
		print('Q is a singular matrix')
		QINV = np.linalg.pinv(Q)

	#np.save(path + 'com_HCmatrix.npy', Q)
	#np.save(path + 'com_HCINVmatrix.npy', QINV)
	#np.savetxt(path + 'com_HC_eig', eiglamO)

	return H3N, np.matmul(covar,np.linalg.inv(H3N)), Q, QINV, eiglamO