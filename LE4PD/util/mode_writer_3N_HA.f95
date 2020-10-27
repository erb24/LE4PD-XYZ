	program inputreader
	IMPLICIT NONE
	integer :: nfrs,n,nmol
	open(unit=5,file="protname.txt",status='old')
	read(5,*)
	read(5,*)n
	read(5,*)nfrs
	close(5)
	!nfrs=5000 !Testing
	write(*,*)n,nfrs
	call project(n,nfrs)

	STOP
	End program inputreader

	subroutine project(n,nfrs)
	IMPLICIT NONE
	integer :: i,j,k,imin,jmin,a
	INTEGER, INTENT(in) :: n,nfrs
	INTEGER :: xicount,xjcount,xcount,ycount,zcount,mineig(1),order(n+1),rcount
	INTEGER :: yicount,yjcount,zicount,zjcount,icount,jcount,xlcount,ylcount,zlcount
	INTEGER :: xmcount,ymcount,zmcount,counter
	DOUBLE PRECISION, dimension(n,nfrs) :: lx,ly,lz,lmag,rx,ry,rz,rvec
	DOUBLE PRECISION :: biglist(3*n,nfrs)
	DOUBLE PRECISION :: Rb,T
	DOUBLE PRECISION :: RO(n,n),cLE(n,n),Q(3*n,3*n),QINV(3*n,3*n)
	DOUBLE PRECISION :: xix(3*n,nfrs),xiy(3*n,nfrs),xiz(3*n,nfrs),xim(3*n,nfrs),xi(3*n,nfrs),xiavg(3*n)
	DOUBLE PRECISION :: traj(3*n,nfrs),mean,var
	character(16)aa,ii,cbins,tipo
	CHARACTER(32) :: protname

	biglist=0.0
	traj=0.0
	rx=0.0
	ry=0.0
	rz=0.0
	

	open(unit=5,file="protname.txt",status='old')
	read(5,'(A)')protname
	close(5)

	open(unit=10,file="com_HCmatrix",status='old')
	open(unit=11,file="com_HCINVmatrix",status='old')
	do i=1,3*n
	  do j=1,3*n
	    READ(10,*)Q(i,j)
	    READ(11,*)QINV(i,j)
	  end do
	end do
	close(10)
	close(11)

	OPEN(unit=24,file='com_biglist.dat',status='old')
	DO k=1,nfrs
	  IF(MOD(k,10000) .EQ. 0)WRITE(*,*)'READING FRAME ',k
	  DO i=1,3*n
	    READ(24,*)traj(i,k)
	  END DO
	END DO
	CLOSE(24)
	!Calculate the trajectories in the coordinates of the modes
	xi=0.0
	xix=0.0
	xiy=0.0
	xiz=0.0
        mean = 0.0
        var = 0.0
	DO a=1,3*n
	  WRITE(aa,*)a
	  aa=adjustl(aa)
	  WRITE(*,*)'CALCULATING THE TRAJECTORY FOR MODE ',a
	  OPEN(unit=25,file='xi_'//TRIM(aa)//'HA.xvg',status='unknown')
	  OPEN(unit=27,file='xi_'//TRIM(aa)//'_HA_aniso.xvg',status='unknown')
	  !OPEN(unit=28,file='diff_mode_'//TRIM(aa)//'.dat',status='unknown')
	  DO k=1,nfrs
	    !jcount=1
	    DO j=1,3*n,3
	      xix(a,k)=QINV(a,j)*traj(j,k)+xix(a,k)
	      xiy(a,k)=QINV(a,j+1)*traj(j+1,k)+xiy(a,k)
	      xiz(a,k)=QINV(a,j+2)*traj(j+2,k)+xiz(a,k)
	      !jcount=jcount+1
	    END DO
            xi(a,k)=xix(a,k)+xiy(a,k)+xiz(a,k)
	    WRITE(25,*)xi(a,k)
            WRITE(27,*)xix(a,k),xiy(a,k),xiz(a,k)
            mean = mean + xi(a,k)
            var = var + xi(a,k)**2
	  END DO
	  CLOSE(25)
          CLOSE(27)
          mean = mean/nfrs
          var = (var/nfrs) - mean**2
          WRITE(*,*)'Mean: ',mean
          WRITE(*,*)'Variance: ',var
          WRITE(*,*)
	END DO
	end subroutine
