	program inputreader
	integer nfrs,n
	open(unit=5,file="protname.txt",status='old')
	read(5,*)
	read(5,*)n
	read(5,*)nfrs
	close(5)
	nbins=60
	write(*,*)n,nfrs
	call umatrix(n,nfrs,nbins)
	End program inputreader

	subroutine umatrix(n,nfrs,nbins)
	integer i,j,k,imin,jmin,a,nfe,ir,nbinsrot,imad,imadmid
	INTEGER :: xcount,ycount,zcount
	REAL,DIMENSION(n) :: avgx,avgy,avgz
	real, dimension(n,nfrs) :: rx,ry,rz,rxref,ryref,rzref
	REAL,DIMENSION(3*n) ::sigfe,fricorr,pvol,
     &fmax,fmin,fmaxp,fminp,fbg,bar
	real, dimension(n) :: lx,ly,lz,lmag
	REAL,DIMENSION(3*n) :: avfe,avfesq,fenorm
	real, dimension(3*n,3*n) :: sigij,rij,qinvm,qm
	real, dimension(3*n,nfrs) :: xix,xiy,xiz,dipcorr,xim,
     &theta,phi,biglist
	real dotpij,um,rrij,bl,hrtheta,hrphi,Rb,T,r,dr,delphi
	integer itheta,iphi
	real hisang(3*n,-nbins:nbins,-nbins:nbins),lambda(3*n)
	character(32)protname
	character(16)aa,ii,cbins
	real hisp,hismax,delha,rdeg,degr,hnorm(3*n),x,y,z,emad
	real feang(3*n,-nbins:nbins,-nbins:nbins),femax,pi,delr
	real testnorm,dc,fmad(nbins*nbins),fmadord(nbins*nbins)
	LOGICAL :: there
c	nfrs=10000 !just for testiing!
	Rb=.00198 !(boltzmanns constant in kcal/mol*K)
	open(unit=10,file='temp')
	read(10,*)T
	close(10)
	felim=0.0
	femin=0.0
	feminp=0.0
	sigfe=0.0
	fricorr=0.0
	avfe=0.0
	avfesq=0.0
	rij=0.0
	hisp=100.0
	hismax=0.0
	xix=0.0
	xiy=0.0
	xiz=0.0
	xim=0.0
	qinvm=0.0
	qm=0.0
	dipcorr=0.0	
	xim=0.0
	theta=0.0
	phi=0.0
	hisang=0.0
	hrtheta=0.0
	hrphi=0.0
	itheta=0
	iphi=0
	pi=3.1415927
	delha=(2.0*360.0)/real(2*nbins)
	degr=((2.0*pi)/360.0) !deg to rad
	rdeg=1.0/degr !rad to deg
	pvol=0.0
	fenorm=0.0
	r=0.0
	dr=5.0/(real(nfrs))
	ir=0
	dc=1.0/real(nfrs)
	fmax=0.0
	fmaxp=0.0
	fmin=0.0
	fminp=200.0
	fbg=0.0
	fmad=0.0
	fmadord=0.0

	delr=delha*degr
	write(*,*)"delr",delr
	hnorm=0.0
	open(unit=21,file='com_RTmatrix',status='old')
	do i=1,3*n
	  do j=1,3*n
	    read(21,*)qinvm(i,j)
	  end do
	end do
	open(unit=5,file="protname.txt",status='old')
	read(5,'(A)')protname
	close(5)
	rx=0.0
	ry=0.0
	rz=0.0
	lx=0.0
	ly=0.0
	lz=0.0
	lmag=0.0
	lavm=0.0
	dotpij=0.0
	sigij=0.0
	um=0.0
	rrij=0.0
	rij=0.0
	bl=0.0
	imin=0
	jmin=0
	bar=0.0
	lambda=0.0

	!read from trajectory
	avgx=0.0
	avgy=0.0
	avgz=0.0
	biglist=0.0
	INQUIRE(FILE='com_biglist.dat',EXIST=there)
	WRITE(*,*)there
	IF ( there ) THEN
	  OPEN(unit=24,file='com_biglist.dat',status='old')
	  DO k=1,nfrs
	    IF(MOD(k,1000) .EQ. 0)WRITE(*,*)'Reading frame ',k
	    DO i=1,3*n
	      READ(24,*)biglist(i,k)
	    END DO
!	calculating instantaneous mode vector
	  do a=1,3*n !mode loop
	    do j=1,3*n !residue loop
	      IF(MOD(j,3) .EQ. 1)xix(a,k)=qinvm(a,j)*biglist(j,k)+xix(a,k)
	      IF(MOD(j,3) .EQ. 2)xiy(a,k)=qinvm(a,j)*biglist(j,k)+xiy(a,k)
	      IF(MOD(j,3) .EQ. 0)xiz(a,k)=qinvm(a,j)*biglist(j,k)+xiz(a,k)
	    end do
	    xim(a,k)=(xix(a,k)**2+xiy(a,k)**2+xiz(a,k)**2)**.5
	  end do
!	calculate theta, phi
	do a=1,3*n
	theta(a,k)=acos(xiz(a,k)/xim(a,k))
	phi(a,k)=atan(xiy(a,k)/xix(a,k))
	if(xix(a,k).lt.0.0)phi(a,k)=phi(a,k)+pi
	theta(a,k)=theta(a,k)*rdeg
	if(phi(a,k).lt.0.0)phi(a,k)=phi(a,k)+2.0*pi
	phi(a,k)=phi(a,k)*rdeg
c	if(a.eq.4)write(*,*)theta(a,k),phi(a,k)
	end do

	!write into histogram: need to break phi up into different number of bins dependent upon sin theta
	do a=1,3*n
	delphi=360./real(nint((nbins*sin(theta(a,k)*degr))))
	hrtheta=theta(a,k)/delha
	hrphi=phi(a,k)/delphi
	itheta=nint(hrtheta)
	iphi=nint(hrphi)
	hisang(a,itheta,iphi)=hisang(a,itheta,iphi)+dc
c	if(a.eq.4)write(*,*)itheta,iphi,hisang(a,itheta,iphi)
	end do	
	  END DO
	  CLOSE(24)
	ELSE
	  open(unit=11,file=trim(protname)//'.g96',status='old')
	  !skip first 7,now read and calculate stuff
	  do i=1,7
	  read(11,*)
	  end do

	  do k=1,nfrs
	    do j=1,n
	      read(11,*)rx(j,k),ry(j,k),rz(j,k)
	      avgx(j)=avgx(j)+rx(j,k)
	      avgy(j)=avgy(j)+ry(j,k)
	      avgz(j)=avgz(j)+rz(j,k)
	    end do
	     !skip 8 lines
	    do j=1,8
	      read(11,*)
	    end do
	  END DO
	  DO j=1,n
	    avgx(j)=avgx(j)/REAL(nfrs)
	    avgy(j)=avgy(j)/REAL(nfrs)
	    avgz(j)=avgz(j)/REAL(nfrs)
	  END DO
          !Calculate reference coordinates
	  DO k=1,nfrs
	    DO i=1,n
	      rxref(i,k)=rx(i,k)-avgx(i)
	      ryref(i,k)=ry(i,k)-avgy(i)
	      rzref(i,k)=rz(i,k)-avgz(i)
	    END DO
  	  END DO
 	  !Make big list of all augmented coordinates
	  xcount=0
	  ycount=0
	  zcount=0
	  WRITE(*,*)'Making the big list of coordinates...'
	  DO k=1,nfrs
	    DO i=1,3*n
	      IF(MOD(i,3) .EQ. 1)THEN
	        xcount=xcount+1
	        biglist(i,k)=rxref(xcount,k)
	      ELSE IF(MOD(i,3) .EQ. 2)THEN
	        ycount=ycount+1
	        biglist(i,k)=ryref(ycount,k)	
	      ELSE IF(MOD(i,3) .EQ. 0)THEN
	        zcount=zcount+1
	        biglist(i,k)=rzref(zcount,k)	 
	      END IF
	    END DO
	    !Reset counters after each time step
	    xcount=0
	    ycount=0
	    zcount=0
!	  calculating instantaneous mode vector
	  do a=1,3*n !mode loop
	    do j=1,3*n !residue loop
	      IF(MOD(j,3) .EQ. 1)xix(a,k)=qinvm(a,j)*biglist(j,k)+xix(a,k)
	      IF(MOD(j,3) .EQ. 2)xiy(a,k)=qinvm(a,j)*biglist(j,k)+xiy(a,k)
	      IF(MOD(j,3) .EQ. 0)xiz(a,k)=qinvm(a,j)*biglist(j,k)+xiz(a,k)
	    end do
	    xim(a,k)=(xix(a,k)**2+xiy(a,k)**2+xiz(a,k)**2)**.5
	  end do
!	calculate theta, phi
	do a=1,3*n
	theta(a,k)=acos(xiz(a,k)/xim(a,k))
	phi(a,k)=atan(xiy(a,k)/xix(a,k))
	if(xix(a,k).lt.0.0)phi(a,k)=phi(a,k)+pi
	theta(a,k)=theta(a,k)*rdeg
	if(phi(a,k).lt.0.0)phi(a,k)=phi(a,k)+2.0*pi
	phi(a,k)=phi(a,k)*rdeg
c	if(a.eq.4)write(*,*)theta(a,k),phi(a,k)
	end do

	!write into histogram: need to break phi up into different number of bins dependent upon sin theta
	do a=1,3*n
	delphi=360./real(nint((nbins*sin(theta(a,k)*degr))))
	hrtheta=theta(a,k)/delha
	hrphi=phi(a,k)/delphi
	itheta=nint(hrtheta)
	iphi=nint(hrphi)
	hisang(a,itheta,iphi)=hisang(a,itheta,iphi)+dc
c	if(a.eq.4)write(*,*)itheta,iphi,hisang(a,itheta,iphi)
	end do
	  END DO
	END IF
	!DO k=1,nfrs   
!	calculating instantaneous mode vector
	!  do a=1,3*n !mode loop
	!    do j=1,3*n !residue loop
	!      IF(MOD(j,3) .EQ. 1)xix(a,k)=qinvm(a,j)*biglist(j,k)+xix(a,k)
	!      IF(MOD(j,3) .EQ. 2)xiy(a,k)=qinvm(a,j)*biglist(j,k)+xiy(a,k)
	!      IF(MOD(j,3) .EQ. 0)xiz(a,k)=qinvm(a,j)*biglist(j,k)+xiz(a,k)
	!    end do
	!    xim(a,k)=(xix(a,k)**2+xiy(a,k)**2+xiz(a,k)**2)**.5
	!  end do
!	calculate theta, phi
	!do a=1,3*n
	!theta(a,k)=acos(xiz(a,k)/xim(a,k))
	!phi(a,k)=atan(xiy(a,k)/xix(a,k))
	!if(xix(a,k).lt.0.0)phi(a,k)=phi(a,k)+pi
	!theta(a,k)=theta(a,k)*rdeg
	!if(phi(a,k).lt.0.0)phi(a,k)=phi(a,k)+2.0*pi
	!phi(a,k)=phi(a,k)*rdeg
c	if(a.eq.4)write(*,*)theta(a,k),phi(a,k)
	!end do

	!write into histogram: need to break phi up into different number of bins dependent upon sin theta
	!do a=1,3*n
	!delphi=360./real(nint((nbins*sin(theta(a,k)*degr))))
	!hrtheta=theta(a,k)/delha
	!hrphi=phi(a,k)/delphi
	!itheta=nint(hrtheta)
	!iphi=nint(hrphi)
	!hisang(a,itheta,iphi)=hisang(a,itheta,iphi)+dc
c	if(a.eq.4)write(*,*)itheta,iphi,hisang(a,itheta,iphi)
	!end do
	

	!come out of time loop
	!end do

	DO a=1,3*n
	  WRITE(aa,*)a
	  aa=adjustl(aa)
	  OPEN(unit=25,file='anly_'//TRIM(aa)//'.dat',status='unknown')
	  WRITE(*,*)'Writing the trajectory for mode ',a
	  DO k=1,nfrs
	    WRITE(25,*)theta(a,k),phi(a,k),xim(a,k)
	  END DO
	  CLOSE(25)
	END DO

	!normalize
c	do i=1,n-1
c	lavm(i)=lavm(i)/(real(nfrs))
c	write(*,*)lavm(i)
c	end do

	!change to probability per solid angle
c	do a=1,n-1
c	do i=1,nbins/2-1
c	iphi=nint(nbins*sin(i*delr))
c	do j=0,iphi
c	hisang(a,i,j)=hisang(a,i,j)/(sin(i*delr)*delr*delr)
c	end do
c	end do
c	end do

	do a=1,3*n
	do i=1,nbins/2-1
	iphi=nint(nbins*sin(i*delr))
	delphi=(2.*pi)/(real(iphi))
c	write(*,*)"theta:",i*delha,"phibins:",iphi
	do j=1,iphi
	hnorm(a)=hnorm(a)+.5*delr*delphi*(hisang(a,i,j)
     &*sin(i*delr)+hisang(a,i-1,j-1)*sin((i-1)*delr))
	end do
	end do
	write(*,*)a,"norm:",hnorm(a)
	end do

	!normalize
	do a=1,3*n
	do i=0,nbins/2-1
	iphi=nint(nbins*sin(i*delr))
c	delphi=(2.*pi)/(real(iphi))
	do j=0,iphi
	hisang(a,i,j)=hisang(a,i,j)/hnorm(a)
	end do
	end do
	end do


	!test normalization
	do a=1,3*n
	testnorm=0.0
	do i=1,nbins/2-1
	iphi=nint(nbins*sin(i*delr))
	delphi=(2.*pi)/(real(iphi))
	do j=1,iphi
	testnorm=testnorm+.5*delr*delphi*(hisang(a,i,j)*
     &sin(i*delr)+hisang(a,i-1,j-1)*sin((i-1)*delr))
	end do
	end do
	write(*,*)a,"testnorm:",testnorm
	end do

	!prob volume
c	open(unit=16,file='pvol.dat')
c	do a=1,n-1
c	do i=2,nbins/2-1
c	do j=1,nbins-1
c	ir=nint(hisang(a,i,j)/dr)
c	if(ir.ne.0)then
c	do k=1,ir
c	pvol(a)=pvol(a)+.5*delr*delr*dr*(((k*dr)**2)
c     &*sin(i*delr)+(((k-1)*dr)**2)*sin((i-1)*delr))
c	end do
c	end if
c	end do
c	end do
c	write(*,*)a,"pvol:",pvol(a)
c	write(16,*)a,pvol(a)
c	end do

	femax=-Rb*T*log(1.0/real(nfrs))
	write(*,*)femax
	do a=1,3*n
	do i=1,nbins/2-1
	iphi=nint(nbins*sin(i*delr))
c	delphi=(2.*pi)/(real(iphi))
	do j=0,iphi
	if(hisang(a,i,j).ne.0.0)then !change to pmf
	feang(a,i,j)=-Rb*T*log(hisang(a,i,j))
	end if
	if(hisang(a,i,j).eq.0.0)feang(a,i,j)=femax
	end do
	end do
	end do

	write(cbins,*)nbins
	cbins=adjustl(cbins)
	open(unit=110,file="fricorr_mp_"//trim(cbins)//"_3N.dat")
	open(unit=111,file="feminmax_mp_"//trim(cbins)//"_3N.dat")
	open(unit=112,file="femin_mp_"//trim(cbins)//"_3N.dat")
	do a=1,3*n
	do i=1,nbins/2-1
	iphi=nint(nbins*sin(i*delr))
c	delphi=(2.*pi)/(real(iphi))
	do j=1,iphi-1
	if(feang(a,i,j).lt..5*femax)then
	 avfe(a)=avfe(a)+feang(a,i,j)
	 avfesq(a)=avfesq(a)+feang(a,i,j)**2
	 fenorm(a)=fenorm(a)+1.0
	  if(feang(a,i,j).le.fminp(a))then
	  fmin(a)=feang(a,i,j)
	  fminp(a)=feang(a,i,j)
	  end if
	  if(feang(a,i,j).gt.fmaxp(a))then
	  fmax(a)=feang(a,i,j)
	  fmaxp(a)=feang(a,i,j)
	  end if
	 end if
	end do
	end do
	avfe(a)=avfe(a)/(fenorm(a))
	avfesq(a)=avfesq(a)/(fenorm(a))
c	write(*,*)"avfe:",avfe(a),"1/4pi",-Rb*T*log(1./(4.*pi))
	write(*,*)"fenorm:",fenorm(a)
	sigfe(a)=(avfesq(a)-avfe(a)**2)**.5
	write(*,*)"?",avfesq(a),avfe(a)**2
	fricorr(a)=exp(sigfe(a)/(Rb*T))
	write(*,*)a,"efluc:",sigfe(a),"kcal/mol"
	write(110,*)sigfe(a)
	write(111,*)fmax(a)-fmin(a)
	write(112,*)fmin(a)
	end do
	close(110)
	close(111)
	close(112)

! typical barrier from ground from median absolute deviation
	open(unit=111,file="fmad_mp_"//trim(cbins)//"_3N.dat")
	do a=1,3*n
	  imad=0
	  fmad=0.0
	  fmadord=0.0
	  do i=1,nbins/2-1
	    iphi=nint(nbins*sin(i*delr))
c	delphi=(2.*pi)/(real(iphi))
	    do j=1,iphi-1
	      if(hisang(a,i,j) .gt. 2.*dc)then
	      imad=imad+1
	      fmad(imad)=abs(feang(a,i,j)-fmin(a))
	      end if
	    end do
	  end do

	imadmid=imad/2
	do i=1,imad
	fmadord(i)=maxval(fmad)
	fmad(maxloc(fmad))=0.0
	end do
	do i=1,imad
	write(*,*)i,"devE:",fmadord(i)
	end do
	emad=fmadord(imadmid)

	write(*,*)a,"mode mad:",emad,"kcal/mol"
	write(111,*)emad
	
	end do
	close(111)


!	open histogram files
	do i=1,3*n
	write(ii,*)i
	ii=adjustl(ii)
	open(unit=100+i,file='fempa_'//trim(ii)//'_3N.dat')
c	open(unit=200+i,file='prob_'//trim(ii)//'.dat')
	end do

!write 2D histograms
	do a=1,3*n
	do j=1,nbins/2-1
	iphi=nint(nbins*sin(j*delr))
	delphi=(360.)/(real(iphi))
	do k=0,iphi
c	if(feang(a,j,k).lt..5*femax)then
	write(100+a,*)j*delha,k*delphi,feang(a,j,k)
c	end if
c	x=hisang(a,j,k)*cos(k*delphi)*sin(j*delr)
c	y=hisang(a,j,k)*sin(k*delphi)*sin(j*delr)
c	z=hisang(a,j,k)*cos(j*delr)
c	write(200+a,*)j*delha,k*delphi,hisang(a,j,k)
c	write(200+a,*)j*delha,k*delha,hisang(a,j,k)
	end do
	write(100+a,*)
c	write(200+a,*)
	end do
	end do

	!open(unit=111,file="fmad_mp_"//trim(cbins)//"_3N.dat")
	!OPEN(unit=24,file='lambda_eig',status='old')
	!OPEN(unit=25,file='covar_eig_rescaled.dat',status='unknown')
	!DO a=1,3*n
	!  READ(24,*)lambda(a)
	!  READ(111,*)bar(a)
	!  WRITE(*,*)bar(a)
	!  bar(a)=exp(bar(a)/(Rb*T))
        !  WRITE(25,*)a,lambda(a)*bar(a)
	!  WRITE(*,*)lambda(a)*bar(a)
	!END DO
	!CLOSE(111)
	!CLOSE(24)
	!CLOSE(25)
	!Write mode trajectories to file
	DO a=1,3*n !Only the first ten for now
	  WRITE(aa,*)a
	  aa=adjustl(aa)
	  OPEN(unit=24,file='xi_'//TRIM(aa)//'_aniso.dat')
	  DO k=1,nfrs
            WRITE(24,*)0.2*k,xix(a,k),xiy(a,k),xiz(a,k),xim(a,k)
	  END DO
	  CLOSE(24)
	END DO
	end subroutine

