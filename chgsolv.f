*
* This subroutine calculates the solvent interaction energy
* at one charge in a system of distributed charges.
*
	SUBROUTINE CHGSOLV(Ntot,Nchg,coi,catnum,pbin,tbin,dbin,abin,
     +   eps,epsp,rmin,rmax,rbin,T,pi,q,xc,yc,zc,x,y,z,rprobe,atrad,
     +   totchg,check_vdw_flag,bphi,brphi)
*
* It loops over all space (r,theta,phi) and at each
* point loops over two angles which define the orientation
* of a solvent dipole in a cavity in the fields of the
* charges. The probability of any given orientation is calculated
* using exp(-beta*V)/sum[exp(-beta*V)] where Vi is the electrostatic
* interaction energy between the charge qi and the hypothetical
* water of dipole moment, mu, at a point r,theta,phi. This probability
* is then used to calculated the average interaction angle between
* the charge of interest and the r,theta,phi point. This angle is
* used to calculate the interaction energy between the point and the
* charge of interest and these energies are summed to give the solvation
* energy with the charge of interest.
*
*

	integer maxrbin,maxabin,maxatom,maxchg
	real eA2Debye,Debye2eA
	real kb,T
	parameter (maxrbin=500,maxabin=100,maxatom=50000)
	parameter (maxchg=2000,eA2Debye=4.796,Debye2eA=0.2085)
c	parameter (kb=1.98717E-03,T=300.)
	parameter (kb=1.98717E-03)

      integer i,j,k,m,n,l
	integer Ntot,Nchg,coi
	integer rstart,rstop,tstart,tstop,pstart,pstop
	integer dstart,dstop,astart,astop

	integer catnum

      real q(maxchg),eps,epsp
      real xc(maxchg),yc(maxchg),zc(maxchg)
	real totchg,qcoi
	real rmin,rmax
	real dx,dy,dz,rd,rdsq
	real r,rx,ry,rz,r2,r2sq,rbin,deltar,rbsf,invr2,invr2sq
	real r1,r1sq,invr1,invr1sq,r3,r3sq,Ubshell,Uborn
	real r2cb,r3cb,drcb,cost2,cost3,dcost
	real theta,theta2,theta3,tbin,deltat,tbsf,sintheta,costheta
	real phi,phi2,phi3,pbin,deltap,pbsf,sinphi,cosphi
	real alpha,alpha2,abin,deltaa,absf,sinalpha,cosalpha
	real delta,delta2,dbin,deltad,dbsf,sindelta,cosdelta
	real u,ux,uy,uz
	real cosgamma(maxchg)
	real V,sumV
	real f(maxabin,maxabin),sumf,sumcgf,P(maxabin,maxabin)
	real dcosgda,dcosgdd,dcosg,avgcosg
	real gr,rho,kconv,Energy,Upt,Usolv,tsf
	real dVol,Vols,Volr,Volt
	real ratio,ratior
	real Potential, Ppt, Pptr, Psolv, Pbshell, Pborn
	real Pratio, Pratior
	real xcoi,ycoi,zcoi,Uptr

	real rprobe,atrad(maxatom)
	real bphi,brphi(maxrbin)
	real x(maxatom),y(maxatom),z(maxatom)
      real xcen(maxatom),ycen(maxatom),zcen(maxatom)

	real pi,twopi,fourpi,negpi,beta,deg2rad

	character*80 nrgfile

	logical overlap

	logical check_vdw_flag(maxatom)


c	print *, 'Ntot=',Ntot,'Nchg=',Nchg,'coi=',coi,'catnum=',catnum
c	print *,'pbin=',pbin,'tbin=',tbin,'dbin=',dbin,'abin=',abin,
c	print *,'eps=',eps,'epsp=',epsp
c	print *,'rmin=',rmin,'rmax=',rmax,'rbin=',rbin,'T=',T,'pi=',pi
*c	print *,'q=',q,'xc=',xc,'yc=',yc,'zc=',zc,'x=',x,'y=',y,'z=',z
*c	print *,'rprobe=',rprobe,'atrad=',atrad,'totchg=',totchg
c	print *,'rprobe=',rprobe,'totchg=',totchg
*c	print *,'check_vdw_flag=',check_vdw_flag,'bphi=',bphi,'brphi=',brphi
c	print *,'bphi=',bphi

*
* Assign constants and initialize variables
*
c	pi=acos(-1.)
	twopi=2.*pi
	fourpi=4.*pi
	negpi=-1.*pi
	beta=1./(kb*T)
	deltar=1./rbin
	rbsf=1.0
	deltat=180./tbin
	tbsf=1.0
	deltap=360./pbin
	pbsf=1.0
	deltaa=180./abin
	absf=1.0
	deltad=360./dbin
	dbsf=1.0
	u=-2.35*Debye2eA 	! dipole that points towards negative end
	V=0.
	sumV=0.
	gr=1.
	rho=0.0334
	kconv=332.
	Uborn=0.
	Ubshell=0.
	ratior=0.
	ratio=0.
	Potential=0.
	Ppt=0.
	Pptr=0.
	Psolv=0.
	Pbshell=0.
	Pborn=0.
	Pratio=0.
	Pratior=0.
	deg2rad=pi/180.
*
* open data files
*
        open(10,form='formatted',access='sequential',
     +  file='escsolv.dat',status='unknown')
        open(11,form='formatted',access='sequential',
     +  file='escsolvp.dat',status='unknown')
        open(12,form='formatted',access='sequential',
     +  file='escsolvv.dat',status='unknown')

*
* read input
*
	qcoi = q(coi)
c	print *, coi,qcoi,xc(coi),yc(coi),zc(coi)
	write(6,5) coi,qcoi,xc(coi),yc(coi),zc(coi)
5	format(5x,
     +   /'The index of the central charge is ',I5,
     +   /'The charge of the central charge is ',F10.5,
     +   /'The coordinates of the central charge are ',3(F10.5,1x))

* Now loop over r from rmin to rmax using bins/steps of
* size 1/2 deltar for the first and last bin and
* deltar for the remaining bins. 

	rstart=nint(rbin*rmin)
	rstop =nint(rbin*rmax)
	
	do 100 i=rstart,rstop
	 r=float(i)*deltar
	 r1=r
	 r1sq=r**2
	 invr1=1./r
	 invr1sq=1./r1sq
	 r2=(float(i)+0.5)*deltar
	 if (r2.gt.rmax) r2=rmax
	 r2sq=r2**2
	 r2cb=r2sq*r2
	 invr2=1./r2
	 invr2sq=1./r2sq
	 r3=(float(i)-0.5)*deltar
	 if (r3.lt.rmin) r3=rmin
	 r3sq=r3**2
	 r3cb=r3sq*r3
	 drcb=(r2cb-r3cb)/3

* Write a . for every r step
	if (mod(i,int(rbin)).ne.1) then
	 if (i.eq.rstart) then
          write(6,'(5x,A,$)') '.'
	 elseif (mod(i,int(rbin)).ne.0) then
	  write(6,'(A,$)') '.'
	 else
	  write(6,'(A,2x,F9.5,A)') '. ',r,' A'
	 endif
	else
	 write(6,'(5x,A,$)') '.'
	endif

* Re-initialize shell dependent variables
	Uptr = 0.		! reinitialize for each r
	Pptr = 0.		! reinitialize for each r
	Vols = 0.		! reinitialize for each r
	Volr = 0.		! reinitialize for each r
	Volt = 0.		! reinitialize for each r

c	write(6,'(//   " *=*=*=* "   //)')
c
c	print *, 'r=',r,'r2=',r2,'r3=',r3,'r scaling factor=',rbsf
c	print *, 'invr2=',invr2,'invr2sq=',invr2sq

* Now loop over theta from 0 to pi and phi from 0 to 2pi
* using bins/steps of size 1/2 deltat and 1/2 deltap for
* the first and last bins and deltat and deltap for the
* remaining bins.

	 pstart=0
c	 pstop =nint(pbin*twopi)-1
	 pstop =nint(pbin)-1

	 do 200 j=pstart,pstop
	  phi=float(j)*deltap
c	  phi2=(float(j)+0.5)*deltap
c	  phi3=(float(j)-0.5)*deltap
c	  if (phi2.gt.twopi) phi2=twopi
	  sinphi=sin(phi*deg2rad)
	  cosphi=cos(phi*deg2rad)

	  tstart=0
c	  tstop =nint(tbin*pi)
	  tstop =nint(tbin)

	  do 300 k=tstart,tstop
	   theta=float(k)*deltat
	   theta2=(float(k)+0.5)*deltat
	   theta3=(float(k)-0.5)*deltat
c	   if (theta2.gt.pi) theta2=pi
	   if (theta2.gt.180.) theta2=180.
	   if (theta3.lt.0.) theta3=0.
	  sintheta=sin(theta*deg2rad)
	  costheta=cos(theta*deg2rad)
	  cost2=cos(theta2*deg2rad)
	  cost3=cos(theta3*deg2rad)
	  dcost=cost3-cost2

* Calculate the xyz components of r
	   rx=r*sintheta*cosphi
	   ry=r*sintheta*sinphi
	   rz=r*costheta

* Initialize the sumf for each r,theta,phi loop
	   sumf=0.
	   sumcgf=0.
	   avgcosg=0.
	   dVol=drcb*dcost*deltap*deg2rad
	   Volt=Volt + dVol

* Check to see if this point overlaps the solute
* by looping over atoms and making sure point is
* greater than rprobe + atrad A away from all atoms.
*

	   overlap=.FALSE.

	   do 350 l=1,Ntot

	    if (check_vdw_flag(l)) then
* Calculate the distance between the dipole and charge l.
	      xcen(l)=x(l)-x(catnum)
	      ycen(l)=y(l)-y(catnum)
	      zcen(l)=z(l)-z(catnum)
	      dx=(rx-xcen(l))
	      dy=(ry-ycen(l))
	      dz=(rz-zcen(l))
	      rdsq=(dx**2) + (dy**2) + (dz**2)
	      rd=sqrt(rdsq) - rprobe - atrad(l)
	      if (rd.lt.0) then
	       overlap=.TRUE.
c	       print *, 'Atom ',l,' overlaps this position.'
c		print *, sqrt(rdsq),rprobe,atrad(l)
	      endif
	     endif
350	    enddo

*
* If this point overlaps any atom then jump to end
* of theta loop
*
	     if (overlap) then
c		print *, 'Skip remainder of loop...'
		goto 300
	     endif

	     Vols=Vols + dVol


* Now loop over alpha from 0 to pi and delta from 0 to 2pi
* using bins/steps of size 1/2 deltaa and 1/2 deltad for
* the first and last bins and deltaa and deltad for the
* remaining bins.

	   dstart=0
c	   dstop =nint(dbin*twopi)-1
	   dstop =nint(dbin)-1

	   do 400 m=dstart,dstop
	    delta=float(m)*deltad
c	    delta2=(float(m)+0.5)*deltad
c	    if (delta2.gt.twopi) delta2=twopi
	    sindelta=sin(delta*deg2rad)
	    cosdelta=cos(delta*deg2rad)

  	    astart=0
c  	    astop =nint(abin*pi)
  	    astop =nint(abin)

	    do 500 n=astart,astop
	     alpha=float(n)*deltaa
c	     alpha2=(float(n)+0.5)*deltaa
c	     if (alpha2.gt.pi) alpha2=pi
c	     if (alpha.gt.pi) alpha=pi
	     if (alpha.gt.180.) alpha=180.
	     sinalpha=sin(alpha*deg2rad)
	     cosalpha=cos(alpha*deg2rad)

* Calculate the xyz components of the dipole moment
	     ux=sinalpha*cosdelta
	     uy=sinalpha*sindelta
	     uz=cosalpha

* Now loop over all charges.

	     sumV=0.

c	     do 600 l=1,Nchg
	     do 600 l=1,coi
	    
* Calculate the angle between the dipole vector and the vector
* between the dipole and lth charge.
	      dx=(rx-xc(l))
	      dy=(ry-yc(l))
	      dz=(rz-zc(l))
	      rdsq=(dx**2) + (dy**2) + (dz**2)
	      rd=sqrt(rdsq)
	      cosgamma(l)=((ux*dx)+(uy*dy)+(uz*dz))/rd
c	print *, 'rdsq=',rdsq,'rd=',rd,'cosgamma=',cosgamma(l)

* Calculate the electrostatic energy between the solvent
* dipole and the current charge as well as the sum of 
* the energy over all charges.

	     if (l.ne.coi) then
	      V=(kconv*u*q(l)*cosgamma(l))/(epsp*rdsq)
	      sumV=sumV + V

c	print *, 'cosgamma=',cosgamma(l),'Energy=',V

c	write(6,150) r,phi,theta,delta,alpha,l,q(l),qcoi
150	format (5(f9.5,1x),i5,1x,3(f9.5,1x))
	     endif
600	     enddo

* Now calculate Probability,P(alpha,delta)
* via P=f/sumf where f=exp(-BV) and B=1./kbT
* If any r,theta,phi point is inside the protein, then
* set f(alpha,delta)=0.

	      f(alpha,delta)=exp(-1.*beta*sumV)

* Now, since avg(M)=SUM_i(Pi*Mi)
* and Pi = fi/SUM_i(f)
* then avgcosg = [SUM_i(fi * cosg(i))]/[SUM_i(fi)]

	     sumf=sumf + f(alpha,delta)
	     sumcgf=sumcgf + (cosgamma(coi)*f(alpha,delta))

c	print *, 'f=',f(alpha,delta),'sumf=',sumf,
c     +    'cosg=',cosgamma(coi),'sumcgf=',sumcgf
	     
500	    enddo
400	   enddo
	
           avgcosg=sumcgf/sumf

c	   print *, 'avgcosg=',avgcosg
c	   write (11,450) r,theta,phi,avgcosg
450        format (4(1x,F12.5))

* Now calculate the contribution this r,theta,phi position
* has to the total solvent-chg_of_interest interaction energy
           Energy= (kconv*u*qcoi*avgcosg)*invr1sq
           Potential= (kconv*u*avgcosg)*invr1sq
c	print *, kconv,u,coi,qcoi,avgcosg,invr2sq
c	   tsf   = rbsf*tbsf*pbsf
c	   dVol  = r2sq*deltar*sintheta*deltat*deltap*tsf
c	   dVol  = r1sq*deltar*sintheta*deltat*deltap*tsf
	   Upt   = rho*gr*Energy*dVol
	   Uptr  = Uptr  + Upt
	   Usolv = Usolv + Upt
	   Ppt   = rho*gr*Potential*dVol
	   Pptr  = Pptr  + Ppt
	   Psolv = Psolv + Ppt
300	  enddo
200      enddo
c*	print *, 'totchg=',totchg
	 Volr=Vols/Volt
c* 	 Ubshell = -1.*kconv*(1.-(1./eps))*(qcoi**2)*((1./r3)-(1./r2))
	 Ubshell = -1.*kconv*Volr*(1.-(1./eps))*(totchg**2)*((1./r3)-
     +		(1./r2))
	 Uborn = Uborn + Ubshell
	 Pbshell = -1.*kconv*Volr*(1.-(1./eps))*(totchg)*((1./r3)-(1./r2))
	 Pborn = Pborn + Pbshell
c	 print *, '----------------------------------------'
c	 print *, 'r=',r,'rin=',r3,'rout=',r2
c	 print *,'Usolv(r)=',Uptr,'Usolv=',Usolv
c	 print *, '----------------------------------------'
	 brphi(i)=Uptr
	 bphi=Usolv
	 ratior =Ubshell/Uptr
	 ratio  =Uborn/Usolv
	 Pratior=Pbshell/Pptr
	 Pratio =Pborn/Psolv
	 write (10,250) r,Uptr,Usolv,Ubshell,Uborn,ratior,ratio,Volr
	 write (11,250) r,Pptr,Psolv,Pbshell,Pborn,Pratior,Pratio,Volr
	 write (12,250) r,Vols,Volt,Volr,fourpi*drcb,Vols/(fourpi*drcb)
250      format (8(1x,F12.5))
100     enddo

	write(6,850) 'The solvent interaction energy ',
     +   ' at charge ',coi,' is ',Usolv
	write(6,850) 'The Born solvent interaction energy ',
     +   ' at charge ',coi,' is ',Uborn
850	format (A,A,I5,A,F10.5)

9999	return
	END
