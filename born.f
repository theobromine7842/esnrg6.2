	SUBROUTINE BORN

* This subroutine calculates the BORN solvation energy
* of a system of charge 'Bchg' beyound a cutoff 'rborn' from the
* central position. In regions containing non-solvent atoms, this
* value is scaled by (1-(Non-solvent volume/Total volume)) for any
* given radial shell of volume.
*
c         integer atmax,nbmax,grpdim
c         parameter (atmax=50000,nbmax=125,grpdim=1000)

	INCLUDE 'dim.inc'

	real ddeg,dangle,drad,dfee,dtheta
	real x(atmax),y(atmax),z(atmax),xc,yc,zc,rc
	real x0,y0,z0,rborn,r,dx,dy,dz,d
	real rmin,rmax,rmax2
	real fee,theta,arange,fmin,fmax,tmin,tmax
	real pi,temp,dvol,pvol(grpdim),rvol(grpdim),avol(atmax)
	real br(grpdim),r1,r2,r3,brphi(grpdim),bphi
	real dielec
	real epsilon(atmax),reps(atmax),atrad(atmax)
	real cenx,ceny,cenz,atpr
	real atrvol(nbmax,grpdim)
	real Bchg,Bchg2
      real div

	integer numatom,n,c,rstrt,rstop,fstrt,fstop
	integer tstrt,tstop,i,j,k,f,t,irborn
	integer atnum,resnum,sresnum,cresnum,catnum
	integer atradtyp(atmax)
      integer IDCLR
      integer idiv

	character*4 resid,atid,segid
	character*4 atname(atmax),resname(atmax)
	character*4 cseg,cresid,catid
	character*1 XYZ_QUEST,AorP
	character*4 segname(atmax),solvseg

	logical atvolchk(nbmax)

	common /COORD/atnum,resnum,sresnum,x,y,z,resid,atid,segid
	common /VDWB/atradtyp,epsilon,reps,atrad
	common /IDB/n,cresnum,catnum,IDCLR,cenx,ceny,cenz,atpr,
     +   atname,resname,cseg,cresid,catid,XYZ_QUEST,AorP,segname
	common /BORNB/irborn,rborn,rmax,pvol,avol,rvol,brphi,bphi,
     +   Bchg,dielec,atrvol,atvolchk,solvseg
	common /CONST/idiv,pi,div

*
* Calculate excluded volume of protein for Born approximation:
*
	write(6,*) 'Calculating excluded volume of protein...'
	numatom = n
	c = catnum
	ddeg = 0.5
	dangle = (180./pi)/ddeg
	drad = 1./div
	dfee = 1./dangle
	dtheta = dfee

c	print *, 'atrad=',atrad(n)

	do 150 n=1,numatom
	 if (segname(n).eq.solvseg) goto 149
*	 write(6,*) n,segname(n),solvseg
	 xc = x(n) - x(c)
	 yc = y(n) - y(c)
	 zc = z(n) - z(c)
	 rc = sqrt(xc**2 + yc**2 + zc**2)
	 if (rc.lt.rborn) goto 149
*	 write(6,*) n,rc
	 rmin = rc - atrad(n)
	 rmax2 = rc + atrad(n)
	 if (rmin.lt.rborn) rmin = rborn
	 if (rmax2.gt.rmax) rmax2 = rmax
*
* Now make sure you don't do the volume
* calculation more once per atom type
*

*!!!! This statement forces redundant checking
	 atvolchk(atradtyp(n))=.FALSE.
*!!!!
	 if (.not. atvolchk(atradtyp(n))) then
	 fee = acos(zc/rc)
	 theta = acos(xc/(rc*sin(fee)))
	 arange = acos(1.-(atrad(n)**2)/(2.*rc**2))
	 fmin = fee - arange
	 fmax = fee + arange
	 tmin = theta - arange
	 tmax = theta + arange
         if ((tmax-tmin).gt.pi) then
          temp = tmin
          tmin = tmax
          tmax = temp + 2.*pi
         endif
	 rstrt = nint(rmin*div)
	 rstop = nint(rmax2*div)
	 fstrt = nint(fmin*dangle)
	 fstop = nint(fmax*dangle)
	 tstrt = nint(tmin*dangle)
	 tstop = nint(tmax*dangle)
*	 write(6,*) n,segname(n),' - ',atname(n),rc
*	 write(6,*) rmin,rmax2
*	 write(6,*) fmin,fmax
*	 write(6,*) tmin,tmax
*	 write(6,*) rstrt,rstop
*	 write(6,*) fstrt,fstop
*	 write(6,*) tstrt,tstop
	 do 151 i=rstrt,rstop
*	  write(6,*) 'i=',i
	  r = real(i)*drad
	  r3 = i-rstrt+1
	  do 152 f=fstrt,fstop
*	   write(6,*) ' f=',f
	   fee = real(f)/dangle
	   z0 = r*cos(fee)
	   do 153 t=tstrt,tstop
*	    write(6,*) '  t=',t
	    theta = real(t)/dangle
	    x0 = r*(cos(theta))*(sin(fee))
	    y0 = r*(sin(theta))*(sin(fee))
            dvol = dtheta*(cos(fee)-cos(fee+dfee))*
     +       ((r+.5*drad)**3-(r-.5*drad)**3)/3.
*	    write(6,*) dtheta,fee,dfee,r,drad 
*	    write(6,*) 'dvol = ',dvol
	    dx = xc - x0
	    dy = yc - y0
	    dz = zc - z0
	    d = sqrt(dx**2 + dy**2 + dz**2)
*	    write(6,*) r,fee,theta
*	    write(6,*) dx,dy,dz,d
*	    write(6,*) xc,yc,zc
*	    write(6,*) x0,y0,z0

c	    write(6,*) 'd=',d,'rad=',atrad(n)

	    if (d.le.(atrad(n))) then
	     pvol(i) = pvol(i) + dvol
	     avol(n) = avol(n) + dvol
	     atrvol(atradtyp(n),r3) =
     +        atrvol(atradtyp(n),r3) + dvol
	     atvolchk(atradtyp(n))=.TRUE.
*	     write(6,*) 'Atom #',n,' excluded.'
	    endif
153	   continue
152	  continue
c	  write(6,*) r,pvol(i),atrvol(atradtyp(n),r3)
151	 continue

	 else

	 rstrt = nint(rmin*div)
	 rstop = nint(rmax2*div)
	 do 1151 i=rstrt,rstop
*	  write(6,*) 'i=',i
	  r = real(i)*drad
	  r3 = i-rstrt+1
	  pvol(i) = pvol(i) + atrvol(atradtyp(n),r3)
	  avol(n) = avol(n) + atrvol(atradtyp(n),r3)
c	  write(6,*) 'Atom #',n,' excluded.'
c	  write(6,*) r,pvol(i),atrvol(atradtyp(n),r3)
1151	 continue

	 endif

*	 write(6,*) 'Atom #',n,' excluded volume = ',avol(n)

149	if (n.eq.1) then
	 write(6,'(5x,A)') 'Analyzing atomic volume ...'
	endif
	if (mod(n,numatom).eq.0) then
	 write(6,'(A,2x,I4,A)') '.',n,' atoms'
	else
	if (mod(n,50).ne.1) then
	 if (mod(n,50).ne.0) then
	  write(6,'(A,$)') '.'
	 else
	  write(6,'(A,2x,I4,A)') '.',n,' atoms'
	 endif
	else
	 write(6,'(5x,A,$)') '.'
	endif
	endif

150	continue
*
* Calculate Born approximation first:
*
	write(6,*) 'Calculating Born energies...'
	write(6,*)
	Bchg2=Bchg*Bchg
c	  write(6,236) 'dist.','non-solv vol','ttl solv vol',
c     +     '1-non/ttl'
        open(4,form='formatted',access='sequential',
     +  file='esborn.dat',status='new')
	do 160 i=1,idiv
	 br(i) = 1.
	 r1 = anint(float(i))/div
	 r2 = anint(float(i+1))/div
	 if (i.ge.irborn) then
	  rvol(i) = 4.*pi*((r1+.5*drad)**3 - (r1-.5*drad)**3)/3.
	  br(i) = 1. - (pvol(i)/rvol(i))
c
c	  write(6,237) r1,pvol(i),rvol(i),br(i)
236	  format (4(A10,1x))
237	  format (4(F10.5,1x))
c
	  brphi(i) = -332.*Bchg2*(1.- (1./dielec))*
     +     (1/r1 - 1/r2)*br(i)
	 else
	  brphi(i) = 0.
	 endif
	 if (i.eq.idiv) brphi(i) = -332.*Bchg2*(1.- (1./dielec))*
     +    (1/r1)*br(i)
	 bphi = bphi + brphi(i)
	 write(4,917) r1,br(i),brphi(i),bphi
917	 format (4(F12.5,1x))
160	continue

*
* Return to main program
*
	 Return
	 End
