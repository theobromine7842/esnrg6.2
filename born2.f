	SUBROUTINE BORN2

* This subroutine calculates the BORN solvation energy
* of a system of charge 'Bchg' beyound a cutoff 'rborn' from the
* central position. In regions containing non-solvent atoms, this
* value is scaled by (1-(Non-solvent volume/Total volume)) for any
* given radial shell of volume.
* BORN2 uses rectangular volume elements
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
	real xmin,ymin,zmin,xmax,ymax,zmax
	real d0,div

	real Bchg,Bchg2

	integer numatom,n,c,rstrt,rstop,fstrt,fstop
	integer tstrt,tstop,i,j,k,f,t,irborn
	integer atnum,resnum,sresnum,cresnum,catnum
	integer atradtyp(atmax)
	integer xstrt,ystrt,zstrt,xstop,ystop,zstop
      integer IDCLR,idiv

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
	 arange = atrad(n)
	 xmin = xc - arange
	 xmax = xc + arange
	 ymin = yc - arange
	 ymax = yc + arange
	 zmin = zc - arange
	 zmax = zc + arange
	 rmin = sqrt(xmin**2 + ymin**2 + zmin**2)
	 rmax2 = sqrt(xmax**2 + ymax**2 + zmax**2)
c	 rmin = rc - arange
c	 rmax2 = rc + arange
c	 if (rmin.lt.rborn) rmin  = rborn
c	 if (rmax2.gt.rmax) rmax2 = rmax
*
* Now make sure you don't do the volume
* calculation more once per atom type
*

*!!!! This statement forces redundant checking
	 atvolchk(atradtyp(n))=.FALSE.
*!!!!
	 if (.not. atvolchk(atradtyp(n))) then
	 xstrt = nint(xmin*div)
	 xstop = nint(xmax*div)
	 ystrt = nint(ymin*div)
	 ystop = nint(ymax*div)
	 zstrt = nint(zmin*div)
	 zstop = nint(zmax*div)

	 do 151 i=xstrt,xstop
*	  write(6,*) 'i=',i
	  x0 = real(i)*drad
	  do 152 f=ystrt,ystop
*	   write(6,*) ' f=',f
	   y0 = real(f)*drad
	   do 153 t=zstrt,zstop
*	    write(6,*) '  t=',t
	    z0 = real(t)*drad
	    dvol = drad**3
	    dx = x0 - xc
	    dy = y0 - yc
	    dz = z0 - zc
	    d  = sqrt(dx**2 + dy**2 + dz**2)
	    d0 = sqrt(x0**2 + y0**2 + z0**2)
	    r  = nint(d0*div)

c	    write(6,*) 'd=',d,'rad=',atrad(n)

	    if (d.le.(atrad(n))) then
	      r3 = nint((d0 - rmin)*div)
	      atrvol(atradtyp(n),r3) =
     +         atrvol(atradtyp(n),r3) + dvol
	      atvolchk(atradtyp(n))=.TRUE.
	     if ((d0 .ge. rborn).and.(d0 .le. rmax)) then
	      pvol(r) = pvol(r) + dvol
	      avol(n) = avol(n) + dvol
*	      write(6,*) 'Atom #',n,' excluded.'
	     endif
	    endif
153	   continue
152	  continue
	  write(6,*) d0,pvol(r),atrvol(atradtyp(n),r3)
151	 continue

	 else

	 rstrt = nint(rmin*div)
	 rstop = nint(rmax2*div)
	 do 1151 i=rstrt,rstop
*	  write(6,*) 'i=',i
	  r = float(i)
	  r3 = r - rstrt
	  d0 = real(r)*drad
	  if ((d0 .ge. rborn).and.(d0 .le. rmax)) then
	   pvol(r) = pvol(r) + atrvol(atradtyp(n),r3)
	   avol(n) = avol(n) + atrvol(atradtyp(n),r3)
	  endif
	  write(6,*) 'Atom #',n,' excluded.'
	  write(6,*) d0,pvol(r),atrvol(atradtyp(n),r3)
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
