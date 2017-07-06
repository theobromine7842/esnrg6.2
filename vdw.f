*
* Subroutine to assign vdw radii.
* (Actually, sigma of LJ expression).
* This sigma obtained from NBOND subroutine.
*
	SUBROUTINE vdw

**
* Set Up Variables
**
*
*
c         integer atmax,nbmax,grpdim,a1dim,a2dim
c         parameter (atmax=50000,nbmax=125,grpdim=1000)
c         parameter (a1dim=200,a2dim=75)

	INCLUDE 'dim.inc'

	real epsilon(atmax),reps(atmax),atrad(atmax)
	real nbterms(nbmax,3)
      real x(atmax),y(atmax),z(atmax)

      real acharge(a1dim,a2dim),q(atmax),gcharge(a1dim,a2dim)

      integer i,j,k,atlen,Pwc,Awc
      integer atnum,resnum,sresnum
	integer numnbtype,nbindex(nbmax,2)
      integer m,n,nres,atomres(a1dim),grp(atmax),numgrp
      integer grpatom(a1dim,a2dim),atomgrp(atmax)
      integer gtype(atmax),grpnum(a1dim),resgrp(atmax)
	integer atradtyp(atmax)

	character*4 nbtype(nbmax)
      character*4 atid,resid,segid
	character*4 akind(a1dim,a2dim),atype(atmax)
      character*4 tatid(a1dim,a2dim),tresid(a1dim)
	character*30 grptitle(12)
	character*80 paramfile

      logical EOF_flag,FIND_flag

c Variables for timing functions
      real*4 timer,crap,agtime,vtime,ttime

	common /GROUP/nres,atomres,grp,numgrp,grpatom,atomgrp,
     +   gtype,grpnum,resgrp,acharge,q,gcharge,grptitle,tatid,
     +   tresid,akind,atype
	common /COORD/atnum,resnum,sresnum,x,y,z,resid,atid,segid
	common /NBONDB/numnbtype,nbindex,nbterms,nbtype,paramfile
	common /VDWB/atradtyp,epsilon,reps,atrad
      common /TIMING/agtime,vtime,ttime

c	print *, 'Starting VDW'
      
**
* Setup any variables necessary
**
	EOF_flag=.FALSE.
	FIND_flag=.FALSE.
	atlen=0
**
* Search nb arrays for atom type
**
	Do 1000 i=1,numnbtype

c	print *, 'Im now in VDW'

	  if (FIND_flag) goto 1000
	  Pwc=nbindex(i,1)
	  Awc=nbindex(i,2)

c	write(6,*) Pwc,Awc,nbtype(i),atype(atnum)

	   if ((Pwc.le.0).and.(Awc.le.0)) then
	    if (nbtype(i).eq.atype(atnum)) FIND_flag=.TRUE.
	   else
	    if (Pwc.gt.0) then
	     atlen=index(atype(atnum),' ') - 1
	     if (atlen.lt.1) atlen=len(atype(atnum))
	     if ((nbtype(i)(1:Pwc-1).eq.atype(atnum)(1:Pwc-1)).and.
     +          (atlen.eq.Pwc)) FIND_flag=.TRUE.
*	print *, 'Pwc=',Pwc
	    endif
	    if (Awc.gt.0) then
	     if (nbtype(i)(1:Awc-1).eq.atype(atnum)(1:Awc-1))
     +        FIND_flag=.TRUE.
	    endif
	   endif
	   if (.not.FIND_flag) goto 1000
	   epsilon(atnum)=nbterms(i,1)
	   reps(atnum)=nbterms(i,2)
	   atrad(atnum)=nbterms(i,3)
	   atradtyp(atnum)=i


c	write(6,*) '  eps=',epsilon(atnum)
c	write(6,*) ' reps=',reps(atnum)
c	write(6,*) 'atrad=',atrad(atnum)
c	write(6,*) 'atradtyp=',atradtyp(atnum)

1000    Continue

*
* Return to main program
*
       vtime=vtime

	 Return
	 End
