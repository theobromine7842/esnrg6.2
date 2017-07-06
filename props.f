        SUBROUTINE PROPS(Cchg)
*
* Subroutine to assign certain properties of the central residue
* or atom:
*

c         integer atmax,grpdim,a1dim,a2dim
c         parameter (atmax=50000,grpdim=1000)
c         parameter (a1dim=200,a2dim=75)

	INCLUDE 'dim.inc'

	real x(atmax),y(atmax),z(atmax),cenx,ceny,cenz,atpr

      real acharge(a1dim,a2dim),q(atmax),gcharge(a1dim,a2dim),Cchg

	integer i,j,k
      integer n,catnum,sresnum,cresnum,resnum,cresat(atmax)
	integer iatpr,atnum
      integer m,nres,atomres(a1dim),grp(atmax),numgrp
      integer grpatom(a1dim,a2dim),atomgrp(atmax)
      integer gtype(atmax),grpnum(a1dim),resgrp(atmax)
      integer IDCLR

      character*4 atid,catid,segid,cseg,resid,cresid
      character*4 atname(atmax),resname(atmax),AorP
	character*4 segname(atmax)
	character*1 XYZ_QUEST
	character*30 grptitle(12)
	character*4 akind(a1dim,a2dim),atype(atmax)
      character*4 tatid(a1dim,a2dim),tresid(a1dim)
	
	common /COORD/atnum,resnum,sresnum,x,y,z,resid,atid,segid
	common /AVG/cresat,iatpr
	common /IDB/n,cresnum,catnum,IDCLR,cenx,ceny,cenz,atpr,
     +   atname,resname,cseg,cresid,catid,XYZ_QUEST,AorP,segname
	common /GROUP/nres,atomres,grp,numgrp,grpatom,atomgrp,
     +   gtype,grpnum,resgrp,acharge,q,gcharge,grptitle,tatid,
     +   tresid,akind,atype


	do 10 i=cresat(1),cresat(iatpr)
c	write(6,9) i,segname(i),resname(i),atname(i)
c9	format (5x,I5,3(1x,A))

c        if (((segid.eq.cseg).and.(sresnum.eq.cresnum)).and.
c     +   ((atid.eq.catid).or.(XYZ_QUEST.eq.'y'))) then
	  Cchg=Cchg + q(i)
c	print *, 'Cchg=',Cchg
c        endif
10	enddo

        Return

        End
