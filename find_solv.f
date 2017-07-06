* Detm. if an atom is a member of a solvent residue outside rborn
	SUBROUTINE FIND_SOLV(i,j,k,segname2)
        
c         integer atmax,nbmax,maxchg,grpdim,a1dim,a2dim
c         parameter (atmax=50000,nbmax=125,maxchg=1000,grpdim=1000)
c         parameter (a1dim=200,a2dim=75)

	INCLUDE 'dim.inc'

	real x(atmax),y(atmax),z(atmax),xc,yc,zc,rc
	real rmin,rmax,rborn
	real pi,dvol,pvol(grpdim),rvol(grpdim),avol(atmax)
	real brphi(grpdim),bphi
	real dielec,dielecp
	real epsilon(atmax),reps(atmax),atrad(atmax)
	real cenx,ceny,cenz,atpr
	real atrvol(nbmax,grpdim)
	real xqc(maxchg),yqc(maxchg),zqc(maxchg)
	real qqc(maxchg),rprobe
	real xsum,ysum,zsum,rsolv
	real Bchg,Bchg2,Cchg
      real acharge(a1dim,a2dim),q(atmax),gcharge(a1dim,a2dim)
      real div

	integer numatom,c
	integer i,j,k,l,f,t,irborn
	integer atnum,resnum,sresnum,cresnum,catnum
	integer atradtyp(atmax)
	integer Nchg,coi
      integer m,n,nres,atomres(a1dim),grp(atmax),numgrp
      integer grpatom(a1dim,a2dim),atomgrp(atmax)
      integer gtype(atmax),grpnum(a1dim),resgrp(atmax)
      integer IDCLR,idiv

	character*4 resid,atid,segid
	character*4 atname(atmax),resname(atmax)
	character*4 cseg,cresid,catid
	character*1 XYZ_QUEST,AorP
	character*30 grptitle(12)
	character*4 segname(atmax),solvseg,segname2(atmax)
      character*4 akind(a1dim,a2dim),atype(atmax)
      character*4 tatid(a1dim,a2dim),tresid(a1dim)

	logical atvolchk(nbmax)
	logical check_vdw_flag(atmax)

	common /COORD/atnum,resnum,sresnum,x,y,z,resid,atid,segid
	common /VDWB/atradtyp,epsilon,reps,atrad
	common /IDB/n,cresnum,catnum,IDCLR,cenx,ceny,cenz,atpr,
     +   atname,resname,cseg,cresid,catid,XYZ_QUEST,AorP,segname
	common /BORNB/irborn,rborn,rmax,pvol,avol,rvol,brphi,bphi,
     +   Bchg,dielec,atrvol,atvolchk,solvseg
	common /CONST/idiv,pi,div
	common /CHGS/Nchg,coi,xqc,yqc,zqc,qqc,rprobe,Cchg,
     +   dielecp,check_vdw_flag
	common /GROUP/nres,atomres,grp,numgrp,grpatom,atomgrp,
     +   gtype,grpnum,resgrp,acharge,q,gcharge,grptitle,tatid,
     +   tresid,akind,atype

c	 print *,k,' ',segname(k),'-',solvseg
c	 print *,k,' ',segname2(k),'-',solvseg
c	if (segname(k).eq.solvseg) then
	if (segname2(k).eq.solvseg) then
	 m=0
	 do 10 while ((m.le.nres).and.(tresid(m).ne.resname(k)))
c	  print *,m,nres,tresid(m),resname(k)
	  m=m+1
10	 continue
	 if (j.eq.1) then
	  xsum=0.0
	  ysum=0.0
	  zsum=0.0
	 endif

	  xsum=xsum + (x(k) - x(catnum))
	  ysum=ysum + (y(k) - y(catnum))
	  zsum=zsum + (z(k) - z(catnum))

	 if (j.eq.atomgrp(i)) then
	 xsum=xsum/atomgrp(i)
	 ysum=ysum/atomgrp(i)
	 zsum=zsum/atomgrp(i)
	 rsolv=sqrt(xsum**2 + ysum**2 + zsum**2)
	 if (rsolv.gt.rborn) then
	  do 30 l=k-(atomgrp(i)-1),k
	   check_vdw_flag(l)=.FALSE.
c	print *, l,rsolv,check_vdw_flag(l)
30	  enddo
	 else
	  do 40 l=k-(atomgrp(i)-1),k
	   check_vdw_flag(l)=.TRUE.
c	print *, l,rsolv,check_vdw_flag(l)
40	  enddo
	 endif
	 endif
	endif
	return
	end
