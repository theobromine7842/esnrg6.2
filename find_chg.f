* Detm. if an atom is a member of a charged group (5 or 6)
	SUBROUTINE FIND_CHG(i,j,k)
        
c 	integer atmax,nbmax,maxchg,grpdim,a1dim,a2dim
c         parameter (atmax=50000,nbmax=125,maxchg=1000,grpdim=1000)
c         parameter (a1dim=200,a2dim=75)

	INCLUDE 'dim.inc'

	real x(atmax),y(atmax),z(atmax)
	real xqc(maxchg),yqc(maxchg),zqc(maxchg)
	real qqc(maxchg),rprobe,dielecp

      real acharge(a1dim,a2dim),q(atmax),gcharge(a1dim,a2dim),Cchg

	integer i,j,k
	integer atnum,resnum,sresnum
	integer Nchg,coi
      integer m,n,nres,atomres(a1dim),grp(atmax),numgrp
      integer grpatom(a1dim,a2dim),atomgrp(atmax)
      integer gtype(atmax),grpnum(a1dim),resgrp(atmax)

	character*30 grptitle(12)
	character*4 resid,atid,segid
      character*4 akind(a1dim,a2dim),atype(atmax)
      character*4 tatid(a1dim,a2dim),tresid(a1dim)

	logical check_vdw_flag(atmax)

	common /COORD/atnum,resnum,sresnum,x,y,z,resid,atid,segid
	common /CHGS/Nchg,coi,xqc,yqc,zqc,qqc,rprobe,Cchg,
     +   dielecp,check_vdw_flag
	common /GROUP/nres,atomres,grp,numgrp,grpatom,atomgrp,
     +   gtype,grpnum,resgrp,acharge,q,gcharge,grptitle,tatid,
     +   tresid,akind,atype

	if ((gtype(i).eq.5) .or.
     +      (gtype(i).eq.6)) then
	 Nchg=Nchg + 1
	 qqc(Nchg) = q(k)
	 xqc(Nchg) = x(k)
	 yqc(Nchg) = y(k)
	 zqc(Nchg) = z(k)
	 write(6,'(4x,2(1x,I5),1x,A,1x,I5,4(1x,F9.5))') 
     +    i,gtype(i),atype(k),Nchg,qqc(Nchg),
     +    xqc(Nchg),yqc(Nchg),zqc(Nchg)
	else
	 write(6,'(4x,2(1x,I5),1x,A)') i,gtype(i),atype(k)
	endif
	return
	end
