* This subroutine calculates the solvent interaction energy
* of a charged system beyound a cutoff 'rborn' from the
* central position. In regions containing non-solvent atoms, this
* value is set to zero for any regions overlapping non-solvent
* atoms.
*
	SUBROUTINE CS_INTERFACE(segname2,Cchg2,dielecp2)
c	subroutine CS_INTERFACE

c         integer atmax,nbmax,maxchg,grpdim,a1dim,a2dim
c         parameter (atmax=50000,nbmax=125,maxchg=1000,grpdim=1000)
c         parameter (a1dim=200,a2dim=75)

	INCLUDE 'dim.inc'

	real x(atmax),y(atmax),z(atmax),xc,yc,zc,rc
	real rmin,rmax,rborn,div
	real pi,dvol,pvol(grpdim),rvol(grpdim),avol(atmax)
	real brphi(grpdim),bphi
	real dielec,dielecp,dielecp2
	real epsilon(atmax),reps(atmax),atrad(atmax)
	real cenx,ceny,cenz,atpr
	real atrvol(nbmax,grpdim)
	real xqc(maxchg),yqc(maxchg),zqc(maxchg)
	real qqc(maxchg),rprobe,temp
	real pbin,tbin,dbin,abin

	real Bchg,Bchg2,Cchg,Cchg2
      real acharge(a1dim,a2dim),q(atmax),gcharge(a1dim,a2dim)

	integer numatom,c
	integer i,j,k,f,t,irborn
	integer atnum,resnum,sresnum,cresnum,catnum
	integer atradtyp(atmax)
	integer Nchg,coi
	integer idiv
      integer m,n,nres,atomres(a1dim),grp(atmax),numgrp
      integer grpatom(a1dim,a2dim),atomgrp(atmax)
      integer gtype(atmax),grpnum(a1dim),resgrp(atmax)
      integer IDCLR

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


*
* Identify charged residue atoms and solvent outside the minimum
* radius (rborn). Then set the charge of interest, and pass the
* variables to the CHGSOLV subroutine.
*
	write(6,*) 'Gathering solvation information...'

c	print *, segname(1),' ',segname(433)
c	print *, segname2(1),' ',segname2(433)

	numatom = n
	c = catnum

	Cchg=Cchg2
	dielecp=dielecp2
	print *,'Cchg=',Cchg,'Cchg2=',Cchg2

	k=0
	do 150 i=1,numgrp
	 do 151 j=1,atomgrp(i)
	  k=k+1
	  call FIND_CHG(i,j,k)
c	print *,'Cchg=',Cchg,'Cchg2=',Cchg2
c	print *, '1=',qqc(1),xqc(1),yqc(1),zqc(1)
	  call FIND_SOLV(i,j,k,segname2)
c	  call FIND_SOLV(i,j,k)
c	print *,'Cchg=',Cchg,'Cchg2=',Cchg2
c	print *, '1=',qqc(1),xqc(1),yqc(1),zqc(1)
151	 enddo
150	enddo

	Cchg2=Cchg
c	print *,'Cchg=',Cchg,'Cchg2=',Cchg2

	call FIND_COI

	write(6,*)
	write(6,'(5x,A,A,$)') 'What temperature (K) do you ',
     +  'wish to use for distributions (REAL): '
	read(5,'(f9.5)') temp
	write(6,'(3x,f9.5)') temp

	write(6,*)
	write(6,'(5x,A,A,$)') 'How many phi bins per 360 degrees ',
     +  'do you wish to use (Real): '
	read(5,'(f9.5)') pbin
	write(6,'(3x,f9.5)') pbin

	write(6,*)
	write(6,'(5x,A,A,$)') 'How many theta bins per 180 degrees ',
     +  'do you wish to use (Real): '
	read(5,'(f9.5)') tbin
	write(6,'(3x,f9.5)') tbin

	write(6,*)
	write(6,'(5x,A,A,$)') 'How many delta bins per 360 degrees ',
     +  'do you wish to use (Real): '
	read(5,'(f9.5)') dbin
	write(6,'(3x,f9.5)') dbin

	write(6,*)
	write(6,'(5x,A,A,$)') 'How many alpa bins per 180 degrees ',
     +  'do you wish to use (Real): '
	read(5,'(f9.5)') abin
	write(6,'(3x,f9.5)') abin

	write(6,*)
	write(6,'(5(5x,A,/))')
     +   'This program excludes the contributions from solvent whose',
     +   'VdW radii would overlap those of other atoms. This is',
     +   'determined by calculating if the distance between a solvent',
     +   'location and an atom is less than the sum of a probe radius',
     +   'and the radius of the atom.'
	write(6,200) 'Please enter a solvent probe radius: '
	read(5,'(f5.2)') rprobe
	write(6,'(2x,F5.2)') rprobe
200	format (5x,A,$)

c	print *, pbin,tbin,dbin,abin,dielecp,Bchg

c* Calc. solvent interaction energy
	call CHGSOLV(n,Nchg,coi,catnum,pbin,tbin,dbin,abin,dielec,
     +   dielecp,rborn,rmax,div,temp,pi,qqc,xqc,yqc,zqc,x,y,z,
     +   rprobe,atrad,Bchg,check_vdw_flag,bphi,brphi)

	return
	end
