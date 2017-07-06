c---------------------------------------
c CMN.INC : Main program variables held
c           in COMMON blocks
c---------------------------------------

c Topology file variables
      integer TOPV,NNTA
      character*80 topfile

c Variables defining GROUP properties
      integer nres,atomres(a1dim),grpatom(a1dim,a2dim),grpnum(a1dim)
      integer grp(atmax),atomgrp(atmax),resgrp(atmax),numgrp
      integer gtype(atmax)
      real q(atmax),acharge(a1dim,a2dim),gcharge(a1dim,a2dim)
      character*30 grptitle(12)
      character*4 tatid(a1dim,a2dim),tresid(a1dim)
      character*4 akind(a1dim,a2dim),atype(atmax)

c Variables describing the coordinate (CRD) file
      integer atnum,resnum,sresnum
      real x(atmax),y(atmax),z(atmax)
      character*4 resid,atid,segid

c Variables used to hold values for the central residue
      integer cresat(atmax),iatpr

c Variables describing NONBONDed parameters
      integer numnbtype,nbindex(nbmax,2)
      real nbterms(nbmax,3)
      character*4 nbtype(nbmax)
      character*80 paramfile

c Variables describing VDW properties
      integer atradtyp(atmax)
      real epsilon(atmax),reps(atmax),atrad(atmax)

c More variables describing the central residue
      integer n,cresnum,catnum,IDCLR
      real cenx,ceny,cenz,atpr
      character*1 XYZ_QUEST
      character*4 atname(atmax),resname(atmax)
      character*4 cseg,cresid,catid
      character*4 AorP,segname(atmax)
        
c Variable used for calculating the BORN solvation contribution
      integer irborn
      real rborn,rmax
      real pvol(grpdim),avol(atmax),rvol(grpdim),brphi(grpdim),bphi
      real Bchg,dielec,atrvol(nbmax,grpdim)
      logical atvolchk(nbmax)
      character*4 solvseg

c Constant maintained throughout the code
      integer idiv
      real pi,div

c Variables describing the CHARGES of various groups
      integer Nchg,coi
      real xqc(maxchg),yqc(maxchg),zqc(maxchg),qqc(maxchg)
      real rprobe,Cchg,dielecp
      logical check_vdw_flag(atmax)

c Variables for timing functions
      real*4 timer,crap,agtime,vtime,ttime

c COMMON BLOCKS using the above variables
      common /TOP/TOPV,NNTA,topfile
      common /GROUP/nres,atomres,grp,numgrp,grpatom,atomgrp,
     +   gtype,grpnum,resgrp,acharge,q,gcharge,grptitle,tatid,
     +   tresid,akind,atype
      common /COORD/atnum,resnum,sresnum,x,y,z,resid,atid,segid
      common /AVG/cresat,iatpr
      common /NBONDB/numnbtype,nbindex,nbterms,nbtype,paramfile
      common /VDWB/atradtyp,epsilon,reps,atrad
      common /IDB/n,cresnum,catnum,IDCLR,cenx,ceny,cenz,atpr,
     +   atname,resname,cseg,cresid,catid,XYZ_QUEST,AorP,segname
      common /BORNB/irborn,rborn,rmax,pvol,avol,rvol,brphi,bphi,
     +   Bchg,dielec,atrvol,atvolchk,solvseg
      common /CONST/idiv,pi,div
      common /CHGS/Nchg,coi,xqc,yqc,zqc,qqc,rprobe,Cchg,
     +   dielecp,check_vdw_flag
      common /TIMING/agtime,vtime,ttime
