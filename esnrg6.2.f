********************************************************
        Program ESNRG
* Version 5.0
* This program is designed to calculate electrostatic  
* energy as a function of distance.                    
* Version 1 based on potential8lb4.f   		       
* Version 2.1 & 2.2 includes minor tuning              
* Version 2.3 Corrects avg. dist for resd files        
* Version 3.0 Adds Born correction via Yelle method    
* Version 3.1 prevents redundant atom looping in BORN  
*	      but this doesn't work yet, so it's been  
*	      forced to hit each atom.		      
* Version 3.2 adds a numerical integration solvation  
*	      energy term as well as cleaning up       
*             other functions.			       
* Version 3.3 fixes a problem with GLY or PRO as the   
*             N-terminal residue                       
* Version 4.0 Allows use of MacKerell all-atom params.
*             This version is also the first version
*             with split code requiring a Makefile  
* Version 4.1 Adds some timing features using the
*             external function dtime
* Version 4.2 Make timing features more generic to OS
*             by adding timer.f function
*******************************************************
* Author: Brian W. Beck & Robert B. Yelle 10/11/95     
*  Last Modified by Brian W. Beck    11/30/98          
*******************************************************
* Version 4.3 Extend the number of atoms and time scale
* Modified by Mingliang Tan	03/03/2004
* Version 5.0 Extend the number of dcd files (77) MLTAN 02/26/05
* Version 5.1 Add time series of ES between the center & specified 
*             residues  MLTAN 9/26/2007 (with MLT07) 
* Version 6.9 See CHANGELOG.md 2017-05-08 Kelly N. Tran

*
*     Set Up Variables 
*
      INCLUDE 'dim.inc'
      INCLUDE 'cmn.inc'

      real rmin,r,num4
      real lenx,leny,lenz,xi,yi,zi,ri
      real xsum,ysum,zsum,xc,yc,zc,gx,gy,gz,gr,rdist
      real sigma,lenb
      real phi,gphi(atmax),rphi(totgrp,grpdim),tphi(totgrp),sphi(totgrp)  
      real srphi(totgrp,grpdim),srphi2(totgrp,grpdim),
     +    sigrphi(totgrp,grpdim)
      real trphi(totgrp,grpdim),trphi2(totgrp,grpdim),tphi2(totgrp)
      real resphi(resdim,totgrp),sigresphi(resdim,totgrp)
      real tresphi(resdim,totgrp),tresphi2(resdim,totgrp)
      real timephi(totgrp,maxtime),arphi(totgrp,grpdim),ts
      real resdist(resdim,totgrp),trmax
      real invatpr,invatgrp,invnfrm,invdiv
      real tvol, Bchg2

      integer NCTL,ctat,dtat
      integer i10,i,j,k,l,si,t,numframe
      integer dist,natfree,natfix
      integer strtframe,endframe,skipframe,totframe,frame
      integer numdcd,dcdnum,du,s,sk,timestep(maxtime),j2
      integer reshit(resdim,totgrp),imdiv,sdiv,count1,count2

      character*80 crdfile,dcdfile,phifile
      character*80 dummy,title(20)
      character*4 psegid,patid,PBC_type
      character*4 filetype
cMLT07i
      integer nres2res,nresspec(20),nopen,hrr
      character*4 tres2res
      real timerrphi(5,totgrp,maxtime),ftimerrphi(5,totgrp)
cMLT07f
      character*4 calc_type,tempres,cut_type
      character*4 SOLV_type
      character*9 word1,word2,word3,word4
      character*1 crdast,BORN_quest,BCHG_quest
      character*1 SOLV_quest,CCHG_quest

      logical EOF_flag,pbc_flag,tx_flag,ty_flag,tz_flag
      logical Err_flag,BORN_flag,BCHG_flag,CHGS_flag,SOLV_flag
      logical CCHG_flag,CUT_flag,qcrys
c KNT
      character*4 hdr
      integer icntrl(20)

*
*   Clear arrays and variables
*
      pi = acos(-1.0)
      grptitle(1)  = 'Backbone Amide'
      grptitle(2)  = 'Backbone Carbonyl'
      grptitle(3)  = 'Polar Side Chains (PSC)'
      grptitle(4)  = 'Solvent (SOLV)'
      grptitle(5)  = 'Charged Side Chains (CSC)'
      grptitle(6)  = 'Counterions (CI)'
      grptitle(7)  = 'Backbone (BB)'
      grptitle(8)  = 'BB + PSC'
      grptitle(9)  = 'BB + PSC + SOLV'
      grptitle(10) = 'BB + PSC + CSC + CI'
      grptitle(11) = 'CSC + CI'
      grptitle(12) = 'TOTAL'
      nres = 0
      lenb = 0.
      NCTL = 0
      atnum = 0
      atpr = 0.0
      iatpr = 0
      IDCLR = 0
      Bchg = 0.
      Cchg = 0.
      cenx = 0.0
      ceny = 0.0
      cenz = 0.0
      rmin = 0.0
      xi = 0.0
      yi = 0.0
      zi = 0.0
      ri = 0.0
      tvol = 0.0
      num4 = 0.0000000000
      word4 = '         '
      cresid = '    '
      CHGS_flag =.FALSE.
      BORN_flag =.FALSE.
      rborn  = 0.
      irborn = 0
      dielec = 1.0
      dielecp = 1.0
      CCHG_quest = 'N'
      CCHG_flag = .FALSE.
      BCHG_quest = 'N'
      BCHG_flag = .FALSE.
      bphi = 0.0
      SOLV_quest='N'
      SOLV_flag=.FALSE.
      SOLV_type='NONE' 
      CUT_flag=.FALSE.
      do 1 j=1,12
      do 2 i=0,grpdim
      rphi(j,i) = 0.0
      trphi(j,i) = 0.0
      srphi(j,i) = 0.0
      trphi2(j,i) = 0.0
      arphi(j,i) = 0.0
      if (j.eq.4) brphi(i)=0.0
*
*	  write(6,*) j,i, arphi(j,i)
*
2     continue
      do 3 i=1,2500
*          j2 = int(j*0.5)
      resphi(i,j) = 0.0
      tresphi(i,j) = 0.0
      tresphi2(i,j) = 0.0
      sigresphi(i,j) = 0.0
3     continue
      do 4 i=1,atmax
      timephi(j,int(i*0.1)) = 0.0
      if (j.eq.1) then
      gphi(i) = 0.0
      x(i) = 0.0
      y(i) = 0.0
      z(i) = 0.0
      q(i) = 0.0000000000
      gtype(i) = 0
      grp(i) = 0
      atomgrp(i) = 0
      resgrp(i) = 0
      atname(i)  = '    '
      atype(i)   = '    '
      resname(i) = '    '
      epsilon(i)=0.0
      reps(i)=0.0
      atrad(i)=0.0
      atradtyp(i)=0
      check_vdw_flag(i)=.TRUE.
      endif 
4     continue
      sphi(j) = 0.0
      tphi(j) = 0.0
      tphi2(j) = 0.0
1     Continue
      do 5 i=1,a1dim
      do 6 j=1,a2dim
      akind(i,j)   = '    '
      acharge(i,j) = 0.0000000000
      gcharge(i,j) = 0.0000000000
      grpatom(i,j) = 0
      tatid(i,j) = '    '
*          write(6,*) i,j,acharge(i,j),gcharge(i,int(j*0.5))
6     continue
      tresid(i) = '    '
      atomres(i) = 0
      grpnum(i) = 0
5     continue  
      do 71 i=1,50
      atvolchk(i)=.FALSE.
      do 81 j=1,grpdim
      atrvol(i,j)=0.0
81    continue
71    continue
      strtframe = 1
      endframe = 1
      totframe = 1
      skipframe = 1
*
*     Ask user for names of input and output files
*
      write(6,'(5x,A)') ' '
      write(6,'(5x,A)') 'Welcome to ESNRG v. 6.1'
      write(6,'(5x,A)') ' '
      write(6,'(5x,A)') 'This program calculates the electrostatic'
      write(6,'(5x,A)') 'energy at a specified atom or residue. '
      write(6,'(5x,A,A)') '[or the Luzzatti error ',
     +    'sum((q1*q2)**2/dist**4) thereof.]'
      write(6,'(5x,A)') '(Experimental version)'
      write(6,'(5x,A)') ' '
      write(6,'(5x,A)') 'Do you wish to calculate electrostatic'
      write(6,'(5x,A,$)') 'energy (ESE) or the Luzzati error (ERR):'
      read(5,'(A)') calc_type
      write(6,'(A)') calc_type
      if ((calc_type(1:3).eq.'ERR').or.(calc_type(1:3).eq.'err')) then
         Err_flag=.TRUE.
      else
         Err_flag=.FALSE.
      endif
      write(6,'(5x,A)') ' '
7     write(6,'(5x,A)') 'Do you want to calculate the energy from'
      write(6,'(5x,A)') 'a coordinate file (type CRD) or '
      write(6,'(5x,A,$)') 'a trajectory file (type DCD):'
      read(5,'(A)') filetype
      if (filetype.eq.'crd') filetype='CRD'
      if (filetype.eq.'dcd') filetype='DCD'
      if ((filetype.eq.'CRD').or.(filetype.eq.'DCD')) then
         write(6,'(A)') filetype
      else
         write(6,'(A,$)') 'Please type in CRD or DCD:'
         goto 7
      endif
      write(6,'(5x,A)') ' '
      write(6,'(5x,A,A,$)') 'Enter the name of the CHARMM ',
     +    'Topology file:'
      read(5,'(A)') topfile
      write(6,'(5x,A)') topfile
      write(6,'(5x,A)') ' '
      write(6,'(5x,A,A,$)') 'Enter the name of the CHARMM ',
     +    'parameter file:'
      read(5,'(A)') paramfile
      write(6,'(5x,A)') paramfile
      write(6,'(5x,A)') ' '
      write(6,'(5x,A,$)') 'Enter the name of the CHARMM CRD file:'
      read(5,'(A)') crdfile
      write(6,'(5x,A)') crdfile
      write(6,'(5x,A)') ' '
      open(2,form='formatted',access='sequential',
     +    file=crdfile,status='old')
      write (6, '(5x, A, $)') 'What type of switching? (none ',
     +    'or afsw=atom force switch)'
      read (5, '(A)') cut_type
      write (6, '(5x, A)') cut_type
      write (6, '(5x, A)') ' '
      if (cut_type.eq.'afsw') then
         CUT_flag=.TRUE.
         write (6,'(5x,A,$)') 'Enter cuton, cutoff'
         read (5, '(2F9.4)') ron,roff
         write (6, '(5x, 2F9.5)') ron,roff 
         gamma = (roff**2 - ron**2)**3
         write (6, *) 'GAMMA',gamma 
         Fafsw = 8.*(((ron*roff)**2)*(roff-ron) -
     +    (roff**5 - ron**5)/5.)/gamma
         write (6, *) 'Fafsw', Fafsw 
         Aafsw = roff**4*(roff**2 - 3*ron**2)/gamma
         write (6, *) 'Aafsw', Aafsw 
         Bafsw = 6.*(ron*roff)**2/gamma
         write (6, *) 'Bafsw', Bafsw 
         Cafsw = - (ron**2 + roff**2)/gamma
         write (6, *) 'Cafsw',Cafsw
         Dafsw = 2./(5.*gamma)
         write (6, *) 'Dafsw', Dafsw 
*        else 
*	  CUT_flag=.FALSE.
*	  write (6, '(5x,A)') 'No cutoffs on energy or force selected' 
      endif

*
* open up data files:
*
      open(13,form='formatted',access='sequential',
     +    file='dat/esgrp.dat',status='unknown')
      open(14,form='formatted',access='sequential',
     +    file='dat/esrgrp.dat',status='unknown')
      open(15,form='formatted',access='sequential',
     +    file='dat/esgrptot.dat',status='unknown')
      open(16,form='formatted',access='sequential',
     +    file='dat/esres.dat',status='unknown')
      open(17,form='formatted',access='sequential',
     +    file='dat/esresd.dat',status='unknown')

      if (filetype.eq.'DCD') then
8     write(6,'(5x,A,$)') 'Enter the number of CHARMM DCD files:'
      read(5,'(I2)') numdcd
      write(6,'(5x,I2)') numdcd
      write(6,*)
      if (numdcd.gt.77) then
      write(6,'(A)') 'No more than 77 dcd files can be read!'
      goto 8
      endif
      do 9 i=1,numdcd
      write(6,'(A,I2,A,$)') 'Enter the name of CHARMM DCD file #',
     +    i,': '
      read(5,'(A)') dcdfile
      write(6,'(5x,A)') dcdfile
      open(i+22,form='unformatted',access='sequential',
     +    file=dcdfile,status='old')
9     continue

      open(18,form='formatted',access='sequential',
     +    file='dat/esgf.dat',status='unknown')
      open(19,form='formatted',access='sequential',
     +    file='dat/esgftot.dat',status='unknown')
      open(12,form='formatted',access='sequential',
     +    file='dat/esrf.dat',status='unknown')
      open(20,form='formatted',access='sequential',
     +    file='dat/estime.dat',status='unknown')
      open(21,form='formatted',access='sequential',
     +    file='dat/estottime.dat',status='unknown')
      endif
*
* Ask user for the box length, radial distance, and the number of 
* divisions per unit length. 
*
      write(6,'(5x,A)') 'What type of periodic boundary conditions'
      write(6,'(5x,A,$)') 'do you want to use (cube, octa, or none):'
      read(5,'(A)') PBC_type
      if (PBC_type.eq.'CUBE') PBC_type='cube'
      if (PBC_type.eq.'OCTA') PBC_type='octa'
      if (PBC_type.eq.'NONE') PBC_type='none'
      write(6,'(5x,A)') PBC_type
      write(6,'(1x,A)') ' '
      if (PBC_type.ne.'none') then
      write(6,'(5x,A,$)') 'Please enter the box dimensions (x,y,z): '
      read(5,'(3(F9.4))') lenx,leny,lenz
*         write(6,'(5x,A)') 'x: '
*         read(5,'(F9.4)') lenx
*         write(6,'(5x,A)') 'y: '
*         read(5,'(F9.4)') leny
*         write(6,'(5x,A)') 'z: '
*         read(5,'(F9.4)') lenz
      write(6,'(3(1x,f9.4))') lenx,leny,lenz
      write(6,'(1x,A)') ' '
      if (PBC_type.eq.'octa') then
          trmax=0.75/sqrt((1./lenx)**2+(1./leny)**2+(1./lenz)**2)
      else
          trmax=0.5*amin1(lenx,leny,lenz)
      endif
      endif
      write(6,'(5x,A)') 'Enter the minimum distance from center'
      write(6,'(5x,A,$)') 'atom to calculate the energy around:'
      read(5,'(F9.4)') rmin
      write(6,*) rmin
      write(6,'(1x,A)') ' '
      write(6,'(5x,A,F9.4,A)') '(Closest Edge of box = ',
     +    trmax,' )'
      write(6,'(1x,A)') ' '
      write(6,'(5x,A)') 'Enter the maximum distance from center'
      write(6,'(5x,A,$)') 'atom to calculate the energy around:'
      read(5,'(F9.4)') rmax
      write(6,*) rmax
      write(6,'(1x,A)') ' '
      write(6,'(5x,2(A),$)') 'Please enter the number of divisions',
     +    ' per unit length (real): '
      read(5,'(F4.1)') div 
      write(6,*) div
      imdiv = int(rmin*div)
      idiv  = int(rmax*div)
      invdiv= 1./div
*
* Determine whether to use solvent correction.
* Don't do corrections if calculating Luzzati error
*
      if (.not. Err_flag) then
          write(6,'(1x,A)') ' '
          write(6,'(5x,A,A,$)') 'Would you like to include a ',
     +    'solvent correction ? (y/n): '
          read(5,'(A1)') SOLV_quest
          write(6,*) SOLV_quest
      else
          SOLV_quest='N'
          SOLV_flag=.FALSE.
      endif
      if ((SOLV_quest.eq.'y') .or. (SOLV_quest.eq.'Y')) then
          SOLV_quest='Y'
          SOLV_flag=.TRUE.
          write(6,'(1x,A)') ' '
          write(6,'(3(5x,A/),5x,A,$)') 'Would you like to use:',
     +    'The NUMERICAL INTEGRATION method (NUMI)',
     +    'The BORN method (BORN)',
     +    'or BOTH (NUMB) : '
          read(5,'(A4)') SOLV_type
          write(6,'(3x,A4)') SOLV_type
      else
          SOLV_quest='N'
          SOLV_flag=.FALSE.
      endif

      if (SOLV_type.eq.'NUMI') then
          CHGS_flag =.TRUE.
      write(6,*)
      write(6,'(2(5x,A,/))')'The NUMI method numerically integrates ',
     +    'the theoretical solvent contribution over r,theta, and phi.'
      write(6,*)
      write(6,'(5x,A,A)') 'The NUMI correction requires a ',
     +    'minimum distance.'
      write(6,'(5x,A,A,$)') 'Enter the start distance for the ',
     +    'NUMI correction: '
      read(5,'(F9.4)') rborn
               write(6,*) rborn
      irborn = int(rborn*div)
      write(6,*)
      write(6,'(5x,A,$)') 'Enter the bulk dielectric constant: '
      read(5,'(F9.4)') dielec
      write(6,*) dielec
      write(6,*)
      write(6,'(5x,A,$)') 'Enter the local dielectric constant: '
      read(5,'(F9.4)') dielecp
      write(6,*) dielecp
      write(6,*)
      write(6,'(5x,A,A,$)') 'Do you wish to use the ',
     +    'net charge of the center group ? (y/n): '
      read(5,'(A1)') CCHG_quest
      write(6,*) CCHG_quest
      if ((CCHG_quest.eq.'y').or.(CCHG_quest.eq.'Y')) then
          CCHG_quest = 'Y'
          CCHG_flag = .TRUE.
      else
          write(6,'(5x,A,$)') 'Please enter the charge to use: '
          read(5,'(F9.4)') Cchg
          write(6,*) Cchg
          CCHG_quest = 'N'
          CCHG_flag = .FALSE.
      endif
      write(6,'(5x,A,A,$)') 'Please enter the solvent segment ',
     +    'name: '
      read(5,'(A4)') solvseg
      write(6,*) solvseg
      write(6,*)
*      else
*      CHGS_flag =.FALSE.
*      rborn  = rmax
*      irborn = idiv
*      dielecp = 1.0
*      CCHG_quest = 'N'
*      CCHG_flag = .FALSE.
*      Cchg = 1.0
*      bphi = 0.0
*      endif

* I included the BCHG_flag=.TRUE. in NUMI so
* could keep the BORN approx. in the CHGSOLV subroutine
      BCHG_flag=.TRUE.

      else if (SOLV_type.eq.'BORN') then
c      BORN_quest = 'Y'
      BORN_flag =.TRUE.
      write(6,*)
      write(6,'(5x,A,A)') 'The BORN expression used is ',
     +    '(1-1/E)kq^2/r .'
      write(6,*)
      write(6,'(5x,A)') 'The Born correction needs a spherical cutoff.'
      write(6,'(5x,A,A,$)') 'Enter the cut-off for the ',
     +    'Born correction: '
      read(5,'(F9.4)') rborn
      write(6,*) rborn
      irborn = int(rborn*div)
      write(6,*)
      write(6,'(5x,A,$)') 'Enter the bulk dielectric constant: '
      read(5,'(F9.4)') dielec
      write(6,*) dielec
      write(6,*)
      write(6,'(5x,A,A,$)') 'Do you wish to use the system ',
     +    'net charge for q ? (y/n): '
      read(5,'(A1)') BCHG_quest
      write(6,*) BCHG_quest
      if ((BCHG_quest.eq.'y').or.
     +    (BCHG_quest.eq.'Y')) then
      BCHG_quest = 'Y'
      BCHG_flag = .TRUE.
      else
      write(6,'(5x,A,$)') 'Please enter the charge to use: '
      read(5,'(F9.4)') Bchg
      write(6,*) Bchg
      BCHG_quest = 'N'
      BCHG_flag = .FALSE.
      endif
      write(6,'(5x,A,A,$)') 'Please enter the solvent segment ',
     +    'name: '
      read(5,'(A4)') solvseg
      write(6,*) solvseg
      write(6,*)
*      else
c      BORN_quest = 'N'
*      BORN_flag =.FALSE.
*      rborn  = rmax
*      irborn = idiv
*      dielec = 1.0
*      BCHG_quest = 'N'
*      BCHG_flag = .FALSE.
*      Bchg = 1.0
*      bphi = 0.0
*      endif

      else if (SOLV_type.eq.'NUMB') then
      CHGS_flag =.TRUE.
      BORN_flag =.TRUE.
      write(6,*)
      write(6,'(5x,A,A,/,A)') 'The NUMB method both numerically ',
     +    'integrates the theoretical solvent contribution over ',
     +    'r,theta, and phi,'
      write(6,'(5x,A,A)') 'and adds a BORN correction past the ',
     +    'edge of the system.'
      write(6,*)
      write(6,'(5x,A,A)') 'The NUMB correction requires a ',
     +    'minimum distance.'
         write(6,'(5x,A,A,$)') 'Enter the start distance for the ',
     +    'NUMB correction: '
      read(5,'(F9.4)') rborn
      write(6,*) rborn
      irborn = int(rborn*div)
      write(6,*)
      write(6,'(5x,A,$)') 'Enter the bulk dielectric constant: '
      read(5,'(F9.4)') dielec
      write(6,*) dielec
      write(6,*)
      write(6,'(5x,A,$)') 'Enter the local dielectric constant: '
      read(5,'(F9.4)') dielecp
      write(6,*) dielecp
      write(6,*)
      write(6,'(5x,A,A,$)') 'Do you wish to use the ',
     +    'net charge of the center group ? (y/n): '
      read(5,'(A1)') CCHG_quest
      write(6,*) CCHG_quest
      if ((CCHG_quest.eq.'y').or.
     +    (CCHG_quest.eq.'Y')) then
          CCHG_quest = 'Y'
          CCHG_flag = .TRUE.
	 else
          write(6,'(5x,A,$)') 'Please enter the charge to use: '
          read(5,'(F9.4)') Cchg
          write(6,*) Cchg
          CCHG_quest = 'N'
          CCHG_flag = .FALSE.
      endif
      write(6,*)
      write(6,'(5x,A,A,$)') 'Do you wish to use the system ',
     +    'net charge for the BORN charge ? (y/n): '
      read(5,'(A1)') BCHG_quest
      write(6,*) BCHG_quest
      if ((BCHG_quest.eq.'y').or.
     +    (BCHG_quest.eq.'Y')) then
          BCHG_quest = 'Y'
          BCHG_flag = .TRUE.
      else
          write(6,'(5x,A,$)') 'Please enter the charge to use: '
          read(5,'(F9.4)') Bchg
          write(6,*) Bchg
          BCHG_quest = 'N'
          BCHG_flag = .FALSE.
      endif
      write(6,'(5x,A,A,$)') 'Please enter the solvent segment ',
     +    'name: '
      read(5,'(A4)') solvseg
      write(6,*) solvseg
      write(6,*)

      else
          CHGS_flag =.FALSE.
          BORN_flag =.FALSE.
          rborn  = rmax
          irborn = idiv
          dielec = 1.0
          dielecp = 1.0
          CCHG_quest = 'N'
          CCHG_flag = .FALSE.
          Cchg = 0.0
          BCHG_quest = 'N'
          BCHG_flag = .FALSE.
          Bchg = 0.0
          bphi = 0.0
      endif

*
* Ask user for atom or coordinate to center around and
* if necessary, find atom in coordinate file.
*
919   write(6,'(5x,A,A)') 'What location would you like to center ',
     +    'the calculation around :'
      write(6,'(5x,A)') '(Please enter ATOM to center around an atom)'
      write(6,'(5x,A)') '(or enter RESI to use a whole residue.)'
      read(5,'(A4)') AorP
      write(6,'(5x,A)') AorP
      write(6,'(5x,A)') 
      if (AorP.eq.'RESI') then
          XYZ_QUEST='y'            
      else
          XYZ_QUEST='n'
          cseg='NONE'
          cresnum=0
          catid='NONE'
      endif
      if (AorP.eq.'ATOM') then
          write(6,'(5x,A,A)') 'Which atom would you like to use ',
     +    'as the center for this calculation : '
          write(6,'(5x,A)') '(Please use charmm CRD format.)'
          write(6,'(5x,A)') '(Note that the contributions from any )'
          write(6,'(5x,A)') '(member of this residue will be ignored.)'
          Write(6,'(1x,A)') 'SEGID:'
          read(5,'(A4)') cseg
          Write(6,'(1x,A)') 'RESNUM:'
          read(5,'(I4)') cresnum
          Write(6,'(1x,A)') 'ATOM:'
          read(5,'(A4)') catid
          Write(6,'(1x,A4,1x,I4,1x,A4)') cseg,cresnum,catid
          write(6,'(5x,A)') ' '
      else if (AorP.eq.'RESI') then
          XYZ_QUEST='y'
          catid='ALL'
          write(6,'(5x,A,A)') 'Which residue would you like to use ',
     +   'as the center for this calculation : '
          write(6,'(5x,A)') '(Please use charmm CRD format.)'
          write(6,'(5x,A)') '(Note that the contributions from any )'
          write(6,'(5x,A)') '(member of this residue will be ignored.)'
          Write(6,'(1x,A)') 'SEGID:'
          read(5,'(A4)') cseg
          Write(6,'(1x,A)') 'RESNUM:'
          read(5,'(I4)') cresnum
          Write(6,'(1x,A4,1x,I4)') cseg,cresnum
          write(6,'(5x,A)') ' '
      else
          write(6,'(5x,A)') 'Please enter ATOM or RESI.'
          write(6,*)
          goto 919
      endif
      EOF_flag = .FALSE.
*
* Read topology file.  Initialize residues, groups, and charges.
*
c      call TOPREAD
      call TOPREAD2
*
* Read Parameter file and retrieve all non-bonded parameters
*
      call NBOND
*
* Make group and charge assignments to the coordinate file.
*
      write(6,*)
      write(6,*) 'Reading coordinate file...'
      Do 100 While (.not. EOF_flag)
          read(2,'(A1)',end=110) crdast
          if (crdast .eq. '*' ) then
              NCTL = NCTL + 1
          else
              EOF_flag =.TRUE.
          endif
100   Continue
110   EOF_flag =.TRUE.
      Continue

      Rewind 2
      Do 120 i=1,NCTL
          read(2,'(A80)',end=180) dummy
          write(6,'(1x,A80)') dummy
120   Continue
      read(2,'(I5)',end=180) ctat
      write(6,'(I5)') ctat

      EOF_flag=.FALSE. 
      n = 0

      Do 130 while (.not. EOF_flag) 
          n = n + 1
          read(2,135,end=180) atnum,resnum,resid,atid,
     +    x(n),y(n),z(n),segid,sresnum
135   format(2I5,2(1x,A4),3(f10.5),1x,A4,1x,I4)
c---------
c           write(6,135) n,resnum,resid,atid,
c      +    x(n),y(n),z(n),segid,sresnum
c---------
*          resname(n) = resid
*          atname(n) = atid
      tempres = resid
c           call id

c	write(6,'(5x,A)') ' '
c	write(6,913) segid,cseg
c	write(6,914) sresnum,cresnum
c	write(6,915) atid,catid
c	write(6,916) XYZ_QUEST
c913	format(5x,'segid=',A4,5x,'cseg=',A4)
c914	format(5x,'sresnum=',I4,5x,'cresnum=',I4)
c915	format(5x,'atid=',A4,5x,'catid=',A4)
c916	format(5x,'XYZ_QUEST=',A1)

*          if (((segid.eq.cseg).and.(sresnum.eq.cresnum)).and.
*     +    ((atid.eq.catid).or.(XYZ_QUEST.eq.'y'))) then
*           catnum = n
*           cresid = resid
*           write(6,*) segid,resnum,atid,catnum,cresid,sresnum
*          endif

      if ((atid.eq.'HT1').or.(atid.eq.'HN1')) then
      psegid=segid
      patid=atid
      if (TOPV.eq.22) then
c-------
c      print *, '*** TOPV=22 ***'
c-------
      if (resid.eq.'PRO') then
         NNTA=8
         resid='PROP'
      elseif (resid.eq.'GLY') then
         NNTA=7
         resid='GLYP'
      else
         NNTA=6
         resid='NTER'
      endif
      n = n - 1
      atnum = atnum - 1
      numgrp = numgrp - 1
      atid='N'
      call assigngrp
      call vdw
      resid = tempres
      IDCLR=1
      call id
      IDCLR=0
*      resname(n) = resid
*      atname(n) = atid
c---------
c      write(6,135) atnum,resnum,resid,atid,
c     +    x(n),y(n),z(n),segid,sresnum
c---------
      n = n + 1
      atnum = atnum + 1
c      numgrp = numgrp - 1
      atid=patid 

      do 140 i=1,NNTA-1
          if (resid.eq.'PRO') then
              NNTA=8
              resid='PROP'
          elseif (resid.eq.'GLY') then
              NNTA=7
              resid='GLYP'
          else
              NNTA=6
              resid='NTER'
          endif
      call assigngrp
      call vdw
      resid = tempres
      call id
*      resname(n) = resid
*      atname(n) = atid
c---------
c      write(6,135) atnum,resnum,resid,atid,
c      +    x(n),y(n),z(n),segid,sresnum
c      print *, 'resid=',resid, 'call assigngrp'
c---------
      n = n + 1
      read(2,135,end=180) atnum,resnum,resid,atid,
     +    x(n),y(n),z(n),segid,sresnum
140   continue

      else
      do 141 i=1,NNTA
          if (resid.eq.'PRO') then
              resid='PROP'
          elseif (resid.eq.'GLY') then
              resid='GLYP'
          else
              resid='NTER'
          endif
          call assigngrp
          call vdw
          resid = tempres
          call id
*          resname(n) = resid
*          atname(n) = atid
c---------
c          write(6,135) atnum,resnum,resid,atid,
c      +    x(n),y(n),z(n),segid,sresnum
c          print *, 'resid=',resid, 'call assigngrp'
c---------
          n = n + 1
          read(2,135,end=180) atnum,resnum,resid,atid,
     +    x(n),y(n),z(n),segid,sresnum
141   continue
      endif
      elseif (atid.eq.'CF') then
          psegid=segid
          resid='FORM'
*          write(6,135) atnum,resnum,resid,atid,
*     +    x(n),y(n),z(n),segid,sresnum
          call assigngrp
          call vdw
          resid = tempres
          call id
          resid = 'FORM'
          n = n + 1
          read(2,135,end=180) atnum,resnum,resid,atid,
     +    x(n),y(n),z(n),segid,sresnum
          resid='FORM'
*          write(6,135) atnum,resnum,resid,atid,
*     +    x(n),y(n),z(n),segid,sresnum
          call assigngrp
          call vdw
          resid = tempres
          call id
          resid = 'FORM'
          n = n + 1
          read(2,135,end=180) atnum,resnum,resid,atid,
     +    x(n),y(n),z(n),segid,sresnum
      elseif (atid.eq.'OT1') then
          n = n - 1
          atnum = atnum - 1
          numgrp = numgrp - 1
          resid='CTER'
          atid='C'
          call assigngrp
          call vdw
          resid = tempres
          IDCLR=1
          call id
          IDCLR=0
*          resname(n) = resid
*          atname(n) = atid
c---------
c          write(6,135) atnum,resnum,resid,atid,
c      +    x(n),y(n),z(n),segid,sresnum
c---------
          resid='CTER'
          n = n + 1
          atnum = atnum + 1
C          numgrp = numgrp - 1
          atid='OT1' 
          call assigngrp
          call vdw
          resid = tempres
          call id
*          resname(n) = resid
*          atname(n) = atid
c---------
c          write(6,135) atnum,resnum,resid,atid,
c      +    x(n),y(n),z(n),segid,sresnum
c---------
          resid='CTER'
          n = n + 1
          read(2,135,end=180) atnum,resnum,dummy,atid,
     +    x(n),y(n),z(n),segid,sresnum
          atid='OT2' 
      endif
      call assigngrp
      call vdw
      resid = tempres
      call id
*      resname(n) = resid
*      atname(n) = atid
c---------
c      write(6,135) atnum,resnum,resid,atid,
c      +    x(n),y(n),z(n),segid,sresnum
c---------
      resid='CTER'

      if (atnum.eq.ctat) EOF_flag=.TRUE.
130   Continue

      if (BCHG_flag) then
      Do 118 j=1,ctat
          Bchg=Bchg + q(j)
c          write(6,117) 'atype(j)=',atype(j),'q(j)=',q(j),
c     +    'Bchg=',Bchg,'atrad(j)=',atrad(j)
117   format (2x,A,A,3(2x,A,F10.5))
118   Continue
c      print *, 'Bchg=',Bchg
      endif

c      print *, 'CHGS_flag=',CHGS_flag
      if (CHGS_flag) then
          call PROPS(Cchg)
      endif

      invatpr=1./atpr

      if (catnum.eq.0) then
          write(6,*) cseg,cresnum,catid,' has not been found.'
          write(6,*) 'Please try again!'
          goto 9999
      endif
      if (AorP.eq.'RESI') then
          catnum=ctat+1
          x(catnum)=cenx*invatpr
          y(catnum)=ceny*invatpr
          z(catnum)=cenz*invatpr
          write(6,*) 'atoms per central residue = ',atpr
          write(6,'(5x,A,A,I5)') 'Geom. Center of residue ',
     +    'stored in index #',catnum
          write(6,*) 'coord=',x(catnum),y(catnum),z(catnum)
      endif

      if (BCHG_flag)
     +    write(6,'(5x,A,F10.5)') 'System net charge = ',Bchg

*
* Now calculate solvation energies, as required by above switches
*
c      if (BORN_flag) then
      if (SOLV_type.eq.'BORN') then
      call BORN
c      call BORN2
      endif

      if (CHGS_flag) then
      q(catnum)=Cchg
      write(6,'(5x,A,F10.5)') 'Charge of central atom or residue = ',
     +    Cchg
c          print *, segname(1),' ',segname(433)
c          print *,'Cchg=',Cchg
c          call CS_INTERFACE
      call CS_INTERFACE(segname,Cchg,dielecp)
c          call CHK_SEG(segname,n)
      endif

      if (SOLV_type.eq.'NUMB') then
          bphi = bphi - brphi(idiv)
          Bchg2=Bchg*Bchg
          brphi(idiv) = brphi(idiv) + (-332.*Bchg2*(1.- (1./dielec))*
     +    (1/rmax))
          bphi = bphi + brphi(i)
      endif

      rewind 2
      if (filetype.eq.'CRD') then 
          write(6,*) 'Calculating elec. energy from coordinate file.'
          goto 206 
      else
          goto 190
      endif

180   write(6,*) 'EOF encountered in CRD file.'
      write(6,*) 'ESNRG program aborting.' 
      EOF_flag=.TRUE.
      Continue
      GOTO 9999
190   Write(6,'(5x,A)') 'Assignment of groups completed.'
      Write(6,'(A)') ' '
      Rewind 2
      EOF_flag =.FALSE.
*
*  Read dcd file type, # frames, total # atoms, etc. 
*
      write(6,'(A)') 'Which frame do you want to start at?'
      read(5,'(I6)') strtframe
      write(6,'(I6)') strtframe
191   write(6,'(A)') 'Which frame do you want to end at?'
      read(5,'(I6)') endframe
      write(6,'(I6)') endframe
      if (endframe.lt.strtframe) then
      write(6,'(A,A)') 'The last frame must be greater ', 
     +    'than the starting frame!'
      goto 191
      endif
      write(6,'(A)') 'How many frames do you want to skip over?'
      read(5,'(I6)') skipframe
      write(6,'(I6)') skipframe
      write(6,'(A,I6,A,I6)') 'Going from ',strtframe,' to ',endframe
      write(6,'(A,I6,A)') 'skipping every ',skipframe,' frames.'

      if (AorP.eq.'RESI') then
          cenx = 0.0
          ceny = 0.0
          cenz = 0.0
          endif

cMLT07i
      write(6,'(5x,A)') 'Do you wish to calculate ES between center '
      write(6,'(5x,A)') 'residue and a specified residue? (yes/no)'
      read(5,'(A)') tres2res
      if (tres2res.eq.'yes') then
      write(6,'(5x,A)') 'How many specified residue(s)?'
      read(5,'(I2)') nres2res
      do i=1,nres2res
      write(6,'(5x,A,x,I2)') 'Residue # for specified residue',i
      read(5,'(I4)') nresspec(i)
      nopen=100-i
      open(nopen,form='formatted',access='sequential',
     +    status='unknown')
      enddo
      endif
cMLT07f

      dcdnum = 0
      numframe = 0.0

195   dcdnum = dcdnum+1     

      du = dcdnum+22
      write(6,*) du,dcdnum
      read(du) hdr,icntrl
      write(6,*) 'Header: ',hdr
      write(6,*) 'ICNTRL: ',icntrl
      if (hdr.ne.'CORD') then
          write(6,'(3x,A)') 'The DCD file you specified is not',
     +    'a coordinate file!'
          write(6,'(3x,A)') 'ESNRG program aborting!'
          goto 9999
      endif

      qcrys=(icntrl(11).eq.1) 
      totframe = icntrl(1)
      write(6,'(A,I8)') 'Number of frames = ',totframe
      write(6,'(A,I9)') 'Initial step = ',icntrl(2)
      write(6,'(A,I6)') 'Write frequency = ',icntrl(3)
      if (strtframe.gt.totframe) then
          strtframe=strtframe-totframe
          dcdnum=dcdnum+1
          goto 195
      endif
      if ((totframe.lt.endframe).and.(dcdnum.eq.numdcd)) then
          write(6,*) 'There are only ',totframe,' remaining frames!'
          endframe=totframe
      endif
      if (totframe.gt.endframe) totframe=endframe
          endframe=endframe-totframe

      read(du) i10,(title(j),j=1,i10)
      Do 200 j=1,i10
          Write(6,'(A)') title(j)
200   Continue
      read(du) dtat
      if (icntrl(9).ne.0) read(du) natfix
      natfix = icntrl(9)      
      natfree = dtat-natfix
      write(6,'(A,I5)') 'The total number of atoms is ',dtat
      write(6,'(A,I5)') 'The number of fixed atoms is ',natfix

      if (dtat.ne.ctat) then
      write(6,*) 'dtat = ',dtat,'	ctat = ',ctat
      write(6,'(3x,A)') 'The number of atoms from the CRD file'
      write(6,'(3x,A)') 'does not match the number of atoms from '
      write(6,'(3x,A)') 'the dynamics file! '
      write(6,'(3x,A)') 'EPOT program aborting!'
      goto 9999
      endif
*
* Read coordinates from dcd file, calculate electrostatic energy.
*
      write(6,'(A,I5,A)') 'Skipping the first ',strtframe-1,' frames.'

      if (strtframe.ne.1) then
      Do 205 j=1,strtframe-1
          n = 0
          if (qcrys) read(du) xtl
          if (j.eq.1) then
              read(du) (x(i),i=1,dtat)
              read(du) (y(i),i=1,dtat)
              read(du) (z(i),i=1,dtat)
          else
              read(du) (x(i),i=natfix+1,dtat)
              read(du) (y(i),i=natfix+1,dtat)
              read(du) (z(i),i=natfix+1,dtat)
          endif           
205   Continue
      endif

206   frame = 0
      numframe = 0
      write(6,*) 'strtframe  ','  totframe  ','  skipframe'
      write(6,*) strtframe,totframe,skipframe

      Do 210 j=strtframe,totframe,skipframe

      n = 0
      numframe = numframe + 1.
      frame = frame + 1
      timestep(int(numframe)) = icntrl(2)+(j-1)*icntrl(3)
c------
c      write(6,*) 'numframe = ',numframe, '  frame = ',frame
c------

      if (filetype.eq.'DCD') then 
      if (qcrys) read(du) xtl 
      if (j.eq.1) then
          read(du) (x(i),i=1,dtat)
          read(du) (y(i),i=1,dtat)
          read(du) (z(i),i=1,dtat)         
      else
          if (j.eq.strtframe) then
          si=skipframe
      else
          si=1
      endif
      do 215 s=si,skipframe
          sk=j-(skipframe-s)
          if (s.ne.skipframe) 
     +    write(6,'(A,I5)') '  Skipping frame #',sk
          read(du) (x(i),i=natfix+1,dtat)
          read(du) (y(i),i=natfix+1,dtat)
          read(du) (z(i),i=natfix+1,dtat)
215   continue
      endif

      if (j.eq.1) then
          write(6,'(5x,A)') 'Analyzing'
      endif
      if (mod(j,totframe).eq.0) then
          write(6,'(A,2x,I4)') '.',j
      else
          if (mod(j,50).ne.1) then
          if (mod(j,50).ne.0) then
              write(6,'(A,$)') '.'
          else
              write(6,'(A,2x,I4)') '.',j
          endif
          else
              write(6,'(5x,A,$)') '.'
          endif
      endif

c      if (AorP.eq.'RESI') then
c      write(6,'(A,I5,$)') 'Frame # ',j
c      else
c      write(6,'(A,I5)') 'Frame # ',j
c      endif

      endif

      if (AorP .eq. 'RESI') then
      if (filetype.eq.'DCD') then
      Do 216 i=1,iatpr
          cenx = cenx + x(cresat(i))
          ceny = ceny + y(cresat(i))
          cenz = cenz + z(cresat(i))
216   continue
      x(catnum)=cenx*invatpr
      y(catnum)=ceny*invatpr
      z(catnum)=cenz*invatpr
      cenx=0.
      ceny=0.
      cenz=0.
c      write(6,'(5x,A,3(F9.5,2x))') 'center xyz= ',
c     +    x(catnum),y(catnum),z(catnum)
      endif
c------
c 	print *, 'iatpr=',iatpr
c------
      Do 912 i=1,iatpr
          n = 0
      Do 2220 k=1,numgrp

      xsum = 0.
      ysum = 0.
      zsum = 0.
      gphi(k) = 0.
c------
c      print *, 'atomgrp(',k,')=',atomgrp(k)
      invatgrp= 1./float(atomgrp(k))
c      print *, 'invatgrp=',invatgrp
c------
      tx_flag = .FALSE.
      ty_flag = .FALSE.
      tz_flag = .FALSE.

      Do 2230 l=1,atomgrp(k)

      n = n + 1

      if (resname(n).eq.cresid) then
c      write(6,*) k,l,n,resname(n),cresid
          phi = 0.
          goto 2230
      endif
c      write(6,*) k,l,n,resname(n),cresid
 
      xc = x(n) - x(catnum)
      yc = y(n) - y(catnum)
      zc = z(n) - z(catnum)

      xi = x(n) - x(cresat(i))
      yi = y(n) - y(cresat(i))
      zi = z(n) - z(cresat(i))
      ri = sqrt(xi**2 + yi**2 + zi**2)

      x0 = xc
      y0 = yc
      z0 = zc
      r0 = sqrt(x0**2 + y0**2 + z0**2)

c
      write(6,2225) n,l,x(n),y(n),z(n)
2225  format(I5,1x,I5,1x,3(f7.2,1x))
c
      resid = resname(n)
c--------
c      print *, 'resname(',n,')=',resname(n)
c--------

      if (PBC_type.ne.'none') then
      if ((resname(n).eq.'TP4E').and.((l.eq.1).or.
     +    (tx_flag).or.(ty_flag).or.(tz_flag))) then 
      call translate(l,xc,yc,zc,lenx,leny,lenz,tx_flag,
     +    ty_flag,tz_flag,PBC_type)
      endif
      if ((resname(n).eq.'POT').or.(resname(n).eq.'CLA')) then
      call translate(l,xc,yc,zc,lenx,leny,lenz,tx_flag,
     +    ty_flag,tz_flag,PBC_type)
      endif
      if ((resid(1:3).eq.'DMF').and.((l.eq.1).or.
     +    (tx_flag).or.(ty_flag).or.(tz_flag))) then
      call translate(l,xc,yc,zc,lenx,leny,lenz,tx_flag,
     +    ty_flag,tz_flag,PBC_type)
      endif
      endif

      r = sqrt(xc**2 + yc**2 + zc**2)

c      if ((j.eq.1).and.((tx_flag).or.(ty_flag).or.
c     +    (tz_flag))) then
c      write(6,2225) n,l,x0,y0,z0,r0
c      write(6,2225) n,l,xc,yc,zc,r
c2225  format(I5,1x,I5,1x,4(f7.2,1x))
c      endif

      if (Err_flag) then
          phi = ((332.*q(n)*q(cresat(i)))**2)/(ri**4)
      elseif (CUT_flag) then
      if (ri.le.ron) then
          phi = 332.*q(n)*q(cresat(i))*(1./ri + Fafsw)
      else if (ri.le.roff) then
          phi=332.*q(n)*q(cresat(i))*(Aafsw*(1./ri-1./roff)
     +    + Bafsw*(roff-ri)+Cafsw*(roff**3-ri**3)
     +    + Dafsw*(roff**5-ri**5))
      else
      phi=0.0
      endif
      else 					
      phi = 332.*q(n)*q(cresat(i))/ri
      endif

      xsum = xsum + xc
      ysum = ysum + yc
      zsum = zsum + zc
      gphi(k) = gphi(k) + phi

*      if ((j.eq.2).and.(resname(n).eq.'TP4E').and.
*     +    (r0.lt.6.5)) then
*      write(6,*) n,l,r0,r,gphi(k)
*      endif

c------
c      write(6,2226) n,resname(n),atname(n),x(n),y(n),z(n),r,
c      +    q(n),phi,gphi(k) 
2226  format(I5,1x,2(A5,1x),8(f10.5,1x))
c------

2230  Continue
 
      gx = xsum*invatgrp
      gy = ysum*invatgrp
      gz = zsum*invatgrp
      gr = sqrt(gx**2 + gy**2 + gz**2)
      dist = nint(gr*div)

      if ((dist.ge.imdiv).and.(dist.le.idiv)) then

* Don't count contributions of solvent outside Born radius:
      if (gtype(k).eq.4) then
      if (gr.gt.rborn) gphi(k) = 0.
      endif

* Place the potential in the appropriate bin:
      rphi(gtype(k),dist) = rphi(gtype(k),dist) + 
     +    gphi(k)
      sphi(gtype(k)) = sphi(gtype(k)) + gphi(k)

*      if (gtype(k).eq.2) then
*      write(6,*) k,gtype(k),rphi(gtype(k),1)
*      write(6,*) k,gtype(k),gr,gphi(k),sphi(k)
*      if ((j.eq.1).and.(gtype(k).eq.4)) then
*      tot=tot+gphi(k)
*      write(6,2235) k,gr,gphi(k),tot
*2235  format(6x,I5,3(5x,f7.3))
*      endif

c------
c      write(6,'(5x,A,F10.5)') 'sphi=',sphi(gtype(k))
c------
      sphi(12) = sphi(12) + gphi(k)
      rphi(12,dist) = rphi(12,dist) + gphi(k) 
      resphi(resgrp(k),gtype(k))=resphi(resgrp(k),gtype(k))+ 
     +    gphi(k)
      resphi(resgrp(k),12)=resphi(resgrp(k),12)+ 
     +    gphi(k)
      resdist(resgrp(k),gtype(k))=resdist(resgrp(k),gtype(k))+
     +    gr
      reshit(resgrp(k),gtype(k))=reshit(resgrp(k),gtype(k))+1
      resdist(resgrp(k),12)=resdist(resgrp(k),12)+
     +    gr
      reshit(resgrp(k),12)=reshit(resgrp(k),12)+1
      if (gtype(k).le.2) then 
      sphi(7) = sphi(7) + gphi(k)
      rphi(7,dist) = rphi(7,dist) + gphi(k)
      resphi(resgrp(k),7)=resphi(resgrp(k),7)+ 
     +    gphi(k)
      resdist(resgrp(k),7)=resdist(resgrp(k),7)+
     +    gr
      reshit(resgrp(k),7)=reshit(resgrp(k),7)+1
*      write(6,*) resgrp(k),resphi(resgrp(k),1),resphi(resgrp(k),2),
*     +    resphi(resgrp(k),7)
      endif

      if (gtype(k).le.3) then
      sphi(8) = sphi(8) + gphi(k)
      rphi(8,dist) = rphi(8,dist) + gphi(k)
      resphi(resgrp(k),8)=resphi(resgrp(k),8)+ 
     +    gphi(k)
      resdist(resgrp(k),8)=resdist(resgrp(k),8)+
     +    gr
      reshit(resgrp(k),8)=reshit(resgrp(k),8)+1
      endif
      if (gtype(k).le.4) then
          sphi(9) = sphi(9) + gphi(k)
          rphi(9,dist) = rphi(9,dist) + gphi(k)
          resphi(resgrp(k),9)=resphi(resgrp(k),9)+ 
     +    gphi(k)
          resdist(resgrp(k),9)=resdist(resgrp(k),9)+
     +    gr
          reshit(resgrp(k),9)=reshit(resgrp(k),9)+1
      endif
      if (gtype(k).ne.4) then
          sphi(10) = sphi(10) + gphi(k)
          rphi(10,dist) = rphi(10,dist) + gphi(k)
          resphi(resgrp(k),10)=resphi(resgrp(k),10)+ 
     +    gphi(k)
          resdist(resgrp(k),10)=resdist(resgrp(k),10)+
     +    gr
          reshit(resgrp(k),10)=reshit(resgrp(k),10)+1
      endif
      if ((gtype(k).eq.5).or.(gtype(k).eq.6)) then
          sphi(11)= sphi(11) + gphi(k)
          rphi(11,dist) = rphi(11,dist) + gphi(k)
          resphi(resgrp(k),11)=resphi(resgrp(k),11)+ 
     +    gphi(k)
          resdist(resgrp(k),11)=resdist(resgrp(k),11)+
     +    gr
          reshit(resgrp(k),11)=reshit(resgrp(k),11)+1
      endif

      endif

      gphi(k) = 0.0

2220  Continue

912   Continue

      else

      Do 220 k=1,numgrp

          xsum = 0.
          ysum = 0.
          zsum = 0.
          gphi(k) = 0.
c------
c      print *, 'atomgrp(',k,')=',atomgrp(k)
          invatgrp= 1./float(atomgrp(k))
c      print *, 'invatgrp=',invatgrp
c------
          tx_flag = .FALSE.
          ty_flag = .FALSE.
          tz_flag = .FALSE.

      Do 230 l=1,atomgrp(k)
          n = n + 1

          if (resname(n).eq.cresid) then
c          write(6,*) k,l,n,resname(n),cresid
              phi = 0.
              goto 230
          endif
c          write(6,*) k,l,n,resname(n),cresid
 
      xc = x(n) - x(catnum)
      yc = y(n) - y(catnum)
      zc = z(n) - z(catnum)

      x0 = xc
      y0 = yc
      z0 = zc
      r0 = sqrt(x0**2 + y0**2 + z0**2)

c
      write(6,2225) n,l,x(n),y(n),z(n)
c
      resid = resname(n)

      if (PBC_type.ne.'none') then
      if ((resname(n).eq.'TP4E').and.((l.eq.1).or.
     +    (tx_flag).or.(ty_flag).or.(tz_flag))) then 
      call translate(l,xc,yc,zc,lenx,leny,lenz,tx_flag,
     +    ty_flag,tz_flag,PBC_type)
      endif
      if ((resname(n).eq.'POT').or.(resname(n).eq.'CLA')) then
      call translate(l,xc,yc,zc,lenx,leny,lenz,tx_flag,
     +    ty_flag,tz_flag,PBC_type)
      endif
      if ((resid(1:3).eq.'DMF').and.((l.eq.1).or.
     +    (tx_flag).or.(ty_flag).or.(tz_flag))) then
      call translate(l,xc,yc,zc,lenx,leny,lenz,tx_flag,
     +    ty_flag,tz_flag,PBC_type)
      endif
      endif

      r = sqrt(xc**2 + yc**2 + zc**2)

      if ((j.eq.1).and.((tx_flag).or.(ty_flag).or.
     +    (tz_flag))) then
      write(6,225) n,l,x0,y0,z0,r0
      write(6,225) n,l,xc,yc,zc,r
225   format(I5,1x,I5,1x,4(f7.2,1x))
      endif

      if (Err_flag) then
          phi = ((332.*q(n)*q(catnum))**2)/(r**4)
      else
          phi = 332.*q(n)*q(catnum)/r
      endif

      xsum = xsum + xc
      ysum = ysum + yc
      zsum = zsum + zc
      gphi(k) = gphi(k) + phi

*      if ((j.eq.2).and.(resname(n).eq.'TP4E').and.
*     +    (r0.lt.6.5)) then
*      write(6,*) n,l,r0,r,gphi(k)
*      endif

*      write(6,226) n,resname(n),atname(n),x(n),y(n),z(n),r
*     +    r,q(n),phi,gphi(k) 
*226   format(I5,1x,2(A5,1x),4(f10.5,1x))
230   Continue
 
      gx = xsum*invatgrp
      gy = ysum*invatgrp
      gz = zsum*invatgrp
      gr = sqrt(gx**2 + gy**2 + gz**2)
      dist = nint(gr*div)

      if ((dist.ge.imdiv).and.(dist.le.idiv)) then

* Don't count contributions of solvent outside Born radius:
      if (gtype(k).eq.4) then
      if (gr.gt.rborn) gphi(k) = 0.
      endif

* Place the potential in the appropriate bin:
      rphi(gtype(k),dist) = rphi(gtype(k),dist) + 
     +    gphi(k)
      sphi(gtype(k)) = sphi(gtype(k)) + gphi(k)

*      if (gtype(k).eq.2) then
*      write(6,*) k,gtype(k),rphi(gtype(k),1)
*      write(6,*) k,gtype(k),gr,gphi(k),sphi(k)
*      if ((j.eq.1).and.(gtype(k).eq.4)) then
*      tot=tot+gphi(k)
*      write(6,235) k,gr,gphi(k),tot
*235   format(6x,I5,3(5x,f7.3))
*      endif

      sphi(12) = sphi(12) + gphi(k)
      rphi(12,dist) = rphi(12,dist) + gphi(k) 
      resphi(resgrp(k),gtype(k))=resphi(resgrp(k),gtype(k))+ 
     +    gphi(k)
*      resphi(resgrp(k),12)=resphi(resgrp(k),12)+ 
*     +    gphi(k)
      resdist(resgrp(k),gtype(k))=resdist(resgrp(k),gtype(k))+
     +    gr
      reshit(resgrp(k),gtype(k))=reshit(resgrp(k),gtype(k))+1
      resdist(resgrp(k),12)=resdist(resgrp(k),12)+
     +    gr
      reshit(resgrp(k),12)=reshit(resgrp(k),12)+1
      if (gtype(k).le.2) then 
      sphi(7) = sphi(7) + gphi(k)
      rphi(7,dist) = rphi(7,dist) + gphi(k)
      resphi(resgrp(k),7)=resphi(resgrp(k),7)+ 
     +    gphi(k)
          resdist(resgrp(k),7)=resdist(resgrp(k),7)+
     +    gr
      reshit(resgrp(k),7)=reshit(resgrp(k),7)+1
*      write(6,*) resgrp(k),resphi(resgrp(k),1),resphi(resgrp(k),2),
*     +    resphi(resgrp(k),7)
      endif
      if (gtype(k).le.3) then
      sphi(8) = sphi(8) + gphi(k)
      rphi(8,dist) = rphi(8,dist) + gphi(k)
*      resphi(resgrp(k),8)=resphi(resgrp(k),8)+ 
*     +    gphi(k)
      resdist(resgrp(k),8)=resdist(resgrp(k),8)+
     +    gr
      reshit(resgrp(k),8)=reshit(resgrp(k),8)+1
      endif 
      if (gtype(k).le.4) then
          sphi(9) = sphi(9) + gphi(k)
          rphi(9,dist) = rphi(9,dist) + gphi(k)
*          resphi(resgrp(k),9)=resphi(resgrp(k),9)+ 
*     +    gphi(k)
          resdist(resgrp(k),9)=resdist(resgrp(k),9)+
     +    gr
          reshit(resgrp(k),9)=reshit(resgrp(k),9)+1
      endif
      if (gtype(k).ne.4) then
          sphi(10) = sphi(10) + gphi(k)
          rphi(10,dist) = rphi(10,dist) + gphi(k)
*          resphi(resgrp(k),10)=resphi(resgrp(k),10)+ 
*     +    gphi(k)
          resdist(resgrp(k),10)=resdist(resgrp(k),10)+
     +    gr
          reshit(resgrp(k),10)=reshit(resgrp(k),10)+1
      endif
      if ((gtype(k).eq.5).or.(gtype(k).eq.6)) then
           sphi(11)= sphi(11) + gphi(k)
           rphi(11,dist) = rphi(11,dist) + gphi(k)
*           resphi(resgrp(k),11)=resphi(resgrp(k),11)+ 
*     +    gphi(k)
           resdist(resgrp(k),11)=resdist(resgrp(k),11)+
     +    gr
           reshit(resgrp(k),11)=reshit(resgrp(k),11)+1
      endif

      endif

      gphi(k) = 0.0

220   Continue

      endif

      sphi(4) = sphi(4) + bphi

      do 240 i=1,12

      tphi(i) = tphi(i) + sphi(i)
      tphi2(i) = tphi2(i) + sphi(i)**2
      timephi(i,int(numframe)) = sphi(i)
*      write(6,*) sphi(i),tphi(i),tphi2(i)

cMLT07i time series of different components
      if (tres2res.eq.'y') then
      do krr=1,nres2res
      timerrphi(nresspec(krr),i,int(numframe))=
     +    resphi(nresspec(krr),i)
      enddo
      endif
cMLT07f

*      if (i.le.7) then
      do 245 k=1,resnum
      tresphi(k,i) = tresphi(k,i) + resphi(k,i)
      tresphi2(k,i) = tresphi2(k,i) + resphi(k,i)**2
*      write(6,*) resphi(k,i),tresphi(k,i),tresphi2(k,i)
      resphi(k,i) = 0.0
c      if (reshit(k,i).lt.1) reshit(k,i)=1
c      resdist(k,i)=resdist(k,i)/float(reshit(k,i))
245   continue
*      endif

      srphi(i,1) = rphi(i,1)

      if (imdiv .ge. 2) then
          sdiv=imdiv
      else
          sdiv=2
      endif

      do 250 k=sdiv,idiv

*
* I have to clear arrays again at this point.
*
      if ((dcdnum .le. 1) .and. (i .le. 1) 
     +    .and. (k .le. sdiv)) then
      do 193 count1=1,12
      do 194 count2=0,grpdim
      arphi(count1,count2) = 0.0
194   continue
193   continue
      endif

      srphi(i,k) = srphi(i,k-1) + rphi(i,k)
      if ((i.eq.4).and.(k.gt.irborn)) srphi(i,k) = srphi(i,k) +
     +    brphi(k)
      if ((i.eq.9).and.(k.gt.irborn)) srphi(i,k) = srphi(i,k) +
     +    brphi(k)
      if ((i.eq.12).and.(k.gt.irborn)) srphi(i,k) = srphi(i,k) +
     +    brphi(k)
*      write(6,255) i,k,rphi(i,k),srphi(i,k),srphi(7,k)
*255   format(2I3,3f10.4)
      trphi(i,k) = trphi(i,k) + srphi(i,k)
      trphi2(i,k) = trphi2(i,k) + srphi(i,k)**2
      arphi(i,k) = arphi(i,k) + rphi(i,k)
      srphi(i,k-1) = 0.0
      rphi(i,k) = 0.0
250   continue
      sdiv=1

      srphi(i,idiv) = 0.0
      sphi(i) = 0.0
*      do 265 k = imdiv,idiv
*      rphi(i,k) = 0.0
*265   continue

240   Continue

210   Continue

      if (filetype.eq.'DCD') then
          close(du)
      write(6,'(A,I2)') 'Finished with DCD file #',dcdnum
      write(6,'(I6,A)') frame,' frames were analyzed.'
      write(6,'(I6,A)') endframe,' remaining frames.'
      strtframe = (frame*skipframe)+strtframe-totframe
      if ((endframe.ne.0).or.(dcdnum.lt.numdcd)) goto 195
      endif 
      write(6,*)
      write(6,'(A,I6,A)') 'A total of ',int(numframe), 
     +    ' frames were analyzed.'
      write(6,*)

      invnfrm = 1./numframe
*
*  Write electrostatic energy data to file.
*
      Do 270 i=1,12

      if (Err_flag) then
      Do 279 j=imdiv,idiv
      arphi(i,j) = sqrt(arphi(i,j)*invnfrm)
      trphi(i,j) = sqrt(trphi(i,j)*invnfrm)
      trphi2(i,j) = sqrt(trphi2(i,j)*invnfrm)
      sigrphi(i,j) = sqrt(trphi2(i,j) - trphi(i,j)**2)
279   continue

      tphi(i) = sqrt(tphi(i)*invnfrm)
      tphi2(i) = sqrt(tphi2(i)*invnfrm)
      sigma = sqrt(tphi2(i) - (tphi(i))**2)
*      write(6,'(I5,2x,f10.4,2x,f10.4)') i,tphi(i),sigma
      write(6,284) grptitle(i),tphi(i),sigma
284   format(A30,2f14.5)

      do 299 k=1,resnum  
      tresphi(k,i) = sqrt(tresphi(k,i)*invnfrm)
      tresphi2(k,i) = sqrt(tresphi2(k,i)*invnfrm)
      sigresphi(k,i) = sqrt(tresphi2(k,i) - 
     +    tresphi(k,i)**2)
      if (reshit(k,i).lt.1) reshit(k,i)=1
      resdist(k,i)=resdist(k,i)/float(reshit(k,i))
c--------
	print *, 'resdist(',k,',',i,')=',resdist(k,i)
c--------
c      resdist(k,i)=resdist(k,i)*invnfrm
299   continue
      else
          Do 280 j=imdiv,idiv
          arphi(i,j) = arphi(i,j)*invnfrm
          trphi(i,j) = trphi(i,j)*invnfrm
          trphi2(i,j) = trphi2(i,j)*invnfrm
          sigrphi(i,j) = sqrt(trphi2(i,j) - trphi(i,j)**2)
280   continue

      tphi(i) = tphi(i)*invnfrm
      tphi2(i) = tphi2(i)*invnfrm
      sigma = sqrt(tphi2(i) - (tphi(i))**2)
*      write(6,'(I5,2x,f10.4,2x,f10.4)') i,tphi(i),sigma
      write(6,285) grptitle(i),tphi(i),sigma
285   format(A30,2f14.5) 

      do 300 k=1,resnum  
      tresphi(k,i) = tresphi(k,i)*invnfrm
      tresphi2(k,i) = tresphi2(k,i)*invnfrm
      sigresphi(k,i) = sqrt(tresphi2(k,i) - 
     +    tresphi(k,i)**2)
      if (reshit(k,i).lt.1) reshit(k,i)=1
      resdist(k,i)=resdist(k,i)/float(reshit(k,i))
c      resdist(k,i)=resdist(k,i)*invnfrm
300   continue
      endif

270   continue

      do 310 j=1,idiv
      rdist = anint(float(j))*invdiv
      write(14,305) rdist,arphi(1,j),arphi(2,j),
     +    arphi(3,j),arphi(4,j),arphi(5,j),arphi(6,j),
     +    arphi(7,j)
305   format(8(f10.4,1x))
      write(13,305) rdist,trphi(1,j),trphi(2,j),
     +    trphi(3,j),trphi(4,j),trphi(5,j),trphi(6,j),
     +    trphi(7,j)
      write(15,305) rdist,trphi(7,j),trphi(8,j),
     +    trphi(9,j),trphi(10,j),trphi(11,j),trphi(12,j),
     +    brphi(j)
*307      format(5(f10.4,1x))
      if (filetype.eq.'DCD') then 
      write(18,305) rdist,sigrphi(1,j),sigrphi(2,j),
     +    sigrphi(3,j),sigrphi(4,j),sigrphi(5,j),sigrphi(6,j),
     +    sigrphi(7,j)
          write(19,305) rdist,sigrphi(7,j),sigrphi(8,j),
     +    sigrphi(9,j),sigrphi(10,j),sigrphi(11,j),sigrphi(12,j)
      endif
310   continue

      do 320 k=1,resnum
      write(16,314) float(k),tresphi(k,1),tresphi(k,2),
     +    tresphi(k,3),tresphi(k,4),tresphi(k,5),tresphi(k,6),
     +    tresphi(k,7),tresphi(k,8),tresphi(k,9),tresphi(k,10),
     +    tresphi(k,11),tresphi(k,12)
          write(17,314) float(k),resdist(k,1),resdist(k,2),
     +    resdist(k,3),resdist(k,4),resdist(k,5),resdist(k,6),
     +    resdist(k,7),resdist(k,8),resdist(k,9),resdist(k,10),
     +    resdist(k,11),resdist(k,12)
          if (filetype.eq.'DCD') then
            write(12,315) float(k),sigresphi(k,1),sigresphi(k,2),
     +    sigresphi(k,3),sigresphi(k,4),sigresphi(k,5),
     +    sigresphi(k,6),sigresphi(k,7)
      endif
314   format(13f9.3)
315   format(8f9.3)
320   continue

      if (filetype.eq.'DCD') then
cMLT07i write title for different files
      if (tres2res.eq.'y') then
      do i=1,nres2res
      nopen=100-i
      write(nopen,*)'ES between center-residue and residue ',
     +    nresspec(i) 
      write(nopen,*)'Timestep, Backbone, BB+PSC, BB+PSC+SOLV', 
     +    'BB+PSC+CSC+CI, CSC+CI, TOTAL'
c      + BB+PSC+SOLV+CI, CSC+CI, TOTAL'!a typo corrected on 9/15/08 MLT
c      fluct
      ftimerrphi(nresspec(i),12)=0.
      enddo
      endif
cMLT07f
      do 330 t=1,int(numframe)
      ts = float(timestep(t))
      write(20,325) ts,timephi(1,t),timephi(2,t),
     +    timephi(3,t),timephi(4,t),timephi(5,t),timephi(6,t),
     +    timephi(7,t)
      write(21,325) ts,timephi(7,t),timephi(8,t),
     +    timephi(9,t),timephi(10,t),timephi(11,t),timephi(12,t)
cMLT07i write time series of different components
      if (tres2res.eq.'y') then
      do i=1,nres2res
      nopen=100-i
c      write(nopen,329) ts,timerrphi(nresspec(i),1,t),
c    +    timerrphi(nresspec(i),2,t), timerrphi(nresspec(i),3,t),
c    +    timerrphi(nresspec(i),4,t), timerrphi(nresspec(i),5,t),
c    +    timerrphi(nresspec(i),6,t), timerrphi(nresspec(i),7,t),
c    +    timerrphi(nresspec(i),8,t), timerrphi(nresspec(i),9,t),
c    +    timerrphi(nresspec(i),10,t),timerrphi(nresspec(i),11,t),
c    +    timerrphi(nresspec(i),12,t)
      write(nopen,327) ts,timerrphi(nresspec(i),7,t),
     +  timerrphi(nresspec(i),8,t), timerrphi(nresspec(i),9,t),
     +  timerrphi(nresspec(i),10,t),timerrphi(nresspec(i),11,t),
     +  timerrphi(nresspec(i),12,t)
c      fluct
      ftimerrphi(nresspec(i),12)=ftimerrphi(nresspec(i),12)+
     +(tresphi(nresspec(i),12)-timerrphi(nresspec(i),12,t))**2
      enddo
      endif
327   format(f10.1,1x,6(f10.3,1x))
329   format(f10.1,1x,12(f10.3,1x))
cMLT07f
325   format(f10.1,1x,7(f10.3,1x))
330   continue
cMLT07i Fluct of total
      if (tres2res.eq.'y') then
      write(6,*) ' '
      write(6,*)'Ave and Fluct of ES between center-res and'
      do krr=1,nres2res
      write(6,*) '  reside ',nresspec(krr),': ',
     +     tresphi(nresspec(krr),12),', ',
     +     sqrt(ftimerrphi(nresspec(krr),12)/real(numframe))
      enddo
      endif
cMLT07f
      endif

      if (SOLV_flag) then
      if (SOLV_type.eq.'BORN') then
      do 289 i=1,ctat
      tvol=tvol+avol(i)
289   Continue
      write(6,287)
      write(6,286) 'BORN ENERGY',bphi
      write(6,286) 'BORN CHARGE',Bchg
      write(6,*)
      write(6,283) 'Non-solvent Volume = ',tvol
283   format (A21,9x,f14.5)
286   format(A11,19x,f14.5) 
287   format(80('_'))

      else if (SOLV_type.eq.'NUMI') then
      write(6,287)
      write(6,286) 'SOLV ENERGY',bphi
      write(6,286) 'CEN. CHARGE',Cchg
      write(6,*)

      else if (SOLV_type.eq.'NUMB') then
      write(6,287)
      write(6,286) 'SOLV ENERGY',bphi
      write(6,286) 'CEN. CHARGE',Cchg
      write(6,286) 'BORN CHARGE',Bchg
      write(6,*)
      endif

      endif

      write(6,*)
      write(6,*) 'ESNRG complete!'
*      write(6,*) 'Total from water = ',tot

      write(6,*)
      write(6,*) ' Assign Group time = ',agtime
      write(6,*) '     Read VdW time = ',vtime
      write(6,*) 'Topology Read time = ',ttime

9999      END
