********************************************************
	Program GOPHERAT
* v.6                                                  *
* This program is designed to calculate g(r) values of *
* one set of atoms (ex. waters) around another atom    *
* (ex. Fe). This program uses CHARMM .CRD and .DCD     *
* files as input.				       *
* Modified to accomodate rectangular PBC's.            *
* Modified again to accept additional dcd files and    *
* constrained atoms.                                   *
********************************************************
* Authors: Brian Beck and Robert Yelle  6/22/95       *
********************************************************
********************
* Set Up Variables *
********************
*
* PARAMETERS
* atmax = maximum number of atoms that can be used in DCD or CRD file
* icdim = max array length for icntrl array. Used in DCD file
* dmax  = max length of radius array * number of divisions per unit length
*
* REAL
* gr    = holds value of radial distribution
* rmax  = half the box length (actual rmax = sqrt(3)/2 * this rmax)
* lenx,leny,lenz  = length of the simulation box
* lene  = length of an edge of the simulation box
* vol   = the spherical volume from r(i-1) to r(i)
* rho   = the number of "hits" per vol
* volmax= total volume out to sqrt(3)/2 * rmax
* rhomax= total number of g(r) atoms in volmax
* div   = the number of divisions per unit length
*         (i.e. div=10 has a spacing a 1/10 between points)
* x, y, z = the coordinates of any atom
* r	= the spherical radius of any atom
* xc, yc, zc, rc = same as x, y, z, r after centering around
*		   specifc atom (ex. Fe)
*
* INTEGER
* i, j, k = loop counter variables
* dist   = radius * number of division per unit length
* oresnum= overall residue number of any atom from CRD file
* resnum = segment residue number of any atom from CRD file
* atnum  = atom number of any atom from CRD file
* cresnum = residue number of center atom from CRD file
* ctat   = Total number of atoms from CRD file
* catnum = atom number of center atom from CRD file
* numgratT = Total number of atoms to use in g(r) calculation
* icntrl = array to hold various simulation parameters from DCD file.
* i10    = number of title lines in DCD file
* dtat   = Total number of atoms from DCD file
* NCTL   = number of title lines from CRD file
* Idiv   = the number of division per unit length (input value here)
* FRAME_NUM = Number of DCD frames to use if not default
*
* CHARACTER
* crdast = hold value of 1st character in each line of CRD file
* FRAME_QUEST = Y/N question variable
* PBC_switch = turns PBC's on or off.
* crdfile  = name of the CHARMM coordinate file
* dcdfile  = name of the CHARMM dnamics coordinate file
* gofrfile = name of the output [r vs. g(r)] file.
* title    = holds titles from DCD files
* dummy    = hold 1 title line from CRD files
* catid    = central atom type (ex. FE )
* segid    = segment name for any atom
* resid    = residue name for any atom
* atid     = atom type for any atom
* cseg     = segment name for central atom
* gratid   = atom type for atoms to be used in g(r) calculation
* ftype    = value from DCD file which tells whether file is DCD or DVL
* cresid   = residue name for central atom
* grseg    = segment name for atoms to be used in g(r) calculation
*
* BOOLEAN
* EOF_flag   = .TRUE. if EOF reached
* skip_flag  = .FALSE. if atom is to be used in g(r) calculation
* PBC_flag   = .TRUE. uses OCTAPBC or CUBEPBC subroutine to reflect
*              atoms across cube or octahedral periodic boundaries.
*
	integer atmax,icdim,dmax
	parameter (atmax=10000,icdim=20,dmax=50000)
	real gr(dmax),rmax,lenx,leny,lenz,lene,vol,rho,pi
	real phi(atmax),theta(atmax),rhomax,volmax,div
	real*4 r(atmax),x(atmax),y(atmax),z(atmax)
	real*4 xc(atmax),yc(atmax),zc(atmax),N,Ntot,d
	integer i,j,k,resnum,oresnum,atnum,cresnum,ctat,catnum,numgratT
	integer icntrl(icdim),i10,dist(dmax),dtat,NCTL,Idiv,FRAME_NUM
        integer natfix,numframe,totframe,endframe,dcdnum,numdcd,du
	integer ttlframes
	character*1 crdast,FRAME_QUEST
	character*3 PBC_switch
	character*80 crdfile,dcdfile,gofrfile,title(20),dummy
	character*4 catid,segid,resid,atid,cseg,gratid,ftype
	character*4 cresid,grseg,PBC_type
	logical EOF_flag,skip_flag(atmax),PBC_flag
*
*	Clear arrarys and variables
*
	pi = acos(-1.0)
	i=0
	j=0
	k=0
	lenb=0.0
	NCTL = 0
	atnum = 0
	oresnum = 0
	numgratT = 0
	rhomax = 0.0
	EOF_flag =.FALSE.
	PBC_flag =.TRUE.
	cresid = '    '
	Do 18 i=0,dmax
	 dist(i) = 0
	 gr(i)   = 0.0
18	Continue
	
*
* Ask user for names of input and output files
*
        write(6,'(5x,A)') ' '
        write(6,'(5x,A)') 'Welcome to GOPHER v. 4'
        write(6,'(5x,A)') ' '
        write(6,'(5x,A)') 'Enter the name of the CHARMM CRD file:'
        read(5,'(A)') crdfile
	write(6,'(5x,A)') crdfile
        open(1,form='formatted',access='sequential',
     +  file=crdfile,status='old')
        write(6,'(5x,A)') ' '
8        write(6,'(5x,A)') 'Enter the number of CHARMM DCD files:'
         read(5,'(I2)') numdcd
         write(6,'(5x,I2)') numdcd
         write(6,'(5x,A)') ' '
         if (numdcd.gt.20) then
          write(6,'(A)') 'No more than 20 dcd files can be read!'
          goto 8
         endif
         do 9 i=1,numdcd
          write(6,'(A,I2)') 'Enter the name of CHARMM DCD file #',i
          read(5,'(A)') dcdfile
          write(6,'(5x,A)') dcdfile
          open(i+10,form='unformatted',access='sequential',
     +     file=dcdfile,status='old')
9        continue
        write(6,'(5x,A)') ' '
        write(6,'(5x,A)') 'Enter the name for the g(r) output file:'
        read(5,'(A)') gofrfile
	write(6,'(5x,A)') gofrfile
	open(2,form='formatted',access='sequential',
     +  file=gofrfile,status='new')
	write(6,'(5x,A)') ' '
*
* Ask user for the maximum box length
*
	write(6,'(1x,A)') ' '
	write(6,'(5x,A)') 'Please enter the box dimensions (x,y,z):'
	read(5,'(3(F9.4))') lenx,leny,lenz
*        read(5,'(f9.4)') lenx
*        read(5,'(f9.4)') leny
*        read(5,'(f9.4)') lenz
	write(6,'(3(F9.4))') lenx,leny,lenz
	write(6,'(1x,A)') ' '
*
* Ask user for and Radius to calculate G(r) for and 
* for the Number of Divisions per unit length
*
        write(6,'(1x,A)') 'Please enter the maximum radius for ',
     +   'calculating G(r) for:'
        read(5,'(f9.4)') rmax
	write(6,'(1x,A)') ' '
	write(6,'(5x,2(A))') 'Please Enter the Number of Divisions',
     + ' per unit length : (REAL)'
        read(5,'(F10.5)') div
        write(6,'(5x,F10.5)') div
*
* Ask user for atom to center around and find atom in coordinate file.
*
	write(6,'(5x,A)') 'What atom would you like to center the g(r)'
	write(6,'(5x,A)') 'calculation around :'
	write(6,'(5x,A)') '(Please use charmm CRD format.)'
    	Write(6,'(1x,A)') 'SEGID:'
        read(5,'(A4)') cseg
    	Write(6,'(1x,A)') 'RESNUM:'
        read(5,'(I4)') cresnum
    	Write(6,'(1x,A)') 'ATOM:'
        read(5,'(A4)') catid
	Write(6,'(1x,A4,1x,I4,1x,A4)') cseg,cresnum,catid
	write(6,'(5x,A)') ' '
	EOF_flag = .FALSE.
	Do 1 While (.not. EOF_flag)
	 read(1,'(A1)',end=2) crdast
	 if (crdast .eq. '*' ) then
	  NCTL = NCTL + 1
	 else
	  EOF_flag =.TRUE.
	 endif
1       Continue
2       EOF_flag =.TRUE.
        Continue
        Rewind 1
        EOF_flag =.FALSE.
        Do 3 While (.not. EOF_flag)
         Do 4 i=1,NCTL
          read(1,'(A)',end=5) dummy
          Write(6,'(1x,A)') dummy
4        Continue
         read(1,'(1x,I4)',end=5) ctat
  	 Write(6,'(1x,I4)') ctat

	 Do 6 i=1,ctat
	  read(1,7,end=80) atnum,oresnum,resid,atid,segid,resnum
7	  format(2(1x,I4),2(1x,A4),31x,A4,1x,I4)
	  if ((segid.eq.cseg).and.(resnum.eq.cresnum).and.
     + (atid.eq.catid)) then
           catnum = atnum
	   EOF_flag=.TRUE.
          endif
6	 Continue
3       Continue
	goto 10
80       write(6,90) 'EOF encountered before ',
     + cseg,cresnum,catid,' found!'
90       format(3x,A,A4,I4,A4,A)
	EOF_flag=.TRUE.
        Continue
	GOTO 9999
5       write(6,'(3x,A)') ' EOF before coordinate data found!'
    	EOF_flag=.TRUE.
        Continue
	GOTO 9999
10	Write(6,'(5x,A4, A)') catid, ' has been found.'
	Write(6,'(A)') ' '
	Rewind 1
	EOF_flag =.FALSE.
*
* Ask user for the type of atom to calculate g(r) for
* and flag them not to be skipped.
*
	write(6,'(5x,A,A)') 'Now enter the atom type you ',
     +  'would like to calculate'
	write(6,'(5x,A)') 'the g(r) for. Please use CHARMM format: '
	write(6,'(5x,A)') 'SEGID: '
	read(5,'(A4)') grseg
	if (grseg.eq.cseg) then
	 PBC_flag=.FALSE.
	endif
	write(6,'(5x,A)') 'ATOM TYPE: '
	read(5,'(A4)') gratid
	write(6,'(2(5x,A4))') grseg,gratid
	if (.not. PBC_flag) then
	 write(6,'(3x,A)') 'PBCs have been turned off because'
	 write(6,'(3x,A4,A,A4)') gratid,' is in the same segment as ',catid
	 write(6,'(3x,A)') 'If you would like them on, enter ON now'
	 write(6,'(3x,A)') ' (All other entries will maintain default)'
	 read(5,'(A2)') PBC_switch
	 write(6,'(3x,A)') PBC_switch
	 if ((PBC_switch.eq.'ON').or.(PBC_switch.eq.'on')) then
	  PBC_flag=.TRUE.
	 else
	  PBC_flag=.FALSE.
	 endif
	else
	 write(6,'(3x,A)') 'PBCs have been turned on because'
	 write(6,'(3x,A4,A,A4)') gratid,' is not in the segment of ',catid
	 write(6,'(3x,A)') 'If you would like them OFF, enter OFF now'
	 write(6,'(3x,A)') ' (All other entries will maintain default)'
	 read(5,'(A3)') PBC_switch
	 write(6,'(3x,A)') PBC_switch
	 if ((PBC_switch.eq.'OFF').or.(PBC_switch.eq.'off')) then
	  PBC_flag=.FALSE.
	 else
	  PBC_flag=.TRUE.
	 endif
	endif

	 write(6,'(5x,A)') 'What is the shape of the box, cubic or'
	 write(6,'(5x,A)') 'truncated octahedral? (CUBE or OCTA)'
	 read(5,'(A4)') PBC_type
	 write(6,'(5x,A4)') PBC_type

       Do 11 While (.not. EOF_flag)
	 Do 12 i=1,NCTL+1
	  read(1,'(A)',end=13) dummy
12	 Continue
	 Do 14 i=1,ctat
	  read(1,7,end=13) atnum,oresnum,resid,atid,segid,resnum
	  if ((atid.eq.gratid).and.(atnum.ne.catnum).and.
     +   (segid.eq.grseg)) then
	   skip_flag(i)=.FALSE.
	   numgratT= numgratT + 1
	  else
	   skip_flag(i)=.TRUE.
	  endif
14       Continue
	 EOF_flag=.TRUE.
11	Continue
        write(6,'(5x,A,I4)') 'Total # of atom used = ',numgratT
	Goto 15
13	write(6,'(A,I4,A)') 'EOF detected before ',ctat,' atoms examined!'
	EOF_flag=.TRUE.
	Continue
	GOTO 9999
15      Write(6,'(1x,A)') ' '
*
* read dcd file type,# frames,total number atoms,etc
* 
*
*     ftype for coordinate files is "CORD"
*     ftype for velocity   files is "VELD"
*     icntrl(1)=NUMBER OF COORDINATE SETS IN FILE
*     icntrl(2)=first step written
*     icntrl(3)=FREQUENCY FOR SAVING COORDINATES
*     icntrl(4)=NUMBER OF STEPS FOR CREATION RUN
*     icntrl(5)=frequency for saving velocities
*     icntrl(6)=0
*     icntrl(7)=0
*     icntrl(8)=number of degrees of freedom
*     icntrl(9)=number of "fixed" atoms=0
*     i10=number of lines in title
*     all the titles start with "*"
*
	do 95 j=1,numdcd
	 du=j+10
	 read(du) ftype,icntrl
	 ttlframes=ttlframes + icntrl(1)
	 rewind du
95	Continue

        Write(6,'(1x,A)') ' '
        Write(6,'(3x,I5,A,A)') ttlframes,' frames detected in all ',
     +   'specified trajectory files.'
        Write(6,'(1x,A)') ' '
	write(6,'(3x,A)') 'How many frames do you wish to use?'
	read(5,'(I6)') endframe
	write(6,'(I6)') endframe
	write(6,'(3x,I6,1x,A)') endframe,'Frames being used.'

        dcdnum = 0

100     dcdnum = dcdnum+1
        du = dcdnum+10
	read(du) ftype,icntrl
	Write(6,'(A)') ftype

	if (ftype.ne.'CORD') then
         write(6,'(3x,A)') 'The DCD file you specified is not
     + a coordinate file!'
	 write(6,'(3x,A)') 'GOPHERAT ending without calculation!'
	 goto 9999
	endif

        totframe=icntrl(1)
        if (totframe.gt.endframe) totframe=endframe
        endframe=endframe-totframe

	read(du) i10,(title(j),j=1,i10)
	Do 17 j=1,i10
 	 Write(6,'(A)') title(j)
17 	Continue
	read(du) dtat
	write(6,'(I6)') dtat

        if (icntrl(9).ne.0) read(du) natfix
        natfix=icntrl(9)
        natfree=dtat-natfix
        write(6,'(A,I5)') 'The total number of atoms is ',dtat
        write(6,'(A,I5)') 'The number of fixed atoms is ',natfix

	if (dtat.ne.ctat) then
         write(6,'(A)') 'The total number of atoms from the CRD file'
         write(6,'(A)') 'does not match the total number of atoms from'
         write(6,'(A)') 'the dynamics file! '
	 write(6,'(A)') 'GOFR ending without calculation!'
	 goto 9999
	endif

	Do 22 j=1,totframe

         if (j.eq.1) then
	  read(du) (x(i),i=1,dtat)
	  read(du) (y(i),i=1,dtat)
	  read(du) (z(i),i=1,dtat)
         else
          read(du) (x(i),i=natfix+1,dtat)
          read(du) (y(i),i=natfix+1,dtat)
          read(du) (z(i),i=natfix+1,dtat)
         endif

	 Do 21 i=1,dtat
	  if (skip_flag(i)) goto 21
	  xc(i) = x(i) - x(catnum)
	  yc(i) = y(i) - y(catnum)
	  zc(i) = z(i) - z(catnum)
	  if (PBC_flag) then
	   if (PBC_type.eq.'OCTA') then
	    call OCTAPBC(i,xc,yc,zc,lenx,leny,lenz)
	   else if (PBC_type.eq.'CUBE') then
	    call CUBEPBC(i,xc,yc,zc,lenx,leny,lenz)
	   else
	    write(6,'(5x,A)') 'Incorrect PBC type entered!'
	    write(6,'(5x,A)') 'Simple CUBIC PBCs being used!'
	    call CUBEPBC(i,xc,yc,zc,lenx,leny,lenz)
	   endif
	  endif
	  r(i) = sqrt(xc(i)**2 + yc(i)**2 + zc(i)**2)
	  dist(nint(r(i)*div)) = dist(nint(r(i)*div)) + 1
21       Continue

	if (j.eq.1) then
	 write(6,'(5x,A)') 'Analyzing'
	endif
	if (mod(j,totframe).eq.0) then
	 write(6,'(A,2x,I5)') '.',j
	else
	if (mod(j,50).ne.1) then
	 if (mod(j,50).ne.0) then
	  write(6,'(A,$)') '.'
	 else
	  write(6,'(A,2x,I5)') '.',j
	 endif
	else
	 write(6,'(5x,A,$)') '.'
	endif
	endif

        numframe = numframe + 1
22	Continue

        close(du)
        write(6,*) 'Finished with DCD file #',dcdnum
        write(6,'(I5,A)') numframe,' frames were analyzed.'
        write(6,'(I5,A)') endframe,' remaining frames.'
        if ((endframe.ne.0).or.(dcdnum.lt.numdcd)) goto 100

        write(6,*)
        write(6,'(A,I5,A)') 'A total of ',numframe, 
     + ' frames were analyzed.'

        if (PBC_type.eq.'OCTA') volmax=0.5*lenx*leny*lenz
        if (PBC_type.eq.'CUBE') volmax=lenx*leny*lenz
        if (PBC_type.eq.'SPHERE') volmax=(4.*pi/3.)*rmax**3
	rhomax = numgratT/volmax

	Do 25 i=1,int(rmax*div)
	 if (i.lt. 1) then
	  vol = (4./3.)*pi*((0.5)**3)/div**3
         else
	  vol = (4./3.)*pi*((float(i)+0.5)**3 - (float(i)-0.5)**3)/div**3
	 endif
	 rho = float(dist(i))/vol
         d = float(i)/div
	 gr(i)=(rho/rhomax)/(float(numframe))
         N = vol*rhomax*gr(i)
         Ntot = Ntot + N
*         write(6,*) vol,N,Ntot
         write(2,26) d,gr(i),N,Ntot   
26	 format(4(1x,F9.4))
25	Continue

*
* Close all open units/files
*
9999    close (1)
        close (2)

*
* Goodbye message
* 
      write(6,'(5x,A)') ' '
      write(6,'(5x,A,f9.2)') 'Volume = ',volmax
      write(6,'(5x,A,I4)') 'Total number of atoms = ',numgratT
      write(6,'(5x,A,f9.7)') 'Density = ',rhomax
      write(6,'(5x,A)') ' '
      write(6,'(5x,A)') 'GOPHER complete'
*
* End of main prgm.
*
        end
*
* Begin Subroutine OCTAPBC
* This subroutine takes around data input centered about
* the atom of interest, and then makes any needed periodic
* translations to get all the g(r) atoms within the box.
* Currently OCTAPBC uses truncated octahedral periodic 
* boundary conditions.
*
	 Subroutine OCTAPBC(i,xc,yc,zc,lenx,leny,lenz)
*
* Set up Variables
*
	 integer dtat,i,atmax,icdim,dmax
	 parameter (atmax=5000,icdim=20,dmax=50000)
	 real lenx,leny,lenz, corr, r75
	 real*4 xc(atmax),yc(atmax), zc(atmax)
	 parameter (r75 = 4.0/3.0)
*
* reduce coords to unitary scale
*
	xc(i) = xc(i)/lenx
	yc(i) = yc(i)/leny
	zc(i) = zc(i)/lenz
*
* first perform a cubic translation
*
	 xc(i)= xc(i) - anint((xc(i)))
	 yc(i)= yc(i) - anint((yc(i)))
	 zc(i)= zc(i) - anint((zc(i)))
*
* Now correct for cube octahedral differences from cube
*
	corr = 0.5 * aint(r75*(abs(xc(i))+abs(yc(i))+abs(zc(i))))
*
	xc(i) = xc(i) - sign(corr, xc(i))
	yc(i) = yc(i) - sign(corr, yc(i))
	zc(i) = zc(i) - sign(corr, zc(i))
*
* return coords to proper length
*
	xc(i) = xc(i)*lenx
	yc(i) = yc(i)*leny
	zc(i) = zc(i)*lenz
*
* Return to main program
*
	 Return
	 End
*
* Begin Subroutine CUBEPBC
* This subroutine takes around data input centered about
* the atom of interest, and then makes any needed periodic
* translations to get all the g(r) atoms within the box.
* Currently CUBEPBC uses simple cubic periodic 
* boundary conditions.
*
	 Subroutine CUBEPBC(i,xc,yc,zc,lenx,leny,lenz)
*
* Set up Variables
*
	 integer dtat,i,atmax,icdim,dmax
	 parameter (atmax=5000,icdim=20,dmax=50000)
	 real volmax,lenx,leny,lenz
	 real*4 xc(atmax),yc(atmax), zc(atmax)
*
* reduce coords to unitary scale
*
	xc(i) = xc(i)/lenx
	yc(i) = yc(i)/leny
	zc(i) = zc(i)/lenz
*
* Now perform a cubic translation
*
	 xc(i)= xc(i) - anint((xc(i)))
	 yc(i)= yc(i) - anint((yc(i)))
	 zc(i)= zc(i) - anint((zc(i)))
*
* return coords to proper length
*
	xc(i) = xc(i)*lenx
	yc(i) = yc(i)*leny
	zc(i) = zc(i)*lenz
*
* Return to main program
*
	 Return
	 End
