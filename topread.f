*
* This subroutine is designed to extract group , atom type, and
* charge information from a CHARMM  parameter file.
*
	 subroutine TOPREAD
**
* Set Up Variables
**
*
*
c         integer atmax,grpdim,a1dim,a2dim
c         parameter (atmax=50000,grpdim=1000)
c         parameter (a1dim=200,a2dim=75)

	INCLUDE 'dim.inc'
        
	real num4

      real acharge(a1dim,a2dim),q(atmax),gcharge(a1dim,a2dim)

      integer i,j,k,m,n
      integer nres,atomres(a1dim),grp(atmax),numgrp
      integer grpatom(a1dim,a2dim),atomgrp(atmax)
      integer gtype(atmax),grpnum(a1dim),resgrp(atmax)

	character*30 grptitle(12)
	character*4 akind(a1dim,a2dim),atype(atmax)
      character*4 tatid(a1dim,a2dim),tresid(a1dim)
      character*9 word1,word2,word3,word4
	character*80 topfile

      logical EOF_flag,FIND_flag

	common /TOP/topfile
	common /GROUP/nres,atomres,grp,numgrp,grpatom,atomgrp,
     +   gtype,grpnum,resgrp,acharge,q,gcharge,grptitle,tatid,
     +   tresid,akind,atype

*
* Read topology file.  Initialize residues, groups, and charges.
*
        write(6,*) 'Reading topology file...'
        open(1,form='formatted',access='sequential',
     +  file=topfile,status='old')
10      read(1,15,end=99) word1,word2 
*        write(6,15) word1,word2
15      format(A4,1x,A4)
        if ((word1.eq.'RESI').or.(word1.eq.'PRES')) then
          nres=nres+1
          tresid(nres)=word2
*          write(6,'(I4,1X,A4)') nres,tresid(nres)
20        read(1,25,end=99) word1,word2,word3,word4
*          write(6,25) word1,word2,word3,word4
25        format(3(A4,1X),2X,A8)  
          if (word1.eq.'GROU') then
            grpnum(nres) = grpnum(nres) + 1
            goto 20
          else if (word1.eq.'ATOM') then 
            atomres(nres)=atomres(nres) + 1
            if (atomres(nres).gt.a2dim) then
              write(6,*) 'Residue ',tresid(nres),' contains more'
              write(6,*) 'than 40 atoms.'
              goto 10
            endif
            grpatom(nres,atomres(nres)) = grpnum(nres)
            if (tresid(nres).eq.'CTER') then
             grpatom(nres,atomres(nres)) = 20
            endif
            tatid(nres,atomres(nres)) = word2
            call numread(word4,num4)
	    akind(nres,atomres(nres)) = word3
            acharge(nres,atomres(nres)) = num4
            gcharge(nres,grpnum(nres)) = gcharge(nres,grpnum(nres)) + 
     +       num4
*            write(6,35) nres,tresid(nres),atomres(nres),
*     +       acharge(nres,atomres(nres)),grpnum(nres),
*     +       gcharge(nres,grpnum(nres))
35          format(I4,1X,A5,I4,1X,F7.4,1x,I4,1x,f7.4)
            goto 20
          else if (word1.eq.'DELE') then
            goto 20
          endif
        endif
        goto 10
99      close(1)
*        write(6,'(I4)') nres
              
*        do 50 i=1,nres
*          do 60 j=1,atomres(i)
*            write(12,*) i,j,tresid(i),tatid(i,j),grpatom(i,j),
*     +       acharge(i,j)
*65          format(i4,1x,i4,1x,a4,1x,a4,1x,i4,1x,f7.4)
*60        continue
*          do 70 k=1, grpnum(i)
*            write(12,*) i,k,gcharge(i,k)
*75          format(i4,1x,i4,1x,f6.3)
*70        continue
*50      continue
*
* Return to main program
*
	 Return
	 End
