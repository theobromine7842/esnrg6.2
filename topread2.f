*
* This subroutine is designed to extract group , atom type, and
* charge information from a CHARMM  parameter file.
*
	 SUBROUTINE TOPREAD2
**
* Set Up Variables
**
*
*
	INCLUDE 'dim.inc'

	real num4
      real acharge(a1dim,a2dim),q(atmax),gcharge(a1dim,a2dim)

      integer i,j,k,m,n
      integer nres,atomres(a1dim),grp(atmax),numgrp
      integer grpatom(a1dim,a2dim),atomgrp(atmax)
      integer gtype(atmax),grpnum(a1dim),resgrp(atmax)
	integer spcindx1,spcindx2,spcindx3,spcindx4
	integer gaplen,NTTL,TOPV,NNTA

	character*1 crdast
	character*4 akind(a1dim,a2dim),atype(atmax)
      character*4 tatid(a1dim,a2dim),tresid(a1dim)
      character*9 word1,word2,word3,word4
	character*30 grptitle(12)
	character*80 topfile,dummy
	character*132 line

      logical EOF_flag,FIND_flag

c Variables for timing functions
      real*4 timer,crap,agtime,vtime,ttime

	common /TOP/TOPV,NNTA,topfile
	common /GROUP/nres,atomres,grp,numgrp,grpatom,atomgrp,
     +   gtype,grpnum,resgrp,acharge,q,gcharge,grptitle,tatid,
     +   tresid,akind,atype
      common /TIMING/agtime,vtime,ttime

      
*
* Clear variables
*
	word1='         '
	word2='         '
	word3='         '
	word4='         '
      NTTL=0
*
* Read topology file.  Initialize residues, groups, and charges.
*
	write(6,*)
        write(6,*) 'Reading topology file...'
        open(1,form='formatted',access='sequential',
     +  file=topfile,status='old')
c--------------
c Add routine to distinguish between All-atom (TOPV=22)
c  and United-Atom (TOPV=21)
c--------------
cori        Do 100 While (.not. EOF_flag)
cori next line added by MLTAN 03-08-05
         do i=1,100
         read(1,'(A1)',end=110) crdast
         if (crdast .eq. '*' ) then
          NTTL = NTTL + 1
         else
          EOF_flag =.TRUE.
         endif
cori next line added by MLTAN 03-08-05
        enddo
100     Continue
	Goto 111
110     EOF_flag =.TRUE.
        Continue

111     Rewind 1
        Do 120 i=1,NTTL
          read(1,'(A80)') dummy
          write(6,'(1x,A80)') dummy
120     Continue
        read(1,'(I6)',end=122) TOPV
        write(6,'(I6)') TOPV
	Goto 125
122     EOF_flag =.TRUE.
        Continue

125     EOF_flag=.FALSE. 

	if (TOPV.eq.22) then
         NNTA=6
        else
         NNTA=5
        endif
c--------------
10	line(1:33)='                                 '
	line=line(1:33)//line(1:33)//line(1:33)//line(1:33)
	read(1,15,end=99) line(1:130)
15      format(A130)
	line(131:132)='! '

	if ((line(1:5).eq.'     ').or.
     +      (line(1:1).eq.'!').or.
     +      (line(1:1).eq.'*')) goto 10

	spcindx1=index(line,' ')
	word1=line(1:spcindx1-1)
	gaplen=1
16	spcindx2=index(line(spcindx1+gaplen:),' ')
c	if (spcindx2.eq.spcindx1+gaplen) then
	if (spcindx2.eq.1) then
	 gaplen=gaplen+1
	 goto 16
	else if (spcindx2.eq.0) then
	 word2='    '
	else
	 spcindx2=spcindx1+gaplen+spcindx2-1
	 word2=line(spcindx1+gaplen:spcindx2)
	endif
c---------
c        write(6,17) word1,word2
17	format(2(A,'~'))
c---------
18      if ((word1.eq.'RESI').or.(word1.eq.'PRES')) then
          nres=nres+1
          tresid(nres)=word2
*          write(6,'(I4,1X,A4)') nres,tresid(nres)
20	   line(1:33)='                                 '
	   line=line(1:33)//line(1:33)//line(1:33)//line(1:33)
	   read(1,15,end=99) line(1:130)
	   line(131:132)='! '

	   if ((line(1:1).eq.' ').or.
     +      (line(1:1).eq.'!').or.
     +      (line(1:1).eq.'*')) goto 20

	   spcindx1=index(line,' ')
	   word1=line(1:spcindx1-1)
	   gaplen=1
26	   spcindx2=index(line(spcindx1+gaplen:),' ')
c	   if (spcindx2.eq.spcindx1+gaplen) then
	   if (spcindx2.eq.1) then
	    gaplen=gaplen+1
	    goto 26
	   else if (spcindx2.eq.0) then
	    word2='    '
	    word3='    '
	    word4='        '
	   else
	    spcindx2=spcindx1+gaplen+spcindx2-1
	    word2=line(spcindx1+gaplen:spcindx2)
	    gaplen=1
27	    spcindx3=index(line(spcindx2+gaplen:),' ')
c	    if (spcindx3.eq.spcindx2+gaplen) then
	    if (spcindx3.eq.1) then
	     gaplen=gaplen+1
	     goto 27
	    else if (spcindx3.eq.0) then
	     word3='    '
	     word4='        '
	    else
	     spcindx3=spcindx2+gaplen+spcindx3-1
	     word3=line(spcindx2+gaplen:spcindx3)
	     gaplen=1
28	     spcindx4=index(line(spcindx3+gaplen:),' ')
c	     if (spcindx4.eq.spcindx3+gaplen) then
	     if (spcindx4.eq.1) then
	      gaplen=gaplen+1
	      goto 28
	     else if (spcindx4.eq.0) then
	      word4='        '
	     else
	      spcindx4=spcindx3+gaplen+spcindx4-1
	      word4=line(spcindx3+gaplen:spcindx4)
	     endif
	    endif
	   endif
c---------
c          write(6,25) word1,word2,word3,word4
25        format(3(A4,'~'),A8,'~')  
c---------
          if (word1(1:4).eq.'GROU') then
            grpnum(nres) = grpnum(nres) + 1
            goto 20
          else if (word1(1:4).eq.'ATOM') then 
            atomres(nres)=atomres(nres) + 1
            if (atomres(nres).gt.a2dim) then
              write(6,'(A,I4,A)') 'Residue ',tresid(nres),
     +          ' contains more'
              write(6,'(A,I4,A)') 'than ',a2dim,' atoms.'
              goto 10
            endif
            grpatom(nres,atomres(nres)) = grpnum(nres)
            if (tresid(nres).eq.'CTER') then
             grpatom(nres,atomres(nres)) = 20
            endif
            tatid(nres,atomres(nres)) = word2(1:4)
	    akind(nres,atomres(nres)) = word3(1:4)
            call numread(word4,num4)
            acharge(nres,atomres(nres)) = num4
            gcharge(nres,grpnum(nres)) = gcharge(nres,grpnum(nres)) + 
     +       num4
c--------
c             write(6,35) nres,tresid(nres),atomres(nres),
c      +       acharge(nres,atomres(nres)),grpnum(nres),
c      +       gcharge(nres,grpnum(nres))
35          format(I4,1X,A5,I4,1X,F7.4,1x,I4,1x,f7.4)
c--------
            goto 20
          else if (word1(1:4).eq.'DELE') then
            goto 20
          else if ((word1(1:4).eq.'RESI').or.
     +             (word1(1:4).eq.'PRES')) then
            goto 18
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
c	print *,'leaving TOPREAD2'

       ttime=ttime

	 Return
	 End
