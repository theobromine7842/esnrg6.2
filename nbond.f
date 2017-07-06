*
* This subroutine is designed to extract Non-bonded information
* from a CHARMM  parameter file.
*
	 SUBROUTINE NBOND
**
* Set Up Variables
**
*
*
	INCLUDE 'dim.inc'

      real Rmin,epsilon,polz,sign,num,tos
	real nbterms(nbmax,3)
      integer i,j,k,Pwc,Awc,dot,numchar,numlen,numread
	integer endnum
	integer numnbtype,nbindex(nbmax,2)
	character*80 paramfile,line
	character*4 nbtype(nbmax)
      logical EOF_flag,NB_flag,VALUE_flag
	logical SIGN_flag
	common /NBONDB/numnbtype,nbindex,nbterms,nbtype,paramfile

**
* Setup any variables necessary
**
	EOF_flag=.FALSE.
	Rmin=0.0
	epsilon=0.0
	NB_flag=.FALSE.
	numnbtype = 0
	tos=1./(2**(1.0/6.0))
	Do 1 i=1,nbmax
	 nbtype(i)='    '
	 Do 2 j=1,3
	  nbterms(i,j)=0.0
	  if (j.lt.3) nbindex(i,j)=-1
2	 Continue
1	Continue
	 
**
* open parameter file
**
	open(3,form='formatted',access='sequential',
     +  file=paramfile,status='old')
**
* Read file
**
	Do 1000 While (.not. EOF_flag)
	 read(3,'(A80)',end=1001) line
	 if ((line(1:4).ne.'NBOND').and.
     +       (line(1:4).ne.'NONB').and.
     +       (.not.NB_flag)) then
	  goto 1000
	 else
*
*	 print *, line(1:4)
*
	  NB_flag=.TRUE.

	 if ((line(1:4).eq.'NBOND').or.
     +       (line(1:4).eq.'NONB')) goto 1000

	  if ((line(1:4).eq.'NBFI').or.
     +        (line(1:4).eq.'HBON').or.
     +        (line(1:4).eq.'BOND').or.
     +        (line(1:4).eq.'THET').or.
     +        (line(1:4).eq.'ANGL').or.
     +        (line(1:3).eq.'PHI') .or.
     +        (line(1:4).eq.'DIHE').or.
     +        (line(1:4).eq.'IMPR').or.
     +        (line(1:4).eq.'IMPH').or.
     +        (line(1:3).eq.'END')) then
	    NB_flag=.FALSE.
	    goto 1000
          endif
	  if ((ichar(line(1:1)).lt.65).or.
     +        (ichar(line(1:1)).gt.90)) goto 1000
	   Pwc=index(line(1:4),'%')
	   Awc=index(line(1:4),'*')
*
*	 print *, line(1:Pwc),line(1:Awc)
*
	  numnbtype=numnbtype+1
	  nbtype(numnbtype) = line(1:4)
	  nbindex(numnbtype,1) = Pwc
	  nbindex(numnbtype,2) = Awc

	   if ((Pwc.eq.0).and.(Awc.eq.0)) then
	    nbindex(numnbtype,1) = -1
	    nbindex(numnbtype,2) = -1
	   endif
*
* Initialize Free Format Read Variables
*
	    dot=0
	    sign=1.0
	    VALUE_flag=.FALSE.
	    SIGN_flag=.FALSE.
	    numchar=0
	    numlen=0
	    numread=0
	    num=0.0
*
	     Do 1002 i=5,len(line)
	      if (numread.eq.3) goto 1002
	      if (dot.lt.1) then
	       dot=index(line(5+numlen:len(line)),'.')+dot
	       endnum=index(line(5+numlen+dot:len(line)),' ')
     +                +4 +numlen +dot
	      endif
c-----
c	print *, i,index(line(5+numlen:len(line)),'.'),
c     +           index(line(5+numlen+dot:len(line)),' ')
c	print *, 'Dot=',dot,'endnum=',endnum
c-----
	      if (i.eq.endnum) then
	       dot=0
	       numlen=i-5
	       numchar=0
	       num=sign*num
	       numread=numread+1
	       if (numread.eq.1) then
		polz=num
	        sign=1.0
*	print *, 'polz=',polz
		num=0.0
	       endif
	       if (numread.eq.2) then
		epsilon=num
		sign=1.0
*	print *, 'epsilon=',epsilon
		num=0.0
	       endif
	       if (numread.eq.3) then
		Rmin=num
	      	sign=1.0
*	print *, 'Rmin=',Rmin
		goto 1002
	       endif
	      endif	       
	      if (line(i:i).eq.' ') then
	       VALUE_flag=.FALSE.
	       dot=dot-1
*	print *, 'Space','  dot=',dot
	      endif
	      if (line(i:i).eq.'-') then
	       VALUE_flag=.TRUE.
	       sign=-1.0
	       SIGN_flag=.TRUE.
	       dot=dot-1
*	print *, 'Minus'
	      endif
	      if (line(i:i).eq.'+') then
	       VALUE_flag=.TRUE.
	       sign=1.0
	       SIGN_flag=.TRUE.
	       dot=dot-1
*	print *, 'Plus'
	      endif
	      if ((ichar(line(i:i)).gt.47).and.
     +            (ichar(line(i:i)).lt.58)) then
	       numchar=numchar+1
	       num=(ichar(line(i:i))-48)*(10.0**(float(dot-numchar-1)))
     +             + num
*	print *, 'numchar=',numchar,'num=',num
	      endif
1002	   Continue
	 nbterms(numnbtype,1) = epsilon
	 nbterms(numnbtype,2) = Rmin
! Convert Rmin to Sigma and store
	 nbterms(numnbtype,3) = Rmin*tos
c------
c 	 write(6,999) nbtype(numnbtype),nbindex(numnbtype,1),
c      +    nbindex(numnbtype,2),nbterms(numnbtype,1),
c      +    nbterms(numnbtype,2),nbterms(numnbtype,3)
999	format(A4,2(1x,I5),3(1x,F10.5))
c------

	 endif
1000    Continue
1001	EOF_flag=.TRUE.
	Continue
*
*
1099      close(3)
*
* Return to main program
*
	 Return
	 End
