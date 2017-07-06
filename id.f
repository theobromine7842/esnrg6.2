*
* Subroutine to assign id:
*
        SUBROUTINE id

	INCLUDE 'dim.inc'

	real x(atmax),y(atmax),z(atmax),cenx,ceny,cenz,atpr
      integer n,catnum,sresnum,cresnum,resnum,cresat(atmax)
	integer iatpr,atnum
      integer IDCLR
      character*4 atid,catid,segid,cseg,resid,cresid
      character*4 atname(atmax),resname(atmax),AorP
	character*4 segname(atmax)
	character*1 XYZ_QUEST
	
	common /COORD/atnum,resnum,sresnum,x,y,z,resid,atid,segid
	common /AVG/cresat,iatpr
	common /IDB/n,cresnum,catnum,IDCLR,cenx,ceny,cenz,atpr,
     +   atname,resname,cseg,cresid,catid,XYZ_QUEST,AorP,segname

c	print *, 'Starting ID'

        resname(n) = resid
        atname(n) = atid
        segname(n) = segid

c	write(6,'(5x,A)') ' '
c	write(6,913) segid,cseg
c	write(6,914) sresnum,cresnum
c	write(6,915) atid,catid
c	write(6,916) XYZ_QUEST
c913	format(5x,'segid=',A4,5x,'cseg=',A4)
c914	format(5x,'sresnum=',I4,5x,'cresnum=',I4)
c915	format(5x,'atid=',A4,5x,'catid=',A4)
c916	format(5x,'XYZ_QUEST=',A1)



        if (((segid.eq.cseg).and.(sresnum.eq.cresnum)).and.
     +   ((atid.eq.catid).or.(XYZ_QUEST.eq.'y'))) then
         catnum = n
         cresid = resid
c         write(6,*) segid,resnum,atid,catnum,cresid,sresnum
	 write(6,111) catnum,resnum,atid,cresid,segid,sresnum
111	 format(2(1x,I5),3(1x,A),1x,I5)

	 if (AorP.eq.'RESI') then
        if (IDCLR.eq.0) then
	   cenx=cenx+x(n)
	   ceny=ceny+y(n)
	   cenz=cenz+z(n)
	   atpr=atpr+1.0
	   iatpr=int(atpr)
	   cresat(iatpr)=n
*	    write(6,*) atpr,x(n),cenx,y(n),ceny,z(n),cenz
        endif
	 else if (AorP.eq.'ATOM') then
	  cenx=x(catnum)
	  ceny=y(catnum)
	  cenz=z(catnum)
	  atpr=1.0
	  iatpr=int(atpr)
	  cresat(iatpr)=catnum
	 endif

        endif

        Return

        End
