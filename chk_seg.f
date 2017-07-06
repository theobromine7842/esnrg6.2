	SUBROUTINE CHK_SEG(segname,n)

        INCLUDE 'dim.inc'

	integer i,n
	character*4 segname(atmax)

	do 10 i=1,n,10
	write(6,20) 'Segment name of atom ',i,' is ',segname(i)
20	format(5x,A,I5,A,A)
10	enddo
	return
	end
