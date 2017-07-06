*
* Assignment subroutine
*
        SUBROUTINE assigngrp

        INCLUDE 'dim.inc'

      integer m,n,nres,atomres(a1dim),grp(atmax),numgrp
      integer grpatom(a1dim,a2dim),atomgrp(atmax),atnum,resnum
      integer gtype(atmax),grpnum(a1dim),resgrp(atmax)
        integer sresnum

        real x(atmax),y(atmax),z(atmax)
      real acharge(a1dim,a2dim),q(atmax),gcharge(a1dim,a2dim)
      
c Variables for timing functions
      real*4 timer,crap,agtime,vtime,ttime

        character*30 grptitle(12)
      character*4 atid,resid,akind(a1dim,a2dim),atype(atmax)
      character*4 tatid(a1dim,a2dim),tresid(a1dim),segid

        common /GROUP/nres,atomres,grp,numgrp,grpatom,atomgrp,
     +   gtype,grpnum,resgrp,acharge,q,gcharge,grptitle,tatid,
     +   tresid,akind,atype
        common /COORD/atnum,resnum,sresnum,x,y,z,resid,atid,segid
      common /TIMING/agtime,vtime,ttime

c        print *, 'Starting ASSIGNGRP'
        
        
        m=1
        do 310 while ((m.le.nres).and.(tresid(m).ne.resid))
c          write(6,'(i4,1x,i4,1x,a4)') m,nres,tresid(m)(1:4)
          m = m + 1
310     continue
        if (m.gt.nres) goto 390 
c        write(6,*)
c        write(6,'(A,A4,A)') 'Residue ',resid,' found!'
        n=1
c        write(6,'(i4,1x,i4)') n,atomres(m)
c        print *, '-',tatid(m,n),'-',atid,'-'
        do 320 while (tatid(m,n).ne.atid)
*          write(6,'(i4,1x,a4)') n,tatid(m,n)
          n = n + 1
320     continue
c         write(6,'(i4,1x,i4,1x,a4,1x,a4)') n,atomres(m),tatid(m,n),
c      +   atid
        if (n.gt.atomres(m)) then
          goto 380
        endif
*        write(6,'(A,A4,A)') 'Atid ',atid,' found!'
        q(atnum) = acharge(m,n)
        atype(atnum) = akind(m,n)
        grp(atnum) = grpatom(m,n)

C-----
C        write(6,*) '===================================='
C        write(6,*) atnum,grp(atnum),grp(atnum-1),resnum,
C     +    resgrp(numgrp)
C-----

        if ((grp(atnum).eq.(grp(atnum-1))).and.
     +   (resnum.eq.resgrp(numgrp))) then
          atomgrp(numgrp) = atomgrp(numgrp) + 1
        else
          numgrp = numgrp + 1
          atomgrp(numgrp) = 1
        endif
        resgrp(numgrp) = resnum
*        write(6,*) m,n,numgrp,gtype(numgrp),grpatom(m,n),
*     +   gcharge(m,grpatom(m,n))
        if (gtype(numgrp).eq.0) then
          if (abs(gcharge(m,grpatom(m,n))).gt.0.0001) then
            gtype(numgrp) = 5
          elseif (atid.eq.'N') then
            gtype(numgrp) = 1
          elseif (atid.eq.'C') then
            gtype(numgrp) = 2
          elseif (abs(q(atnum)).gt.0.0001) then
            gtype(numgrp) = 3
          endif
        endif
        if ((resid.eq.'NTER').or.(resid.eq.'PROP').or.
     +      (resid.eq.'GLYP')) gtype(numgrp) = 5
        if (resid.eq.'CTER') gtype(numgrp) = 5
        if (resid(1:3).eq.'TP4') gtype(numgrp) = 4
        if (atid.eq.'OH2') gtype(numgrp) = 4
        if (resid(1:3).eq.'DMF') gtype(numgrp) = 4
        if (atid.eq.'NA') gtype(numgrp) = 6
        if (atid.eq.'MNA') gtype(numgrp) = 6
        if (atid.eq.'CL') gtype(numgrp) = 6
        if (atid.eq.'XCL') gtype(numgrp) = 6
        if (atid.eq.'POT') gtype(numgrp) = 6
        if (atid.eq.'CLA') gtype(numgrp) = 6

c-----
c        write(6,324) 'resi','atom','rgrp','ngrp','agrp','gtyp',
c     +   'grptitle                 ','charge','grpchg'
324     format(6(A4,1x),A25,1x,2(A6,1x))
c        write(6,325) resid,atid,resgrp(numgrp),numgrp,
c     +   atomgrp(numgrp),gtype(numgrp),grptitle(gtype(numgrp))(1:25),
c     +   q(atnum),gcharge(m,grpatom(m,n))
325     format(A4,1x,A4,1x,I4,1x,I4,1x,I4,1x,I4,1x,A,1x,f6.3,1x,f6.3)
c------
        agtime=agtime 
        goto 399

380     write(6,'(A,A4,A)') 'Atom name ',atid,' not found in topology'
        write(6,'(A,A4,A)') 'file for residue ',resid
        goto 399
390     write(6,'(A,A4,A)') 'Residue name ',resid,' not found in',
     + ' topology file.'

399     return
        end
