*
* Begin Subroutine TRANSLATE
* This subroutine takes around data input centered about
* the atom of interest, and then makes any needed periodic
* boundary translation necessary for EFIELD calculations.
* Currently TRANSLATE uses truncated octahedral periodic 
* boundary conditions.
*
         subroutine translate(l,xc,yc,zc,lenx,leny,lenz,
     +     tx_flag,ty_flag,tz_flag,PBC_type)
*
* Set up Variables
*
         real lenx,leny,lenz,corr,r75
         real xc,yc,zc,x0,y0,z0
         integer l,dx,dy,dz
         logical tx_flag,ty_flag,tz_flag,txyz_flag
         character*4 PBC_type
         parameter (r75 = 4.0/3.0)
*
* Apply truncated octahedral translations by the RB method,
* if necessary.  The purpose of the logical variables is to
* ensure that the water hydrogen atoms undergo the same
* translation as the water oxygen (l=1).
*
        txyz_flag=.FALSE.
        x0 = xc
        y0 = yc
        z0 = zc

        if ((tx_flag).and.(ty_flag).and.(tz_flag)) txyz_flag=.TRUE.

c        if (txyz_flag) write(6,*) 'XYZ trans'
c	change from sign to anint in case more than 2 boxl?
        if ((abs(xc).gt.lenx*0.5).or.(tx_flag))
c     +    xc = xc - sign(lenx,xc)
     +    xc = x0 - lenx*anint(x0/lenx)
        if ((abs(yc).gt.leny*0.5).or.(ty_flag))
c     +    yc = yc - sign(leny,yc)
     +    yc = y0 - leny*anint(y0/leny)
        if ((abs(zc).gt.lenz*0.5).or.(tz_flag))
c     +    zc = zc - sign(lenz,zc)
     +    zc = z0 - lenz*anint(z0/lenz)

c        if ((PBC_type.eq.'OCTA').or.(PBC_type.eq.'octa')) then
c        if ((((abs(xc)/lenx)+(abs(yc)/leny)+(abs(zc)/lenz)).gt.
c     +     0.75).and.(l.eq.1)) txyz_flag=.TRUE.
c        if (txyz_flag) then 
c        xc = xc - sign(lenx*0.5,xc)
c        yc = yc - sign(leny*0.5,yc)
c        zc = zc - sign(lenz*0.5,zc)
c        endif
c        endif

        dx = nint(xc-x0)
        dy = nint(yc-y0)
        dz = nint(zc-z0)

        if (l.eq.1) then
         if (dx.ne.0) tx_flag = .TRUE.
         if (dy.ne.0) ty_flag = .TRUE.
         if (dz.ne.0) tz_flag = .TRUE.
        endif
*
* Return to main program
*
         Return

         End
