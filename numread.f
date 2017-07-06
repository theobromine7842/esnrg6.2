*
* Number read subroutine, free format
*
        SUBROUTINE numread(word,num)

        integer d,l,n,e,inum
        real num,rnum
        character*9 word
        character*1 w
        logical s

        s=.FALSE.
        e=0
        d=0
        num=0.0000000000
        rnum=0.0000000000
        do 400 l=1,9
          w=word(l:l)
          n=ichar(w)-48
          if ((d.eq.0).and.(w.eq.'-')) s=.TRUE.
          if (w.eq.'.') d=l
          if ((d.ne.0).and.(n.ge.0).and.(n.le.9)) then
            e=l
          endif
*          write(6,*) l,w,n,num
400     continue
        n=0
        do 410 l=1,e-d
          n = ichar(word(d+l:d+l))-48
*          rnum = real(n)/(real(10.0**real(l)))
          rnum = real(n)*(0.1**l)
          num = num + rnum
*          write(6,*) l,n,rnum,num
410     continue
        n=ichar(word(d-1:d-1))-48
        if ((n.ge.0).and.(n.le.9)) num = num + n
        if (s.eqv..TRUE.) num = 0.0000000000 - num
        rnum = anint(num*10000.0)
*        num = rnum/10000.0
        num = rnum*0.0001
*        write(6,*) n,rnum,num
        return
        end
