      subroutine snapout(dir,mdir)
      implicit none
C
      include 'pom2k.c'
C-----------------------------------------
      integer mend,mdir,iday, ihour, imin
      character dir*80, day*3, hour*3, mint*3
      character dir2*120
c----------------------------------------
      mend=20190414
c
      day='000'
      iday=int(time+0.01)
      write(day(1:3),'(i3.3)')iday
c
      hour='_00'
      ihour=int((time-iday)*24+0.01)
      write(hour(2:3),'(i2.2)')ihour
c
      mint='_00'
      imin=int((time-iday)*24*60-ihour*60+0.01)
      write(mint(2:3),'(i2.2)')imin
c            
c      open(85,file='log.txt')
c      write(85,*)dir,mdir,mode
c      close(85)
c---------------------------------------
      if(mode.eq.2) then
         open(70,file=dir(1:mdir)//'/UAVAE_'//day//hour//mint,
     1        form='unformatted')
            write(70)ua,va,el,aam2d,wubot,wvbot,mend
         close(70)
c
      else if (mode.eq.3) then
c
c      dir2=trim(dir(1:mdir)//'/UVW_'//day//hour//mint)
c      write(6,*)dir2,mdir
c      pause
         open(71,file=dir(1:mdir)//'/UVW_'//day//hour//mint,
     1        form='unformatted',STATUS='NEW')
          write(71)u,v,w,mend
         close(71)
c
         open(72,file=dir(1:mdir)//'/TSR_'//day//hour//mint,
     1        form='unformatted',STATUS='NEW')
          write(72)t,s,rho,mend
         close(72)
c
         open(73,file=dir(1:mdir)//'/EtUaVa_'//day//hour//mint,
     1        form='unformatted',STATUS='NEW')
          write(73)et,ua,va,mend
         close(73)
c
         open(74,file=dir(1:mdir)//'/AamKhKm_'//day//hour//mint,
     1        form='unformatted',STATUS='NEW')
          write(74)aam,kh,km,mend
         close(74)
c
c         write(70)u,v,t,s,et,ua,va,kh,km,aam,mend
c
c         open(70,file=dir(1:mdir)//'/u/U_'//day//hour//min,
c     1        form='unformatted')
c         write(70)u,mend
c         close(70)
c         
c         open(70,file=dir(1:mdir)//'/v/V_'//day//hour//min,
c     1        form='unformatted')
c         write(70)v,mend
c         close(70)
c
c         open(70,file=dir(1:mdir)//'/w/W_'//day//hour//min,
c     1        form='unformatted')
c         write(70)w,mend
c         close(70)
c         
c         open(70,file=dir(1:mdir)//'/t/T_'//day//hour//min,
c     1        form='unformatted')
c         write(70)t,mend
c        close(70)
c         
c         open(70,file=dir(1:mdir)//'/s/S_'//day//hour//min,
c     1        form='unformatted')
c         write(70)s,mend
c         close(70)
c         
c         open(70,file=dir(1:mdir)//'/et/ET_'//day//hour//min,
c     1        form='unformatted')
c         write(70)et,mend
c         close(70)
c         
c         open(70,file=dir(1:mdir)//'/density/rho_'//day//hour//min,
c     1        form='unformatted')
c         write(70)rho,rmean,mend
c         close(70)
c         changed by gaojie 20190320         

      else if(mode.eq.4) then
         open(70,file=dir(1:mdir)//'/UVE_'//day//hour//mint,
     1        form='unformatted')
         write(70)u,v,et,kh,aam,mend
         close(70)
      end if

      return
c
      end
c
c----------------------------------------------------------------
      subroutine chop(buf,ic)
c
      character*80 buf
c
      do i = 80,1,-1
         if (buf(i:i).eq.'/') then
c            ic = i - 1
            ic=len_trim(buf)-1
            goto 9
         endif
         if (buf(i:i).ne.' ') then
c            ic=i
            ic = len_trim(buf)
            goto 9
         endif
      enddo
c
 9    return
      end
                     
