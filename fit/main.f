      program main
      
      include 'parameter.f'
      
      character kname*100
      character outname*100
      character stemp*100
      integer nfit
      real period
      real dataarray(npmax)
      real beg,dt
      real ao, so, wo, tp, tg, wow
      real chi
      integer nlen, nerr, ierr
      
c Read in sac file name
      call getarg(1, kname)
c Read in center period
      call getarg(2, stemp)
      read(stemp,*) period
      wow = 6.283/period
c     print *,wow
c Read in nfit
      call getarg(3, stemp)
      read(stemp,*) nfit

c Read in output file name
      call getarg(4, outname)

      call rsac1(kname, dataarray, nlen, beg, del, npmax, nerr)

      if(nerr .NE. 0) then
          write(*,*)'Error reading in file: ',kname
      call exit(-1)
      endif


      do ii= 1, nlen
        dataarray(ii)=dataarray(ii)*1e9
      end do

      call fit(nlen,del,dataarray,wow,nfit,ao,so,wo,tp,tg,chi,ierr)

      open (unit=10,file=outname,status="unknown")

      ao = ao/1e9

      write (10,*) ao,so,wo,tp,tg,chi,ierr

      print *, ao,so,wo,tp,tg,chi,ierr

      close (10)


      end
