	subroutine fit(npts,dt,c,wow,nfit,ao,so,wo,tp,tg,chi,ierr)
c
c     subroutine which performs the waveform fitting of the Gaussian wavelet
c       for the estimation of the generalized seismological data functionals
c
c	input:  c     - correlation function
c             npts  - number of points in the correlation function
c             dt    - spacing of time points in sec
c             wow   - center frequency for determination of the fitting window
c             nfit  - parameter which controls the size of the fitting window
c	output: ao    - amplitude factor (always positive)
c             so    - frequency half-width
c             wo    - frequency of the cosine carrier function
c             tp    - phase shift of cosine carrier function
c             tg    - group shift of Gaussian
c             chi   - final value the chisq evaluation
c             ierr  - error code
c
c     we are fitting a five-parameter expression of the form:
c
c     cest = ao * exp{-0.5*[so(t - tg)]**2} * cos[wo(t - tp)]
c
c     02/27/91 - first pass, passed on the oringial subroutine wit2
c
      include 'parameter.f'
c
      parameter (ninter = 30)
      parameter (ncol = 5)
c
      character*1 ans
c
      real*4 c(npmax), cest(npmax), d(npmax)
      real*4 wt(npmax), tt(npmax)
      real*4 m(ncol), a(npmax,5), ml(ncol)
      real*4 w(ncol), u(npmax,ncol), v(ncol,ncol), x(ncol)
      real*4 cov(ncol,ncol)
      real*4 chisq(ninter), sigs(ncol)
      real*4 ao, so, wo, tp, tg, wow
      real*4 wtg,wtp,wao,wso,wwo
      real*4 a0, a1, a3
      real*4 b1, b2, b3
      real*4 tpeak, tenvel
      real*4 mhz, gamma
c
      integer*4 npts, ierr, ndum2
      integer*4 j0, j1, j2, j3
c
c     commons
c
      common /c_fit/ tpeak, tenvel, gamma, fit_len
      common /c_error/ sigs
c
c     data statements
c
      data tol1, tol2, tol3 /0.0000001, 0.0000001, 300.0/
      data itmax /30/
      data tpi /6.2831853071796/
      data index, index2 /1, 0/
c
      ierr = 1
      alpha = 0.0
c
      mhz = 1000./tpi
c
c     make a starting estimate of these five parameters to begin the inversion
c
c     location the largest peak and its amplitude
c
      call maxsp(c, 1, npts, a0, j0, a1, j1, a3, j3, 1)
      ndum2 = npts/2
      tpeak = float(j0 - 1 - ndum2)*dt
      t1 = float(j1 - 1 - ndum2)*dt
      t3 = float(j3 - 1 - ndum2)*dt
c
c     location of the envelope peak and its amplitude
c
      call envelop(npts, c, cest)
      call maxsp(cest, 1, npts, b2, j2, b1, j1, b3, j3, 1)
      do ii = 1, npmax
        cest(ii) = 0.0
      end do
      tenvel = float(j2 - 1 - ndum2)*dt
c
c     ao
c
      ao = a0
c
c     wo and so
c
      so = wow * 0.1

      if (so .ne. 0.) then
        so2 = so**2
      else
        so2 = -(alog(abs(a1/a0))/(t1**2) + alog(abs(a3/a0))/(t3**2))
        so = sqrt(so2)
      endif
c      if (wo .eq. 0.) then
        wo = tpi/(t3 - t1)
c      endif
c     print *, 'so=',so
c     print *, 'wo=',wo
c
c     tp and tg
c
      tp = tpeak
      tg = tenvel
c
c     now form the fitting window - lots of bookkeeping required
c
c     nfil is the desired length of the Hanning window, but correponds to 1/2 the
c       total lenth of the taper
c
      nt = nfit*tpi/wow/dt
      nfil = 2*nt + 1
      wn = tpi/(2*nfil*dt)
      fit_len = nfil*dt
      print*,' length of fitting window: ', fit_len
c
c     one-sided offset
c
      npad = j0 - nt - 1
c
c     check for a window starting before or after the cross-correlation
c
      if (npad .lt. 0) then
        print*,' window starts before correlation '
        ierr = -1
        return
      endif
c
c     compute the time offset and store it
c
      ipad = ndum2 + 1
      do ii = 1, nfil
        tt(ii) = float(ii - ipad + npad)*dt
        wt(ii) = cos(wn*float(ii - nt - 1)*dt)**2
      end do
c
c     fill the array for the estimated values
c
      m(1) = ao
      m(2) = wo
      m(3) = so2
      m(4) = tp
      m(5) = tg

      wao=10000
      wwo=1
      wso=1
      wtp=1/wo
      wtg=1/so2
c
c     loop of iterations
c
      iter = 0
c
100   continue
c
      iter = iter + 1
      do ii = 1, 5
        ml(ii) = m(ii)
      end do
c
      ao = m(1)
      wo = m(2)
      so2 = abs(m(3))
      tp = m(4)
      tg = m(5)
      so = sqrt(so2)
c
c     compute differential data vector and the partial derivatives
c
c     print *, 'begin'
      do ii = 1, nfil
        t = tt(ii)
        ek = exp(-0.5*so2*(t - tg)**2)
        ck = cos(wo*(t - tp))
        sk = sin(wo*(t - tp))
        wnd = wt(ii)
c
c       analytic cross-correlation
c
        cest(ii) = ao*ek*ck
c
c       differential data vector
c
        d(ii) = (c(npad + ii) - cest(ii))*wnd
c
c       partials - dc/dao
c
        a(ii,1) = ek*ck*wnd*wao
c
c       partials - dc/dwo
c
        a(ii,2) = -ao*ek*sk*(t - tp)*wnd*wwo
c
c       partials - dc/dso2
c
        a(ii,3) = -0.5*ao*ek*ck*((t - tg)**2)*wnd*wso
c
c       partials - dc/dtp
c
        a(ii,4) = ao*ek*sk*wo*wnd*wtp
c
c       partials - dc/dtg
c
        a(ii,5) = ao*ek*ck*so2*(t - tg)*wnd*wtg

c        print *, a(ii,1),a(ii,2),a(ii,3),a(ii,4),a(ii,5)
c
        end do
c       print *, 'end'
    
c
c     calculate chisqr fit and test for convergence
c
      chi1 = 0.0
      chi2 = 0.0
      do ii = 1, nfil
        chi1 = chi1 + d(ii)**2
        chi2 = chi2 + (cest(ii)*wt(ii))**2
      end do
      chisq(iter) = chi1/chi2
c
c      write(6,'(a,i3,5f10.5,f10.5)') 
c     & '  estimates: #, ao, wo, so, tp, tg, chisq ',
c     & iter, ao, wo*mhz, so*mhz, tp, tg, chisq(iter)
c
c     check on the number of iterations
c
      if (iter .gt. itmax) then
        print*,' '
        print*,' maximum number of iterations exceeded'
        go to 200
        print*,' continue iterating?'
        read(*,'(a)') ans
        if (ans .eq. 'y') then
          iter = 0
          go to 100
        else
          go to 200
        endif
      elseif (iter .ne. 1) then
        dch = chisq(iter) - chisq(iter-1)
c
c       absolute value
c
        if (chisq(iter) .lt. tol1) then
          print*, ' good enough'
          go to 200
c
c       convergence
c
        elseif (abs(dch) .lt. tol2) then
          print*,' inversion converged'
          go to 200
c
c       all hell breaking loose
c
        elseif (chisq(iter) .gt. tol3) then
          print*,' inversion blowing up'
c
c         un-update estimates and do not invert again
c
          do j = 1, ninv
            m(j) = m(j) - x(j)
          end do
          ierr = -1
          go to 300
        endif
      end if
c
c     am = d
c
      ninv = ncol
c
c     form svd of the partial derivative matrix 
c
      call svd1(npmax,ncol,nfil,ninv,a,u,v,w,index)
c
c     preform backsubstitution
c
      call inv(npmax,ncol,nfil,ninv,u,v,w,d,x,alpha,cov)

c
c     update estimates
c
        x(1)=x(1)*wao
        x(2)=x(2)*wwo
        x(3)=x(3)*wso
        x(4)=x(4)*wtp
        x(5)=x(5)*wtg
      do j = 1, ninv
        m(j) = m(j) + x(j)
      end do
c
c     check for a bad problem
c
      if (m(1) .lt. 0.0) then
        print*,' negative amplitude! - help'
        do ii = 1, 5
          m(ii) = ml(ii)
        end do
      end if
c
c     loop back for another test
c
      go to 100
c
c     passed convergence test
c
200   continue
      print*,' '
300   continue
c
      ao = m(1)
      wo = m(2)
      so2 = m(3)
      tp = m(4)
      tg = m(5)
      so = sqrt(so2)
      chi = chisq(iter)
c
c     estimate errors on these parameters
c
      chi1 = 0.0
      do ii = 1, nfil
        chi1 = chi1 + (d(ii) - cest(ii)*wt(ii))**2
      end do
      do ii = 1, ninv
        sigs(ii) = sqrt(cov(ii,ii))*sqrt(chi1/(nfil-2))
      end do
c
c     do ii = 1, 5
c       write(6,'(3f15.7)') m(ii), sigs(ii), sqrt(cov(ii,ii))
c     end do
c
c     all done
c
      return
      end
