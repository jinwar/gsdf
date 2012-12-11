c
c
c
      subroutine coolb(nn,datai,signi)
c
c     fast Fourier transform of complex input datai
c     (used as real array to speed things)
c     from Caltech initially; to MIT by C Jones, Sept 1986
c
      dimension datai(*)
      n=2**(nn+1)
      j=1
      do 5 i=1,n,2
        if(i-j)1,2,2
    1   tempr=datai(j)
        tempi=datai(j+1)
        datai(j)=datai(i)
        datai(j+1)=datai(i+1)
        datai(i)=tempr
        datai(i+1)=tempi
    2   m=n/2
    3   if(j-m)5,5,4
    4   j=j-m
        m=m/2
        if(m-2)5,3,3
    5   j=j+m
      mmax=2
    6 if(mmax-n)7,10,10
    7 istep=2*mmax
      theta=signi*6.28318531/float(mmax)
      sinth=sin(theta/2.)
      wstpr=-2.0  *sinth*sinth
      wstpi= sin(theta)
      wr=1.
      wi=0.
      do 9 m=1,mmax,2
        do 8 i=m,n,istep
          j=i+mmax
          tempr=wr*datai(j)-wi*datai(j+1)
          tempi=wr*datai(j+1)+wi*datai(j)
          datai(j)=datai(i)-tempr
          datai(j+1)=datai(i+1)-tempi
          datai(i)=datai(i)+tempr
    8     datai(i+1)=datai(i+1)+tempi
        tempr=wr
        wr=wr*wstpr-wi*wstpi+wr
    9   wi=wi*wstpr+tempr*wstpi+wi
      mmax=istep
      go to 6
   10 return
      end
