c
c
c
      subroutine envelop(nnww,jnww,jenvl)
c
c     subroutine to create an envelope of trace in jnww(1:nnww)
c     stores envelope in jenvl
c     stolen from numerous sstudy programs of J.B. Minster
c     to MIT via Univ. of Utah (J. Pechmann)
c     modified slightly by C Jones, Sept 1986
c
      include 'parameter.f'
c
      real jnww(*),jenvl(*)
      complex carray(npmax)
c
c     find power of 2
c
      nl=1
      itmp=2
10    nl=nl+1
      itmp=itmp*2
      if(itmp.gt.nnww)go to 20
      go to 10
20    np=itmp
      if(np .gt. npmax) then
        print*,' parameter npmax exceeded'
        return
      endif
25    np2=itmp/2+1
      do 30 i=1,nnww
30      carray(i)=jnww(i)
      nwp1=nnww+1
      do 31 i=nwp1,np
31       carray(i)=(0.,0.)
      call coolb(nl,carray,-1.)
      carray(1)=carray(1)*0.5
      do 40 i=np2,np
40      carray(i)=(0.,0.)
      call coolb(nl,carray,1.)
      xx=2./float(np)
      do 50 i=1,np
        x=cabs(carray(i))*xx
50      jenvl(i)=x
      return
      end
