

      subroutine kessler(dt,tauto,taccr,tevar,ruh,rvh,rmh,pi0,th0,tmp,   &
                         rho,rr,pp3d,th3d,prs,                        &
                         qv3d,qc3d,qr3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real :: dt
      real*8 :: tauto,taccr,tevar
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh,pi0,th0
      real, dimension(ib:ie,jb:je,kb:ke) :: tmp,rho,rr,pp3d,    &
                                            th3d,prs,qv3d,qc3d,qr3d

      integer :: i,j,k
      real :: qvnew,qcnew,qrnew
      real :: ar,cr,qvs,er,term1,cpml,cvml,rm,tem
      real*8 dum
      real*8, dimension(nk) :: bud1,bud2,bud3

      real rslf

!-------------------------------------------------------------------

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        tmp(i,j,k)=(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
      enddo
      enddo
      enddo

!$omp parallel do default(shared)  &
!$omp private(k)
      do k=1,nk
        bud1(k)=0.0d0
        bud2(k)=0.0d0
        bud3(k)=0.0d0
      enddo

    IF(isnd.ne.4)THEN

!$omp parallel do default(shared)              &
!$omp private(i,j,k,qvnew,qcnew,qrnew,ar,cr,qvs,er,term1,cpml,cvml,rm,tem)
      do k=1,nk
      do j=1,nj
      do i=1,ni

        qvnew=qv3d(i,j,k)
        qcnew=qc3d(i,j,k)
        qrnew=qr3d(i,j,k)

        ar=0.0
        cr=0.0
        er=0.0

        !! autoconversion of cloud to rain
        if(qcnew.gt.0.0)then
          ar=max(0.001*(qcnew-0.001),0.0)
          ar=ar*dt
        endif

        !! accretion of cloud by rain
        if(qcnew.gt.0.0 .and. qrnew.gt.0.0)then
          cr=2.2*qcnew*(qrnew**0.875)
          cr=cr*dt
        endif

        !! evap of rain to vapor
        if(qrnew.gt.0.0)then
          qvs=rslf(prs(i,j,k),tmp(i,j,k))
          if(qvnew.lt.qvs)then
            er=(1.6+30.3922*((rho(i,j,k)*qrnew)**0.2046))*          &
                (1.0-(qvnew/qvs))*                                  &
                ((rho(i,j,k)*qrnew)**0.525)/                        &
               ( (2.03e4+9.584e6/(qvs*prs(i,j,k))) * rho(i,j,k) )
            er=min(er*dt,qrnew)
            if( (qvnew+er).gt.qvs )then
              er=qvs-qvnew
            endif
          endif
        endif

        if((ar+cr).gt.qcnew)then
          term1=ar+cr
          ar=qcnew*ar/term1
          cr=qcnew*cr/term1
        endif

        qvnew=qvnew+er
        qcnew=qcnew-(ar+cr)
        qrnew=qrnew+(ar+cr-er)

      if(er.gt.1.0e-7)then
        if(neweqts.ge.1)then
          cpml=cp+cpv*qvnew+cpl*(qcnew+qrnew)
          cvml=cv+cvv*qvnew+cpl*(qcnew+qrnew)
          rm=rd+rv*qvnew
          th3d(i,j,k)=th3d(i,j,k)-er*(                                     &
                  (lv1-lv2*tmp(i,j,k))/(cpdcv*cvml*(pi0(i,j,k)+pp3d(i,j,k)))   &
                - (th0(i,j,k)+th3d(i,j,k))*(rv/cvml)*(1.0-rovcp*cpml/rm) )
          pp3d(i,j,k)=((rho(i,j,k)*rm*(th0(i,j,k)+th3d(i,j,k))*rp00)**rddcv)-pi0(i,j,k)
          prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
        else
          th3d(i,j,k)=th3d(i,j,k)-er*( (lv1-lv2*tmp(i,j,k))         &
                                      /(cp*(pi0(i,j,k)+pp3d(i,j,k))) )
          rho(i,j,k)=prs(i,j,k)   &
             /(rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))*(1.0+qvnew*reps))
        endif
      endif

        tem=ruh(i)*rvh(j)*rmh(i,j,k)

        bud1(k)=bud1(k)+rr(i,j,k)*ar*tem
        bud2(k)=bud2(k)+rr(i,j,k)*cr*tem
        bud3(k)=bud3(k)+rr(i,j,k)*er*tem

        qv3d(i,j,k)=qvnew
        qc3d(i,j,k)=qcnew
        qr3d(i,j,k)=qrnew

      enddo
      enddo
      enddo

    ENDIF

      dum=dx*dy*dz

      do k=1,nk
        tauto=tauto+bud1(k)*dum
      enddo

      do k=1,nk
        taccr=taccr+bud2(k)*dum
      enddo

      do k=1,nk
        tevar=tevar+bud3(k)*dum
      enddo

      if(timestats.ge.1) time_microphy=time_microphy+mytime()
 
      RETURN
      END


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine satadj(nrk,tcond,tevac,ruh,rvh,rmh,pi0,th0,   &
                        rho,rr,pp3d,prs,th3d,q3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      integer nrk
      real*8 :: tcond,tevac
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh,pi0,th0
      real, dimension(ib:ie,jb:je,kb:ke) :: rho,rr,pp3d,prs,th3d
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d

      integer :: i,j,k,n,nmax,iflag
      real :: tnew,qvs,qvnew,qcnew,thn,pin,cpml,cvml,rm,lhv,thlast,dqv
      real :: converge,t1,p1,d1,tem,ql,qi
      real*8 dum
      real*8, dimension(nk) :: bud1,bud2
      real rslf
      logical :: doit

!--------------------------------------------------------------------
!  iterative sat adj.

    nmax=0
    iflag=0

    IF(nrk.ge.3)THEN
!$omp parallel do default(shared)  &
!$omp private(k)
      do k=1,nk
        bud1(k)=0.0d0
        bud2(k)=0.0d0
      enddo
    ENDIF

    IF(neweqts.ge.1)THEN

      if(nrk.eq.4)then
!!!        converge=0.0005
        converge=5.0*tsmall
      else
!!!        converge=0.01
        converge=100.0*tsmall
      endif

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,tnew,qvs,qvnew,qcnew,thn,pin,cpml,cvml,rm,lhv,   &
!$omp thlast,dqv,t1,p1,d1,tem,ql,qi,doit)
      do k=1,nk
      do j=1,nj
      do i=1,ni

        tnew=(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
        qvs=rslf(prs(i,j,k),tnew)
              
        IF(q3d(i,j,k,2).gt.0.0 .or. q3d(i,j,k,nqv).gt.qvs)THEN

          qvnew=q3d(i,j,k,nqv)
          qcnew=q3d(i,j,k,2)
          thn=th3d(i,j,k)
          pin=pp3d(i,j,k)
          ql=0.0
          qi=0.0
          do n=nql1,nql2
            ql=ql+q3d(i,j,k,n)
          enddo
          if(iice.eq.1)then
            do n=nqs1,nqs2
              qi=qi+q3d(i,j,k,n)
            enddo
          endif
          cpml=cp+cpv*qvnew+cpl*ql+cpi*qi
          cvml=cv+cvv*qvnew+cpl*ql+cpi*qi
          rm=rd+rv*q3d(i,j,k,nqv)
          lhv=lv1-lv2*tnew

          t1=lhv/(cpdcv*cvml*(pi0(i,j,k)+pp3d(i,j,k)))   &
            -(th0(i,j,k)+th3d(i,j,k))*(rv/cvml)*(1.0-rovcp*cpml/rm)
          p1=rovcp*( lhv/(cvml*(th0(i,j,k)+th3d(i,j,k)))   &
            -(pi0(i,j,k)+pp3d(i,j,k))*rv*cpml/(rm*cvml) )
!!!          d1=t1*(pi0(i,j,k)+pp3d(i,j,k))*4097.8531
          d1=t1*(pi0(i,j,k)+pp3d(i,j,k))*17.67*243.5

          n=0
          thlast=thn
          doit=.true.

          do while( doit )
            n=n+1
!!!            dqv=(qvs-qvnew)/(1.0+d1*qvs/((tnew-35.86)**2) )
            dqv=(qvs-qvnew)/(1.0+d1*qvs/((tnew-29.65)**2) )
            dqv=min(dqv,qcnew)
            if(  (qvnew+dqv).lt.1.0e-20 ) dqv=1.0e-20-qvnew

            qvnew=qvnew+dqv
            qcnew=qcnew-dqv
            thn=thn-dqv*t1
            pin=pin-dqv*p1
            prs(i,j,k)=p00*((pi0(i,j,k)+pin)**cpdrd)

            doit = .false.
            if( abs(thn-thlast).gt.converge )then
              thlast=thn
              tnew=(th0(i,j,k)+thn)*(pi0(i,j,k)+pin)
              qvs=rslf(prs(i,j,k),tnew)
              doit = .true.
            endif

            if(n.gt.50.and.n.lt.100)then
              print *,n,thn,pin
            elseif(n.eq.100)then
              print *,'  infinite loop!'
              print *,'  i,j,k=',i,j,k
              iflag=1
              doit=.false.
            endif

          enddo

          tem=ruh(i)*rvh(j)*rmh(i,j,k)

          bud1(k)=bud1(k)+rr(i,j,k)*max(qcnew-q3d(i,j,k,2),0.0)*tem
          bud2(k)=bud2(k)+rr(i,j,k)*max(qvnew-q3d(i,j,k,nqv),0.0)*tem

          th3d(i,j,k)=thn
          q3d(i,j,k,2)=qcnew
          q3d(i,j,k,nqv)=qvnew
          pp3d(i,j,k)=pin

          nmax=max(n,nmax)

        ENDIF

      enddo
      enddo
      enddo

    ELSE

      nmax=1

!$omp parallel do default(shared)  &
!$omp private(i,j,k,qvnew,qcnew,tnew,qvs,lhv,dqv,tem)
      do k=1,nk
      do j=1,nj
      do i=1,ni

        qvnew=q3d(i,j,k,nqv)
        qcnew=q3d(i,j,k,2)
        tnew=(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
        qvs=rslf(prs(i,j,k),tnew)
        lhv=lv1-lv2*tnew
!!!        dqv=(qvs-qvnew)/(1.0+lhv*qvs*4097.8531/(cp*((tnew-35.86)**2)))
        dqv=(qvs-qvnew)/(1.0+lhv*qvs*17.67*243.5/(cp*((tnew-29.65)**2)))
        if(  (qvnew+dqv).lt.1.0e-20 ) dqv=1.0e-20-qvnew
        dqv=min(dqv,max(0.0,qcnew))

        qvnew=qvnew+dqv
        qcnew=qcnew-dqv
        tem=ruh(i)*rvh(j)*rmh(i,j,k)
        bud1(k)=bud1(k)+rr(i,j,k)*max(qcnew-q3d(i,j,k,2),0.0)*tem
        bud2(k)=bud2(k)+rr(i,j,k)*max(qvnew-q3d(i,j,k,nqv),0.0)*tem

        th3d(i,j,k)=th3d(i,j,k)-dqv*( lhv/(cp*(pi0(i,j,k)+pp3d(i,j,k))) )
        q3d(i,j,k,2)=qcnew
        q3d(i,j,k,nqv)=qvnew
        rho(i,j,k)=prs(i,j,k)   &
             /(rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))*(1.0+qvnew*reps))

      enddo
      enddo
      enddo

    ENDIF

    IF(nrk.ge.3)THEN
      dum=dx*dy*dz
      do k=1,nk
        tcond=tcond+bud1(k)*dum
      enddo

      do k=1,nk
        tevac=tevac+bud2(k)*dum
      enddo
    ENDIF

!!!      print *,'  nmax=',nmax

      if(iflag.ne.0)then
        print *
        print *,' Convergence cannot be reached in satadj subroutine.'
        print *
        print *,' This may be a problem with the algorithm in satadj.'
        print *,' However, the model may have became unstable somewhere'
        print *,' else and the symptoms first appeared here.'
        print *
        print *,' Try decreasing the timestep (dtl and/or nsound).'
        print *
        print *,'  ... stopping cm1 ... '
        print *
        call stopcm1
      endif

      if(timestats.ge.1) time_satadj=time_satadj+mytime()

      RETURN
      END


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine k_fallout(rho,qr3d,vr)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, dimension(ib:ie,jb:je,kb:ke) :: rho,qr3d,vr
 
      integer i,j,k

!--------------------------------------------------------------------
!  Get fall velocities

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        vr(i,j,k)=0.0
      enddo
      enddo
      enddo

      IF(isnd.ne.4)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          if(qr3d(i,j,k).gt.1.0e-12)then
            vr(i,j,k)=14.34*((rho(i,j,k)*qr3d(i,j,k))**0.1346)    &
                           *sqrt(1.15/rho(i,j,k))
          endif
        enddo
        enddo
        enddo
      ENDIF

      if(timestats.ge.1) time_fall=time_fall+mytime()
 
      RETURN
      END


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine fallout(dt,train,ruh,rvh,zh,mh,mf,rain,rr,rho,   &
                         q3d,vq)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real :: dt
      real*8 :: train
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,nrain) :: rain
      real, dimension(ib:ie,jb:je,kb:ke) :: rr,rho,q3d,vq

      integer :: i,j,k,n,nr,nrk
      real :: qmax,f1,f2,f3,b1,b2,b3,w1,w2,w3,dttmp
      real, dimension(-1:nk+3):: rq
      real, dimension(0:nk) :: ffk
      real, dimension(nk) :: qtmp

      integer, dimension(ni,nj) :: nfall
      real :: crmax,dtfall

      real epsilon,onedsix,thdtw
      parameter(epsilon=1.0e-12)
      parameter(onedsix=1.0/6.0)
      parameter(thdtw=13.0/12.0)
      real*8 tem
      real*8 bud(nj)
      real*8 t1,t2,t3

!--------------------------------------------------------------------

!$omp parallel do default(shared)  &
!$omp private(j)
      do j=1,nj
        bud(j)=0.0d0
      enddo

!$omp parallel do default(shared)            &
!$omp private(i,j,k,crmax)
      do j=1,nj
      do i=1,ni
        crmax = 0.0
        do k=1,nk
          crmax = max( crmax , vq(i,j,k)*dt*rdz*mh(i,j,k) )
        enddo
        nfall(i,j) = max( 1 , int(crmax+1.0) )
      enddo
      enddo

!--------------------------------------------------------------------

  IF(ifall.eq.1)THEN

!$omp parallel do default(shared)            &
!$omp private(i,j,k,rq,n,dtfall,nr)
    do j=1,nj
    do i=1,ni
      dtfall=dt/nfall(i,j)
      do n=1,nfall(i,j)
        do k=1,nk
          rq(k)=rho(i,j,k)*vq(i,j,k)*max(0.0,q3d(i,j,k))
        enddo
        rq(nk+1)=0.0
        do k=1,nk
          q3d(i,j,k)=q3d(i,j,k)+dtfall*(rq(k+1)-rq(k))   &
                                     *rdz*mf(i,j,k+1)/rho(i,j,k)
        enddo
        do nr=1,nrain
          rain(i,j,nr)=rain(i,j,nr)+0.1*dtfall*rq(1)
        enddo
        bud(j)=bud(j)+dtfall*rr(i,j,1)*vq(i,j,1)*max(0.0,q3d(i,j,1))*ruh(i)*rvh(j)
      enddo
    enddo
    enddo

!--------------------------------------------------------------------

  ELSEIF(ifall.eq.2)THEN

!$omp parallel do default(shared)            &
!$omp private(i,j,k,rq,n,dtfall,nr)
    do j=1,nj
    do i=1,ni
      dtfall=dt/nfall(i,j)
      do n=1,nfall(i,j)
        do k=1,nk
          rq(k)=rho(i,j,k)*vq(i,j,k)*q3d(i,j,k)
        enddo
        rq(nk+1)=0.0
        do k=2,nk
          q3d(i,j,k)=q3d(i,j,k)+dtfall*(rq(k+1)-rq(k-1))   &
                                     /(zh(i,j,k+1)-zh(i,j,k-1))/rho(i,j,k)
        enddo
        q3d(i,j,1)=q3d(i,j,1)+dtfall*(rq(2)-rq(1))*rdz*mf(i,j,2)/rho(i,j,1)
        do k=1,nk
          q3d(i,j,k)=max(q3d(i,j,k),0.0)
        enddo
        do nr=1,nrain
          rain(i,j,nr)=rain(i,j,nr)+0.1*dtfall*rq(1)
        enddo
        bud(j)=bud(j)+dtfall*rr(i,j,1)*vq(i,j,1)*max(0.0,q3d(i,j,1))*ruh(i)*rvh(j)
      enddo
    enddo
    enddo

!--------------------------------------------------------------------

  ELSEIF(ifall.eq.3)THEN

!$omp parallel do default(shared)            &
!$omp private(i,j,k,nrk,qmax,dttmp,rq,ffk,qtmp,    &
!$omp f1,f2,f3,b1,b2,b3,w1,w2,w3,t1,t2,t3,n,dtfall,nr)
    do j=1,nj
    do i=1,ni
      qmax=0.
      do k=1,nk
        qmax=max(qmax,q3d(i,j,k))
      enddo
      IF(qmax.gt.1.0e-8)THEN
        dtfall=dt/nfall(i,j)
        do n=1,nfall(i,j)
          do k=1,nk
            qtmp(k)=q3d(i,j,k)
          enddo
          do nrk=1,3
            do k=1,nk
              rq(k)=rho(i,j,k)*vq(i,j,k)*max(0.0,qtmp(k))
            enddo
            rq(-1)=rq(1)
            rq( 0)=rq(1)
            rq(nk+1)=0.0
            rq(nk+2)=0.0
            rq(nk+3)=0.0
            do k=0,nk
              f1=(  2*rq(k+3)-7*rq(k+2)+11*rq(k+1) )*onedsix
              f2=( -rq(k+2)+5*rq(k+1)+2*rq(k  ) )*onedsix
              f3=(  2*rq(k+1)+5*rq(k  )-rq(k-1) )*onedsix

              b1=thdtw*(rq(k+3)-2*rq(k+2)+rq(k+1))**2    &
                +0.25*(rq(k+3)-4*rq(k+2)+3*rq(k+1))**2
              b2=thdtw*(rq(k+2)-2*rq(k+1)+rq(k  ))**2    &
                +0.25*(rq(k+2)-rq(k  ))**2
              b3=thdtw*(rq(k+1)-2*rq(k  )+rq(k-1))**2    &
                +0.25*(3*rq(k+1)-4*rq(k  )+rq(k-1))**2

              w1=0.10/((epsilon+b1)**2)
              w2=0.60/((epsilon+b2)**2)
              w3=0.30/((epsilon+b3)**2)

              ffk(k)=((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
            enddo
            dttmp=dtfall/(4-nrk)
            do k=1,nk
              qtmp(k)=q3d(i,j,k)+dttmp*(ffk(k)-ffk(k-1))   &
                                    *rdz*mh(i,j,k)/rho(i,j,k)
            enddo
          enddo
          t1=0.0d0
          t2=0.0d0
          do k=1,nk
            t1=t1+rho(i,j,k)*qtmp(k)
            q3d(i,j,k)=max(0.0,qtmp(k))
            t2=t2+rho(i,j,k)*q3d(i,j,k)
          enddo
          t3=(t1+1.0d-20)/(t2+1.0d-20)
          if(t3.lt.0.0) t3=1.0d0
          do k=1,nk
            q3d(i,j,k)=t3*q3d(i,j,k)
          enddo
          do nr=1,nrain
            rain(i,j,nr)=rain(i,j,nr)+0.1*dtfall*max(rq(1),0.0)
          enddo
          bud(j)=bud(j)+dtfall*rr(i,j,1)*vq(i,j,1)*max(0.0,q3d(i,j,1))*ruh(i)*rvh(j)
        enddo
      ENDIF
    enddo
    enddo

  ENDIF

!--------------------------------------------------------------------

      tem=dx*dy
      do j=1,nj
        train=train+bud(j)*tem
      enddo

      if(timestats.ge.1) time_fall=time_fall+mytime()

      RETURN
      END


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getefall(dt,cpx,ruh,rvh,mf,pi0,th0,t,cvm,rr,   &
                          pp3d,th3d,q3d,vr)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real :: dt,cpx
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,th0
      real, dimension(ib:ie,jb:je,kb:ke) :: t,cvm,rr
      real, dimension(ib:ie,jb:je,kb:ke) :: pp3d,th3d,q3d,vr

      integer :: i,j,k
      real :: ften

!$omp parallel do default(shared)  &
!$omp private(i,j,k,ften)
      do k=1,nk-1
      do j=1,nj
      do i=1,ni
        ften=q3d(i,j,k)*vr(i,j,k)*(                                  &
                     cpx*(t(i,j,k+1)-t(i,j,k))*rdz*mf(i,j,k+1) + g   &
                                  )/cvm(i,j,k)
        th3d(i,j,k)=th3d(i,j,k)+dt*ften/( cpdcv*(pi0(i,j,k)+pp3d(i,j,k)) )
        pp3d(i,j,k)=pp3d(i,j,k)+dt*ften*rovcp/(th0(i,j,k)+th3d(i,j,k))
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_fall=time_fall+mytime()

      RETURN
      END


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine geterain(dt,cpx,lx1,erain,ruh,rvh,t,rr,q3d,vr)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real :: dt,cpx,lx1
      real*8 :: erain
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: t,rr
      real, dimension(ib:ie,jb:je,kb:ke) :: q3d,vr

      integer :: i,j
      real*8 :: bud(nj)

!$omp parallel do default(shared)  &
!$omp private(j)
      do j=1,nj
        bud(j)=0.0d0
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        bud(j)=bud(j)+dt*vr(i,j,1)*rr(i,j,1)*q3d(i,j,1)*(cpx*t(i,j,1)-lx1)*ruh(i)*rvh(j)
      enddo
      enddo

      do j=1,nj
        erain=erain+bud(j)*dx*dy
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_fall=time_fall+mytime()

      RETURN
      END


