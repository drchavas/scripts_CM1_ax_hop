

      subroutine sound(dt,xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,zh,mh,rmh,mf,rmf,  &
                        pi0,thv0,rho0,rr0,rf0,th0,dzdx,dzdy,            &
                        radbcw,radbce,radbcs,radbcn,                    &
                        fwk,ppd,dpdzx,dpdzy,pk,mm,qk,fpk,               &
                        gx,u0,ua,u3d,uten,gy,v0,va,v3d,vten,wa,w3d,wten,      &
                        ppi,pp3d,ppten,tha,th3d,thten,thterm,tk,        &
                        thv,ppterm,ftk,bk,sk,nrk,rtime,                 &
                        reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,             &
                        uw31,uw32,ue31,ue32,us31,us32,un31,un32,        &
                        vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,        &
                        ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,        &
                        sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,        &
                        pw31,pw32,pe31,pe32,ps31,ps32,pn31,pn32,        &
                        pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,rxh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,rho0,rr0,rf0,th0
      real, intent(in), dimension(itb:ite,jtb:jte) :: dzdx,dzdy
      real, intent(in), dimension(itb:ite+1,jtb:jte,ktb:kte) :: gx
      real, intent(in), dimension(itb:ite,jtb:jte+1,ktb:kte) :: gy
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: fwk,ppd,dpdzx,dpdzy,pk,mm,qk,fpk
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua,u3d,uten
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va,v3d,vten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa,w3d,wten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tha,th3d,thten,thterm,tk
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: thv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppterm,ftk,bk,sk
      integer, intent(in) :: nrk
      real, intent(in)  :: rtime
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p
      real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, intent(inout), dimension(imp,cmp,kmp)   :: vs31,vs32,vn31,vn32
      real, intent(inout), dimension(cmp,jmp,kmp-1) :: ww31,ww32,we31,we32
      real, intent(inout), dimension(imp,cmp,kmp-1) :: ws31,ws32,wn31,wn32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: pw31,pw32,pe31,pe32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ps31,ps32,pn31,pn32
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2

!-----

      integer :: i,j,k,n,nloop
      real :: tem,r1,dts
      real :: aa,bb,cc,dd

!!!      real alph,beta,asq
!!!      parameter(alph=0.60)
!!!      parameter(beta=1.0-alph)
!!!      parameter(asq=alph*alph)

      real :: beta,asq

#ifdef MPI
      include 'mpif.h'
#endif

!---------------------------------------------------------------------

      beta=1.0-alph
      asq=alph*alph

!---------------------------------------------------------------------
!  Prepare for acoustic steps

      if(nrk.ge.2)then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          u3d(i,j,k)=ua(i,j,k)
        enddo
        enddo
        enddo
 
      IF(axisymm.eq.0)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          v3d(i,j,k)=va(i,j,k)
        enddo
        enddo
        enddo
      ENDIF
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          w3d(i,j,k)=wa(i,j,k)
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          pp3d(i,j,k)=ppi(i,j,k)
          ppd(i,j,k)=pp3d(i,j,k)
        enddo
        enddo
        enddo

      else

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          ppd(i,j,k)=pp3d(i,j,k)
        enddo
        enddo
        enddo

      endif

        IF(thsmall.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            th3d(i,j,k)=tha(i,j,k)
          enddo
          enddo
          enddo
        ENDIF

!---------------------------------------------------------------------
 
      if(nrk.eq.1)then
!!!        nloop=1
!!!        dts=dt/3.
        nloop=nint(float(nsound)/3.0)
        dts=dt/(nloop*3.0)
        if( dts.gt.(dt/nsound) )then
          nloop=nloop+1
          dts=dt/(nloop*3.0)
        endif
      elseif(nrk.eq.2)then
        nloop=0.5*nsound
        dts=dt/nsound
      elseif(nrk.eq.3)then
        nloop=nsound
        dts=dt/nsound
      endif

!!!      print *,'  nloop,dts,dttmp = ',nloop,dts,nloop*dts

!-----------------------------------------------------------------

      tem = dts*g*0.5/cp
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        pk(i,j,k)=tem/thv0(i,j,k)
      enddo
      enddo
      enddo

      tem = dts*cp*0.5*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        mm(i,j,k)=tem*mf(i,j,k)*( thv(i,j,k-1)+thv(i,j,k) )
      enddo
      enddo
      enddo

      IF(thsmall.eq.1)THEN
        tem = dts*g
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          bk(i,j,k)=tem/(th0(i,j,k-1)+th0(i,j,k))
        enddo
        enddo
        enddo
      ENDIF

        IF(thsmall.eq.1.and.neweqts.ge.1.and.imoist.eq.1)THEN
          tem = dts*rdz
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            sk(i,j,k)=tem*thterm(i,j,k)*(th0(i,j,k)+th3d(i,j,k))*mh(i,j,k)
          enddo
          enddo
          enddo
        ENDIF

      if(timestats.ge.1) time_sound=time_sound+mytime()

!---------------------------------------------------------------------

      DO N=1,NLOOP

!-----

        if(irbc.eq.2)then
 
          if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,u3d)
 
          if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,v3d)
 
        endif

!-----

        if(wbc.eq.2.and.ibw.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            u3d(1,j,k)=u3d(1,j,k)+dts*( -radbcw(j,k)          &
                      *(u3d(2,j,k)-u3d(1,j,k))*rdx*uh(1)   &
                         +uten(1,j,k) )
          enddo
          enddo
        endif

        if(ebc.eq.2.and.ibe.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            u3d(ni+1,j,k)=u3d(ni+1,j,k)+dts*( -radbce(j,k)              &
                         *(u3d(ni+1,j,k)-u3d(ni  ,j,k))*rdx*uh(ni)   &
                         +uten(ni+1,j,k) )
          enddo
          enddo
        endif

        if(timestats.ge.1) time_sound=time_sound+mytime()

!-----

        IF(roflux.eq.1)THEN
          call restrict_openbc_we(rvh,rmh,rho0,u3d)
        ENDIF

!-----
#ifdef MPI
        IF(n.gt.1)THEN
          call comm_1s_end_we(ppd,pw1,pw2,pe1,pe2,   &
                               ps1,ps2,pn1,pn2,reqs_p)
        ENDIF
#endif
!-----

        tem=rdx*cp*0.5

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1+ibw,ni+1-ibe
          u3d(i,j,k)=u3d(i,j,k)+dts*( uten(i,j,k)                &
                  -tem*(ppd(i,j,k)-ppd(i-1,j,k))*uf(i)*   &
                       (thv(i,j,k)+thv(i-1,j,k)) )
        enddo
        enddo
        enddo

        IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            fwk(i,j,k)=0.5*(ppd(i-1,j,k)+ppd(i,j,k))
          enddo
          enddo
          enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=2,nk-1
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            dpdzx(i,j,k)=( fwk(i,j,k+1)-fwk(i,j,k-1) )*rdz2
          enddo
          enddo
          enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            dpdzx(i,j,1 )=( fwk(i,j,2 )-fwk(i,j,1   ) )*rdz
            dpdzx(i,j,nk)=( fwk(i,j,nk)-fwk(i,j,nk-1) )*rdz
          enddo
          enddo

          tem=dts*cp*0.5

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)                                  &
                  -tem*(thv(i,j,k)+thv(i-1,j,k))*gx(i,j,k)*dpdzx(i,j,k)
          enddo
          enddo
          enddo

        ENDIF

!----------------------------------------------
!  convergence forcing:

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. nx.gt.1 )THEN
            call convinitu(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibw,ibe,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xf,yh,zh,u0,u3d)
          ENDIF
        ENDIF

!----------------------------------------------

        if(timestats.ge.1) time_sound=time_sound+mytime()

        if(n.eq.nloop)then
          call bcu(u3d)
        endif

#ifdef MPI
        if(n.eq.nloop)then
          call comm_3u_start(u3d,uw31,uw32,ue31,ue32,   &
                                 us31,us32,un31,un32,reqs_u)
        endif
#endif

!-----
#ifdef MPI
        IF(n.gt.1)THEN
          call comm_1s_end_sn(ppd,pw1,pw2,pe1,pe2,   &
                               ps1,ps2,pn1,pn2,reqs_p)
        ENDIF
#endif
!-----

      IF(axisymm.eq.0)THEN

        if(sbc.eq.2.and.ibs.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            v3d(i,1,k)=v3d(i,1,k)+dts*( -radbcs(i,k)          &
                      *(v3d(i,2,k)-v3d(i,1,k))*rdy*vh(1)   &
                      +vten(i,1,k) )
          enddo
          enddo
        endif
 
        if(nbc.eq.2.and.ibn.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            v3d(i,nj+1,k)=v3d(i,nj+1,k)+dts*( -radbcn(i,k)              &
                         *(v3d(i,nj+1,k)-v3d(i,nj  ,k))*rdy*vh(nj)   &
                         +vten(i,nj+1,k) )
          enddo
          enddo
        endif

        if(timestats.ge.1) time_sound=time_sound+mytime()

!-----

        IF(roflux.eq.1)THEN
          call restrict_openbc_sn(ruh,rmh,rho0,v3d)
        ENDIF

!-----

        tem=rdy*cp*0.5

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1+ibs,nj+1-ibn
        do i=1,ni
          v3d(i,j,k)=v3d(i,j,k)+dts*( vten(i,j,k)                &
                  -tem*(ppd(i,j,k)-ppd(i,j-1,k))*vf(j)*   &
                       (thv(i,j,k)+thv(i,j-1,k)) )
        enddo
        enddo
        enddo

        IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            fwk(i,j,k)=0.5*(ppd(i,j-1,k)+ppd(i,j,k))
          enddo
          enddo
          enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=2,nk-1
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            dpdzy(i,j,k)=( fwk(i,j,k+1)-fwk(i,j,k-1) )*rdz2
          enddo
          enddo
          enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            dpdzy(i,j,1 )=( fwk(i,j,2   )-fwk(i,j,1   ) )*rdz
            dpdzy(i,j,nk)=( fwk(i,j,nk  )-fwk(i,j,nk-1) )*rdz
          enddo
          enddo

          tem=dts*cp*0.5

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            v3d(i,j,k)=v3d(i,j,k)                                  &
                  -tem*(thv(i,j,k)+thv(i,j-1,k))*gy(i,j,k)*dpdzy(i,j,k)
          enddo
          enddo
          enddo

        ENDIF

!----------------------------------------------
!  convergence forcing:

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. ny.gt.1 )THEN
            call convinitv(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibs,ibn,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xh,yf,zh,v0,v3d)
          ENDIF
        ENDIF

!----------------------------------------------

        if(timestats.ge.1) time_sound=time_sound+mytime()

        if(n.eq.nloop)then
          call bcv(v3d)
        endif

#ifdef MPI
        if(n.eq.nloop)then
          call comm_3v_start(v3d,vw31,vw32,ve31,ve32,   &
                                 vs31,vs32,vn31,vn32,reqs_v)
        endif
#endif

      ENDIF

!---------------------------------------------------
!  update implicit solver arrays

        tem = dts*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          qk(i,j,k)=tem*ppterm(i,j,k)*(pi0(i,j,k)+pp3d(i,j,k))*mh(i,j,k)
        enddo
        enddo
        enddo

      IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          fpk(i,j,k)=                 beta*(                                 &
            w3d(i,j,k)*(pk(i,j,k)+qk(i,j,k))                                 &
           +w3d(i,j,k+1)*(pk(i,j,k)-qk(i,j,k)) )                             &
           +dts*( ppten(i,j,k)   &
           -ppterm(i,j,k)*(pi0(i,j,k)+pp3d(i,j,k))                       &
                             *( (u3d(i+1,j,k)-u3d(i,j,k))*rdx*uh(i)  &
                               +(v3d(i,j+1,k)-v3d(i,j,k))*rdy*vh(j) ) )
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          fpk(i,j,k)=                 beta*(                                 &
            w3d(i,j,k)*(pk(i,j,k)+qk(i,j,k))                                 &
           +w3d(i,j,k+1)*(pk(i,j,k)-qk(i,j,k)) )                             &
           +dts*( ppten(i,j,k)   &
           -ppterm(i,j,k)*(pi0(i,j,k)+pp3d(i,j,k))                       &
                             *( (xf(i+1)*u3d(i+1,j,k)-xf(i)*u3d(i,j,k))*rdx*uh(i)*rxh(i) ) )
        enddo
        enddo
        enddo

      ENDIF

        IF(thsmall.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ftk(i,j,k)=dts*( thten(i,j,k)+beta*rr0(i,j,k)*(    &
            w3d(i,j,k)*tk(i,j,k)+w3d(i,j,k+1)*tk(i,j,k+1)  ) )
        enddo
        enddo
        enddo

        IF(neweqts.ge.1 .and. imoist.eq.1)THEN

        IF(axisymm.eq.0)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ftk(i,j,k)=ftk(i,j,k)+beta*sk(i,j,k)*( w3d(i,j,k)-w3d(i,j,k+1) )   &
             -dts*thterm(i,j,k)*(th0(i,j,k)+th3d(i,j,k))                       &
                               *( (u3d(i+1,j,k)-u3d(i,j,k))*rdx*uh(i)  &
                                 +(v3d(i,j+1,k)-v3d(i,j,k))*rdy*vh(j) )
          enddo
          enddo
          enddo
        ELSE
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ftk(i,j,k)=ftk(i,j,k)+beta*sk(i,j,k)*( w3d(i,j,k)-w3d(i,j,k+1) )   &
             -dts*thterm(i,j,k)*(th0(i,j,k)+th3d(i,j,k))                       &
                               *( (xf(i+1)*u3d(i+1,j,k)-xf(i)*u3d(i,j,k))*rdx*uh(i)*rxh(i) )
          enddo
          enddo
          enddo
        ENDIF   ! for axisymm

        ENDIF   ! for neweqts

        ENDIF   ! for thsmall

        IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=2,nk-1
          do j=1,nj
          do i=1,ni+1
            dpdzx(i,j,k)=gx(i,j,k)*(u3d(i,j,k+1)-u3d(i,j,k-1))*rdz2
          enddo
          enddo
          enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni+1
            dpdzx(i,j, 1)=gx(i,j,1 )*(u3d(i,j,2 )-u3d(i,j,1   ))*rdz
            dpdzx(i,j,nk)=gx(i,j,nk)*(u3d(i,j,nk)-u3d(i,j,nk-1))*rdz
          enddo
          enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=2,nk-1
          do j=1,nj+1
          do i=1,ni
            dpdzy(i,j,k)=gy(i,j,k)*(v3d(i,j,k+1)-v3d(i,j,k-1))*rdz2
          enddo
          enddo
          enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj+1
          do i=1,ni
            dpdzy(i,j, 1)=gy(i,j,1 )*(v3d(i,j,2 )-v3d(i,j,1   ))*rdz
            dpdzy(i,j,nk)=gy(i,j,nk)*(v3d(i,j,nk)-v3d(i,j,nk-1))*rdz
          enddo
          enddo

          tem = dts*0.5

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            fpk(i,j,k)=fpk(i,j,k)-tem*ppterm(i,j,k)*(pi0(i,j,k)+pp3d(i,j,k))    &
                       *( (dpdzx(i,j,k)+dpdzx(i+1,j,k))             &
                         +(dpdzy(i,j,k)+dpdzy(i,j+1,k)) )
          enddo
          enddo
          enddo

          IF(thsmall.eq.1.and.neweqts.ge.1.and.imoist.eq.1)THEN

          tem = dts*0.5

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ftk(i,j,k)=ftk(i,j,k)-tem*thterm(i,j,k)*(th0(i,j,k)+th3d(i,j,k))    &
                       *( (dpdzx(i,j,k)+dpdzx(i+1,j,k))             &
                         +(dpdzy(i,j,k)+dpdzy(i,j+1,k)) )
          enddo
          enddo
          enddo

          ENDIF

        ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          fwk(i,j,k)=dts*wten(i,j,k)                            &
             -beta*mm(i,j,k)*(pp3d(i,j,k)-pp3d(i,j,k-1))
        enddo
        enddo
        enddo

        IF(thsmall.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            fwk(i,j,k)=fwk(i,j,k)                             &
               +beta*bk(i,j,k)*(th3d(i,j,k)+th3d(i,j,k-1))
          enddo
          enddo
          enddo

        ENDIF

!--------------------------------------
!  Get boundary conditions for w

        if(terrain_flag)then
          if(timestats.ge.1) time_sound=time_sound+mytime()
          call bcwsfc(dzdx,dzdy,u3d,v3d,w3d)
        else
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            w3d(i,j,1)=0.0
          enddo
          enddo
        endif

!---------------------------------------------
!  the implicit solver ....
 
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dpdzx(i,j,1)=0.0
          dpdzy(i,j,1)=w3d(i,j,1)
        enddo
        enddo

      IF(thsmall.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,aa,bb,cc,dd,r1)
        do j=1,nj
        do k=2,nk
        do i=1,ni
          aa=asq*( mm(i,j,k)*(qk(i,j,k)-pk(i,j,k))   &
                  +bk(i,j,k)*(dts*tk(i,j,k+1)*rr0(i,j,k)-sk(i,j,k)) )
          bb=1.0+asq*( mm(i,j,k)*( qk(i,j,k)+qk(i,j,k-1)                    &
                                  +pk(i,j,k)-pk(i,j,k-1) )                  &
                      -bk(i,j,k)*(dts*tk(i,j,k)*(rr0(i,j,k)+rr0(i,j,k-1))   &
                                  +sk(i,j,k)-sk(i,j,k-1)) )
          cc=asq*( mm(i,j,k)*(pk(i,j,k-1)+qk(i,j,k-1))   &
                  +bk(i,j,k)*(dts*tk(i,j,k-1)*rr0(i,j,k-1)+sk(i,j,k-1)) )
          dd=w3d(i,j,k)+fwk(i,j,k)                                                    &
              -alph*( mm(i,j,k)*(pp3d(i,j,k)+fpk(i,j,k)-pp3d(i,j,k-1)-fpk(i,j,k-1))   &
                     -bk(i,j,k)*(th3d(i,j,k)+ftk(i,j,k)+th3d(i,j,k-1)+ftk(i,j,k-1)) )
          r1 = 1.0/(bb-cc*dpdzx(i,j,k-1))
          dpdzx(i,j,k)=aa*r1
          dpdzy(i,j,k)=(dd+cc*dpdzy(i,j,k-1))*r1
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k,aa,bb,cc,dd,r1)
        do j=1,nj
        do k=2,nk
        do i=1,ni
          aa=asq*mm(i,j,k)*(qk(i,j,k)-pk(i,j,k))
          bb=1.0+asq*mm(i,j,k)*( qk(i,j,k)+qk(i,j,k-1)    &
                            +pk(i,j,k)-pk(i,j,k-1) )
          cc=asq*mm(i,j,k)*(pk(i,j,k-1)+qk(i,j,k-1))
          dd=w3d(i,j,k)+fwk(i,j,k)    &
              -alph*mm(i,j,k)*(pp3d(i,j,k)+fpk(i,j,k)-pp3d(i,j,k-1)-fpk(i,j,k-1))
          r1 = 1.0/(bb-cc*dpdzx(i,j,k-1))
          dpdzx(i,j,k)=aa*r1
          dpdzy(i,j,k)=(dd+cc*dpdzy(i,j,k-1))*r1
        enddo
        enddo
        enddo

      ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          w3d(i,j,nk+1)=0.0
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do j=1,nj
        do k=nk,2,-1
        do i=1,ni
          w3d(i,j,k)=dpdzx(i,j,k)*w3d(i,j,k+1)+dpdzy(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_sound=time_sound+mytime()

        if(n.eq.nloop) call bcw(w3d,1)

        if(terrain_flag) call bcwsfc(dzdx,dzdy,u3d,v3d,w3d)

#ifdef MPI
        if(n.eq.nloop)then
          call comm_3w_start(w3d,ww31,ww32,we31,we32,   &
                                 ws31,ws32,wn31,wn32,reqs_w)
        endif
#endif

!-------------------------
#ifdef MPI
        IF(n.gt.1)THEN
          call comm_1s_end_misc(ppd,pw1,pw2,pe1,pe2,   &
                               ps1,ps2,pn1,pn2,reqs_p)
        ENDIF
#endif
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ppd(i,j,k)=pp3d(i,j,k)
          pp3d(i,j,k)=pp3d(i,j,k)+fpk(i,j,k)+alph*(    &
               w3d(i,j,k  )*(pk(i,j,k)+qk(i,j,k))          &
              +w3d(i,j,k+1)*(pk(i,j,k)-qk(i,j,k)) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          ppd(i,j,k)=pp3d(i,j,k)-ppd(i,j,k)-dts*ppten(i,j,k)
          ppd(i,j,k)=pp3d(i,j,k)+kdiv*ppd(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_sound=time_sound+mytime()

        if(n.eq.nloop .and. nrk.eq.3 .and. imoist.eq.1)then
          tem=0
        else
          if(n.eq.nloop) call bcs(pp3d)
        endif

#ifdef MPI
        if(n.eq.nloop)then
          if(nrk.lt.3.or.imoist.eq.0)then
            call comm_3s_start(pp3d,pw31,pw32,pe31,pe32,   &
                                    ps31,ps32,pn31,pn32,reqs_p)
          endif
        endif
#endif

        IF(n.lt.nloop)THEN
          call bcs(ppd)
#ifdef MPI
          call comm_1s_start(ppd,pw1,pw2,pe1,pe2,   &
                                 ps1,ps2,pn1,pn2,reqs_p)
#endif
        ENDIF

!--------------------------------------------------------------------

      IF(thsmall.eq.1)THEN

        tem = alph*dts

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          th3d(i,j,k)=th3d(i,j,k)+ftk(i,j,k)+tem*rr0(i,j,k)*(    &
            w3d(i,j,k)*tk(i,j,k)+w3d(i,j,k+1)*tk(i,j,k+1)  )
        enddo
        enddo
        enddo

        IF(neweqts.ge.1 .and. imoist.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            th3d(i,j,k)=th3d(i,j,k)+alph*sk(i,j,k)*( w3d(i,j,k)-w3d(i,j,k+1) )
          enddo
          enddo
          enddo
        ENDIF
        if(timestats.ge.1) time_sound=time_sound+mytime()

        if(n.eq.nloop)then
          if(nrk.lt.3.or.imoist.eq.0)then
            call bcs(th3d)
#ifdef MPI
            call comm_3s_start(th3d,sw31,sw32,se31,se32,   &
                                    ss31,ss32,sn31,sn32,reqs_s)
#endif
          endif
        endif

      ENDIF

!--------------------------------------------------------------------

      ENDDO

!!!      if(terrain_flag)then
!!!        if(timestats.ge.1) time_sound=time_sound+mytime()
!!!        call bcwsfc(dzdx,dzdy,u3d,v3d,w3d)
!!!        call bc2d(w3d(ib,jb,1))
!!!      endif

!--------------------------------------------------------------------

      return
      end


