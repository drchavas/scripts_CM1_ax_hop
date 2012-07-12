



      subroutine adv6s(bflag,bsq,xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,      &
                       rho0,rr0,rf0,rrf0,advx,advy,advz,dum,divx,mass,   &
                       rru,rrv,rrw,s1,s,sten,pdef,dt)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      integer bflag
      real*8 bsq
      real, dimension(ib:ie) :: xh,rxh,uh,ruh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je) :: vh,rvh
      real, dimension(itb:ite,jtb:jte) :: gz
      real, dimension(ib:ie,jb:je,kb:ke) :: mh,rmh,rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: advx,advy,advz,dum,divx,mass
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, dimension(ib:ie,jb:je,kb:ke) :: s1,s,sten
      integer pdef
      real dt
 
      integer i,j,k,i1,i2,j1,j2
      real, parameter :: tem = 1.0/60.0

!----------------------------------------------------------------
! Advection in x-direction

      i1 = 1
      i2 = ni+1

      if(wbc.eq.2 .and. ibw.eq.1)then
        i1=4
      endif

      if(ebc.eq.2 .and. ibe.eq.1)then
        i2=ni-2
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=i1,i2
        dum(i,j,k)=rru(i,j,k)*( 37.0*(s(i  ,j,k)+s(i-1,j,k))     &
                               -8.0*(s(i+1,j,k)+s(i-2,j,k))     &
                                   +(s(i+2,j,k)+s(i-3,j,k)) )*tem
      enddo
      enddo
      enddo

      IF(wbc.eq.2 .and. ibw.eq.1) call advbcsw(dum,rru,s)

      IF(ebc.eq.2 .and. ibe.eq.1) call advbcse(xf,dum,rru,s)

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        advx(i,j,k)=-(dum(i+1,j,k)-dum(i,j,k))*rdx*uh(i)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        advx(i,j,k)=-(xf(i+1)*dum(i+1,j,k)-xf(i)*dum(i,j,k))*rdx*uh(i)*rxh(i)
      enddo
      enddo
      enddo

    ENDIF

      IF(stat_qsrc.eq.1.and.(wbc.eq.2.or.ebc.eq.2).and.bflag.eq.1)THEN
        if(timestats.ge.1) time_advs=time_advs+mytime()
        call bsx(dt,bsq,rvh,rmh,dum)
      ENDIF

      IF(pdscheme.eq.1 .and. pdef.eq.1)then
        if(timestats.ge.1) time_advs=time_advs+mytime()
        call pdefx(xh,rho0,advx,dum,mass,s1,dt)
      ENDIF

      IF(wbc.eq.2 .and. ibw.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
          i=1
          if(rru(i,j,k).ge.0.0)then
            advx(i,j,k)=advx(i,j,k)-s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)
          endif
        enddo
        enddo
      ENDIF

      IF(ebc.eq.2 .and. ibe.eq.1)THEN

      IF(axisymm.eq.0)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
          i=ni+1
          if(rru(i,j,k).lt.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)
          endif
        enddo
        enddo
      ELSE
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
          i=ni+1
          if(rru(i,j,k).lt.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-s(i,j,k)*(xf(i+1)*rru(i+1,j,k)-xf(i)*rru(i,j,k))*rdx*uh(i)*rxh(i)
          endif
        enddo
        enddo
      ENDIF

      ENDIF

!----------------------------------------------------------------
! Advection in y-direction

    IF(axisymm.eq.0)THEN

      j1 = 1
      j2 = nj+1

      if(sbc.eq.2 .and. ibs.eq.1)then
        j1=4
      endif

      if(nbc.eq.2 .and. ibn.eq.1)then
        j2=nj-2
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
      do i=1,ni
        dum(i,j,k)=rrv(i,j,k)*( 37.0*(s(i,j  ,k)+s(i,j-1,k))     &
                               -8.0*(s(i,j+1,k)+s(i,j-2,k))     &
                                   +(s(i,j+2,k)+s(i,j-3,k)) )*tem
      enddo
      enddo
      enddo

      IF(sbc.eq.2 .and. ibs.eq.1) call advbcss(dum,rrv,s)

      IF(nbc.eq.2 .and. ibn.eq.1) call advbcsn(dum,rrv,s)

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        advy(i,j,k)=-(dum(i,j+1,k)-dum(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo

      IF(stat_qsrc.eq.1.and.(sbc.eq.2.or.nbc.eq.2).and.bflag.eq.1)THEN
        if(timestats.ge.1) time_advs=time_advs+mytime()
        call bsy(dt,bsq,ruh,rmh,dum)
      ENDIF

      IF(pdscheme.eq.1 .and. pdef.eq.1)THEN
        if(timestats.ge.1) time_advs=time_advs+mytime()
        call pdefy(rho0,advy,dum,mass,s1,dt)
      ENDIF

      IF(sbc.eq.2 .and. ibs.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do i=1,ni
          j=1
          if(rrv(i,j,k).ge.0.0)then
            advy(i,j,k)=advy(i,j,k)-s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)
          endif
        enddo
        enddo
      ENDIF

      IF(nbc.eq.2 .and. ibn.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do i=1,ni
          j=nj+1
          if(rrv(i,j,k).lt.0.0)then
            j=nj
            advy(i,j,k)=advy(i,j,k)-s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)
          endif
        enddo
        enddo
      ENDIF

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        advy(i,j,k)=0.0
      enddo
      enddo
      enddo

    ENDIF

!----------------------------------------------------------------

      if(vadvorder.eq.5)then
        call vadv5s(gz,mh,rho0,advz,dum,mass,rrw,s1,s,pdef,dt)
      elseif(vadvorder.eq.6)then
        call vadv6s(gz,mh,rho0,advz,dum,mass,rrw,s1,s,pdef,dt)
      endif

!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=sten(i,j,k)+( advx(i,j,k)+advy(i,j,k)+advz(i,j,k)    &
                                 +s(i,j,k)*divx(i,j,k) )*rr0(i,j,k)
      enddo
      enddo
      enddo

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advs=time_advs+mytime()
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine vadv6s(gz,mh,rho0,advz,dum,mass,rrw,s1,s,pdef,dt)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(itb:ite,jtb:jte) :: gz
      real, dimension(ib:ie,jb:je,kb:ke) :: mh,rho0
      real, dimension(ib:ie,jb:je,kb:ke) :: advz,dum,mass
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, dimension(ib:ie,jb:je,kb:ke) :: s1,s
      integer pdef
      real dt

      integer i,j,k
      real, parameter :: tem1 = 1.0/60.0
      real, parameter :: tem2 = 1.0/12.0

!----------------------------------------------------------------
! Advection in z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=4,nk-2
      do j=1,nj
      do i=1,ni
        dum(i,j,k)=rrw(i,j,k)                                 &
                            *( 37.0*(s(i,j,k  )+s(i,j,k-1))          &
                               -8.0*(s(i,j,k+1)+s(i,j,k-2))          &
                                   +(s(i,j,k+2)+s(i,j,k-3)) )*tem1
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=3,(nk-1),(nk-4)
      do j=1,nj
      do i=1,ni
        dum(i,j,k)=rrw(i,j,k)                                &
                            *( 7.0*(s(i,j,k  )+s(i,j,k-1))          &
                                  -(s(i,j,k+1)+s(i,j,k-2)) )*tem2
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk,(nk-2)
      do j=1,nj
      do i=1,ni
        dum(i,j,k)=rrw(i,j,k)*0.5*(s(i,j,k-1)+s(i,j,k))
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,(nk+1),nk
      do j=1,nj
      do i=1,ni
        dum(i,j,k)=0.0
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        advz(i,j,k)=-(dum(i,j,k+1)-dum(i,j,k))*rdz*mh(i,j,k)
      enddo
      enddo
      enddo

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          advz(i,j,k)=advz(i,j,k)/gz(i,j)
        enddo
        enddo
        enddo

      ENDIF

      IF(pdscheme.eq.1 .and. pdef.eq.1)THEN
        if(timestats.ge.1) time_advs=time_advs+mytime()
        call pdefz(rho0,advz,dum,mass,s1,dt)
      ENDIF

!----------------------------------------------------------------
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine adv6u(xf,rxf,uf,vh,gz,mh,rho0,rr0,rf0,rrf0,rru0,dum,advx,advy,advz,divx, &
                       rru,u3d,uten,rrv,rrw)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, dimension(ib:ie+1) :: xf,rxf,uf
      real, dimension(jb:je) :: vh
      real, dimension(itb:ite,jtb:jte) :: gz
      real, dimension(ib:ie,jb:je,kb:ke) :: mh,rho0,rr0,rf0,rrf0,rru0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advx,advy,advz,divx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru,u3d,uten
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw
 
      integer i,j,k,i1,i2,j1,j2,id1,id2
      real, parameter :: tem = 1.0/120.0

!------------------------------------------------------------

      if(ibw.eq.1)then
        i1=2
      else
        i1=1
      endif
 
      if(ibe.eq.1)then
        i2=ni+1-1
      else
        i2=ni+1
      endif
 
!----------------------------------------------------------------
! Advection in x-direction

      id1 = i1-1
      id2 = i2

      if(wbc.eq.2 .and. ibw.eq.1)then
        id1 = 3
      endif

      if(ebc.eq.2 .and. ibe.eq.1)then
        id2 = ni-2
      endif

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=id1,id2
        dum(i,j,k)=(rru(i,j,k)+rru(i+1,j,k))*             &
                   ( 37.0*(u3d(i+1,j,k)+u3d(i  ,j,k))     &
                     -8.0*(u3d(i+2,j,k)+u3d(i-1,j,k))     &
                         +(u3d(i+3,j,k)+u3d(i-2,j,k)) )*tem
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=id1,id2
        dum(i,j,k)=(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))*             &
                   ( 37.0*(u3d(i+1,j,k)+u3d(i  ,j,k))     &
                     -8.0*(u3d(i+2,j,k)+u3d(i-1,j,k))     &
                         +(u3d(i+3,j,k)+u3d(i-2,j,k)) )*tem
      enddo
      enddo
      enddo

    ENDIF

      IF(wbc.eq.2 .and. ibw.eq.1) call advbcuw(dum,rru,u3d)

      IF(ebc.eq.2 .and. ibe.eq.1) call advbcue(xf,dum,rru,u3d)

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=i1,i2
        advx(i,j,k)=-(dum(i,j,k)-dum(i-1,j,k))*rdx*uf(i)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=i1,i2
        advx(i,j,k)=-(dum(i,j,k)-dum(i-1,j,k))*rdx*uf(i)*rxf(i)
      enddo
      enddo
      enddo

    ENDIF

!----------------------------------------------------------------
! Advection in y-direction

    IF(axisymm.eq.0)THEN

      j1 = 1
      j2 = nj+1

      if(sbc.eq.2 .and. ibs.eq.1)then
        j1=4
      endif

      if(nbc.eq.2 .and. ibn.eq.1)then
        j2=nj-2
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
      do i=i1,i2
        dum(i,j,k)=(rrv(i,j,k)+rrv(i-1,j,k))*             &
                   ( 37.0*(u3d(i,j  ,k)+u3d(i,j-1,k))     &
                     -8.0*(u3d(i,j+1,k)+u3d(i,j-2,k))     &
                         +(u3d(i,j+2,k)+u3d(i,j-3,k)) )*tem
      enddo
      enddo
      enddo

      IF(sbc.eq.2 .and. ibs.eq.1) call advbcus(dum,u3d,rrv,i1,i2)

      IF(nbc.eq.2 .and. ibn.eq.1) call advbcun(dum,u3d,rrv,i1,i2)

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=i1,i2
        advy(i,j,k)=-(dum(i,j+1,k)-dum(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo

      IF(sbc.eq.2 .and. ibs.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do i=i1,i2
          j=1
          if((rrv(i,j,k)+rrv(i-1,j,k)).ge.0.0)then
            advy(i,j,k)=advy(i,j,k)-u3d(i,j,k)*0.5*(                    &
                            (rrv(i-1,j+1,k)-rrv(i-1,j,k))               &
                           +(rrv(i  ,j+1,k)-rrv(i  ,j,k)) )*rdy*vh(j)
          endif
        enddo
        enddo
      ENDIF

      IF(nbc.eq.2 .and. ibn.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do i=i1,i2
          j=nj+1
          if((rrv(i,j,k)+rrv(i-1,j,k)).lt.0.0)then
            j=nj
            advy(i,j,k)=advy(i,j,k)-u3d(i,j,k)*0.5*(                    &
                            (rrv(i-1,j+1,k)-rrv(i-1,j,k))               &
                           +(rrv(i  ,j+1,k)-rrv(i  ,j,k)) )*rdy*vh(j)
          endif
        enddo
        enddo
      ENDIF

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=i1,i2
        advy(i,j,k)=0.0
      enddo
      enddo
      enddo

    ENDIF

!----------------------------------------------------------------

      if(vadvorder.eq.5)then
        call vadv5u(gz,mh,dum,advz,u3d,rrw,i1,i2)
      elseif(vadvorder.eq.6)then
        call vadv6u(gz,mh,dum,advz,u3d,rrw,i1,i2)
      endif

!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=i1,i2
        uten(i,j,k)=uten(i,j,k)+( advx(i,j,k)+advy(i,j,k)+advz(i,j,k)    &
                   +u3d(i,j,k)*0.5*(divx(i,j,k)+divx(i-1,j,k)) )*rru0(i,j,k)
      enddo
      enddo
      enddo

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advu=time_advu+mytime()
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine vadv6u(gz,mh,dum,advz,u3d,rrw,i1,i2)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(itb:ite,jtb:jte) :: gz
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advz
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      integer :: i1,i2

      integer i,j,k
      real, parameter :: tem1 = 1.0/120.0
      real, parameter :: tem2 = 1.0/24.0

!----------------------------------------------------------------
! Advection in z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=4,nk-2
      do j=1,nj
      do i=i1,i2
        dum(i,j,k)=(rrw(i,j,k)+rrw(i-1,j,k))*           &
                   ( 37.0*(u3d(i,j,k  )+u3d(i,j,k-1))          &
                     -8.0*(u3d(i,j,k+1)+u3d(i,j,k-2))          &
                         +(u3d(i,j,k+2)+u3d(i,j,k-3)) )*tem1
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=3,(nk-1),(nk-4)
      do j=1,nj
      do i=i1,i2
        dum(i,j,k)=(rrw(i,j,k)+rrw(i-1,j,k))*          &
                   ( 7.0*(u3d(i,j,k  )+u3d(i,j,k-1))          &
                        -(u3d(i,j,k+1)+u3d(i,j,k-2)) )*tem2
      enddo
      enddo
      enddo
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk,(nk-2)
      do j=1,nj
      do i=i1,i2
        dum(i,j,k)=0.25*(rrw(i,j,k)+rrw(i-1,j,k))    &
                       *(u3d(i,j,k-1)+u3d(i,j,k))
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,(nk+1),nk
      do j=1,nj
      do i=i1,i2
        dum(i,j,k)=0.0
      enddo
      enddo
      enddo

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=i1,i2
          advz(i,j,k)=-(dum(i,j,k+1)-dum(i,j,k))*rdz*0.5*(mh(i-1,j,k)/gz(i-1,j)+mh(i,j,k)/gz(i,j))
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=i1,i2
          advz(i,j,k)=-(dum(i,j,k+1)-dum(i,j,k))*rdz*0.5*(mh(i-1,j,k)+mh(i,j,k))
        enddo
        enddo
        enddo

      ENDIF

!----------------------------------------------------------------
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine adv6v(xh,rxh,uh,xf,vf,gz,mh,rho0,rr0,rf0,rrf0,rrv0,dum,advx,advy,advz,divx, &
                       rru,rrv,v3d,vten,rrw)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, dimension(ib:ie) :: xh,rxh,uh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je+1) :: vf
      real, dimension(itb:ite,jtb:jte) :: gz
      real, dimension(ib:ie,jb:je,kb:ke) :: mh,rho0,rr0,rf0,rrf0,rrv0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advx,advy,advz,divx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv,v3d,vten
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw
 
      integer i,j,k,i1,i2,j1,j2,jd1,jd2
      real, parameter :: tem = 1.0/120.0

!------------------------------------------------------------

      if(ibs.eq.1)then
        j1=2
      else
        j1=1
      endif
 
      if(ibn.eq.1)then
        j2=nj+1-1
      else
        j2=nj+1
      endif
 
!----------------------------------------------------------------
! Advection in x-direction

      i1 = 1
      i2 = ni+1

      if(wbc.eq.2 .and. ibw.eq.1)then
        i1=4
      endif

      if(ebc.eq.2 .and. ibe.eq.1)then
        i2=ni-2
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
      do i=i1,i2
        dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*             &
                   ( 37.0*(v3d(i  ,j,k)+v3d(i-1,j,k))     &
                     -8.0*(v3d(i+1,j,k)+v3d(i-2,j,k))     &
                         +(v3d(i+2,j,k)+v3d(i-3,j,k)) )*tem
      enddo
      enddo
      enddo

      IF(wbc.eq.2 .and. ibw.eq.1) call advbcvw(dum,rru,v3d,j1,j2)

      IF(ebc.eq.2 .and. ibe.eq.1) call advbcve(xf,dum,rru,v3d,j1,j2)

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
      do i=1,ni
        advx(i,j,k)=-(dum(i+1,j,k)-dum(i,j,k))*rdx*uh(i)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
      do i=1,ni
        advx(i,j,k)=-(xf(i+1)*dum(i+1,j,k)-xf(i)*dum(i,j,k))*rdx*uh(i)*rxh(i)
      enddo
      enddo
      enddo

    ENDIF

      IF(wbc.eq.2 .and. ibw.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=j1,j2
          i=1
          if((rru(i,j,k)+rru(i,j-1,k)).ge.0.0)then
            advx(i,j,k)=advx(i,j,k)-v3d(i,j,k)*0.5*(            &
                    (rru(i+1,j-1,k)-rru(i,j-1,k))               &
                   +(rru(i+1,j  ,k)-rru(i,j  ,k)) )*rdx*uh(i)
          endif
        enddo
        enddo
      ENDIF

      IF(ebc.eq.2 .and. ibe.eq.1)THEN

      IF(axisymm.eq.0)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=j1,j2
          i=ni+1
          if((rru(i,j,k)+rru(i,j-1,k)).lt.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-v3d(i,j,k)*0.5*(            &
                    (rru(i+1,j-1,k)-rru(i,j-1,k))               &
                   +(rru(i+1,j  ,k)-rru(i,j  ,k)) )*rdx*uh(i)
          endif
        enddo
        enddo
      ELSE
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=j1,j2
          i=ni+1
          if((rru(i,j,k)+rru(i,j-1,k)).lt.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-v3d(i,j,k)*0.5*(            &
                    (xf(i+1)*rru(i+1,j-1,k)-xf(i)*rru(i,j-1,k))               &
                   +(xf(i+1)*rru(i+1,j  ,k)-xf(i)*rru(i,j  ,k)) )*rdx*uh(i)*rxh(i)
          endif
        enddo
        enddo
      ENDIF

      ENDIF

!----------------------------------------------------------------
! Advection in y-direction

    IF(axisymm.eq.0)THEN

      jd1 = j1-1
      jd2 = j2

      if(sbc.eq.2 .and. ibs.eq.1)then
        jd1 = 3
      endif

      if(nbc.eq.2 .and. ibn.eq.1)then
        jd2 = nj-2
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=jd1,jd2
      do i=1,ni
        dum(i,j,k)=(rrv(i,j,k)+rrv(i,j+1,k))*             &
                   ( 37.0*(v3d(i,j+1,k)+v3d(i,j  ,k))     &
                     -8.0*(v3d(i,j+2,k)+v3d(i,j-1,k))     &
                         +(v3d(i,j+3,k)+v3d(i,j-2,k)) )*tem
      enddo
      enddo
      enddo

      IF(sbc.eq.2 .and. ibs.eq.1) call advbcvs(dum,rrv,v3d)

      IF(nbc.eq.2 .and. ibn.eq.1) call advbcvn(dum,rrv,v3d)

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
      do i=1,ni
        advy(i,j,k)=-(dum(i,j,k)-dum(i,j-1,k))*rdy*vf(j)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
      do i=1,ni
        advy(i,j,k)=0.0
      enddo
      enddo
      enddo

    ENDIF

!----------------------------------------------------------------

      if(vadvorder.eq.5)then
        call vadv5v(gz,mh,dum,advz,v3d,rrw,j1,j2)
      elseif(vadvorder.eq.6)then
        call vadv6v(gz,mh,dum,advz,v3d,rrw,j1,j2)
      endif

!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
      do i=1,ni
        vten(i,j,k)=vten(i,j,k)+( advx(i,j,k)+advy(i,j,k)+advz(i,j,k)    &
                   +v3d(i,j,k)*0.5*(divx(i,j,k)+divx(i,j-1,k)) )*rrv0(i,j,k)
      enddo
      enddo
      enddo

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advv=time_advv+mytime()
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine vadv6v(gz,mh,dum,advz,v3d,rrw,j1,j2)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(itb:ite,jtb:jte) :: gz
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advz
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      integer :: j1,j2

      integer i,j,k
      real, parameter :: tem1 = 1.0/120.0
      real, parameter :: tem2 = 1.0/24.0

!----------------------------------------------------------------
! Advection in z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=4,nk-2
      do j=j1,j2
      do i=1,ni
        dum(i,j,k)=(rrw(i,j,k)+rrw(i,j-1,k))*           &
                   ( 37.0*(v3d(i,j,k  )+v3d(i,j,k-1))          &
                     -8.0*(v3d(i,j,k+1)+v3d(i,j,k-2))          &
                         +(v3d(i,j,k+2)+v3d(i,j,k-3)) )*tem1
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=3,(nk-1),(nk-4)
      do j=j1,j2
      do i=1,ni
        dum(i,j,k)=(rrw(i,j,k)+rrw(i,j-1,k))*          &
                   ( 7.0*(v3d(i,j,k  )+v3d(i,j,k-1))          &
                        -(v3d(i,j,k+1)+v3d(i,j,k-2)) )*tem2
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk,(nk-2)
      do j=j1,j2
      do i=1,ni
        dum(i,j,k)=0.25*(rrw(i,j,k)+rrw(i,j-1,k))    &
                       *(v3d(i,j,k-1)+v3d(i,j,k))
      enddo
      enddo
      enddo
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,(nk+1),nk
      do j=j1,j2
      do i=1,ni
        dum(i,j,k)=0.0
      enddo
      enddo
      enddo

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=j1,j2
        do i=1,ni
          advz(i,j,k)=-(dum(i,j,k+1)-dum(i,j,k))*rdz*0.5*(mh(i,j-1,k)/gz(i,j-1)+mh(i,j,k)/gz(i,j))
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=j1,j2
        do i=1,ni
          advz(i,j,k)=-(dum(i,j,k+1)-dum(i,j,k))*rdz*0.5*(mh(i,j-1,k)+mh(i,j,k))
        enddo
        enddo
        enddo

      ENDIF

!----------------------------------------------------------------
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine adv6w(xh,rxh,uh,xf,vh,gz,mf,rho0,rr0,rf0,rrf0,dum,advx,advy,advz,divx, &
                       rru,rrv,rrw,w3d,wten)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, dimension(ib:ie) :: xh,rxh,uh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je) :: vh
      real, dimension(itb:ite,jtb:jte) :: gz
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advx,advy,advz,divx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw,w3d,wten
 
      integer i,j,k,i1,i2,j1,j2
      real, parameter :: tem = 1.0/120.0

!----------------------------------------------------------------
! Advection in x-direction

      i1 = 1
      i2 = ni+1

      if(wbc.eq.2 .and. ibw.eq.1)then
        i1=4
      endif

      if(ebc.eq.2 .and. ibe.eq.1)then
        i2=ni-2
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=i1,i2
        dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*             &
                   ( 37.0*(w3d(i  ,j,k)+w3d(i-1,j,k))     &
                     -8.0*(w3d(i+1,j,k)+w3d(i-2,j,k))     &
                         +(w3d(i+2,j,k)+w3d(i-3,j,k)) )*tem
      enddo
      enddo
      enddo

      IF(wbc.eq.2 .and. ibw.eq.1) call advbcww(dum,rru,w3d)

      IF(ebc.eq.2 .and. ibe.eq.1) call advbcwe(xf,dum,rru,w3d)

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        advx(i,j,k)=-(dum(i+1,j,k)-dum(i,j,k))*rdx*uh(i)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        advx(i,j,k)=-(xf(i+1)*dum(i+1,j,k)-xf(i)*dum(i,j,k))*rdx*uh(i)*rxh(i)
      enddo
      enddo
      enddo

    ENDIF

      IF(wbc.eq.2 .and. ibw.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
          i=1
          if((rru(i,j,k)+rru(i,j,k-1)).ge.0.0)then
            advx(i,j,k)=advx(i,j,k)-w3d(i,j,k)*0.5*(      &
                    (rru(i+1,j,k-1)-rru(i,j,k-1))         &
                   +(rru(i+1,j,k  )-rru(i,j,k  )) )*rdx*uh(i)
          endif
        enddo
        enddo
      ENDIF

      IF(ebc.eq.2 .and. ibe.eq.1)THEN

      IF(axisymm.eq.0)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
          i=ni+1
          if((rru(i,j,k)+rru(i,j,k-1)).lt.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-w3d(i,j,k)*0.5*(      &
                    (rru(i+1,j,k-1)-rru(i,j,k-1))         &
                   +(rru(i+1,j,k  )-rru(i,j,k  )) )*rdx*uh(i)
          endif
        enddo
        enddo
      ELSE
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
          i=ni+1
          if((rru(i,j,k)+rru(i,j,k-1)).lt.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-w3d(i,j,k)*0.5*(      &
                    (xf(i+1)*rru(i+1,j,k-1)-xf(i)*rru(i,j,k-1))         &
                   +(xf(i+1)*rru(i+1,j,k  )-xf(i)*rru(i,j,k  )) )*rdx*uh(i)*rxh(i)
          endif
        enddo
        enddo
      ENDIF

      ENDIF

!----------------------------------------------------------------
! Advection in y-direction

    IF(axisymm.eq.0)THEN

      j1 = 1
      j2 = nj+1

      if(sbc.eq.2 .and. ibs.eq.1)then
        j1=4
      endif

      if(nbc.eq.2 .and. ibn.eq.1)then
        j2=nj-2
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=j1,j2
      do i=1,ni
        dum(i,j,k)=(rrv(i,j,k)+rrv(i,j,k-1))*             &
                   ( 37.0*(w3d(i,j  ,k)+w3d(i,j-1,k))     &
                     -8.0*(w3d(i,j+1,k)+w3d(i,j-2,k))     &
                         +(w3d(i,j+2,k)+w3d(i,j-3,k)) )*tem
      enddo
      enddo
      enddo

      IF(sbc.eq.2 .and. ibs.eq.1) call advbcws(dum,rrv,w3d)

      IF(nbc.eq.2 .and. ibn.eq.1) call advbcwn(dum,rrv,w3d)

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        advy(i,j,k)=-(dum(i,j+1,k)-dum(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo

      IF(sbc.eq.2 .and. ibs.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do i=1,ni
          j=1
          if((rrv(i,j,k)+rrv(i,j,k-1)).ge.0.0)then
            advy(i,j,k)=advy(i,j,k)-w3d(i,j,k)*0.5*(       &
                           (rrv(i,j+1,k-1)-rrv(i,j,k-1))   &
                          +(rrv(i,j+1,k  )-rrv(i,j,k  )) )*rdy*vh(j)
          endif
        enddo
        enddo
      ENDIF

      IF(nbc.eq.2 .and. ibn.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do i=1,ni
          j=nj+1
          if((rrv(i,j,k)+rrv(i,j,k-1)).lt.0.0)then
            j=nj
            advy(i,j,k)=advy(i,j,k)-w3d(i,j,k)*0.5*(       &
                           (rrv(i,j+1,k-1)-rrv(i,j,k-1))   &
                          +(rrv(i,j+1,k  )-rrv(i,j,k  )) )*rdy*vh(j)
          endif
        enddo
        enddo
      ENDIF

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        advy(i,j,k)=0.0
      enddo
      enddo
      enddo

    ENDIF

!----------------------------------------------------------------

      if(vadvorder.eq.5)then
        call vadv5w(gz,mf,dum,advz,rrw,w3d)
      elseif(vadvorder.eq.6)then
        call vadv6w(gz,mf,dum,advz,rrw,w3d)
      endif
 
!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        wten(i,j,k)=wten(i,j,k)+( advx(i,j,k)+advy(i,j,k)+advz(i,j,k)    &
                   +w3d(i,j,k)*0.5*(divx(i,j,k)+divx(i,j,k-1)) )*rrf0(i,j,k)
      enddo
      enddo
      enddo

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advw=time_advw+mytime()
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine vadv6w(gz,mf,dum,advz,rrw,w3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(itb:ite,jtb:jte) :: gz
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advz
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw,w3d

      integer i,j,k
      real, parameter :: tem1 = 1.0/120.0
      real, parameter :: tem2 = 1.0/24.0

!----------------------------------------------------------------
! Advection in z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=3,nk-2
      do j=1,nj
      do i=1,ni
        dum(i,j,k)=(rrw(i,j,k)+rrw(i,j,k+1))*  &
                   ( 37.0*(w3d(i,j,k+1)+w3d(i,j,k  ))          &
                     -8.0*(w3d(i,j,k+2)+w3d(i,j,k-1))          &
                         +(w3d(i,j,k+3)+w3d(i,j,k-2)) )*tem1
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,(nk-1),(nk-3)
      do j=1,nj
      do i=1,ni
        dum(i,j,k)=(rrw(i,j,k)+rrw(i,j,k+1))*         &
                   ( 7.0*(w3d(i,j,k+1)+w3d(i,j,k  ))                  &
                        -(w3d(i,j,k+2)+w3d(i,j,k-1)) )*tem2
      enddo
      enddo
      enddo
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk,(nk-1)
      do j=1,nj
      do i=1,ni
        dum(i,j,k)=0.25*(rrw(i,j,k)+rrw(i,j,k+1))    &
                       *(w3d(i,j,k)+w3d(i,j,k+1))
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        advz(i,j,k)=-(dum(i,j,k)-dum(i,j,k-1))*rdz*mf(i,j,k)
      enddo
      enddo
      enddo

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          advz(i,j,k)=advz(i,j,k)/gz(i,j)
        enddo
        enddo
        enddo

      ENDIF

!----------------------------------------------------------------
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcsw(dum,rru,s)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je,kb:ke) :: s

      integer i,j,k
      real, parameter :: tem = 1.0/6.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        i=3
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*(-s(i-2,j,k)+5.*s(i-1,j,k)+2.*s(i,j,k))*tem
        else
          dum(i,j,k)=rru(i,j,k)*(-s(i+1,j,k)+5.*s(i,j,k)+2.*s(i-1,j,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        i=2
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*s(i-1,j,k)
        else
          dum(i,j,k)=rru(i,j,k)*s(i,j,k)
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        i=1
        if(rru(i,j,k).lt.0.0)then
          dum(i,j,k)=rru(i,j,k)*s(i,j,k)
        else
          dum(i,j,k)=dum(i+1,j,k)
        endif
      enddo
      enddo

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcse(xf,dum,rru,s)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie+1) :: xf
      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je,kb:ke) :: s

      integer i,j,k
      real, parameter :: tem = 1.0/6.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        i=ni-1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*(-s(i-2,j,k)+5.*s(i-1,j,k)+2.*s(i,j,k))*tem
        else
          dum(i,j,k)=rru(i,j,k)*(-s(i+1,j,k)+5.*s(i,j,k)+2.*s(i-1,j,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        i=ni
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*s(i-1,j,k)
        else
          dum(i,j,k)=rru(i,j,k)*s(i,j,k)
        endif
      enddo
      enddo

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        i=ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*s(i-1,j,k)
        else
          dum(i,j,k)=dum(i-1,j,k)
        endif
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        i=ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*s(i-1,j,k)
        else
          dum(i,j,k)=dum(i-1,j,k)*xf(ni)/xf(ni+1)
        endif
      enddo
      enddo

    ENDIF

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcss(dum,rrv,s)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, dimension(ib:ie,jb:je,kb:ke) :: s

      integer i,j,k
      real, parameter :: tem = 1.0/6.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=1,ni
        j=3
        if(rrv(i,j,k).ge.0.0)then
          dum(i,j,k)=rrv(i,j,k)*(-s(i,j-2,k)+5.*s(i,j-1,k)+2.*s(i,j,k))*tem
        else
          dum(i,j,k)=rrv(i,j,k)*(-s(i,j+1,k)+5.*s(i,j,k)+2.*s(i,j-1,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=1,ni
        j=2
        if(rrv(i,j,k).ge.0.0)then
          dum(i,j,k)=rrv(i,j,k)*s(i,j-1,k)
        else
          dum(i,j,k)=rrv(i,j,k)*s(i,j,k)
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=1,ni
        j=1
        if(rrv(i,j,k).lt.0.0)then
          dum(i,j,k)=rrv(i,j,k)*s(i,j,k)
        else
          dum(i,j,k)=dum(i,j+1,k)
        endif
      enddo
      enddo

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcsn(dum,rrv,s)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, dimension(ib:ie,jb:je,kb:ke) :: s

      integer i,j,k
      real, parameter :: tem = 1.0/6.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=1,ni
        j=nj-1
        if(rrv(i,j,k).ge.0.0)then
          dum(i,j,k)=rrv(i,j,k)*(-s(i,j-2,k)+5.*s(i,j-1,k)+2.*s(i,j,k))*tem
        else
          dum(i,j,k)=rrv(i,j,k)*(-s(i,j+1,k)+5.*s(i,j,k)+2.*s(i,j-1,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=1,ni
        j=nj
        if(rrv(i,j,k).ge.0.0)then
          dum(i,j,k)=rrv(i,j,k)*s(i,j-1,k)
        else
          dum(i,j,k)=rrv(i,j,k)*s(i,j,k)
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=1,ni
        j=nj+1
        if(rrv(i,j,k).ge.0.0)then
          dum(i,j,k)=rrv(i,j,k)*s(i,j-1,k)
        else
          dum(i,j,k)=dum(i,j-1,k)
        endif
      enddo
      enddo

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcuw(dum,rru,u3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru,u3d

      integer i,j,k
      real, parameter :: tem = 1.0/12.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        i=2
        if((rru(i,j,k)+rru(i+1,j,k)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i+1,j,k))*                     &
                 (-u3d(i-1,j,k)+5.*u3d(i  ,j,k)+2.*u3d(i+1,j,k))*tem
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i+1,j,k))*                       &
                 (-u3d(i+2,j,k)+5.*u3d(i+1,j,k)+2.*u3d(i  ,j,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        i=1
        if((rru(i,j,k)+rru(i+1,j,k)).ge.0.)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i+1,j,k))*u3d(i  ,j,k)
        else
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i+1,j,k))*u3d(i+1,j,k)
        endif
      enddo
      enddo

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcue(xf,dum,rru,u3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie+1) :: xf
      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru,u3d

      integer i,j,k
      real, parameter :: tem = 1.0/12.0

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        i=ni-1
        if((rru(i,j,k)+rru(i+1,j,k)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i+1,j,k))*                     &
                 (-u3d(i-1,j,k)+5.*u3d(i  ,j,k)+2.*u3d(i+1,j,k))*tem
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i+1,j,k))*                       &
                 (-u3d(i+2,j,k)+5.*u3d(i+1,j,k)+2.*u3d(i  ,j,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        i=ni
        if((rru(i,j,k)+rru(i+1,j,k)).ge.0.)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i+1,j,k))*u3d(i  ,j,k)
        else
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i+1,j,k))*u3d(i+1,j,k)
        endif
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        i=ni-1
        if((rru(i,j,k)+rru(i+1,j,k)).ge.0.)then
          dum(i,j,k)=(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))*                     &
                 (-u3d(i-1,j,k)+5.*u3d(i  ,j,k)+2.*u3d(i+1,j,k))*tem
        else
          dum(i,j,k)=(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))*                       &
                 (-u3d(i+2,j,k)+5.*u3d(i+1,j,k)+2.*u3d(i  ,j,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        i=ni
        if((rru(i,j,k)+rru(i+1,j,k)).ge.0.)then
          dum(i,j,k)=0.5*(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))*u3d(i  ,j,k)
        else
          dum(i,j,k)=0.5*(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))*u3d(i+1,j,k)
        endif
      enddo
      enddo

    ENDIF

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcus(dum,u3d,rrv,i1,i2)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      integer i1,i2

      integer i,j,k
      real, parameter :: tem = 1.0/12.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=i1,i2
        j=3
        if((rrv(i,j,k)+rrv(i-1,j,k)).ge.0.0)then
          dum(i,j,k)=(rrv(i,j,k)+rrv(i-1,j,k))    &
                    *(-u3d(i,j-2,k)+5.*u3d(i,j-1,k)+2.*u3d(i,j,k))*tem
        else
          dum(i,j,k)=(rrv(i,j,k)+rrv(i-1,j,k))    &
                    *(-u3d(i,j+1,k)+5.*u3d(i,j,k)+2.*u3d(i,j-1,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=i1,i2
        j=2
        if((rrv(i,j,k)+rrv(i-1,j,k)).ge.0.0)then
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i-1,j,k))*u3d(i,j-1,k)
        else
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i-1,j,k))*u3d(i,j,k)
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=i1,i2
        j=1
        if((rrv(i,j,k)+rrv(i-1,j,k)).lt.0.0)then
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i-1,j,k))*u3d(i,j,k)
        else
          dum(i,j,k)=dum(i,j+1,k)
        endif
      enddo
      enddo

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcun(dum,u3d,rrv,i1,i2)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      integer i1,i2

      integer i,j,k
      real, parameter :: tem = 1.0/12.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=i1,i2
        j=nj-1
        if((rrv(i,j,k)+rrv(i-1,j,k)).ge.0.0)then
          dum(i,j,k)=(rrv(i,j,k)+rrv(i-1,j,k))    &
                    *(-u3d(i,j-2,k)+5.*u3d(i,j-1,k)+2.*u3d(i,j,k))*tem
        else
          dum(i,j,k)=(rrv(i,j,k)+rrv(i-1,j,k))    &
                    *(-u3d(i,j+1,k)+5.*u3d(i,j,k)+2.*u3d(i,j-1,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=i1,i2
        j=nj
        if((rrv(i,j,k)+rrv(i-1,j,k)).ge.0.0)then
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i-1,j,k))*u3d(i,j-1,k)
        else
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i-1,j,k))*u3d(i,j,k)
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=i1,i2
        j=nj+1
        if((rrv(i,j,k)+rrv(i-1,j,k)).ge.0.0)then
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i-1,j,k))*u3d(i,j-1,k)
        else
          dum(i,j,k)=dum(i,j-1,k)
        endif
      enddo
      enddo

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcvw(dum,rru,v3d,j1,j2)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      integer j1,j2

      integer i,j,k
      real, parameter :: tem = 1.0/12.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
        i=3
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*    &
                 (-v3d(i-2,j,k)+5.*v3d(i-1,j,k)+2.*v3d(i  ,j,k))*tem
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*    &
                 (-v3d(i+1,j,k)+5.*v3d(i  ,j,k)+2.*v3d(i-1,j,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
        i=2
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j-1,k))*v3d(i-1,j,k)
        else
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j-1,k))*v3d(i  ,j,k)
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
        i=1
        if((rru(i,j,k)+rru(i,j-1,k)).lt.0.0)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j-1,k))*v3d(i,j,k)
        else
          dum(i,j,k)=dum(i+1,j,k)
        endif
      enddo
      enddo

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcve(xf,dum,rru,v3d,j1,j2)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie+1) :: xf
      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      integer j1,j2

      integer i,j,k
      real, parameter :: tem = 1.0/12.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
        i=ni-1
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*    &
                 (-v3d(i-2,j,k)+5.*v3d(i-1,j,k)+2.*v3d(i  ,j,k))*tem
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*    &
                 (-v3d(i+1,j,k)+5.*v3d(i  ,j,k)+2.*v3d(i-1,j,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
        i=ni
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j-1,k))*v3d(i-1,j,k)
        else
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j-1,k))*v3d(i  ,j,k)
        endif
      enddo
      enddo

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
        i=ni+1
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.0)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j-1,k))*v3d(i-1,j,k)
        else
          dum(i,j,k)=dum(i-1,j,k)
        endif
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
        i=ni+1
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.0)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j-1,k))*v3d(i-1,j,k)
        else
          dum(i,j,k)=dum(i-1,j,k)*xf(ni)/xf(ni+1)
        endif
      enddo
      enddo

    ENDIF

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcvs(dum,rrv,v3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv,v3d

      integer i,j,k
      real, parameter :: tem = 1.0/12.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=1,ni
        j=2
        if((rrv(i,j,k)+rrv(i,j+1,k)).ge.0.)then
          dum(i,j,k)=(rrv(i,j,k)+rrv(i,j+1,k))*                     &
                 (-v3d(i,j-1,k)+5.*v3d(i,j  ,k)+2.*v3d(i,j+1,k))*tem
        else
          dum(i,j,k)=(rrv(i,j,k)+rrv(i,j+1,k))*                       &
                 (-v3d(i,j+2,k)+5.*v3d(i,j+1,k)+2.*v3d(i,j  ,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=1,ni
        j=1
        if((rrv(i,j,k)+rrv(i,j+1,k)).ge.0.)then
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i,j+1,k))*v3d(i,j  ,k)
        else
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i,j+1,k))*v3d(i,j+1,k)
        endif
      enddo
      enddo

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcvn(dum,rrv,v3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv,v3d

      integer i,j,k
      real, parameter :: tem = 1.0/12.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=1,ni
        j=nj-1
        if((rrv(i,j,k)+rrv(i,j+1,k)).ge.0.)then
          dum(i,j,k)=(rrv(i,j,k)+rrv(i,j+1,k))*                     &
                 (-v3d(i,j-1,k)+5.*v3d(i,j  ,k)+2.*v3d(i,j+1,k))*tem
        else
          dum(i,j,k)=(rrv(i,j,k)+rrv(i,j+1,k))*                       &
                 (-v3d(i,j+2,k)+5.*v3d(i,j+1,k)+2.*v3d(i,j  ,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do i=1,ni
        j=nj
        if((rrv(i,j,k)+rrv(i,j+1,k)).ge.0.)then
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i,j+1,k))*v3d(i,j  ,k)
        else
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i,j+1,k))*v3d(i,j+1,k)
        endif
      enddo
      enddo

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcww(dum,rru,w3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w3d

      integer i,j,k
      real, parameter :: tem = 1.0/12.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
        i=3
        if((rru(i,j,k)+rru(i,j,k-1)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*     &
                 (-w3d(i-2,j,k)+5.*w3d(i-1,j,k)+2.*w3d(i  ,j,k))*tem
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*     &
                 (-w3d(i+1,j,k)+5.*w3d(i  ,j,k)+2.*w3d(i-1,j,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
        i=2
        if((rru(i,j,k)+rru(i,j,k-1)).ge.0.)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j,k-1))*w3d(i-1,j,k)
        else
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j,k-1))*w3d(i  ,j,k)
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
        i=1
        if((rru(i,j,k)+rru(i,j,k-1)).lt.0.0)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j,k-1))*w3d(i,j,k)
        else
          dum(i,j,k)=dum(i+1,j,k)
        endif
      enddo
      enddo

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcwe(xf,dum,rru,w3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie+1) :: xf
      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w3d

      integer i,j,k
      real, parameter :: tem = 1.0/12.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
        i=ni-1
        if((rru(i,j,k)+rru(i,j,k-1)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*     &
                 (-w3d(i-2,j,k)+5.*w3d(i-1,j,k)+2.*w3d(i  ,j,k))*tem
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*     &
                 (-w3d(i+1,j,k)+5.*w3d(i  ,j,k)+2.*w3d(i-1,j,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
        i=ni
        if((rru(i,j,k)+rru(i,j,k-1)).ge.0.)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j,k-1))*w3d(i-1,j,k)
        else
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j,k-1))*w3d(i  ,j,k)
        endif
      enddo
      enddo

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
        i=ni+1
        if((rru(i,j,k)+rru(i,j,k-1)).ge.0.0)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j,k-1))*w3d(i-1,j,k)
        else
          dum(i,j,k)=dum(i-1,j,k)
        endif
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
        i=ni+1
        if((rru(i,j,k)+rru(i,j,k-1)).ge.0.0)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j,k-1))*w3d(i-1,j,k)
        else
          dum(i,j,k)=dum(i-1,j,k)*xf(ni)/xf(ni+1)
        endif
      enddo
      enddo

    ENDIF

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcws(dum,rrv,w3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w3d

      integer i,j,k
      real, parameter :: tem = 1.0/12.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do i=1,ni
        j=3
        if((rrv(i,j,k)+rrv(i,j,k-1)).ge.0.0)then
          dum(i,j,k)=(rrv(i,j,k)+rrv(i,j,k-1))    &
                    *(-w3d(i,j-2,k)+5.*w3d(i,j-1,k)+2.*w3d(i,j,k))*tem
        else
          dum(i,j,k)=(rrv(i,j,k)+rrv(i,j,k-1))    &
                    *(-w3d(i,j+1,k)+5.*w3d(i,j,k)+2.*w3d(i,j-1,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do i=1,ni
        j=2
        if((rrv(i,j,k)+rrv(i,j,k-1)).ge.0.0)then
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i,j,k-1))*w3d(i,j-1,k)
        else
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i,j,k-1))*w3d(i,j,k)
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do i=1,ni
        j=1
        if((rrv(i,j,k)+rrv(i,j,k-1)).lt.0.0)then
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i,j,k-1))*w3d(i,j,k)
        else
          dum(i,j,k)=dum(i,j+1,k)
        endif
      enddo
      enddo

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advbcwn(dum,rrv,w3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: dum
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w3d

      integer i,j,k
      real, parameter :: tem = 1.0/12.0

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do i=1,ni
        j=nj-1
        if((rrv(i,j,k)+rrv(i,j,k-1)).ge.0.0)then
          dum(i,j,k)=(rrv(i,j,k)+rrv(i,j,k-1))    &
                    *(-w3d(i,j-2,k)+5.*w3d(i,j-1,k)+2.*w3d(i,j,k))*tem
        else
          dum(i,j,k)=(rrv(i,j,k)+rrv(i,j,k-1))    &
                    *(-w3d(i,j+1,k)+5.*w3d(i,j,k)+2.*w3d(i,j-1,k))*tem
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do i=1,ni
        j=nj
        if((rrv(i,j,k)+rrv(i,j,k-1)).ge.0.0)then
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i,j,k-1))*w3d(i,j-1,k)
        else
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i,j,k-1))*w3d(i,j,k)
        endif
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do i=1,ni
        j=nj+1
        if((rrv(i,j,k)+rrv(i,j,k-1)).ge.0.0)then
          dum(i,j,k)=0.5*(rrv(i,j,k)+rrv(i,j,k-1))*w3d(i,j-1,k)
        else
          dum(i,j,k)=dum(i,j-1,k)
        endif
      enddo
      enddo

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefx(xh,rho0,advx,dum,mass,s1,dt)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie) :: xh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0
      real, dimension(ib:ie,jb:je,kb:ke) :: advx,dum,mass,s1
      real :: dt

      integer i,j,k
      real foo1,foo2,foo3,rdt
      logical, dimension(-1:ni+2) :: flag








!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=rho0(i,j,k)*s1(i,j,k)+dt*advx(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pdef=time_pdef+mytime()

        call bcs(dum)




        rdt=1.0/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=-1,ni+2
          mass(i,j,k)=0.0
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pdef=time_pdef+mytime()





        if(ibw.eq.1.and.wbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            dum(-1,j,k)=0.0
            dum( 0,j,k)=0.0
            dum( 1,j,k)=0.0
          enddo
          enddo
        endif

        if(ibe.eq.1.and.ebc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            dum(ni  ,j,k)=0.0
            dum(ni+1,j,k)=0.0
            dum(ni+2,j,k)=0.0
          enddo
          enddo
        endif

      IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,foo1,foo2,foo3,flag)
        do k=1,nk
        do j=1,nj
        do i=-1,ni+2
          flag(i)=.false.
        enddo
        do i=0,ni+1
          if(dum(i,j,k).lt.0.0)then
            foo1=max(0.0,dum(i-1,j,k))
            foo2=max(0.0,dum(i+1,j,k))
            if(foo1+foo2.gt.1.0e-30)then
              foo3=max(dum(i,j,k),-(foo1+foo2))/(foo1+foo2)
              mass(i-1,j,k)=mass(i-1,j,k)+foo1*foo3
              mass(i  ,j,k)=mass(i  ,j,k)-(foo1+foo2)*foo3
              mass(i+1,j,k)=mass(i+1,j,k)+foo2*foo3
              if(dum(i-1,j,k).gt.1.0e-30) flag(i-1)=.true.
                                          flag(i  )=.true.
              if(dum(i+1,j,k).gt.1.0e-30) flag(i+1)=.true.
            endif
          endif
        enddo
        do i=1,ni
        if(flag(i))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advx(i,j,k)=(dum(i,j,k)-rho0(i,j,k)*s1(i,j,k))*rdt
        endif
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k,foo1,foo2,foo3,flag)
        do k=1,nk
        do j=1,nj
        do i=-1,ni+2
          flag(i)=.false.
        enddo
        do i=0,ni+1
          if(dum(i,j,k).lt.0.0)then
            foo1=max(0.0,dum(i-1,j,k))
            foo2=max(0.0,dum(i+1,j,k))
            if(foo1+foo2.gt.1.0e-30)then
              foo3=max(xh(i)*dum(i,j,k),-(xh(i-1)*foo1+xh(i+1)*foo2))   &
                                        /(xh(i-1)*foo1+xh(i+1)*foo2)
              mass(i-1,j,k)=mass(i-1,j,k)+foo1*foo3
              mass(i  ,j,k)=mass(i  ,j,k)-(foo1+foo2)*foo3
              mass(i+1,j,k)=mass(i+1,j,k)+foo2*foo3
              if(dum(i-1,j,k).gt.1.0e-30) flag(i-1)=.true.
                                          flag(i  )=.true.
              if(dum(i+1,j,k).gt.1.0e-30) flag(i+1)=.true.
            endif
          endif
        enddo
        do i=1,ni
        if(flag(i))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advx(i,j,k)=(dum(i,j,k)-rho0(i,j,k)*s1(i,j,k))*rdt
        endif
        enddo
        enddo
        enddo

      ENDIF

!----------------------------------------------------------------

      if(timestats.ge.1) time_pdef=time_pdef+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefy(rho0,advy,dum,mass,s1,dt)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: rho0
      real, dimension(ib:ie,jb:je,kb:ke) :: advy,dum,mass,s1
      real :: dt

      integer i,j,k
      real foo1,foo2,foo3,rdt
      logical, dimension(-1:nj+2) :: flag








!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=rho0(i,j,k)*s1(i,j,k)+dt*advy(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pdef=time_pdef+mytime()

        call bcs(dum)




        rdt=1.0/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=-1,nj+2
        do i=1,ni
          mass(i,j,k)=0.0
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pdef=time_pdef+mytime()





        if(ibs.eq.1.and.sbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            dum(i,-1,k)=0.0
            dum(i, 0,k)=0.0
            dum(i, 1,k)=0.0
          enddo
          enddo
        endif

        if(ibn.eq.1.and.nbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            dum(i,nj  ,k)=0.0
            dum(i,nj+1,k)=0.0
            dum(i,nj+2,k)=0.0
          enddo
          enddo
        elseif(ibn.eq.1.and.nbc.eq.3)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            dum(i,nj+1,k)=dum(i,nj,k)
            dum(i,nj+2,k)=dum(i,nj,k)
          enddo
          enddo
        endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k,foo1,foo2,foo3,flag)
        do k=1,nk
        do i=1,ni
        do j=-1,nj+2
          flag(j)=.false.
        enddo
        do j=0,nj+1
          if(dum(i,j,k).lt.0.0)then
            foo1=max(0.0,dum(i,j-1,k))
            foo2=max(0.0,dum(i,j+1,k))
            if(foo1+foo2.gt.1.0e-30)then
              foo3=max(dum(i,j,k),-(foo1+foo2))/(foo1+foo2)
              mass(i,j-1,k)=mass(i,j-1,k)+foo1*foo3
              mass(i,j  ,k)=mass(i,j  ,k)-(foo1+foo2)*foo3
              mass(i,j+1,k)=mass(i,j+1,k)+foo2*foo3
              if(dum(i,j-1,k).gt.1.0e-30) flag(j-1)=.true.
                                          flag(j  )=.true.
              if(dum(i,j+1,k).gt.1.0e-30) flag(j+1)=.true.
            endif
          endif
        enddo
        do j=1,nj
        if(flag(j))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advy(i,j,k)=(dum(i,j,k)-rho0(i,j,k)*s1(i,j,k))*rdt
        endif
        enddo
        enddo
        enddo

!----------------------------------------------------------------

      if(timestats.ge.1) time_pdef=time_pdef+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefz(rho0,advz,dum,mass,s1,dt)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: rho0
      real, dimension(ib:ie,jb:je,kb:ke) :: advz,dum,mass,s1
      real :: dt

      integer i,j,k
      real foo1,foo2,foo3,rdt
      logical, dimension(0:nk+1) :: flag

!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=rho0(i,j,k)*s1(i,j,k)+dt*advz(i,j,k)
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum(i,j, 0)=0.0
          dum(i,j,nk+1)=0.0
        enddo
        enddo

        rdt=1.0/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
          mass(i,j,k)=0.0
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k,foo1,foo2,foo3,flag)
      do j=1,nj
      do i=1,ni
        do k=0,nk+1
          flag(k)=.false.
        enddo
        do k=1,nk
          if(dum(i,j,k).lt.0.0)then
            foo1=max(0.0,dum(i,j,k-1))
            foo2=max(0.0,dum(i,j,k+1))
            if(foo1+foo2.gt.1.0e-30)then
              foo3=max(dum(i,j,k),-(foo1+foo2))/(foo1+foo2)
              mass(i,j,k-1)=mass(i,j,k-1)+foo1*foo3
              mass(i,j,k  )=mass(i,j,k  )-(foo1+foo2)*foo3
              mass(i,j,k+1)=mass(i,j,k+1)+foo2*foo3
              if(dum(i,j,k-1).gt.1.0e-30) flag(k-1)=.true.
                                          flag(k  )=.true.
              if(dum(i,j,k+1).gt.1.0e-30) flag(k+1)=.true.
            endif
          endif
        enddo
        do k=1,nk
        if(flag(k))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advz(i,j,k)=(dum(i,j,k)-rho0(i,j,k)*s1(i,j,k))*rdt
        endif
        enddo
      enddo
      enddo

!----------------------------------------------------------------

      if(timestats.ge.1) time_pdef=time_pdef+mytime()

      return
      end


