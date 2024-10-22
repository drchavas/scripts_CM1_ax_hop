



      subroutine adv5s(bflag,bsq,xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,          &
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
      real, parameter :: tem2 = 1.0/6.0

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
        ! this seems to be faster on most other platforms:
        if(rru(i,j,k).ge.0.)then
          dum(i,j,k)=rru(i,j,k)*( 2.*s(i-3,j,k)-13.*s(i-2,j,k)   &
                +47.*s(i-1,j,k)+27.*s(i,j,k)-3.*s(i+1,j,k) )*tem
        else
          dum(i,j,k)=rru(i,j,k)*( 2.*s(i+2,j,k)-13.*s(i+1,j,k)   &
                +47.*s(i,j,k)+27.*s(i-1,j,k)-3.*s(i-2,j,k) )*tem
        endif
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

      IF(pdscheme.eq.1 .and. pdef.eq.1)THEN
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
        ! this seems to be faster on most other platforms:
        if(rrv(i,j,k).ge.0.)then
          dum(i,j,k)=rrv(i,j,k)*( 2.*s(i,j-3,k)-13.*s(i,j-2,k)   &
                +47.*s(i,j-1,k)+27.*s(i,j,k)-3.*s(i,j+1,k) )*tem
        else
          dum(i,j,k)=rrv(i,j,k)*( 2.*s(i,j+2,k)-13.*s(i,j+1,k)   &
                +47.*s(i,j,k)+27.*s(i,j-1,k)-3.*s(i,j-2,k) )*tem
        endif
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


      subroutine vadv5s(gz,mh,rho0,advz,dum,mass,rrw,s1,s,pdef,dt)
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
      real, parameter :: tem2 = 1.0/6.0
      real, parameter :: tem3 = 1.0/12.0

!----------------------------------------------------------------
! Advection in z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=4,nk-2
      do j=1,nj
      do i=1,ni
        ! this seems to be faster on most other platforms:
        if(rrw(i,j,k).ge.0.)then
          dum(i,j,k)=rrw(i,j,k)*( 2.*s(i,j,k-3)-13.*s(i,j,k-2)      &
                +47.*s(i,j,k-1)+27.*s(i,j,k)-3.*s(i,j,k+1) )*tem1
        else
          dum(i,j,k)=rrw(i,j,k)*( 2.*s(i,j,k+2)-13.*s(i,j,k+1)      &
                +47.*s(i,j,k)+27.*s(i,j,k-1)-3.*s(i,j,k-2) )*tem1
        endif
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=3,(nk-1),(nk-4)
      do j=1,nj
      do i=1,ni
        ! this seems to be faster on most other platforms:
        if(rrw(i,j,k).ge.0.)then
          dum(i,j,k)=rrw(i,j,k)    &
                    *(-s(i,j,k-2)+5.*s(i,j,k-1)+2.*s(i,j,k))*tem2
        else
          dum(i,j,k)=rrw(i,j,k)    &
                    *(-s(i,j,k+1)+5.*s(i,j,k)+2.*s(i,j,k-1))*tem2
        endif
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


      subroutine adv5u(xf,rxf,uf,vh,gz,mh,rho0,rr0,rf0,rrf0,rru0,dum,advx,advy,advz,divx, &
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
        ! this seems to be faster on most other platforms:
        if((rru(i,j,k)+rru(i+1,j,k)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i+1,j,k))*                     &
                 ( 2.*u3d(i-2,j,k)-13.*u3d(i-1,j,k)+47.*u3d(i,j,k)   &
                  +27.*u3d(i+1,j,k)-3.*u3d(i+2,j,k) )*tem
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i+1,j,k))*                       &
                 ( 2.*u3d(i+3,j,k)-13.*u3d(i+2,j,k)+47.*u3d(i+1,j,k)   &
                  +27.*u3d(i,j,k)-3.*u3d(i-1,j,k) )*tem
        endif
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=id1,id2
        ! this seems to be faster on most other platforms:
        if((rru(i,j,k)+rru(i+1,j,k)).ge.0.)then
          dum(i,j,k)=(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))*                     &
                 ( 2.*u3d(i-2,j,k)-13.*u3d(i-1,j,k)+47.*u3d(i,j,k)   &
                  +27.*u3d(i+1,j,k)-3.*u3d(i+2,j,k) )*tem
        else
          dum(i,j,k)=(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))*                       &
                 ( 2.*u3d(i+3,j,k)-13.*u3d(i+2,j,k)+47.*u3d(i+1,j,k)   &
                  +27.*u3d(i,j,k)-3.*u3d(i-1,j,k) )*tem
        endif
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
        ! this seems to be faster on most other platforms:
        if((rrv(i,j,k)+rrv(i-1,j,k)).ge.0.)then
          dum(i,j,k)=(rrv(i,j,k)+rrv(i-1,j,k))*                       &
                 ( 2.*u3d(i,j-3,k)-13.*u3d(i,j-2,k)+47.*u3d(i,j-1,k)   &
                  +27.*u3d(i,j,k)-3.*u3d(i,j+1,k) )*tem
        else
          dum(i,j,k)=(rrv(i,j,k)+rrv(i-1,j,k))*                     &
                 ( 2.*u3d(i,j+2,k)-13.*u3d(i,j+1,k)+47.*u3d(i,j,k)   &
                  +27.*u3d(i,j-1,k)-3.*u3d(i,j-2,k) )*tem
        endif
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


      subroutine vadv5u(gz,mh,dum,advz,u3d,rrw,i1,i2)
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
      real, parameter :: tem2 = 1.0/12.0
      real, parameter :: tem3 = 1.0/24.0

!----------------------------------------------------------------
! Advection in z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=4,nk-2
      do j=1,nj
      do i=i1,i2
        ! this seems to be faster on most other platforms:
        if((rrw(i,j,k)+rrw(i-1,j,k)).ge.0.)then
          dum(i,j,k)=(rrw(i,j,k)+rrw(i-1,j,k))*                       &
                 ( 2.*u3d(i,j,k-3)-13.*u3d(i,j,k-2)+47.*u3d(i,j,k-1)   &
                  +27.*u3d(i,j,k)-3.*u3d(i,j,k+1) )*tem1
        else
          dum(i,j,k)=(rrw(i,j,k)+rrw(i-1,j,k))*                     &
                 ( 2.*u3d(i,j,k+2)-13.*u3d(i,j,k+1)+47.*u3d(i,j,k)   &
                  +27.*u3d(i,j,k-1)-3.*u3d(i,j,k-2) )*tem1
        endif
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=3,(nk-1),(nk-4)
      do j=1,nj
      do i=i1,i2
        ! this seems to be faster on most other platforms:
        if((rrw(i,j,k)+rrw(i-1,j,k)).ge.0.)then
          dum(i,j,k)=(rrw(i,j,k)+rrw(i-1,j,k))*   &
                 (-u3d(i,j,k-2)+5.*u3d(i,j,k-1)+2.*u3d(i,j,k  ))*tem2
        else
          dum(i,j,k)=(rrw(i,j,k)+rrw(i-1,j,k))*   &
                 (-u3d(i,j,k+1)+5.*u3d(i,j,k  )+2.*u3d(i,j,k-1))*tem2
        endif
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


      subroutine adv5v(xh,rxh,uh,xf,vf,gz,mh,rho0,rr0,rf0,rrf0,rrv0,dum,advx,advy,advz,divx, &
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
        ! this seems to be faster on most other platforms:
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*                        &
                 ( 2.*v3d(i-3,j,k)-13.*v3d(i-2,j,k)+47.*v3d(i-1,j,k)    &
                  +27.*v3d(i,j,k)-3.*v3d(i+1,j,k) )*tem
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*                      &
                 ( 2.*v3d(i+2,j,k)-13.*v3d(i+1,j,k)+47.*v3d(i,j,k)    &
                  +27.*v3d(i-1,j,k)-3.*v3d(i-2,j,k) )*tem
        endif
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
        ! this seems to be faster on most other platforms:
        if((rrv(i,j,k)+rrv(i,j+1,k)).ge.0.)then
          dum(i,j,k)=(rrv(i,j,k)+rrv(i,j+1,k))*                      &
                 ( 2.*v3d(i,j-2,k)-13.*v3d(i,j-1,k)+47.*v3d(i,j,k)    &
                  +27.*v3d(i,j+1,k)-3.*v3d(i,j+2,k) )*tem
        else
          dum(i,j,k)=(rrv(i,j,k)+rrv(i,j+1,k))*                        &
                 ( 2.*v3d(i,j+3,k)-13.*v3d(i,j+2,k)+47.*v3d(i,j+1,k)    &
                  +27.*v3d(i,j,k)-3.*v3d(i,j-1,k) )*tem
        endif
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


      subroutine vadv5v(gz,mh,dum,advz,v3d,rrw,j1,j2)
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
      real, parameter :: tem2 = 1.0/12.0
      real, parameter :: tem3 = 1.0/24.0

!----------------------------------------------------------------
! Advection in z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=4,nk-2
      do j=j1,j2
      do i=1,ni
        ! this seems to be faster on most other platforms:
        if((rrw(i,j,k)+rrw(i,j-1,k)).ge.0.)then
          dum(i,j,k)=(rrw(i,j,k)+rrw(i,j-1,k))*                        &
                 ( 2.*v3d(i,j,k-3)-13.*v3d(i,j,k-2)+47.*v3d(i,j,k-1)    &
                  +27.*v3d(i,j,k)-3.*v3d(i,j,k+1) )*tem1
        else
          dum(i,j,k)=(rrw(i,j,k)+rrw(i,j-1,k))*                      &
                 ( 2.*v3d(i,j,k+2)-13.*v3d(i,j,k+1)+47.*v3d(i,j,k)    &
                  +27.*v3d(i,j,k-1)-3.*v3d(i,j,k-2) )*tem1
        endif
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=3,(nk-1),(nk-4)
      do j=j1,j2
      do i=1,ni
        ! this seems to be faster on most other platforms:
        if((rrw(i,j,k)+rrw(i,j-1,k)).ge.0.)then
          dum(i,j,k)=(rrw(i,j,k)+rrw(i,j-1,k))*    &
                 (-v3d(i,j,k-2)+5.*v3d(i,j,k-1)+2.*v3d(i,j,k  ))*tem2
        else
          dum(i,j,k)=(rrw(i,j,k)+rrw(i,j-1,k))*    &
                 (-v3d(i,j,k+1)+5.*v3d(i,j,k  )+2.*v3d(i,j,k-1))*tem2
        endif
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


      subroutine adv5w(xh,rxh,uh,xf,vh,gz,mf,rho0,rr0,rf0,rrf0,dum,advx,advy,advz,divx, &
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
        ! this seems to be faster on most other platforms:
        if((rru(i,j,k)+rru(i,j,k-1)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*                        &
                 ( 2.*w3d(i-3,j,k)-13.*w3d(i-2,j,k)+47.*w3d(i-1,j,k)    &
                  +27.*w3d(i,j,k)-3.*w3d(i+1,j,k) )*tem
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*                      &
                 ( 2.*w3d(i+2,j,k)-13.*w3d(i+1,j,k)+47.*w3d(i,j,k)    &
                  +27.*w3d(i-1,j,k)-3.*w3d(i-2,j,k) )*tem
        endif
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
        ! this seems to be faster on most other platforms:
        if((rrv(i,j,k)+rrv(i,j,k-1)).ge.0.)then
          dum(i,j,k)=(rrv(i,j,k)+rrv(i,j,k-1))*                        &
                 ( 2.*w3d(i,j-3,k)-13.*w3d(i,j-2,k)+47.*w3d(i,j-1,k)    &
                  +27.*w3d(i,j,k)-3.*w3d(i,j+1,k) )*tem
        else
          dum(i,j,k)=(rrv(i,j,k)+rrv(i,j,k-1))*                      &
                 ( 2.*w3d(i,j+2,k)-13.*w3d(i,j+1,k)+47.*w3d(i,j,k)    &
                  +27.*w3d(i,j-1,k)-3.*w3d(i,j-2,k) )*tem
        endif
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


      subroutine vadv5w(gz,mf,dum,advz,rrw,w3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, dimension(itb:ite,jtb:jte) :: gz
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advz
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw,w3d

      integer i,j,k
      real, parameter :: tem1 = 1.0/120.0
      real, parameter :: tem2 = 1.0/12.0
      real, parameter :: tem3 = 1.0/24.0

!----------------------------------------------------------------
! Advection in z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=3,nk-2
      do j=1,nj
      do i=1,ni
        ! this seems to be faster on most other platforms:
        if((rrw(i,j,k)+rrw(i,j,k+1)).ge.0.)then
          dum(i,j,k)=(rrw(i,j,k)+rrw(i,j,k+1))*       &
                 ( 2.*w3d(i,j,k-2)-13.*w3d(i,j,k-1)+47.*w3d(i,j,k)    &
                  +27.*w3d(i,j,k+1)-3.*w3d(i,j,k+2) )*tem1
        else
          dum(i,j,k)=(rrw(i,j,k)+rrw(i,j,k+1))*       &
                 ( 2.*w3d(i,j,k+3)-13.*w3d(i,j,k+2)+47.*w3d(i,j,k+1)  &
                  +27.*w3d(i,j,k)-3.*w3d(i,j,k-1) )*tem1
        endif
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,(nk-1),(nk-3)
      do j=1,nj
      do i=1,ni
        ! this seems to be faster on most other platforms:
        if((rrw(i,j,k)+rrw(i,j,k+1)).ge.0.)then
          dum(i,j,k)=(rrw(i,j,k)+rrw(i,j,k+1))*    &
                 (-w3d(i,j,k-1)+5.*w3d(i,j,k  )+2.*w3d(i,j,k+1))*tem2
        else
          dum(i,j,k)=(rrw(i,j,k)+rrw(i,j,k+1))*    &
                 (-w3d(i,j,k+2)+5.*w3d(i,j,k+1)+2.*w3d(i,j,k  ))*tem2
        endif
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


