



      subroutine turbtke(dt,dodrag,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,    &
                         nm,defsq,defh,tk,lenscl,lenh,grdscl,rgrdscl,   &
                         kmh,kmv,khh,khv,tkea,tketen,t13,t23,ua,va)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real :: dt
      logical, intent(in) :: dodrag,dosfcflx
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, dimension(ib:ie,jb:je,kb:ke) :: th0
      real, dimension(ib:ie,jb:je) :: thflux,qvflux,rth0s
      real, dimension(ib:ie,jb:je,kb:ke) :: nm,defsq,defh,tk,   &
                                            lenscl,lenh,grdscl,rgrdscl
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tketen
      real, dimension(ib:ie,jb:je,kb:ke) :: t13,t23
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va

!----------------------------------------

      integer i,j,k
      real prinv,tem,tem1,tem2,tem3


!------------------------------------------------------------------
!  get grid scale

    IF(tconfig.eq.1)THEN
      ! single length scale:  appropriate if dx,dy are nearly the same as dz

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk+1
      do j=1,nj
      do i=1,ni
        grdscl(i,j,k)=(dx*ruh(i)*dy*rvh(j)*dz*rmf(i,j,k))**0.33333333
        rgrdscl(i,j,k)=1.0/grdscl(i,j,k)
      enddo
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN
      ! two length scales:  one for horizontal, one for vertical

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem)
      do j=1,nj
      do i=1,ni
        tem=sqrt(dx*ruh(i)*dy*rvh(j))
        do k=1,nkt
          lenh(i,j,k)=tem
        enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk+1
      do j=1,nj
      do i=1,ni
        grdscl(i,j,k)=dz*rmf(i,j,k)
        rgrdscl(i,j,k)=1.0/grdscl(i,j,k)
      enddo
      enddo
      enddo

    ENDIF

!------------------------------------------------------------------
!  Get turbulence length scale

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nkt
      do j=1,nj
      do i=1,ni
        tk(i,j,k)=max(tkea(i,j,k),1.0e-6)
        lenscl(i,j,k)=grdscl(i,j,k)
        if(nm(i,j,k).gt.1.0e-6)then
          lenscl(i,j,k)=0.8165*sqrt(tk(i,j,k)/nm(i,j,k))
          lenscl(i,j,k)=min(lenscl(i,j,k),grdscl(i,j,k))
          lenscl(i,j,k)=max(lenscl(i,j,k),1.0e-6*grdscl(i,j,k))
        endif 
      enddo
      enddo
      enddo

!----------------------------------------------
!  Get km, kh

    IF(tconfig.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,prinv)
      do k=1,nkt
      do j=1,nj
      do i=1,ni
        kmh(i,j,k)=0.10*sqrt(tk(i,j,k))*lenscl(i,j,k)
        kmv(i,j,k)=kmh(i,j,k)
        prinv=3.00
        if(nm(i,j,k).gt.1.0e-6)then
          prinv=min(1.0+2.00*lenscl(i,j,k)*rgrdscl(i,j,k),3.00)
        endif
        khh(i,j,k)=kmh(i,j,k)*prinv
        khv(i,j,k)=khh(i,j,k)
      enddo
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,prinv)
      do k=1,nkt
      do j=1,nj
      do i=1,ni
        kmh(i,j,k)=0.10*sqrt(tk(i,j,k))*lenh(i,j,k)
        kmv(i,j,k)=0.10*sqrt(tk(i,j,k))*lenscl(i,j,k)
        prinv=3.00
        if(nm(i,j,k).gt.1.0e-6)then
          prinv=min(1.0+2.00*lenscl(i,j,k)*rgrdscl(i,j,k),3.00)
        endif
        khh(i,j,k)=kmh(i,j,k)*prinv
        khv(i,j,k)=kmv(i,j,k)*prinv
      enddo
      enddo
      enddo

    ENDIF

!------------------------------------------------------------
!  Buoyancy Term

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nkt
      do j=1,nj
      do i=1,ni
        tketen(i,j,k)=tketen(i,j,k)-khv(i,j,k)*nm(i,j,k)
      enddo
      enddo
      enddo

      IF(dosfcflx)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          tketen(i,j,1)=tketen(i,j,1)   &
                       +g*( thflux(i,j)*rth0s(i,j)+repsm1*qvflux(i,j) )
        enddo
        enddo

      ENDIF

!------------------------------------------------------------
! Shear term 

    IF(tconfig.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        tketen(i,j,k)=tketen(i,j,k)+kmv(i,j,k)*(defsq(i,j,k)+defh(i,j,k))
      enddo
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        tketen(i,j,k)=tketen(i,j,k)+kmv(i,j,k)*defsq(i,j,k)   &
                                   +kmh(i,j,k)*defh(i,j,k)
      enddo
      enddo
      enddo

    ENDIF

    IF(bcturbu.eq.3)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        tketen(i,j,1)=tketen(i,j,1)+(              &
                      ( t13(i  ,j,1)*ua(i  ,j,1)     &
                       +t13(i+1,j,1)*ua(i+1,j,1) )   &
                    + ( t23(i,j  ,1)*va(i,j  ,1)     &
                       +t23(i,j+1,1)*va(i,j+1,1) )   &
                                    )*rdz*mf(i,j,1)
      enddo
      enddo

    ENDIF

!------------------------------------------------------------
!  limit for numerical stability:

      tem1 = 0.125*dx*dx/dt
      tem2 = 0.125*dy*dy/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nkt
      do j=1,nj
      do i=1,ni
        kmh(i,j,k) = min( kmh(i,j,k) , tem1*ruh(i)*ruh(i) , tem2*rvh(j)*rvh(j) )
        khh(i,j,k) = min( khh(i,j,k) , tem1*ruh(i)*ruh(i) , tem2*rvh(j)*rvh(j) )
      enddo
      enddo
      enddo

!------------------------------------------------------------
! Set values at boundaries

      if(timestats.ge.1) time_turb=time_turb+mytime()
      call bct(kmh)
      call bct(kmv)
      call bct(khh)
      call bct(khv)

!------------------------------------------------------------
!  Dissipation Term

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nkt
      do j=1,nj
      do i=1,ni
        tketen(i,j,k)=tketen(i,j,k)                           &
                     -(0.191+0.796*lenscl(i,j,k)*rgrdscl(i,j,k))*   &
                      tk(i,j,k)*sqrt(tk(i,j,k))/lenscl(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_turb=time_turb+mytime()

!--------------------------------------------------------------
!  finished
      
      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbsmag(dt,dodrag,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,  &
                          nm,defsq,defh,lenscl,grdscl,lenh,            &
                          kmh,kmv,khh,khv,t13,t23,ua,va)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real :: dt
      logical, intent(in) :: dodrag,dosfcflx
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, dimension(ib:ie,jb:je,kb:ke) :: th0
      real, dimension(ib:ie,jb:je) :: thflux,qvflux,rth0s
      real, dimension(ib:ie,jb:je,kb:ke) :: nm,defsq,defh,lenscl,grdscl,lenh
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, dimension(ib:ie,jb:je,kb:ke) :: t13,t23
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va

      integer i,j,k
      real :: tem,tem1,tem2,tem3


      real, parameter :: cs      = 0.18
      real, parameter :: prandtl = 1.0/3.00
      real, parameter :: prinv   = 1.0/prandtl

!!!      real, parameter :: c_m = 0.0856
!!!      real, parameter :: c_h = 0.214
!!!      real, parameter :: c_l = 0.816
!!!      real, parameter :: ce1 = 0.191
!!!      real, parameter :: ce2 = 0.654
!!!      real, parameter :: ric = 0.23

!-----------------------------------------------------------------------

    IF(tconfig.eq.1)THEN
      ! single length scale:  appropriate if dx,dy are nearly the same as dz

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk+1
      do j=1,nj
      do i=1,ni
        grdscl(i,j,k)=(dx*ruh(i)*dy*rvh(j)*dz*rmf(i,j,k))**0.33333333
      enddo
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN
      ! two length scales:  one for horizontal, one for vertical

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem)
      do j=1,nj
      do i=1,ni
        tem=sqrt(dx*ruh(i)*dy*rvh(j))
        do k=1,nkt
          lenh(i,j,k)=tem
        enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk+1
      do j=1,nj
      do i=1,ni
        grdscl(i,j,k)=dz*rmf(i,j,k)
      enddo
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------------
!  Interior points:

    IF(tconfig.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nkt-1
      do j=1,nj
      do i=1,ni
        kmh(i,j,k)=((cs*grdscl(i,j,k))**2)     &
                 *sqrt( max(defsq(i,j,k)+defh(i,j,k)-nm(i,j,k)*prinv,0.0) )
        kmv(i,j,k)=kmh(i,j,k)
      enddo
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nkt-1
      do j=1,nj
      do i=1,ni
        kmh(i,j,k)=((cs*lenh(i,j,k))**2)     &
                 *sqrt( max(defh(i,j,k),0.0) )
        kmv(i,j,k)=((cs*grdscl(i,j,k))**2)     &
                 *sqrt( max(defsq(i,j,k)-nm(i,j,k)*prinv,0.0) )
      enddo
      enddo
      enddo

    ENDIF

!------------------------------------------------------------
!  Surface:

    IF(bcturbu.eq.1.or.bcturbu.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,tem1,tem2)
      do j=1,nj
      do i=1,ni
        kmh(i,j,1) = kmh(i,j,2)
        tem1 = 0.0
      IF(dosfcflx)THEN
        tem2 = g*( thflux(i,j)*rth0s(i,j)+repsm1*qvflux(i,j) )
      ELSE
        tem2 = 0.0
      ENDIF
        kmv(i,j,1) = ( ((cs*grdscl(i,j,1))**4)*max(tem1+tem2,0.0) )**0.33333333
      enddo
      enddo

    ELSEIF(bcturbu.eq.3)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,tem1,tem2)
      do j=1,nj
      do i=1,ni
        kmh(i,j,1) = 0.0
        tem1 = ( ( t13(i  ,j,1)*ua(i  ,j,1)     &
                  +t13(i+1,j,1)*ua(i+1,j,1) )   &
               + ( t23(i,j  ,1)*va(i,j  ,1)     &
                  +t23(i,j+1,1)*va(i,j+1,1) )   &
               )*rdz*mf(i,j,1)
      IF(dosfcflx)THEN
        tem2 = g*( thflux(i,j)*rth0s(i,j)+repsm1*qvflux(i,j) )
      ELSE
        tem2 = 0.0
      ENDIF
        kmv(i,j,1) = ( ((cs*grdscl(i,j,1))**4)*max(tem1+tem2,0.0) )**0.33333333
      enddo
      enddo

    ENDIF

!------------------------------------------------------------
!  Top of domain:

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        kmh(i,j,nkt) = 0.0
        kmv(i,j,nkt) = 0.0
      enddo
      enddo

!--------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()
      call bct(kmh)
      call bct(kmv)

!--------------------------------------------------------------

      tem1 = 0.125*dx*dx/dt
      tem2 = 0.125*dy*dy/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nkt+1
      do j=0,nj+1
      do i=0,ni+1
        khh(i,j,k)=kmh(i,j,k)*prinv
        ! limit for numerical stability:
        khh(i,j,k) = min( khh(i,j,k) , tem1*ruh(i)*ruh(i)   &
                                     , tem2*rvh(j)*rvh(j) )
        kmh(i,j,k) = min( kmh(i,j,k) , tem1*ruh(i)*ruh(i)   &
                                     , tem2*rvh(j)*rvh(j) )
        khv(i,j,k)=kmv(i,j,k)*prinv
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=0,nj+1
      do i=0,ni+1
        khh(i,j,0)=kmh(i,j,0)
        khh(i,j,1)=kmh(i,j,1)
        khv(i,j,0)=kmv(i,j,0)
        khv(i,j,1)=kmv(i,j,1)
      enddo
      enddo

!--------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbparam(nstep,zf,dt,dodrag,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,  &
                             nm,defsq,defh,kmh,kmv,khh,khv,t13,t23,ua,va)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      integer, intent(in) :: nstep
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real :: dt
      logical, intent(in) :: dodrag,dosfcflx
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, dimension(ib:ie,jb:je,kb:ke) :: th0
      real, dimension(ib:ie,jb:je) :: thflux,qvflux,rth0s
      real, dimension(ib:ie,jb:je,kb:ke) :: nm,defsq,defh
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, dimension(ib:ie,jb:je,kb:ke) :: t13,t23
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va

      integer i,j,k
      real :: tem1,tem2,tem3


      real, parameter :: prandtl = 1.0

      real, parameter :: prinv   = 1.0/prandtl

!--------------------------------------------------------------
!  Smagorinsky-type scheme for parameterized turbulence:
!--------------------------------------------------------------
!  Interior:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      kmh(i,j,k)=(l_h**2)*sqrt( defh(i,j,k) )
      kmv(i,j,k)=(l_v**2)*sqrt( max(defsq(i,j,k)-nm(i,j,k)*prinv,0.0) )
    enddo
    enddo
    enddo

!--------------------------------------------------------------
!  Surface:

    IF(bcturbu.eq.1.or.bcturbu.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,tem1,tem2)
      do j=1,nj
      do i=1,ni
        kmh(i,j,1) = kmh(i,j,2)
        tem1 = 0.0
      IF(dosfcflx)THEN
        tem2 = g*( thflux(i,j)*rth0s(i,j)+repsm1*qvflux(i,j) )
      ELSE
        tem2 = 0.0
      ENDIF
        kmv(i,j,1) = ( (l_v**4)*max(tem1+tem2,0.0) )**0.33333333
!!!        kmv(i,j,1) = ( (0.0**4)*max(tem1+tem2,0.0) )**0.33333333
      enddo
      enddo

    ELSEIF(bcturbu.eq.3)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,tem1,tem2)
      do j=1,nj
      do i=1,ni
        kmh(i,j,1) = 0.0
        tem1 = ( ( t13(i  ,j,1)*ua(i  ,j,1)     &
                  +t13(i+1,j,1)*ua(i+1,j,1) )   &
               + ( t23(i,j  ,1)*va(i,j  ,1)     &
                  +t23(i,j+1,1)*va(i,j+1,1) )   &
               )*rdz*mf(i,j,1)
      IF(dosfcflx)THEN
        tem2 = g*( thflux(i,j)*rth0s(i,j)+repsm1*qvflux(i,j) )
      ELSE
        tem2 = 0.0
      ENDIF
        kmv(i,j,1) = ( (l_v**4)*max(tem1+tem2,0.0) )**0.33333333
!!!        kmv(i,j,1) = ( (0.0**4)*max(tem1+tem2,0.0) )**0.33333333
      enddo
      enddo

    ENDIF

!------------------------------------------------------------
!  Top of domain:
!  something simple, for now:

!$omp parallel do default(shared)   &
!$omp private(i,j,tem1,tem2)
      do j=1,nj
      do i=1,ni
        kmh(i,j,  1) = 0.0
        kmv(i,j,  1) = 0.0
        kmh(i,j,nkt) = 0.0
        kmv(i,j,nkt) = 0.0
      enddo
      enddo

!--------------------------------------------------------------
! boundary conditions:

      if(timestats.ge.1) time_turb=time_turb+mytime()
      call bct(kmh)
      call bct(kmv)

!--------------------------------------------------------------
!  calculate Kh
!  and also limit horizontal coeffs for numerical stability:

      tem1 = 0.125*dx*dx/dt
      tem2 = 0.125*dy*dy/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=0,nkt+1
      do j=0,nj+1
      do i=0,ni+1
        khh(i,j,k)=kmh(i,j,k)*prinv
        khh(i,j,k) = min( khh(i,j,k) , tem1*ruh(i)*ruh(i) , tem2*rvh(j)*rvh(j) )
        kmh(i,j,k) = min( kmh(i,j,k) , tem1*ruh(i)*ruh(i) , tem2*rvh(j)*rvh(j) )
        khv(i,j,k)=kmv(i,j,k)*prinv
      enddo
      enddo
      enddo

!!$omp parallel do default(shared)   &
!!$omp private(i,j)
!      do j=0,nj+1
!      do i=0,ni+1
!        khh(i,j,0)=kmh(i,j,0)
!        khh(i,j,1)=kmh(i,j,1)
!        khv(i,j,0)=kmv(i,j,0)
!        khv(i,j,1)=kmv(i,j,1)
!      enddo
!      enddo

!--------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine gettau(dodrag,xf,rxf,rho0,rf0,kmh,kmv,t11,t12,t13,t22,t23,t33,ua)
      implicit none
      
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      logical, intent(in) :: dodrag
      real, dimension(ib:ie+1) :: xf,rxf
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rf0
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv
      real, dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
        
      integer i,j,k

!-----------------------------------------------------------------------
! Note:  turb coefficients are now defined on w points

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=0,nj+1
      do i=0,ni+1
        t11(i,j,k)=t11(i,j,k)*(kmh(i,j,k)+kmh(i,j,k+1))*rho0(i,j,k)
        t22(i,j,k)=t22(i,j,k)*(kmh(i,j,k)+kmh(i,j,k+1))*rho0(i,j,k)
        t33(i,j,k)=t33(i,j,k)*(kmv(i,j,k)+kmv(i,j,k+1))*rho0(i,j,k)
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni+1
        t12(i,j,k)=t12(i,j,k)*0.03125                                         &
     *( ( (kmh(i-1,j-1,k  )+kmh(i,j,k  ))+(kmh(i-1,j,k  )+kmh(i,j-1,k  )) )   &
       +( (kmh(i-1,j-1,k+1)+kmh(i,j,k+1))+(kmh(i-1,j,k+1)+kmh(i,j-1,k+1)) ) ) &
           *( (rho0(i-1,j-1,k)+rho0(i,j,k))+(rho0(i-1,j,k)+rho0(i,j-1,k)) )
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=0,ni+1
        t11(i,j,k)=t11(i,j,k)*(kmh(i,j,k)+kmh(i,j,k+1))*rho0(i,j,k)
        t33(i,j,k)=t33(i,j,k)*(kmv(i,j,k)+kmv(i,j,k+1))*rho0(i,j,k)
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=2,ni+1
        t22(i,j,k)=2.0*rho0(1,1,k)     &
                  *0.25*(kmh(i-1,j,k)+kmh(i,j,k)+kmh(i-1,j,k+1)+kmh(i,j,k+1))  &
                  *ua(i,j,k)*rxf(i)
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=2,ni+1
        t12(i,j,k)=t12(i,j,k)*rho0(1,1,k)   &
                  *0.25*(kmh(i,j,k+1)+kmh(i,j,k)+kmh(i-1,j,k+1)+kmh(i-1,j,k))
      enddo
      enddo
      enddo

    ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj+1
      do i=1,ni+1
        t13(i,j,k)=t13(i,j,k)*0.25        &
           *( kmv(i-1,j,k)+kmv(i,j,k) )   &
           *( rf0(i-1,j,k)+rf0(i,j,k) )
        t23(i,j,k)=t23(i,j,k)*0.25        &
           *( kmv(i,j-1,k)+kmv(i,j,k) )   &
           *( rf0(i,j-1,k)+rf0(i,j,k) )
      enddo
      enddo
      enddo

!--------------------------------------------------------------
!  lateral boundary conditions for axisymmetric simulations

    IF(axisymm.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(k)
      do k=1,nk
        t22(1,1,k)=0.0
      enddo

!$omp parallel do default(shared)   &
!$omp private(k)
      do k=1,nk
        t12(1,1,k)=0.0
      enddo

!$omp parallel do default(shared)   &
!$omp private(k)
      do k=1,nk+1
        t13(1,1,k)=0.0
      enddo

    ENDIF

!--------------------------------------------------------------
!  lower boundary conditions

    IF(dodrag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1
        t13(i,j,1)=t13(i,j,1)*0.5*(rf0(i-1,j,1)+rf0(i,j,1))
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        t23(i,j,1)=t23(i,j,1)*0.5*(rf0(i,j-1,1)+rf0(i,j,1))
      enddo
      enddo

    ELSE

      IF(bcturbu.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni+1
          t13(i,j,1)=0.0
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni
          t23(i,j,1)=0.0
        enddo
        enddo

      ELSEIF(bcturbu.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni+1
          t13(i,j,1)=t13(i,j,2)
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni
          t23(i,j,1)=t23(i,j,2)
        enddo
        enddo

      ENDIF

    ENDIF

!--------------------------------------------------------------
!  upper boundary conditions

      IF(bcturbu.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni+1
          t13(i,j,nk+1)=0.0
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni
          t23(i,j,nk+1)=0.0
        enddo
        enddo

      ELSEIF(bcturbu.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni+1
          t13(i,j,nk+1)=t13(i,j,nk)
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni
          t23(i,j,nk+1)=t23(i,j,nk)
        enddo
        enddo

      ENDIF

!--------------------------------------------------------------
!  finished

      if(timestats.ge.1) time_turb=time_turb+mytime()
 
      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine calcdef(dodrag,xh,rxh,uh,xf,rxf,uf,vh,vf,mh,mf,defsq,defh,   &
                         dum3,dum4,ua,va,wa,t11,t12,t13,t22,t23,t33,gx,gy)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      logical, intent(in) :: dodrag
      real, dimension(ib:ie) :: xh,rxh,uh
      real, dimension(ib:ie+1) :: xf,rxf,uf
      real, dimension(jb:je) :: vh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: defsq,defh,dum3,dum4
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, dimension(itb:ite+1,jtb:jte,ktb:kte) :: gx
      real, dimension(itb:ite,jtb:jte+1,ktb:kte) :: gy
        
      integer i,j,k
      real :: tmp11,tmp22,tmp33,tmp12,tmp13,tmp23

!----------------------------------------------------------------------
!
!  Reference:  Mason, 1989, JAS, p. 1497
!
!----------------------------------------------------------------------
!  First, calculations assuming no terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=0,ni+1 
        t11(i,j,k)=(ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)
      enddo
      enddo
      enddo

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1 
      do i=1,ni+1
        t12(i,j,k)=(ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)   &
                  +(va(i,j,k)-va(i-1,j,k))*rdx*uf(i)
      enddo
      enddo
      enddo       

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=2,ni+1
        t12(i,j,k)=xf(i)*(va(i,j,k)*rxh(i)-va(i-1,j,k)*rxh(i-1))*rdx*uf(i)
      enddo
      enddo
      enddo

    ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni+1
        t13(i,j,k)=(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)   &
                  +(ua(i,j,k)-ua(i,j,k-1))*rdz*0.5*(mf(i-1,j,k)+mf(i,j,k))
      enddo
      enddo
      enddo

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=0,nj+1
      do i=1,ni
        t22(i,j,k)=(va(i,j+1,k)-va(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=2,ni+1
        t22(i,j,k)=ua(i,j,k)*rxf(i)
      enddo
      enddo
      enddo

    ENDIF

    IF(axisymm.eq.0)THEN
      
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj+1   
      do i=1,ni
        t23(i,j,k)=(wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j)   &
                  +(va(i,j,k)-va(i,j,k-1))*rdz*0.5*(mf(i,j-1,k)+mf(i,j,k))
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj   
      do i=1,ni
        t23(i,j,k)=(va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)
      enddo
      enddo
      enddo

    ENDIF
      
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        t33(i,j,k)=(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(i,j,k)
      enddo
      enddo
      enddo

!------------------------------------------------------------------
!  lateral boundary conditions for axisymmetric simulations

    IF(axisymm.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(k)
      do k=1,nk
        t22(1,1,k)=0.0
      enddo

!$omp parallel do default(shared)   &
!$omp private(k)
      do k=1,nk
        t12(1,1,k)=0.0
      enddo

!$omp parallel do default(shared)   &
!$omp private(k)
      do k=1,nk+1
        t13(1,1,k)=0.0
      enddo

    ENDIF

!------------------------------------------------------------------
!  lower boundary conditions

  IF(.not.dodrag)THEN

    IF(bcturbu.eq.1)THEN
      
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1 
        t13(i,j,1)=0.0
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        t23(i,j,1)=0.0
      enddo
      enddo

    ELSEIF(bcturbu.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1 
        t13(i,j,1)=t13(i,j,2)
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        t23(i,j,1)=t23(i,j,2)
      enddo
      enddo

    ELSEIF(bcturbu.eq.3)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1
        t13(i,j,1   )= 2.0*ua(i,j,1 )*rdz*0.5*(mf(i-1,j,1   )+mf(i,j,1   ))
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1   
      do i=1,ni
        t23(i,j,1   )= 2.0*va(i,j,1 )*rdz*0.5*(mf(i,j-1,1   )+mf(i,j,1   ))
      enddo
      enddo

    ENDIF

  ENDIF

!------------------------------------------------------------------
!  upper boundary conditions

    IF(bcturbu.eq.1)THEN
      
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1 
        t13(i,j,nk+1)=0.0
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        t23(i,j,nk+1)=0.0
      enddo
      enddo

    ELSEIF(bcturbu.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1 
        t13(i,j,nk+1)=t13(i,j,nk)
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        t23(i,j,nk+1)=t23(i,j,nk)
      enddo
      enddo

    ELSEIF(bcturbu.eq.3)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1
        t13(i,j,nk+1)=-2.0*ua(i,j,nk)*rdz*0.5*(mf(i-1,j,nk+1)+mf(i,j,nk+1))
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1   
      do i=1,ni
        t23(i,j,nk+1)=-2.0*va(i,j,nk)*rdz*0.5*(mf(i,j-1,nk+1)+mf(i,j,nk+1))
      enddo
      enddo

    ENDIF


!------------------------------------------------------------------
!  now, add the "correction" terms for terrain

      IF(terrain_flag)THEN

!-------- t11 --------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk-1
        do j=1,nj
        do i=0,ni+2
          dum3(i,j,k)=gx(i,j,k)*(ua(i,j,k+1)-ua(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=0,ni+2
          dum3(i,j, 1)=gx(i,j,1 )*(ua(i,j,2 )-ua(i,j,1   ))*rdz
          dum3(i,j,nk)=gx(i,j,nk)*(ua(i,j,nk)-ua(i,j,nk-1))*rdz
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=0,ni+1
          t11(i,j,k)=t11(i,j,k)+0.5*(dum3(i,j,k)+dum3(i+1,j,k))
        enddo
        enddo
        enddo

!-------- t12 --------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk-1
        do j=0,nj+1
        do i=1,ni+1
          dum3(i,j,k)=(ua(i,j,k+1)-ua(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=0,nj+1
        do i=1,ni+1
          dum3(i,j,1 )=(ua(i,j,2 )-ua(i,j,1   ))*rdz
          dum3(i,j,nk)=(ua(i,j,nk)-ua(i,j,nk-1))*rdz
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=1,ni+1
          dum3(i,j,k)=dum3(i,j,k)*0.25*( gy(i-1,j  ,k)+gy(i,j  ,k)    &
                                        +gy(i-1,j+1,k)+gy(i,j+1,k) )
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk-1
        do j=1,nj+1
        do i=0,ni+1
          dum4(i,j,k)=(va(i,j,k+1)-va(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=0,ni+1
          dum4(i,j,1 )=(va(i,j,2 )-va(i,j,1   ))*rdz
          dum4(i,j,nk)=(va(i,j,nk)-va(i,j,nk-1))*rdz
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=0,ni+1
          dum4(i,j,k)=dum4(i,j,k)*0.25*( gx(i  ,j-1,k)+gx(i  ,j,k)    &
                                        +gx(i+1,j-1,k)+gx(i+1,j,k) )
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          t12(i,j,k)=t12(i,j,k)+0.5*( (dum3(i,j-1,k)+dum3(i,j,k))  &
                                     +(dum4(i-1,j,k)+dum4(i,j,k)) )
        enddo
        enddo
        enddo

!-------- t13 --------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=0,ni+1
          dum3(i,j,k)=0.25*( (gx(i,j,k-1)+gx(i+1,j,k-1))     &
                            +(gx(i,j,k  )+gx(i+1,j,k  )) )   &
                     *(wa(i,j,k+1)-wa(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni+1
          t13(i,j,k)=t13(i,j,k)+0.5*(dum3(i,j,k)+dum3(i-1,j,k))
        enddo
        enddo
        enddo

!-------- t22 --------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk-1
        do j=0,nj+2
        do i=1,ni
          dum3(i,j,k)=gy(i,j,k)*(va(i,j,k+1)-va(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=0,nj+2
        do i=1,ni
          dum3(i,j, 1)=gy(i,j,1 )*(va(i,j,2 )-va(i,j,1   ))*rdz
          dum3(i,j,nk)=gy(i,j,nk)*(va(i,j,nk)-va(i,j,nk-1))*rdz
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=1,ni
          t22(i,j,k)=t22(i,j,k)+0.5*(dum3(i,j,k)+dum3(i,j+1,k))
        enddo
        enddo
        enddo

!-------- t23 --------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=0,nj+1
        do i=1,ni
          dum3(i,j,k)=0.25*( (gy(i,j,k-1)+gy(i,j+1,k-1))     &
                            +(gy(i,j,k  )+gy(i,j+1,k  )) )   &
                     *(wa(i,j,k+1)-wa(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj+1
        do i=1,ni
          t23(i,j,k)=t23(i,j,k)+0.5*(dum3(i,j,k)+dum3(i,j-1,k))
        enddo
        enddo
        enddo

      ENDIF

!  end of terrain calculations
!----------------------------------------------------------------------
!  if l_h or l_v is zero, set appropriate terms to zero:
!    (just to be sure)

    IF( iturb.eq.3 .and. l_h.lt.tsmall )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        t11(i,j,k) = 0.0
        t22(i,j,k) = 0.0
        t12(i,j,k) = 0.0
      enddo
      enddo
      enddo
    ENDIF

    IF( iturb.eq.3 .and. l_v.lt.tsmall )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        t13(i,j,k) = 0.0
        t23(i,j,k) = 0.0
        t33(i,j,k) = 0.0
      enddo
      enddo
      enddo
    ENDIF

!----------------------------------------------------------------------
!  calculate D term:

    IF(axisymm.eq.0)THEN
      ! Cartesian domain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23)
      do k=2,nk
      do j=1,nj
      do i=1,ni

        tmp11=0.5*( t11(i,j,k-1)**2 + t11(i,j,k)**2 )
        tmp22=0.5*( t22(i,j,k-1)**2 + t22(i,j,k)**2 )
        tmp33=0.5*( t33(i,j,k-1)**2 + t33(i,j,k)**2 )

        tmp12=0.125*( ( ( t12(i,j  ,k-1)**2 + t12(i+1,j+1,k-1)**2 )     &
                      + ( t12(i,j+1,k-1)**2 + t12(i+1,j  ,k-1)**2 ) )   &
                     +( ( t12(i,j  ,k  )**2 + t12(i+1,j+1,k  )**2 )     &
                      + ( t12(i,j+1,k  )**2 + t12(i+1,j  ,k  )**2 ) ) )

        tmp13=0.5*( t13(i,j,k)**2 + t13(i+1,j,k)**2 )

        tmp23=0.5*( t23(i,j,k)**2 + t23(i,j+1,k)**2 )

        defsq(i,j,k)= 2.0*( tmp33 ) + ( tmp13 + tmp23 )

        defh(i,j,k) = 2.0*( tmp11 + tmp22 ) + tmp12

      enddo
      enddo
      enddo

!--------------------------------------------
    ELSE
      ! axisymmetric domain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23)
      do k=2,nk
      do j=1,nj
      do i=1,ni

        tmp11=0.5*( t11(i,j,k-1)**2 + t11(i,j,k)**2 )
        tmp22=0.25*( (t22(i,j,k-1)**2 + t22(i+1,j,k-1)**2)  &
                    +(t22(i,j,k  )**2 + t22(i+1,j,k  )**2) )
        tmp33=0.5*( t33(i,j,k-1)**2 + t33(i,j,k)**2 )

        tmp12=0.25*(  ( t12(i,j  ,k-1)**2 + t12(i+1,j  ,k-1)**2 )     &
                    + ( t12(i,j  ,k  )**2 + t12(i+1,j  ,k  )**2 ) )

        tmp13=0.5*( t13(i,j,k)**2 + t13(i+1,j,k)**2 )

        tmp23=      t23(i,j,k)**2

        defsq(i,j,k)= 2.0*( tmp33 ) + tmp13 + tmp23

        defh(i,j,k) = 2.0*( tmp11 + tmp22 ) + tmp12

      enddo
      enddo
      enddo

    ENDIF  ! endif for axisymm

!--------------------------------------------------------------
!  finished

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine calcnm(mf,pi0,thv0,th0,cloudvar,nm,t,qt,thv,cloud,   &
                        prs,pp,th,qa)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'goddard.incl'

      logical, dimension(maxq) :: cloudvar
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,th0
      real, dimension(ib:ie,jb:je,kb:ke) :: nm,t,qt,thv,cloud,prs
      real, dimension(ib:ie,jb:je,kb:ke) :: pp,th
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa

      integer i,j,k,n
      real pavg,tavg,qtavg,qvs,lhv,cpml,gamma,qiavg,qsavg,qgavg,drdt
      real qlavg,qvl,qvi,fliq,fice
      real rslf,rsif

!----------------------------------------------------------------------
!  Dry nm

    IF(imoist.eq.0)then

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        nm(i,j,k)=alog( (th0(i,j,k)+th(i,j,k))/(th0(i,j,k-1)+th(i,j,k-1)) ) &
                    *g*rdz*mf(i,j,k)
      enddo
      enddo
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        nm(i,j,   1)=0.0
        nm(i,j,nk+1)=0.0
      enddo
      enddo

!-----------------------------------------------------------------------
!  Moist nm

    ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        t(i,j,k)=(th0(i,j,k)+th(i,j,k))*(pi0(i,j,k)+pp(i,j,k))
      enddo
      enddo
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        qt(i,j,k)=0.0
      enddo
      enddo
      enddo

      DO n=1,numq
        IF( (n.eq.nqv) .or.                                 &
            (n.ge.nql1.and.n.le.nql2) .or.                  &
            (n.ge.nqs1.and.n.le.nqs2.and.iice.eq.1) )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qt(i,j,k)=qt(i,j,k)+qa(i,j,k,n)
          enddo
          enddo
          enddo
        ENDIF
      ENDDO

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        thv(i,j,k)=(th0(i,j,k)+th(i,j,k))*(1.0+reps*qa(i,j,k,nqv))   &
                                         /(1.0+qt(i,j,k))
      enddo
      enddo
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        nm(i,j,k)=g*alog(thv(i,j,k)/thv(i,j,k-1))*rdz*mf(i,j,k)
      enddo
      enddo
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        nm(i,j,   1)=0.0
        nm(i,j,nk+1)=0.0
      enddo
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        cloud(i,j,k)=0.0
      enddo
      enddo
      enddo
      do n=1,numq
        if(cloudvar(n))then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            cloud(i,j,k)=cloud(i,j,k)+qa(i,j,k,n)
          enddo
          enddo
          enddo
        endif
      enddo

!    IF(iice.eq.0)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,pavg,tavg,qtavg,qvs,lhv,cpml,drdt,gamma)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        IF( (cloud(i,j,k).ge.clwsat) .and. (cloud(i,j,k-1).ge.clwsat) )THEN
          pavg =0.5*(prs(i,j,k-1)+prs(i,j,k))
          tavg =0.5*(  t(i,j,k-1)+  t(i,j,k))
          qtavg=0.5*( qt(i,j,k-1)+ qt(i,j,k))
          qvs=rslf(pavg,tavg)
          lhv=lv1-lv2*tavg
          cpml=cp+cpv*qvs+cpl*(qtavg-qvs)

          drdt=17.67*(273.15-29.65)*qvs/((tavg-29.65)**2)
          gamma=g*(1.0+qtavg)*(1.0+lhv*qvs/(rd*tavg))/(cpml+lhv*drdt)
          nm(i,j,k)=g*( ( alog(t(i,j,k)/t(i,j,k-1))*rdz*mf(i,j,k)      &
                            +gamma/tavg )*(1.0+tavg*drdt/(eps+qvs))   &
                         -alog((1.0+qt(i,j,k))/(1.0+qt(i,j,k-1)))*rdz*mf(i,j,k) )
        ENDIF
      enddo
      enddo
      enddo

!    ELSE
!
!!$omp parallel do default(shared)  &
!!$omp private(i,j,k,pavg,tavg,qtavg,qlavg,qvl,qvi,fliq,fice,qvs,cpml,   &
!!$omp lhv,drdt,gamma,n)
!      do k=2,nk
!      do j=1,nj
!      do i=1,ni
!        IF( (cloud(i,j,k).ge.clwsat) .and. (cloud(i,j,k-1).ge.clwsat) )THEN
!          pavg =0.5*(prs(i,j,k-1)+prs(i,j,k))
!          tavg =0.5*(  t(i,j,k-1)+  t(i,j,k))
!          qtavg=0.5*( qt(i,j,k-1)+ qt(i,j,k))
!          qlavg=0.0
!          do n=nql1,nql2
!            qlavg=qlavg+0.5*( qa(i,j,k-1,n)+ qa(i,j,k,n) )
!          enddo
!          qvl=rslf(pavg,tavg)
!          qvi=rsif(pavg,tavg)
!          fliq=max(min((tavg-t00k)*rt0,1.0),0.0)
!          fice=1.0-fliq
!          qvs=fliq*qvl+fice*qvi
!          cpml=cp+cpv*qvs+cpl*qlavg+cpi*(qtavg-qlavg-qvs)
!          lhv=fliq*(lv1-lv2*tavg)+fice*(ls1-ls2*tavg)
!
!          drdt=fliq*17.67*(273.15-29.65)*qvl/((tavg-29.65)**2)    &
!              +fice*21.8745584*(273.15-7.66)*qvi/((tavg-7.66)**2)
!          if(tavg.gt.t00k.and.tavg.lt.t0k)then
!            drdt=drdt+(qvl-qvi)*rt0
!          endif
!          gamma=g*(1.0+qtavg)*(1.0+lhv*qvs/(rd*tavg))/(cpml+lhv*drdt)
!          nm(i,j,k)=g*( ( alog(t(i,j,k)/t(i,j,k-1))*rdz*mf(i,j,k)      &
!                            +gamma/tavg )*(1.0+tavg*drdt/(eps+qvs))   &
!                         -alog((1.0+qt(i,j,k))/(1.0+qt(i,j,k-1)))*rdz*mf(i,j,k) )
!        ENDIF
!      enddo
!      enddo
!      enddo
!
!    ENDIF    ! endif for iice

    ENDIF    ! endif for imoist

!----------------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbs(iflux,dt,dosfcflx,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho0,rr0,rf0,sflux,   &
                       turbx,turby,turbz,dum,dsdz,s,sten,khh,khv,gx,gy,doimpl)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      integer iflux
      real :: dt
      logical, intent(in) :: dosfcflx
      real, dimension(ib:ie) :: xh,rxh,uh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je) :: vh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0
      real, dimension(ib:ie,jb:je) :: sflux
      real, dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum,dsdz,s,sten
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: khh,khv
      real, dimension(itb:ite+1,jtb:jte,ktb:kte) :: gx
      real, dimension(itb:ite,jtb:jte+1,ktb:kte) :: gy
      logical, dimension(ib:ie,jb:je) :: doimpl

      integer :: i,j,k
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2
      real, dimension(kb:ke) :: dumz
      real, dimension(nk) :: cfa,cfb,cfc,cfd,s2
      real, dimension(nk) :: lgbth,lgbph

!---------------------------------------------------------------
!  calculate and store dsdz

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk-1
        do j=0,nj+1
        do i=0,ni+1
          dsdz(i,j,k)=(s(i,j,k+1)-s(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          dsdz(i,j,1 )=(s(i,j,2 )-s(i,j,1   ))*rdz
          dsdz(i,j,nk)=(s(i,j,nk)-s(i,j,nk-1))*rdz
        enddo
        enddo

      ENDIF

!---------------------------------------------------------------
!  x-direction

      IF(.not.terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          dum(i,j,k)= -0.125*( rho0(i,j,k)+rho0(i-1,j,k) )         &
                            *(  (khh(i,j,k  )+ khh(i-1,j,k  ))     &
                               +(khh(i,j,k+1)+ khh(i-1,j,k+1)) )   &
                            *(    s(i,j,k)-   s(i-1,j,k) )*rdx*uf(i)
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          dum(i,j,k)= -0.125*( rho0(i,j,k)+rho0(i-1,j,k) )         &
                            *(  (khh(i,j,k  )+ khh(i-1,j,k  ))     &
                               +(khh(i,j,k+1)+ khh(i-1,j,k+1)) )   &
                           *( (s(i,j,k)-s(i-1,j,k))*rdx*uf(i)      &
                             +gx(i,j,k)*0.5*(dsdz(i-1,j,k)+dsdz(i,j,k)) )
        enddo
        enddo
        enddo

      ENDIF

    IF(axisymm.eq.0)THEN
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        turbx(i,j,k)=-(dum(i+1,j,k)-dum(i,j,k))*rdx*uh(i)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        turbx(i,j,k)=-(xf(i+1)*dum(i+1,j,k)-xf(i)*dum(i,j,k))*rdx*uh(i)*rxh(i)
      enddo
      enddo
      enddo

    ENDIF

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=-khh(i,j,k)*0.5*(rho0(i,j,k-1)+rho0(i,j,k))*(  &
                      0.5*( (s(i+1,j,k)+s(i+1,j,k-1))               &
                           -(s(i-1,j,k)+s(i-1,j,k-1)) )*rdx2        &
                    +0.25*( (gx(i,j,k-1)+gx(i+1,j,k-1))             &
                           +(gx(i,j,k  )+gx(i+1,j,k  )) )           &
                         *( s(i,j,k)-s(i,j,k-1) )*rdz    )
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum(i,j,1)=0.0
          dum(i,j,nk+1)=0.0
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=turbx(i,j,k)-0.5*(gx(i,j,k)+gx(i+1,j,k))   &
                                       *(dum(i,j,k+1)-dum(i,j,k))*rdz
        enddo
        enddo
        enddo

      ENDIF

!---------------------------------------------------------------
!  y-direction

    IF(axisymm.eq.0)THEN

      IF(.not.terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          dum(i,j,k)= -0.125*( rho0(i,j,k)+rho0(i,j-1,k) )         &
                            *(  (khh(i,j,k  )+ khh(i,j-1,k  ))     &
                               +(khh(i,j,k+1)+ khh(i,j-1,k+1)) )   &
                           *(    s(i,j,k)-   s(i,j-1,k) )*rdy*vf(j)
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          dum(i,j,k)= -0.125*( rho0(i,j,k)+rho0(i,j-1,k) )         &
                            *(  (khh(i,j,k  )+ khh(i,j-1,k  ))     &
                               +(khh(i,j,k+1)+ khh(i,j-1,k+1)) )   &
                           *( (s(i,j,k)-s(i,j-1,k))*rdy*vf(j)      &
                             +gy(i,j,k)*0.5*(dsdz(i,j-1,k)+dsdz(i,j,k)) )
        enddo
        enddo
        enddo

      ENDIF
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        turby(i,j,k)=-(dum(i,j+1,k)-dum(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=-khh(i,j,k)*0.5*(rho0(i,j,k-1)+rho0(i,j,k))*(  &
                      0.5*( (s(i,j+1,k)+s(i,j+1,k-1))               &
                           -(s(i,j-1,k)+s(i,j-1,k-1)) )*rdx2        &
                    +0.25*( (gy(i,j,k-1)+gy(i,j+1,k-1))             &
                           +(gy(i,j,k  )+gy(i,j+1,k  )) )           &
                         *( s(i,j,k)-s(i,j,k-1) )*rdz    )
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum(i,j,1)=0.0
          dum(i,j,nk+1)=0.0
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          turby(i,j,k)=turby(i,j,k)-0.5*(gy(i,j,k)+gy(i,j+1,k))   &
                                       *(dum(i,j,k+1)-dum(i,j,k))*rdz
        enddo
        enddo
        enddo

      ENDIF

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        turby(i,j,k)=0.0
      enddo
      enddo
      enddo

    ENDIF

!---------------------------------------------------------------------
!  z-direction

    IF( iturb.eq.3 .and. l_v.lt.tsmall )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

    ELSE

      rdt = 1.0/dt
      tema = -1.0*dt*vialpha*rdz*rdz
      temb = dt*vibeta*rdz*rdz
      temc = dt*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k,dumz,r1,r2,cfa,cfb,cfc,cfd,lgbth,lgbph,tem,s2)
    DO j=1,nj
    DO i=1,ni

      IF(doimpl(i,j))THEN
        ! implicit calculation:

        do k=1,nk
          r1 = khv(i,j,k  )*mf(i,j,k  )*rf0(i,j,k  )*mh(i,j,k)*rr0(i,j,k)
          r2 = khv(i,j,k+1)*mf(i,j,k+1)*rf0(i,j,k+1)*mh(i,j,k)*rr0(i,j,k)
          cfa(k) = tema*r1
          cfc(k) = tema*r2
          cfd(k) = s(i,j,k)   &
                 + temb*(r1*s(i,j,k-1)-(r1+r2)*s(i,j,k)+r2*s(i,j,k+1) )
        enddo
        IF(bcturbs.eq.1)THEN
          cfa( 1) = 0.0
          cfc(nk) = 0.0
        ELSEIF(bcturbs.eq.2)THEN
          cfa( 1) = 0.0
          cfc( 1) = 0.0
          cfa(nk) = 0.0
          cfc(nk) = 0.0
        ENDIF
        do k=1,nk
          cfb(k) = 1.0 - cfa(k) - cfc(k)
        enddo
        if(iflux.eq.1 .and. dosfcflx)then
          cfd(1) = cfd(1) + temc*sflux(i,j)*rf0(i,j,1)*mh(i,j,1)*rr0(i,j,1)
        endif

        lgbth(1)=-cfc(1)/cfb(1)
        lgbph(1)= cfd(1)/cfb(1)
        do k=2,nk
          tem = 1.0/(cfa(k)*lgbth(k-1)+cfb(k))
          lgbth(k)=-cfc(k)*tem
          lgbph(k)=(cfd(k)-cfa(k)*lgbph(k-1))*tem
        enddo
        s2(nk)=lgbph(nk)
        do k=nk-1,1,-1
          s2(k)=lgbth(k)*s2(k+1)+lgbph(k)
        enddo

        do k=1,nk
          turbz(i,j,k) = rho0(i,j,k)*(s2(k)-s(i,j,k))*rdt
        enddo

      ELSE
        ! explicit calculation:

        do k=2,nk
          dumz(k)=-khv(i,j,k)*(s(i,j,k)-s(i,j,k-1))*rdz*mf(i,j,k)*rf0(i,j,k)
        enddo

        IF(bcturbs.eq.1)THEN
 
          dumz(1)=0.0
          dumz(nk+1)=0.0

        ELSEIF(bcturbs.eq.2)THEN
 
          dumz(1)=dumz(2)
          dumz(nk+1)=dumz(nk)

        ENDIF

        if(iflux.eq.1 .and. dosfcflx)then
          dumz(1)=sflux(i,j)*rf0(i,j,1)
        endif

        do k=1,nk
          turbz(i,j,k)=-(dumz(k+1)-dumz(k))*rdz*mh(i,j,k)
        enddo

      ENDIF

    ENDDO
    ENDDO

  ENDIF

!---------------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=sten(i,j,k)+(turbx(i,j,k)+turby(i,j,k)+turbz(i,j,k))*rr0(i,j,k)
      enddo
      enddo
      enddo

!---------------------------------------------------------------------

      if(timestats.ge.1) time_tmix=time_tmix+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbt(dt,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho0,rr0,rf0,rrf0,   &
                       turbx,turby,turbz,dum,dtdz,t,tten,kmh,kmv,gx,gy,doimpl)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real :: dt
      real, dimension(ib:ie) :: xh,rxh,uh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je) :: vh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum,dtdz
      real, dimension(ibt:iet,jbt:jet,kbt:ket) :: t,tten
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv
      real, dimension(itb:ite+1,jtb:jte,ktb:kte) :: gx
      real, dimension(itb:ite,jtb:jte+1,ktb:kte) :: gy
      logical, dimension(ib:ie,jb:je) :: doimpl

      integer :: i,j,k
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2
      real, dimension(nk+1) :: cfa,cfb,cfc,cfd,cfe,cff,t2
      real, dimension(kb:ke) :: dumz

!---------------------------------------------------------------
!  calculate and store dtdz

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=0,nj+1
        do i=0,ni+1
          dtdz(i,j,k)=(t(i,j,k+1)-t(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

      ENDIF

!---------------------------------------------------------------
!  x-direction

      IF(.not.terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni+1
          dum(i,j,k)= -0.25*( rf0(i,j,k)+rf0(i-1,j,k) )   &
                           *( kmh(i,j,k)+kmh(i-1,j,k) )   &
                           *(   t(i,j,k)-  t(i-1,j,k) )*rdx*uf(i)
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni+1
          dum(i,j,k)= -0.25*( rf0(i,j,k)+rf0(i-1,j,k) )         &
                           *( kmh(i,j,k)+kmh(i-1,j,k) )         &
                           *( (t(i,j,k)-t(i-1,j,k))*rdx*uf(i)   &
                             +0.25*(gx(i,j,k-1)+gx(i,j,k))*(dtdz(i-1,j,k)+dtdz(i,j,k)) )
        enddo
        enddo
        enddo

      ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        turbx(i,j,k)=-(dum(i+1,j,k)-dum(i,j,k))*rdx*uh(i)
      enddo
      enddo
      enddo

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=-0.5*(kmh(i,j,k)+kmh(i,j,k+1))*rho0(i,j,k)*(  &
                      0.5*( (t(i+1,j,k)+t(i+1,j,k+1))              &
                           -(t(i-1,j,k)+t(i-1,j,k+1)) )*rdx2       &
                     +0.5*( gx(i,j,k  )+gx(i+1,j,k  ) )            &
                         *( t(i,j,k+1)-t(i,j,k) )*rdz    )
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=turbx(i,j,k)-0.25*( (gx(i,j,k-1)+gx(i+1,j,k-1))     &
                                          +(gx(i,j,k  )+gx(i+1,j,k  )) )   &
                                        *(dum(i,j,k)-dum(i,j,k-1))*rdz
        enddo
        enddo
        enddo

      ENDIF

!---------------------------------------------------------------
!  y-direction

      IF(.not.terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj+1
        do i=1,ni
          dum(i,j,k)= -0.25*( rf0(i,j,k)+rf0(i,j-1,k) )   &
                           *( kmh(i,j,k)+kmh(i,j-1,k) )   &
                           *(   t(i,j,k)-  t(i,j-1,k) )*rdy*vf(j)
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj+1
        do i=1,ni
          dum(i,j,k)= -0.25*( rf0(i,j,k)+rf0(i,j-1,k) )         &
                           *( kmh(i,j,k)+kmh(i,j-1,k) )         &
                           *( (t(i,j,k)-t(i,j-1,k))*rdy*vf(j)   &
                             +0.25*(gy(i,j,k-1)+gy(i,j,k))*(dtdz(i,j-1,k)+dtdz(i,j,k)) )
        enddo
        enddo
        enddo

      ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        turby(i,j,k)=-(dum(i,j+1,k)-dum(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=-0.5*(kmh(i,j,k)+kmh(i,j,k+1))*rho0(i,j,k)*(  &
                      0.5*( (t(i,j+1,k)+t(i,j+1,k+1))              &
                           -(t(i,j-1,k)+t(i,j-1,k+1)) )*rdy2       &
                     +0.5*( gy(i,j,k  )+gy(i,j+1,k  ) )            &
                         *( t(i,j,k+1)-t(i,j,k) )*rdz    )
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          turby(i,j,k)=turby(i,j,k)-0.25*( (gy(i,j,k-1)+gy(i,j+1,k-1))     &
                                          +(gy(i,j,k  )+gy(i,j+1,k  )) )   &
                                        *(dum(i,j,k)-dum(i,j,k-1))*rdz
        enddo
        enddo
        enddo

      ENDIF

!---------------------------------------------------------------------
!  z-direction

      rdt = 1.0/dt
      tema = dt*vialpha*rdz*rdz
      temb = dt*vibeta*rdz*rdz
      temc = dt*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,cfe,cff,tem,t2,dumz)
    do j=1,nj
    do i=1,ni

    IF(doimpl(i,j))THEN

        do k=2,nk
          r1 = 0.5*(kmv(i,j,k-1)+kmv(i,j,k  ))*mh(i,j,k-1)*rho0(i,j,k-1)*mf(i,j,k)*rrf0(i,j,k)
          r2 = 0.5*(kmv(i,j,k  )+kmv(i,j,k+1))*mh(i,j,k  )*rho0(i,j,k  )*mf(i,j,k)*rrf0(i,j,k)
          cfa(k) = tema*r1
          cfc(k) = tema*r2
          cfb(k) = 1.0 + cfa(k) + cfc(k)
          cfd(k) = t(i,j,k)    &
                 + temb*(r2*t(i,j,k+1)-(r1+r2)*t(i,j,k)+r1*t(i,j,k-1))
        enddo

        cfe(1)=0.0
        cff(1)=t(i,j,1)
        do k=2,nk
          tem = 1.0/(cfb(k)-cfc(k)*cfe(k-1))
          cfe(k)=cfa(k)*tem
          cff(k)=(cfd(k)+cfc(k)*cff(k-1))*tem
        enddo
        t2(nk+1)=0.0
        do k=nk,2,-1
          t2(k)=cfe(k)*t2(k+1)+cff(k)
        enddo
        do k=2,nk
          turbz(i,j,k) = rf0(i,j,k)*(t2(k)-t(i,j,k))*rdt
        enddo

    ELSE

      do k=1,nk
        dumz(k)=-0.5*(kmv(i,j,k)+kmv(i,j,k+1))*(t(i,j,k+1)-t(i,j,k))*rdz*mh(i,j,k)*rho0(i,j,k)
      enddo
      do k=2,nk
        turbz(i,j,k)=-(dumz(k)-dumz(k-1))*rdz*mf(i,j,k)
      enddo

    ENDIF

    enddo
    enddo

!---------------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        tten(i,j,k)=tten(i,j,k)+2.0*(turbx(i,j,k)+turby(i,j,k)+turbz(i,j,k))*rrf0(i,j,k)
      enddo
      enddo
      enddo

!---------------------------------------------------------------------

      if(timestats.ge.1) time_tmix=time_tmix+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbu(dt,dodrag,xh,ruh,xf,rxf,uf,vh,mh,mf,rmf,rho0,rf0,rru0,   &
           turbx,turby,turbz,dum,u,uten,w,t11,t12,t13,t22,kmv,gx,gy,doimpl)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real :: dt
      logical, intent(in) :: dodrag
      real, dimension(ib:ie) :: xh,ruh
      real, dimension(ib:ie+1) :: xf,rxf,uf
      real, dimension(jb:je) :: vh
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rf0,rru0
      real, dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u,uten
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv
      real, dimension(itb:ite+1,jtb:jte,ktb:kte) :: gx
      real, dimension(itb:ite,jtb:jte+1,ktb:kte) :: gy
      logical, dimension(ib:ie,jb:je) :: doimpl

      integer :: i,j,k,ip
      real :: rdt,tema,temb,temc,temd,teme
      real :: tem,r1,r2,dwdx1,dwdx2
      real, dimension(nk) :: cfa,cfb,cfc,cfd,u2
      real, dimension(nk) :: lgbth,lgbph
      logical :: dow1,dow2

!---------------------------------------------------------------
!  x-direction

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        turbx(i,j,k)=(t11(i,j,k)-t11(i-1,j,k))*rdx*uf(i)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(j,k)
      do k=1,nk
      do j=1,nj
        turbx(1,j,k)=0.0
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=2,ni+1
        turbx(i,j,k)=( (xh(i)*t11(i,j,k)-xh(i-1)*t11(i-1,j,k))*rdx*uf(i)  &
                      -t22(i,j,k) )*rxf(i)
      enddo
      enddo
      enddo

    ENDIF

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk-1
        do j=1,nj
        do i=1,ni+1
          dum(i,j,k)=( (t11(i-1,j,k+1)+t11(i,j,k+1))    &
                      -(t11(i-1,j,k-1)+t11(i,j,k-1)) )*rdz4
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni+1
          dum(i,j,1 )=( (t11(i-1,j,2   )+t11(i,j,2   ))    &
                       -(t11(i-1,j,1   )+t11(i,j,1   )) )*rdz2
          dum(i,j,nk)=( (t11(i-1,j,nk  )+t11(i,j,nk  ))    &
                       -(t11(i-1,j,nk-1)+t11(i,j,nk-1)) )*rdz2
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          turbx(i,j,k)=turbx(i,j,k)+gx(i,j,k)*dum(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

!---------------------------------------------------------------
!  y-direction

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        turby(i,j,k)=(t12(i,j+1,k)-t12(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk-1
        do j=1,nj+1
        do i=1,ni+1
          dum(i,j,k)=0.5*(gy(i-1,j,k)+gy(i,j,k))    &
                        *(t12(i,j,k+1)-t12(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni+1
          dum(i,j,1 )=0.5*(gy(i-1,j,1 )+gy(i,j,1 ))    &
                         *(t12(i,j,2 )-t12(i,j,1))*rdz
          dum(i,j,nk)=0.5*(gy(i-1,j,nk)+gy(i,j,nk  ))    &
                         *(t12(i,j,nk)-t12(i,j,nk-1))*rdz
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          turby(i,j,k)=turby(i,j,k)+0.5*(dum(i,j,k)+dum(i,j+1,k))
        enddo
        enddo
        enddo

      ENDIF

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        turby(i,j,k)=0.0
      enddo
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------
!  z-direction

    IF( iturb.eq.3 .and. l_v.lt.tsmall )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

    ELSE

      rdt = 0.5/dt
      tema = -0.0625*dt*vialpha*rdz*rdz
      temb =  0.0625*dt*vibeta*rdz*rdz
      temd =  0.25*dt*rdz
      temc =  0.5*dt*rdz
      teme = 0.25*dx*dx/dt

      ip = 0
      if( axisymm.eq.1 ) ip = 1

      dow1 = .false.
      dow2 = .false.

      if( iturb.eq.1 .or. iturb.eq.2 ) dow1 = .true.
      if( iturb.eq.3 .and. l_v.gt.tsmall ) dow1 = .true.
      if( dow1 .and. terrain_flag ) dow2 = .true.

      IF( dow2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=0,ni+1
          dum(i,j,k)=0.25*( (gx(i,j,k-1)+gx(i+1,j,k-1))     &
                           +(gx(i,j,k  )+gx(i+1,j,k  )) )   &
                     *(w(i,j,k+1)-w(i,j,k-1))*rdz2
        enddo
        enddo
        enddo
      ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,lgbth,lgbph,tem,u2,dwdx1,dwdx2)
    DO j=1,nj
    DO i=1+ip,ni+1

    IF(doimpl(i-1,j).or.doimpl(i,j))THEN

        do k=1,nk
          tem = (mh(i-1,j,k)+mh(i,j,k))*rru0(i,j,k)
          r1 = (kmv(i-1,j,k  )+kmv(i,j,k  ))*(mf(i-1,j,k  )+mf(i,j,k  ))   &
              *(rf0(i-1,j,k  )+rf0(i,j,k  ))*tem
          r2 = (kmv(i-1,j,k+1)+kmv(i,j,k+1))*(mf(i-1,j,k+1)+mf(i,j,k+1))   &
              *(rf0(i-1,j,k+1)+rf0(i,j,k+1))*tem
          cfa(k) = tema*r1
          cfc(k) = tema*r2
          cfd(k) = u(i,j,k)    &
                 + temb*( r2*u(i,j,k+1)-(r1+r2)*u(i,j,k)+r1*u(i,j,k-1) )
          if( dow1 )then
            dwdx2 = (w(i,j,k+1)-w(i-1,j,k+1))*rdx*uf(i)
            dwdx1 = (w(i,j,k  )-w(i-1,j,k  ))*rdx*uf(i)
          endif
          if( dow2 )then
            if(k.ne.nk) dwdx2 = dwdx2 + 0.5*(dum(i,j,k+1)+dum(i-1,j,k+1))
            if(k.ne. 1) dwdx1 = dwdx1 + 0.5*(dum(i,j,k  )+dum(i-1,j,k  ))
          endif
          if( dow1 )then
            r1 = min( 0.5*(kmv(i-1,j,k  )+kmv(i,j,k  )) , teme*ruh(i)*ruh(i) ) &
                         *(rf0(i-1,j,k  )+rf0(i,j,k  ))*tem
            r2 = min( 0.5*(kmv(i-1,j,k+1)+kmv(i,j,k+1)) , teme*ruh(i)*ruh(i) ) &
                         *(rf0(i-1,j,k+1)+rf0(i,j,k+1))*tem
            cfd(k) = cfd(k) +    &
                   temd*( r2*dwdx2   &
                        - r1*dwdx1 )
          endif
        enddo
        IF(bcturbu.eq.1)THEN
          cfa( 1) = 0.0
          cfc(nk) = 0.0
        ELSEIF(bcturbu.eq.2)THEN
          cfa( 1) = 0.0
          cfc( 1) = 0.0
          cfa(nk) = 0.0
          cfc(nk) = 0.0
        ELSEIF(bcturbu.eq.3)THEN
          cfa( 1) = 0.0
          cfc(nk) = 0.0
        ENDIF
        do k=1,nk
          cfb(k) = 1.0 - cfa(k) - cfc(k)
        enddo
        if(dodrag)then
          tem = temc*t13(i,j,1)*(mh(i-1,j,1)+mh(i,j,1))*rru0(i,j,1)
!!!          cfb(1) = cfb(1) + vialpha*tem/(1.0e-20+u(i,j,1))
!!!          cfd(1) = cfd(1) - vibeta*tem
          ! set back to simpler version:  now usable if pertflx = 1
          cfd(1) = cfd(1) - tem
        endif

        lgbth(1)=-cfc(1)/cfb(1)
        lgbph(1)= cfd(1)/cfb(1)
        do k=2,nk
          tem = 1.0/(cfa(k)*lgbth(k-1)+cfb(k))
          lgbth(k)=-cfc(k)*tem
          lgbph(k)=(cfd(k)-cfa(k)*lgbph(k-1))*tem
        enddo
        u2(nk)=lgbph(nk)
        do k=nk-1,1,-1
          u2(k)=lgbth(k)*u2(k+1)+lgbph(k)
        enddo

        do k=1,nk
          turbz(i,j,k) = (rho0(i-1,j,k)+rho0(i,j,k))*(u2(k)-u(i,j,k))*rdt
        enddo

    ELSE

      do k=1,nk
        turbz(i,j,k)=(t13(i,j,k+1)-t13(i,j,k))*rdz*0.5*(mh(i-1,j,k)+mh(i,j,k))
      enddo

    ENDIF

    ENDDO
    ENDDO

      IF(axisymm.eq.1)THEN
        DO k=1,nk
          turbz(1,1,k) = 0.0
        ENDDO
      ENDIF

  ENDIF

!-----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        uten(i,j,k)=uten(i,j,k)+(turbx(i,j,k)+turby(i,j,k)+turbz(i,j,k))*rru0(i,j,k)
      enddo
      enddo
      enddo

!-------------------------------------------------------------------
!  All done

      if(timestats.ge.1) time_tmix=time_tmix+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbv(dt,dodrag,xh,rxh,uh,xf,rvh,vf,mh,mf,rho0,rf0,rrv0,   &
            turbx,turby,turbz,dum,v,vten,w,t12,t22,t23,kmv,gx,gy,doimpl)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real :: dt
      logical, intent(in) :: dodrag
      real, dimension(ib:ie) :: xh,rxh,uh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je) :: rvh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rf0,rrv0
      real, dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v,vten
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, dimension(ib:ie,jb:je,kb:ke) :: t12,t22,t23
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv
      real, dimension(itb:ite+1,jtb:jte,ktb:kte) :: gx
      real, dimension(itb:ite,jtb:jte+1,ktb:kte) :: gy
      logical, dimension(ib:ie,jb:je) :: doimpl
 
      integer :: i,j,k,ip
      real :: rdt,tema,temb,temc,temd,teme
      real :: tem,r1,r2,dwdy1,dwdy2
      real, dimension(nk) :: cfa,cfb,cfc,cfd,v2
      real, dimension(nk) :: lgbth,lgbph
      logical :: dow1,dow2

!---------------------------------------------------------------
!  x-direction

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        turbx(i,j,k)=(t12(i+1,j,k)-t12(i,j,k))*rdx*uh(i)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        turbx(i,j,k)=(xf(i+1)*xf(i+1)*t12(i+1,j,k)-xf(i)*xf(i)*t12(i,j,k))*rdx*uh(i)*rxh(i)*rxh(i)
      enddo
      enddo
      enddo

    ENDIF

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk-1
        do j=1,nj+1
        do i=1,ni+1
          dum(i,j,k)=0.5*(gx(i,j-1,k)+gx(i,j,k))    &
                        *(t12(i,j,k+1)-t12(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni+1
          dum(i,j,1 )=0.5*(gx(i,j-1,1 )+gx(i,j,1 ))    &
                         *(t12(i,j,2 )-t12(i,j,1   ))*rdz
          dum(i,j,nk)=0.5*(gx(i,j-1,nk)+gx(i,j,nk))    &
                         *(t12(i,j,nk)-t12(i,j,nk-1))*rdz
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          turbx(i,j,k)=turbx(i,j,k)+0.5*(dum(i,j,k)+dum(i+1,j,k))
        enddo
        enddo
        enddo

      ENDIF

!---------------------------------------------------------------
!  y-direction

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        turby(i,j,k)=(t22(i,j,k)-t22(i,j-1,k))*rdy*vf(j)
      enddo
      enddo
      enddo

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk-1
        do j=1,nj+1
        do i=1,ni
          dum(i,j,k)=( (t22(i,j-1,k+1)+t22(i,j,k+1))    &
                      -(t22(i,j-1,k-1)+t22(i,j,k-1)) )*rdz4
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni
          dum(i,j,1 )=( (t22(i,j-1,2   )+t22(i,j,2   ))    &
                       -(t22(i,j-1,1   )+t22(i,j,1   )) )*rdz2
          dum(i,j,nk)=( (t22(i,j-1,nk  )+t22(i,j,nk  ))    &
                       -(t22(i,j-1,nk-1)+t22(i,j,nk-1)) )*rdz2
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          turby(i,j,k)=turby(i,j,k)+gy(i,j,k)*dum(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        turby(i,j,k)=0.0
      enddo
      enddo
      enddo

    ENDIF
 
!-----------------------------------------------------------------
!  z-direction

    IF( iturb.eq.3 .and. l_v.lt.tsmall )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

    ELSE

      rdt = 0.5/dt
      tema = -0.0625*dt*vialpha*rdz*rdz
      temb =  0.0625*dt*vibeta*rdz*rdz
      temd =  0.25*dt*rdz
      temc =  0.5*dt*rdz
      teme = 0.25*dy*dy/dt

      ip = 1
      if( axisymm.eq.1 ) ip = 0

      dow1 = .false.
      dow2 = .false.

      if( iturb.eq.1 .or. iturb.eq.2 ) dow1 = .true.
      if( iturb.eq.3 .and. l_v.gt.tsmall ) dow1 = .true.
      if( axisymm.eq.1 ) dow1 = .false.
      if( dow1 .and. terrain_flag ) dow2 = .true.

      IF( dow2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=0,nj+1
        do i=1,ni
          dum(i,j,k)=0.25*( (gy(i,j,k-1)+gy(i,j+1,k-1))     &
                           +(gy(i,j,k  )+gy(i,j+1,k  )) )   &
                     *(w(i,j,k+1)-w(i,j,k-1))*rdz2
        enddo
        enddo
        enddo
      ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,lgbth,lgbph,tem,v2,dwdy1,dwdy2)
    do j=1,nj+ip
    do i=1,ni

    IF(doimpl(i,j-1).or.doimpl(i,j))THEN

        do k=1,nk
          tem = (mh(i,j-1,k)+mh(i,j,k))*rrv0(i,j,k)
          r1 = (kmv(i,j-1,k  )+kmv(i,j,k  ))*(mf(i,j-1,k  )+mf(i,j,k  ))   &
              *(rf0(i,j-1,k  )+rf0(i,j,k  ))*tem
          r2 = (kmv(i,j-1,k+1)+kmv(i,j,k+1))*(mf(i,j-1,k+1)+mf(i,j,k+1))   &
              *(rf0(i,j-1,k+1)+rf0(i,j,k+1))*tem
          cfa(k) = tema*r1
          cfc(k) = tema*r2
          cfd(k) = v(i,j,k)    &
                 + temb*( r2*v(i,j,k+1)-(r1+r2)*v(i,j,k)+r1*v(i,j,k-1) )
          if( dow1 )then
            dwdy2 = (w(i,j,k+1)-w(i,j-1,k+1))*rdy*vf(j)
            dwdy1 = (w(i,j,k  )-w(i,j-1,k  ))*rdy*vf(j)
          endif
          if( dow2 )then
            if(k.ne.nk) dwdy2 = dwdy2 + 0.5*(dum(i,j,k+1)+dum(i,j-1,k+1))
            if(k.ne. 1) dwdy1 = dwdy1 + 0.5*(dum(i,j,k  )+dum(i,j-1,k  ))
          endif
          if( dow1 )then
            r1 = min( 0.5*(kmv(i,j-1,k  )+kmv(i,j,k  )) , teme*rvh(j)*rvh(j) ) &
                         *(rf0(i,j-1,k  )+rf0(i,j,k  ))*tem
            r2 = min( 0.5*(kmv(i,j-1,k+1)+kmv(i,j,k+1)) , teme*rvh(j)*rvh(j) ) &
                         *(rf0(i,j-1,k+1)+rf0(i,j,k+1))*tem
            cfd(k) = cfd(k) +    &
                   temd*( r2*dwdy2   &
                        - r1*dwdy1 )
          endif
        enddo
        IF(bcturbu.eq.1)THEN
          cfa( 1) = 0.0
          cfc(nk) = 0.0
        ELSEIF(bcturbu.eq.2)THEN
          cfa( 1) = 0.0
          cfc( 1) = 0.0
          cfa(nk) = 0.0
          cfc(nk) = 0.0
        ELSEIF(bcturbu.eq.3)THEN
          cfa( 1) = 0.0
          cfc(nk) = 0.0
        ENDIF
        do k=1,nk
          cfb(k) = 1.0 - cfa(k) - cfc(k)
        enddo
        if(dodrag)then
          tem = temc*t23(i,j,1)*(mh(i,j-1,1)+mh(i,j,1))*rrv0(i,j,1)
          ! set back to simpler version:  now usable if pertflx = 1
!!!          cfb(1) = cfb(1) + vialpha*tem/(1.0e-20+v(i,j,1))
!!!          cfd(1) = cfd(1) - vibeta*tem
          cfd(1) = cfd(1) - tem
        endif

        lgbth(1)=-cfc(1)/cfb(1)
        lgbph(1)= cfd(1)/cfb(1)
        do k=2,nk
          tem = 1.0/(cfa(k)*lgbth(k-1)+cfb(k))
          lgbth(k)=-cfc(k)*tem
          lgbph(k)=(cfd(k)-cfa(k)*lgbph(k-1))*tem
        enddo
        v2(nk)=lgbph(nk)
        do k=nk-1,1,-1
          v2(k)=lgbth(k)*v2(k+1)+lgbph(k)
        enddo

        do k=1,nk
          turbz(i,j,k) = (rho0(i,j-1,k)+rho0(i,j,k))*(v2(k)-v(i,j,k))*rdt
        enddo

    ELSE

      do k=1,nk
        turbz(i,j,k)=(t23(i,j,k+1)-t23(i,j,k))*rdz*0.5*(mh(i,j-1,k)+mh(i,j,k))
      enddo

    ENDIF

    enddo
    enddo

  ENDIF

!-----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        vten(i,j,k)=vten(i,j,k)+(turbx(i,j,k)+turby(i,j,k)+turbz(i,j,k))*rrv0(i,j,k)
      enddo
      enddo
      enddo

!-------------------------------------------------------------------
!  All done
 
      if(timestats.ge.1) time_tmix=time_tmix+mytime()
 
      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
      subroutine turbw(dt,xh,rxh,uh,xf,vh,mh,mf,rho0,rf0,rrf0,   &
           turbx,turby,turbz,dum,w,wten,t13,t23,t33,kmv,gx,gy,doimpl)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real :: dt
      real, dimension(ib:ie) :: xh,rxh,uh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je) :: vh
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w,wten
      real, dimension(ib:ie,jb:je,kb:ke) :: t13,t23,t33
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv
      real, dimension(itb:ite+1,jtb:jte,ktb:kte) :: gx
      real, dimension(itb:ite,jtb:jte+1,ktb:kte) :: gy
      logical, dimension(ib:ie,jb:je) :: doimpl
 
      integer :: i,j,k
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2
      real, dimension(nk+1) :: cfa,cfb,cfc,cfd,cfe,cff,w2

!----------------------------------------------------------------
!  x-direction

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        turbx(i,j,k)=(t13(i+1,j,k)-t13(i,j,k))*rdx*uh(i)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        turbx(i,j,k)=(xf(i+1)*t13(i+1,j,k)-xf(i)*t13(i,j,k))*rdx*uh(i)*rxh(i)
      enddo
      enddo
      enddo

    ENDIF

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          dum(i,j,k)=(t13(i,j,k+1)-t13(i,j,k))*rdz
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=turbx(i,j,k)+0.0625*( (gx(i,j,k-1)+gx(i+1,j,k-1))       &
                                            +(gx(i,j,k  )+gx(i+1,j,k  )) )     &
                                          *( (dum(i,j,k-1)+dum(i+1,j,k-1))     &
                                            +(dum(i,j,k  )+dum(i+1,j,k  )) )
        enddo
        enddo
        enddo

      ENDIF

!----------------------------------------------------------------
!  y-direction

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        turby(i,j,k)=(t23(i,j+1,k)-t23(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          dum(i,j,k)=(t23(i,j,k+1)-t23(i,j,k))*rdz
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          turby(i,j,k)=turby(i,j,k)+0.0625*( (gy(i,j,k-1)+gy(i,j+1,k-1))       &
                                            +(gy(i,j,k  )+gy(i,j+1,k  )) )     &
                                          *( (dum(i,j,k-1)+dum(i,j+1,k-1))     &
                                            +(dum(i,j,k  )+dum(i,j+1,k  )) )
        enddo
        enddo
        enddo

      ENDIF

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        turby(i,j,k)=0.0
      enddo
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------
!  z-direction

    IF( iturb.eq.3 .and. l_v.lt.tsmall )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

    ELSE

      rdt = 1.0/dt
      tema = dt*vialpha*rdz*rdz
      temb = dt*vibeta*rdz*rdz
      temc = dt*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,cfe,cff,tem,w2)
    do j=1,nj
    do i=1,ni

    IF(doimpl(i,j))THEN

        do k=2,nk
          r1 = (kmv(i,j,k-1)+kmv(i,j,k  ))*mh(i,j,k-1)*rho0(i,j,k-1)*mf(i,j,k)*rrf0(i,j,k)
          r2 = (kmv(i,j,k  )+kmv(i,j,k+1))*mh(i,j,k  )*rho0(i,j,k  )*mf(i,j,k)*rrf0(i,j,k)
          cfa(k) = tema*r1
          cfc(k) = tema*r2
          cfb(k) = 1.0 + cfa(k) + cfc(k)
          cfd(k) = w(i,j,k)    &
                 + temb*(r2*w(i,j,k+1)-(r1+r2)*w(i,j,k)+r1*w(i,j,k-1))
        enddo

        cfe(1)=0.0
        cff(1)=w(i,j,1)
        do k=2,nk
          tem = 1.0/(cfb(k)-cfc(k)*cfe(k-1))
          cfe(k)=cfa(k)*tem
          cff(k)=(cfd(k)+cfc(k)*cff(k-1))*tem
        enddo
        w2(nk+1)=0.0
        do k=nk,2,-1
          w2(k)=cfe(k)*w2(k+1)+cff(k)
        enddo
        do k=2,nk
          turbz(i,j,k) = rf0(i,j,k)*(w2(k)-w(i,j,k))*rdt
        enddo

    ELSE

      do k=2,nk
        turbz(i,j,k)=(t33(i,j,k)-t33(i,j,k-1))*rdz*mf(i,j,k)
      enddo

    ENDIF

    enddo
    enddo

  ENDIF

!-----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        wten(i,j,k)=wten(i,j,k)+(turbx(i,j,k)+turby(i,j,k)+turbz(i,j,k))*rrf0(i,j,k)
      enddo
      enddo
      enddo

!-------------------------------------------------------------------
!  All done

      if(timestats.ge.1) time_tmix=time_tmix+mytime()
 
      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine getepsilon(rxh,uh,xf,rxf,uf,yh,vh,yf,vf,mh,mf,rr0,rrf0,   &
                            dum1,dum2,dum3,dum4,dum5,tem1,tem2,epsd,       &
                            t11,t12,t13,t22,t23,t33,ua,va,wa,gx,gy)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie) :: rxh,uh
      real, dimension(ib:ie+1) :: xf,rxf,uf
      real, dimension(jb:je) :: yh,vh
      real, dimension(jb:je+1) :: yf,vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: rr0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,tem1,tem2,epsd
      real, dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, dimension(itb:ite+1,jtb:jte,ktb:kte) :: gx
      real, dimension(itb:ite,jtb:jte+1,ktb:kte) :: gy

      integer :: i,j,k
      real :: tmp11,tmp22,tmp33,tmp12,tmp13,tmp23

!-----------------------------------------------------------------------
!  Cartesian grid:

      IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          dum1(i,j,k)=(ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)   &
                     +(va(i,j,k)-va(i-1,j,k))*rdx*uf(i)
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni+1
          dum2(i,j,k)=(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)   &
                     +(ua(i,j,k)-ua(i,j,k-1))*rdz*0.5*(mf(i-1,j,k)+mf(i,j,k))
        enddo
        enddo
        enddo

        IF(bcturbu.eq.1.or.bcturbu.eq.2)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni+1
            dum2(i,j,1)=0.0
            dum2(i,j,nk+1)=0.0
          enddo
          enddo
        ELSEIF(bcturbu.eq.3)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni+1
            dum2(i,j,1)=2.0*ua(i,j,1)*rdz*0.5*(mf(i-1,j,1)+mf(i,j,1))
            dum2(i,j,nk+1)=-2.0*ua(i,j,nk)*rdz*0.5*(mf(i-1,j,nk+1)+mf(i,j,nk+1))
          enddo
          enddo
        ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj+1
        do i=1,ni
          dum3(i,j,k)=(wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j)   &
                     +(va(i,j,k)-va(i,j,k-1))*rdz*0.5*(mf(i,j-1,k)+mf(i,j,k))
        enddo
        enddo
        enddo

        IF(bcturbu.eq.1.or.bcturbu.eq.2)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj+1
          do i=1,ni
            dum3(i,j,1)=0.0
            dum3(i,j,nk+1)=0.0
          enddo
          enddo
        ELSEIF(bcturbu.eq.3)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj+1
          do i=1,ni
            dum3(i,j,1)=2.0*va(i,j,1)*rdz*0.5*(mf(i,j-1,1)+mf(i,j,1))
            dum3(i,j,nk+1)=-2.0*va(i,j,nk)*rdz*0.5*(mf(i,j-1,nk+1)+mf(i,j,nk+1))
          enddo
          enddo
        ENDIF

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum4(i,j,k)=(ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)
        enddo
        enddo
        enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum5(i,j,k)=(va(i,j+1,k)-va(i,j,k))*rdy*vh(j)
        enddo
        enddo
        enddo

!------------------------------------------------------------------
!  now, add the "correction" terms for terrain

      IF(terrain_flag)THEN

!-------- dum1 --------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk-1
        do j=0,nj+1
        do i=1,ni+1
          tem1(i,j,k)=(ua(i,j,k+1)-ua(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=0,nj+1
        do i=1,ni+1
          tem1(i,j,1 )=(ua(i,j,2 )-ua(i,j,1   ))*rdz
          tem1(i,j,nk)=(ua(i,j,nk)-ua(i,j,nk-1))*rdz
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=1,ni+1
          tem1(i,j,k)=tem1(i,j,k)*0.25*( gy(i-1,j  ,k)+gy(i,j  ,k)    &
                                        +gy(i-1,j+1,k)+gy(i,j+1,k) )
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk-1
        do j=1,nj+1
        do i=0,ni+1
          tem2(i,j,k)=(va(i,j,k+1)-va(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=0,ni+1
          tem2(i,j,1 )=(va(i,j,2 )-va(i,j,1   ))*rdz
          tem2(i,j,nk)=(va(i,j,nk)-va(i,j,nk-1))*rdz
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=0,ni+1
          tem2(i,j,k)=tem2(i,j,k)*0.25*( gx(i  ,j-1,k)+gx(i  ,j,k)    &
                                        +gx(i+1,j-1,k)+gx(i+1,j,k) )
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          dum1(i,j,k)=dum1(i,j,k)+0.5*( (tem1(i,j-1,k)+tem1(i,j,k))  &
                                       +(tem2(i-1,j,k)+tem2(i,j,k)) )
        enddo
        enddo
        enddo

!-------- dum2 --------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=0,ni+1
          tem1(i,j,k)=0.25*( (gx(i,j,k-1)+gx(i+1,j,k-1))     &
                            +(gx(i,j,k  )+gx(i+1,j,k  )) )   &
                     *(wa(i,j,k+1)-wa(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni+1
          dum2(i,j,k)=dum2(i,j,k)+0.5*(tem1(i,j,k)+tem1(i-1,j,k))
        enddo
        enddo
        enddo

!-------- dum3 --------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=0,nj+1
        do i=1,ni
          tem1(i,j,k)=0.25*( (gy(i,j,k-1)+gy(i,j+1,k-1))     &
                            +(gy(i,j,k  )+gy(i,j+1,k  )) )   &
                     *(wa(i,j,k+1)-wa(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj+1
        do i=1,ni
          dum3(i,j,k)=dum3(i,j,k)+0.5*(tem1(i,j,k)+tem1(i,j-1,k))
        enddo
        enddo
        enddo

!-------- dum4 --------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk-1
        do j=1,nj
        do i=0,ni+2
          tem1(i,j,k)=gx(i,j,k)*(ua(i,j,k+1)-ua(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=0,ni+2
          tem1(i,j, 1)=gx(i,j,1 )*(ua(i,j,2 )-ua(i,j,1   ))*rdz
          tem1(i,j,nk)=gx(i,j,nk)*(ua(i,j,nk)-ua(i,j,nk-1))*rdz
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=0,ni+1
          dum4(i,j,k)=dum4(i,j,k)+0.5*(tem1(i,j,k)+tem1(i+1,j,k))
        enddo
        enddo
        enddo

!-------- dum5 --------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk-1
        do j=0,nj+2
        do i=1,ni
          tem1(i,j,k)=gy(i,j,k)*(va(i,j,k+1)-va(i,j,k-1))*rdz2
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=0,nj+2
        do i=1,ni
          tem1(i,j, 1)=gy(i,j,1 )*(va(i,j,2 )-va(i,j,1   ))*rdz
          tem1(i,j,nk)=gy(i,j,nk)*(va(i,j,nk)-va(i,j,nk-1))*rdz
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=1,ni
          dum5(i,j,k)=dum5(i,j,k)+0.5*(tem1(i,j,k)+tem1(i,j+1,k))
        enddo
        enddo
        enddo

      ENDIF

!------------------------------------------------------------------

!$omp parallel do default(shared)  &
!$omp private(i,j,k,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          tmp11 = t11(i,j,k)*dum4(i,j,k)
          tmp22 = t22(i,j,k)*dum5(i,j,k)
          tmp33 = t33(i,j,k)*(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(i,j,k)
          tmp12 = 0.25*( ( t12(i  ,j  ,k)*dum1(i  ,j  ,k)   &
                          +t12(i+1,j+1,k)*dum1(i+1,j+1,k) ) &
                        +( t12(i+1,j  ,k)*dum1(i+1,j  ,k)   &
                          +t12(i  ,j+1,k)*dum1(i  ,j+1,k) ) )
          tmp13 = 0.25*( ( t13(i  ,j,k  )*dum2(i  ,j,k  )   &
                          +t13(i+1,j,k  )*dum2(i+1,j,k  ) ) &
                        +( t13(i  ,j,k+1)*dum2(i  ,j,k+1)   &
                          +t13(i+1,j,k+1)*dum2(i+1,j,k+1) ) )
          tmp23 = 0.25*( ( t23(i,j  ,k  )*dum3(i,j  ,k  )   &
                          +t23(i,j+1,k  )*dum3(i,j+1,k  ) ) &
                        +( t23(i,j  ,k+1)*dum3(i,j  ,k+1)   &
                          +t23(i,j+1,k+1)*dum3(i,j+1,k+1) ) )
          epsd(i,j,k) = rr0(i,j,k)*((tmp11+tmp22+tmp33)+tmp12+(tmp13+tmp23))
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  Axisymmetric grid:

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=2,ni+1
          dum1(i,j,k)=xf(i)*(va(i,j,k)*rxh(i)-va(i-1,j,k)*rxh(i-1))*rdx*uf(i)
        enddo
        enddo
        enddo
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          dum1(1,j,k)=0.0
        enddo
        enddo
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni+1
          dum2(i,j,k)=(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)   &
                     +(ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k)
        enddo
        enddo
        enddo
        IF(bcturbu.eq.1.or.bcturbu.eq.2)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni+1
            dum2(i,j,1)=0.0
            dum2(i,j,nk+1)=0.0
          enddo
          enddo
        ELSEIF(bcturbu.eq.3)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni+1
            dum2(i,j,1)=2.0*ua(i,j,1)*rdz*mf(1,1,1)
            dum2(i,j,nk+1)=-2.0*ua(i,j,nk)*rdz*mf(1,1,nk+1)
          enddo
          enddo
        ENDIF
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj+1
        do i=1,ni
          dum3(i,j,k)=(va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)
        enddo
        enddo
        enddo
        IF(bcturbu.eq.1.or.bcturbu.eq.2)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj+1
          do i=1,ni
            dum3(i,j,1)=0.0
            dum3(i,j,nk+1)=0.0
          enddo
          enddo
        ELSEIF(bcturbu.eq.3)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj+1
          do i=1,ni
            dum3(i,j,1)=2.0*va(i,j,1)*rdz*mf(1,1,1)
            dum3(i,j,nk+1)=-2.0*va(i,j,nk)*rdz*mf(1,1,nk+1)
          enddo
          enddo
        ENDIF
!$omp parallel do default(shared)  &
!$omp private(i,j,k,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23)
        do k=1,nk
        do j=1,nj
        do i=1,ni  
          tmp11 = t11(i,j,k)*(ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)
          tmp33 = t33(i,j,k)*(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(i,j,k)
          tmp22 = 0.5*( t22(i  ,j,k)*ua(i  ,j,k)*rxf(i  )   &
                       +t22(i+1,j,k)*ua(i+1,j,k)*rxf(i+1) )
          tmp12 = 0.5*( t12(i  ,j,k)*dum1(i  ,j,k) &
                       +t12(i+1,j,k)*dum1(i+1,j,k) )
          tmp23 = 0.5*( t23(i,j,k  )*dum3(i,j,k  ) &
                       +t23(i,j,k+1)*dum3(i,j,k+1) )
          tmp13 = 0.25*( t13(i  ,j,k  )*dum2(i  ,j,k  ) &
                        +t13(i+1,j,k  )*dum2(i+1,j,k  ) &
                        +t13(i  ,j,k+1)*dum2(i  ,j,k+1) &
                        +t13(i+1,j,k+1)*dum2(i+1,j,k+1) )
          epsd(i,j,k) = rr0(1,1,k)*(tmp11+tmp22+tmp33+tmp12+tmp13+tmp23)
        enddo
        enddo
        enddo

      ENDIF

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end

