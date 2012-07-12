

      subroutine diff6s(dt,s0,dumx,dumy,p,s,sten)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real :: dt
      real, dimension(ib:ie,jb:je,kb:ke) :: s0
      real, dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy,p
      real, dimension(ib:ie,jb:je,kb:ke) :: s,sten

      integer i,j,k

      real coef
      coef=(kdiff6/64.0/dt)

!------------------------------------------------------------
!  x-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        p(i,j,k)=( 10.0*(s(i  ,j,k)-s(i-1,j,k))     &
                   -5.0*(s(i+1,j,k)-s(i-2,j,k))     &
                       +(s(i+2,j,k)-s(i-3,j,k)) )
      enddo
      enddo
      enddo

      if(mdiff.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          if( p(i,j,k)*(s(i,j,k)-s(i-1,j,k)).le.0.0 )then
            p(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dumx(i,j,k)=coef*(p(i+1,j,k)-p(i,j,k))
      enddo
      enddo
      enddo

!------------------------------------------------------------
!  y-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        p(i,j,k)=( 10.0*(s(i,j  ,k)-s(i,j-1,k))     &
                   -5.0*(s(i,j+1,k)-s(i,j-2,k))     &
                       +(s(i,j+2,k)-s(i,j-3,k)) )
      enddo
      enddo
      enddo

      if(mdiff.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          if( p(i,j,k)*(s(i,j,k)-s(i,j-1,k)).le.0.0 )then
            p(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dumy(i,j,k)=coef*(p(i,j+1,k)-p(i,j,k))
      enddo
      enddo
      enddo

!------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=sten(i,j,k)+(dumx(i,j,k)+dumy(i,j,k))
      enddo
      enddo
      enddo

!------------------------------------------------------------
!  z-direction

    IF(vdiff.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dumx(i,j,k)=s(i,j,k)-s0(i,j,k)
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=4,nk-2
      do j=1,nj
      do i=1,ni
        p(i,j,k)=( 10.0*(dumx(i,j,k  )-dumx(i,j,k-1))     &
                   -5.0*(dumx(i,j,k+1)-dumx(i,j,k-2))     &
                       +(dumx(i,j,k+2)-dumx(i,j,k-3)) )
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        p(i,j,3)   =-p(i,j,5)
        p(i,j,2)   = p(i,j,4)
        p(i,j,nk  )= p(i,j,nk-2)
        p(i,j,nk-1)=-p(i,j,nk-3)
      enddo
      enddo

      if(mdiff.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          if( p(i,j,k)*(dumx(i,j,k)-dumx(i,j,k-1)).le.0.0 )then
            p(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo
      endif

      IF(bcturbs.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          p(i,j,1)    = 0.0
          p(i,j,nk+1) = 0.0
        enddo
        enddo

      ELSEIF(bcturbs.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          p(i,j,1)    = p(i,j,2)
          p(i,j,nk+1) = p(i,j,nk)
        enddo
        enddo

      ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=sten(i,j,k)+coef*(p(i,j,k+1)-p(i,j,k))
      enddo
      enddo
      enddo

    ENDIF

!------------------------------------------------------------

      if(timestats.ge.1) time_diffu=time_diffu+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine diff6u(dt,u0,dumx,dumy,p,u,uten)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real :: dt
      real, dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy,p
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u0,u,uten
 
      integer i,j,k
 
      real coef
      coef=(kdiff6/64.0/dt)

!------------------------------------------------------------
!  x-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+2
        p(i,j,k)=( 10.0*(u(i  ,j,k)-u(i-1,j,k))     &
                   -5.0*(u(i+1,j,k)-u(i-2,j,k))     &
                       +(u(i+2,j,k)-u(i-3,j,k)) )
      enddo
      enddo
      enddo

      if(mdiff.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+2
          if( p(i,j,k)*(u(i,j,k)-u(i-1,j,k)).le.0.0 )then
            p(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        dumx(i,j,k)=coef*(p(i+1,j,k)-p(i,j,k))
      enddo
      enddo
      enddo

!------------------------------------------------------------
!  y-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni+1
        p(i,j,k)=( 10.0*(u(i,j  ,k)-u(i,j-1,k))     &
                   -5.0*(u(i,j+1,k)-u(i,j-2,k))     &
                       +(u(i,j+2,k)-u(i,j-3,k)) )
      enddo
      enddo
      enddo

      if(mdiff.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          if( p(i,j,k)*(u(i,j,k)-u(i,j-1,k)).le.0.0 )then
            p(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        dumy(i,j,k)=coef*(p(i,j+1,k)-p(i,j,k))
      enddo
      enddo
      enddo

!------------------------------------------------------------

!$omp parallel do default(shared)    &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        uten(i,j,k)=uten(i,j,k)+(dumx(i,j,k)+dumy(i,j,k))
      enddo
      enddo
      enddo

!------------------------------------------------------------
!  z-direction

    IF(vdiff.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        dumx(i,j,k)=u(i,j,k)-u0(i,j,k)
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=4,nk-2
      do j=1,nj
      do i=1,ni+1
        p(i,j,k)=( 10.0*(dumx(i,j,k  )-dumx(i,j,k-1))     &
                   -5.0*(dumx(i,j,k+1)-dumx(i,j,k-2))     &
                       +(dumx(i,j,k+2)-dumx(i,j,k-3)) )
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1
        p(i,j,3)   =-p(i,j,5)
        p(i,j,2)   = p(i,j,4)
        p(i,j,nk  )= p(i,j,nk-2)
        p(i,j,nk-1)=-p(i,j,nk-3)
      enddo
      enddo

      if(mdiff.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni+1
          if( p(i,j,k)*(dumx(i,j,k)-dumx(i,j,k-1)).le.0.0 )then
            p(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo
      endif

      IF(bcturbu.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni+1
          p(i,j,1)    = 0.0
          p(i,j,nk+1) = 0.0
        enddo
        enddo

      ELSEIF(bcturbu.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni+1
          p(i,j,1)    = p(i,j,2)
          p(i,j,nk+1) = p(i,j,nk)
        enddo
        enddo

      ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        uten(i,j,k)=uten(i,j,k)+coef*(p(i,j,k+1)-p(i,j,k))
      enddo
      enddo
      enddo

    ENDIF

!------------------------------------------------------------
 
      if(timestats.ge.1) time_diffu=time_diffu+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine diff6v(dt,v0,dumx,dumy,p,v,vten)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real :: dt
      real, dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy,p
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v0,v,vten

      integer i,j,k
 
      real coef
      coef=(kdiff6/64.0/dt)

!------------------------------------------------------------
!  x-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni+1
        p(i,j,k)=( 10.0*(v(i  ,j,k)-v(i-1,j,k))     &
                   -5.0*(v(i+1,j,k)-v(i-2,j,k))     &
                       +(v(i+2,j,k)-v(i-3,j,k)) )
      enddo
      enddo
      enddo

      if(mdiff.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          if( p(i,j,k)*(v(i,j,k)-v(i-1,j,k)).le.0.0 )then
            p(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        dumx(i,j,k)=coef*(p(i+1,j,k)-p(i,j,k))
      enddo
      enddo
      enddo

!------------------------------------------------------------
!  y-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+2
      do i=1,ni
        p(i,j,k)=( 10.0*(v(i,j  ,k)-v(i,j-1,k))     &
                   -5.0*(v(i,j+1,k)-v(i,j-2,k))     &
                       +(v(i,j+2,k)-v(i,j-3,k)) )
      enddo
      enddo
      enddo

      if(mdiff.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+2
        do i=1,ni
          if( p(i,j,k)*(v(i,j,k)-v(i,j-1,k)).le.0.0 )then
            p(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        dumy(i,j,k)=coef*(p(i,j+1,k)-p(i,j,k))
      enddo
      enddo
      enddo

!------------------------------------------------------------

!$omp parallel do default(shared)    &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        vten(i,j,k)=vten(i,j,k)+(dumx(i,j,k)+dumy(i,j,k))
      enddo
      enddo
      enddo

!------------------------------------------------------------
!  z-direction

    IF(vdiff.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        dumx(i,j,k)=v(i,j,k)-v0(i,j,k)
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=4,nk-2
      do j=1,nj+1
      do i=1,ni
        p(i,j,k)=( 10.0*(dumx(i,j,k  )-dumx(i,j,k-1))     &
                   -5.0*(dumx(i,j,k+1)-dumx(i,j,k-2))     &
                       +(dumx(i,j,k+2)-dumx(i,j,k-3)) )
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        p(i,j,3)   =-p(i,j,5)
        p(i,j,2)   = p(i,j,4)
        p(i,j,nk  )= p(i,j,nk-2)
        p(i,j,nk-1)=-p(i,j,nk-3)
      enddo
      enddo

      if(mdiff.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj+1
        do i=1,ni
          if( p(i,j,k)*(dumx(i,j,k)-dumx(i,j,k-1)).le.0.0 )then
            p(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo
      endif

      IF(bcturbu.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni
          p(i,j,1)    = 0.0
          p(i,j,nk+1) = 0.0
        enddo
        enddo

      ELSEIF(bcturbu.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni
          p(i,j,1)    = p(i,j,2)
          p(i,j,nk+1) = p(i,j,nk)
        enddo
        enddo

      ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        vten(i,j,k)=vten(i,j,k)+coef*(p(i,j,k+1)-p(i,j,k))
      enddo
      enddo
      enddo

    ENDIF

!------------------------------------------------------------
 
      if(timestats.ge.1) time_diffu=time_diffu+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine diff6w(dt,dumx,dumy,p,w,wten)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real :: dt
      real, dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy,p
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w,wten

      integer i,j,k
 
      real coef
      coef=(kdiff6/64.0/dt)

!------------------------------------------------------------
!  x-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni+1
        p(i,j,k)=( 10.0*(w(i  ,j,k)-w(i-1,j,k))     &
                   -5.0*(w(i+1,j,k)-w(i-2,j,k))     &
                       +(w(i+2,j,k)-w(i-3,j,k)) )
      enddo
      enddo
      enddo

      if(mdiff.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni+1
          if( p(i,j,k)*(w(i,j,k)-w(i-1,j,k)).le.0.0 )then
            p(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        dumx(i,j,k)=coef*(p(i+1,j,k)-p(i,j,k))
      enddo
      enddo
      enddo

!------------------------------------------------------------
!  y-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj+1
      do i=1,ni
        p(i,j,k)=( 10.0*(w(i,j  ,k)-w(i,j-1,k))     &
                   -5.0*(w(i,j+1,k)-w(i,j-2,k))     &
                       +(w(i,j+2,k)-w(i,j-3,k)) )
      enddo
      enddo
      enddo

      if(mdiff.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj+1
        do i=1,ni
          if( p(i,j,k)*(w(i,j,k)-w(i,j-1,k)).le.0.0 )then
            p(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        dumy(i,j,k)=coef*(p(i,j+1,k)-p(i,j,k))
      enddo
      enddo
      enddo

!------------------------------------------------------------

!$omp parallel do default(shared)    &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        wten(i,j,k)=wten(i,j,k)+(dumx(i,j,k)+dumy(i,j,k))
      enddo
      enddo
      enddo

!------------------------------------------------------------
!  z-direction

    IF(vdiff.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=4,nk-1
      do j=1,nj
      do i=1,ni
        p(i,j,k)=( 10.0*(w(i,j,k  )-w(i,j,k-1))     &
                   -5.0*(w(i,j,k+1)-w(i,j,k-2))     &
                       +(w(i,j,k+2)-w(i,j,k-3)) )
      enddo
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        p(i,j,3)   =-p(i,j,5)
        p(i,j,2)   = p(i,j,4)
        p(i,j,nk+1)= p(i,j,nk-1)
        p(i,j,nk  )=-p(i,j,nk-2)
      enddo
      enddo

      if(mdiff.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk+1
        do j=1,nj
        do i=1,ni
          if( p(i,j,k)*(w(i,j,k)-w(i,j,k-1)).le.0.0 )then
            p(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        wten(i,j,k)=wten(i,j,k)+coef*(p(i,j,k+1)-p(i,j,k))
      enddo
      enddo
      enddo

    ENDIF

!------------------------------------------------------------
 
      if(timestats.ge.1) time_diffu=time_diffu+mytime()
 
      return
      end


