

      subroutine diff2u(flag,rxh,uh,xf,uf,vh,vf,mh,mf,dumx,dumy,diffx,diffy,epsd,ua,uten)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      integer flag
      real, dimension(ib:ie) :: rxh,uh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je) :: vh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy,diffx,diffy,epsd
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua,uten

      integer i,j,k
      real coef

!--------------------------
!
!  flag = 1 does 2nd order artificial diffusion
!  flag = 2 does dns viscosity term
!
!-----------------------------------------------------------------------
!  x-direction

      if(flag.eq.1)then
        coef=kdiff2
      elseif(flag.eq.2)then
        coef=viscosity
      endif

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=0,ni+1
        dumx(i,j,k)=(ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=0,ni+1
        dumx(i,j,k)=(xf(i+1)*ua(i+1,j,k)-xf(i)*ua(i,j,k))*rdx*uh(i)*rxh(i)
      enddo
      enddo
      enddo

    ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        diffx(i,j,k)=coef*(dumx(i,j,k)-dumx(i-1,j,k))*rdx*uf(i)
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------
!  y-direction

      if(flag.eq.1)then
        coef=kdiff2
      elseif(flag.eq.2)then
        coef=viscosity
      endif

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni+1
        dumy(i,j,k)=(ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)
      enddo
      enddo
      enddo
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        diffy(i,j,k)=coef*(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo

    ELSE
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni+1
        dumy(i,j,k)=0.0
        diffy(i,j,k)=0.0
      enddo
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        uten(i,j,k)=uten(i,j,k)+diffx(i,j,k)+diffy(i,j,k)
      enddo
      enddo
      enddo

      IF(idiss.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          epsd(i,j,k)=epsd(i,j,k)+coef*( (dumx(i,j,k)**2)       &
                +0.25*( (dumy(i  ,j,k)**2+dumy(i  ,j+1,k)**2)   &
                       +(dumy(i+1,j,k)**2+dumy(i+1,j+1,k)**2) ) )
        enddo
        enddo
        enddo
      ENDIF

!-----------------------------------------------------------------------
!  z-direction

      IF(vdiff.eq.1.or.flag.eq.2)THEN

        if(flag.eq.1)then
          coef=kdiff2
        elseif(flag.eq.2)then
          coef=viscosity
        endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni+1
          dumx(i,j,k)=(ua(i,j,k)-ua(i,j,k-1))*rdz*0.5*(mf(i-1,j,k)+mf(i,j,k))
        enddo
        enddo
        enddo

        IF(flag.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni+1
            dumx(i,j,1)=0.0
            dumx(i,j,nk+1)=0.0
          enddo
          enddo

        ELSEIF(flag.eq.2)THEN

          if(bc_wind.eq.1)then      ! free slip b.c.

!$omp parallel do default(shared)   &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni+1
              dumx(i,j,1)=0.0
              dumx(i,j,nk+1)=0.0
            enddo
            enddo

          elseif(bc_wind.eq.2)then      ! no slip b.c.

!$omp parallel do default(shared)   &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni+1
              dumx(i,j,1)=2.0*ua(i,j,1)*rdz*0.5*(mf(i-1,j,1)+mf(i,j,1))
            enddo
            enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni+1
              dumx(i,j,nk+1)=-2.0*ua(i,j,nk)*rdz*0.5*(mf(i-1,j,nk+1)+mf(i,j,nk+1))
            enddo
            enddo

          endif

        ENDIF
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          uten(i,j,k)=uten(i,j,k)+coef*(dumx(i,j,k+1)-dumx(i,j,k))*rdz*0.5*(mh(i-1,j,k)+mh(i,j,k))
        enddo
        enddo
        enddo

        IF(idiss.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            epsd(i,j,k)=epsd(i,j,k)+0.25*coef*(               &
                          (dumx(i,j,k  )**2+dumx(i+1,j,k  )**2)   &
                         +(dumx(i,j,k+1)**2+dumx(i+1,j,k+1)**2)   &
                                              )
          enddo
          enddo
          enddo
        ENDIF

      ENDIF

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_diffu=time_diffu+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine diff2v(flag,xh,uh,rxf,uf,vh,vf,mh,mf,dumx,dumy,diffx,diffy,epsd,va,vten)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      integer flag
      real, dimension(ib:ie) :: xh,uh
      real, dimension(ib:ie+1) :: rxf,uf
      real, dimension(jb:je) :: vh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy,diffx,diffy,epsd
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va,vten
 
      integer i,j,k
      real coef
 
!--------------------------
!
!  flag = 1 does 2nd order artificial diffusion
!  flag = 2 does dns viscosity term
!
!-----------------------------------------------------------------------
!  x-direction

      if(flag.eq.1)then
        coef=kdiff2
      elseif(flag.eq.2)then
        coef=viscosity
      endif

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni+1
        dumx(i,j,k)=(va(i,j,k)-va(i-1,j,k))*rdx*uf(i)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(j,k)
      do k=1,nk
      do j=1,nj+1
        dumx(1,j,k)=0.0
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=2,ni+1
        dumx(i,j,k)=(xh(i)*va(i,j,k)-xh(i-1)*va(i-1,j,k))*rdx*uf(i)*rxf(i)
      enddo
      enddo
      enddo

    ENDIF
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        diffx(i,j,k)=coef*(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------
!  y-direction

      if(flag.eq.1)then
        coef=kdiff2
      elseif(flag.eq.2)then
        coef=viscosity
      endif

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=0,nj+1
      do i=1,ni
        dumy(i,j,k)=(va(i,j+1,k)-va(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        diffy(i,j,k)=coef*(dumy(i,j,k)-dumy(i,j-1,k))*rdy*vf(j)
      enddo
      enddo
      enddo

    ELSE
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=0,nj+1
      do i=1,ni
        dumy(i,j,k)=0.0
        diffy(i,j,k)=0.0
      enddo
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        vten(i,j,k)=vten(i,j,k)+diffx(i,j,k)+diffy(i,j,k)
      enddo
      enddo
      enddo

      IF(idiss.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          epsd(i,j,k)=epsd(i,j,k)+coef*( (dumy(i,j,k)**2)       &
                +0.25*( (dumx(i,j  ,k)**2+dumx(i+1,j  ,k)**2)   &
                       +(dumx(i,j+1,k)**2+dumx(i+1,j+1,k)**2) ) )
        enddo
        enddo
        enddo
      ENDIF

!-----------------------------------------------------------------------
!  z-direction

      IF(vdiff.eq.1.or.flag.eq.2)THEN

        if(flag.eq.1)then
          coef=kdiff2
        elseif(flag.eq.2)then
          coef=viscosity
        endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj+1
        do i=1,ni
          dumx(i,j,k)=(va(i,j,k)-va(i,j,k-1))*rdz*0.5*(mf(i,j-1,k)+mf(i,j,k))
        enddo
        enddo
        enddo

        IF(flag.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj+1
          do i=1,ni
            dumx(i,j,1)=0.0
            dumx(i,j,nk+1)=0.0
          enddo
          enddo

        ELSEIF(flag.eq.2)THEN

          if(bc_wind.eq.1)then      ! free slip b.c.

!$omp parallel do default(shared)   &
!$omp private(i,j)
            do j=1,nj+1
            do i=1,ni
              dumx(i,j,1)=0.0
              dumx(i,j,nk+1)=0.0
            enddo
            enddo

          elseif(bc_wind.eq.2)then      ! no slip b.c.

!$omp parallel do default(shared)   &
!$omp private(i,j)
            do j=1,nj+1
            do i=1,ni
              dumx(i,j,1)=2.0*va(i,j,1)*rdz*0.5*(mf(i,j-1,1)+mf(i,j,1))
            enddo
            enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
            do j=1,nj+1
            do i=1,ni
              dumx(i,j,nk+1)=-2.0*va(i,j,nk)*rdz*0.5*(mf(i,j-1,nk+1)+mf(i,j,nk+1))
            enddo
            enddo

          endif

        ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vten(i,j,k)=vten(i,j,k)+coef*(dumx(i,j,k+1)-dumx(i,j,k))*rdz*0.5*(mh(i,j-1,k)+mh(i,j,k))
        enddo
        enddo
        enddo

        IF(idiss.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            epsd(i,j,k)=epsd(i,j,k)+0.25*coef*(               &
                          (dumx(i,j,k  )**2+dumx(i,j+1,k  )**2)   &
                         +(dumx(i,j,k+1)**2+dumx(i,j+1,k+1)**2)   &
                                              )
          enddo
          enddo
          enddo
        ENDIF

      ENDIF
 
!-----------------------------------------------------------------------

      if(timestats.ge.1) time_diffu=time_diffu+mytime()
 
      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine diff2w(flag,dissw,rxh,uh,xf,uf,vh,vf,mh,mf,dumx,dumy,diffx,diffy,epsd,wa,wten)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      integer flag,dissw
      real, dimension(ib:ie) :: rxh,uh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je) :: vh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy,diffx,diffy,epsd
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa,wten
 
      integer i,j,k
      real coef

!--------------------------
!
!  flag = 1 does 2nd order artificial diffusion
!  flag = 2 does dns viscosity term
!
!-----------------------------------------------------------------------
!  x-direction

      if(flag.eq.1)then
        coef=kdiff2
      elseif(flag.eq.2)then
        coef=viscosity
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk+1
      do j=1,nj
      do i=1,ni+1
        dumx(i,j,k)=(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)
      enddo
      enddo
      enddo

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk+1
      do j=1,nj
      do i=1,ni
        diffx(i,j,k)=coef*(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk+1
      do j=1,nj
      do i=1,ni
        diffx(i,j,k)=coef*(xf(i+1)*dumx(i+1,j,k)-xf(i)*dumx(i,j,k))*rdx*uh(i)*rxh(i)
      enddo
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------------
!  y-direction

      if(flag.eq.1)then
        coef=kdiff2
      elseif(flag.eq.2)then
        coef=viscosity
      endif

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk+1
      do j=1,nj+1
      do i=1,ni
        dumy(i,j,k)=(wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j)
      enddo
      enddo
      enddo
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk+1
      do j=1,nj
      do i=1,ni
        diffy(i,j,k)=coef*(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk+1
      do j=1,nj+1
      do i=1,ni
        dumy(i,j,k)=0.0
        diffy(i,j,k)=0.0
      enddo
      enddo
      enddo

    ENDIF
 
!-----------------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk+1
      do j=1,nj
      do i=1,ni
        wten(i,j,k)=wten(i,j,k)+diffx(i,j,k)+diffy(i,j,k)
      enddo
      enddo
      enddo

      IF(idiss.eq.1.and.dissw.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          epsd(i,j,k)=epsd(i,j,k)+0.25*coef*(                     &
                      ( (dumx(i,j,k  )**2+dumx(i+1,j,k  )**2)     &
                       +(dumx(i,j,k+1)**2+dumx(i+1,j,k+1)**2) )   &
                     +( (dumy(i,j,k  )**2+dumy(i,j+1,k  )**2)     &
                       +(dumy(i,j,k+1)**2+dumy(i,j+1,k+1)**2) ) )
        enddo
        enddo
        enddo
      ENDIF

!-----------------------------------------------------------------------
!  z-direction

      IF(vdiff.eq.1.or.flag.eq.2)THEN

        if(flag.eq.1)then
          coef=kdiff2
        elseif(flag.eq.2)then
          coef=viscosity
        endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dumx(i,j,k)=(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(i,j,k)
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          wten(i,j,k)=wten(i,j,k)+coef*(dumx(i,j,k)-dumx(i,j,k-1))*rdz*mf(i,j,k)
        enddo
        enddo
        enddo

        IF(idiss.eq.1.and.dissw.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            epsd(i,j,k)=epsd(i,j,k)+coef*(dumx(i,j,k)**2)
          enddo
          enddo
          enddo
        ENDIF

      ENDIF
 
!-----------------------------------------------------------------------

      if(timestats.ge.1) time_diffu=time_diffu+mytime()
 
      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
 
      subroutine diff2s(flag,rxh,uh,xf,uf,vh,vf,mh,mf,dumx,dumy,dum,s,sten)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      integer flag
      real, dimension(ib:ie) :: rxh,uh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je) :: vh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy,dum
      real, dimension(ib:ie,jb:je,kb:ke) :: s,sten
 
      integer i,j,k
      real coef

!--------------------------
!
!  flag = 1 does 2nd order artificial diffusion
!  flag = 2 does dns conduction term
!
!-----------------------------------------------------------------------
!  x-direction

      if(flag.eq.1)then
        coef=kdiff2
      elseif(flag.eq.2)then
        coef=viscosity/pr_num
      endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        dum(i,j,k)=(s(i,j,k)-s(i-1,j,k))*rdx*uf(i)
      enddo
      enddo
      enddo

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dumx(i,j,k)=coef*(dum(i+1,j,k)-dum(i,j,k))*rdx*uh(i)
      enddo
      enddo
      enddo

    ELSEIF(axisymm.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dumx(i,j,k)=coef*(xf(i+1)*dum(i+1,j,k)-xf(i)*dum(i,j,k))*rdx*uh(i)*rxh(i)
      enddo
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------------
!  y-direction

      if(flag.eq.1)then
        coef=kdiff2
      elseif(flag.eq.2)then
        coef=viscosity/pr_num
      endif

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        dum(i,j,k)=(s(i,j,k)-s(i,j-1,k))*rdy*vf(j)
      enddo
      enddo
      enddo
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dumy(i,j,k)=coef*(dum(i,j+1,k)-dum(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dumy(i,j,k)=0.0
      enddo
      enddo
      enddo

    ENDIF
 
!-----------------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=sten(i,j,k)+dumx(i,j,k)+dumy(i,j,k)
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------
!  z-direction

      IF(vdiff.eq.1.or.flag.eq.2)THEN

        if(flag.eq.1)then
          coef=kdiff2
        elseif(flag.eq.2)then
          coef=viscosity/pr_num
        endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=(s(i,j,k)-s(i,j,k-1))*rdz*mf(i,j,k)
        enddo
        enddo
        enddo

        IF(flag.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum(i,j,1)=0.0
            dum(i,j,nk+1)=0.0
          enddo
          enddo

        ELSEIF(flag.eq.2)THEN

          if(bc_temp.eq.1)then      ! constant theta at boundary

!$omp parallel do default(shared)   &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              dum(i,j,1)=2.0*(s(i,j,1)-ptc_bot)*rdz*mf(i,j,1)
            enddo
            enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              dum(i,j,nk+1)=2.0*(ptc_top-s(i,j,nk))*rdz*mf(i,j,nk+1)
            enddo
            enddo

          elseif(bc_temp.eq.2)then      ! constant flux at boundary

!$omp parallel do default(shared)   &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              dum(i,j,1)=ptc_bot
            enddo
            enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              dum(i,j,nk+1)=ptc_top
            enddo
            enddo

          endif

        ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          sten(i,j,k)=sten(i,j,k)+coef*(dum(i,j,k+1)-dum(i,j,k))*rdz*mh(i,j,k)
        enddo
        enddo
        enddo

      ENDIF
 
!-----------------------------------------------------------------------

      if(timestats.ge.1) time_diffu=time_diffu+mytime()
 
      return
      end


