



      subroutine wenos(bflag,bsq,xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,      &
                       rho0,rr0,rf0,rrf0,advx,advy,advz,dum,   &
                       divx,rru,rrv,rrw,s,sten,dt)
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
      real, dimension(ib:ie,jb:je,kb:ke) :: advx,advy,advz,dum,divx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, dimension(ib:ie,jb:je,kb:ke) :: s,sten
      real :: dt

      integer i,j,k

      real s1,s2,s3,s4,s5
      real f1,f2,f3
      real b1,b2,b3
      real w1,w2,w3
      real tem

      real epsilon,onedsix,thdtw
      parameter(epsilon=1.0e-08)
      parameter(onedsix=1.0/6.0)
      parameter(thdtw=13.0/12.0)

      real, parameter :: tem2 = 1.0/6.0
      real, parameter :: tem3 = 1.0/12.0
 
!----------------------------------------------------------------
! Advection in x-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k,s1,s2,s3,s4,s5,f1,f2,f3,b1,b2,b3,w1,w2,w3)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        if(rru(i,j,k).ge.0.0)then
          s1=s(i-3,j,k)
          s2=s(i-2,j,k)
          s3=s(i-1,j,k)
          s4=s(i  ,j,k)
          s5=s(i+1,j,k)
        else
          s1=s(i+2,j,k)
          s2=s(i+1,j,k)
          s3=s(i  ,j,k)
          s4=s(i-1,j,k)
          s5=s(i-2,j,k)
        endif

        f1=( 2.0*s1 -7.0*s2 +11.0*s3 )*onedsix
        f2=(    -s2 +5.0*s3  +2.0*s4 )*onedsix
        f3=( 2.0*s3 +5.0*s4      -s5 )*onedsix

        b1=thdtw*( s1 -2.0*s2 +s3 )**2 + 0.25*(     s1 -4.0*s2 +3.0*s3 )**2
        b2=thdtw*( s2 -2.0*s3 +s4 )**2 + 0.25*(     s2             -s4 )**2
        b3=thdtw*( s3 -2.0*s4 +s5 )**2 + 0.25*( 3.0*s3 -4.0*s4     +s5 )**2
 
        w1=0.10/(epsilon+b1)**2
        w2=0.60/(epsilon+b2)**2
        w3=0.30/(epsilon+b3)**2
 
        dum(i,j,k)=rru(i,j,k)*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo
      enddo

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

!----------------------------------------------------------------
! Advection in y-direction

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,s1,s2,s3,s4,s5,f1,f2,f3,b1,b2,b3,w1,w2,w3)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        if(rrv(i,j,k).ge.0.0)then
          s1=s(i,j-3,k)
          s2=s(i,j-2,k)
          s3=s(i,j-1,k)
          s4=s(i,j  ,k)
          s5=s(i,j+1,k)
        else
          s1=s(i,j+2,k)
          s2=s(i,j+1,k)
          s3=s(i,j  ,k)
          s4=s(i,j-1,k)
          s5=s(i,j-2,k)
        endif

        f1=( 2.0*s1 -7.0*s2 +11.0*s3 )*onedsix
        f2=(    -s2 +5.0*s3  +2.0*s4 )*onedsix
        f3=( 2.0*s3 +5.0*s4      -s5 )*onedsix

        b1=thdtw*( s1 -2.0*s2 +s3 )**2 + 0.25*(     s1 -4.0*s2 +3.0*s3 )**2
        b2=thdtw*( s2 -2.0*s3 +s4 )**2 + 0.25*(     s2             -s4 )**2
        b3=thdtw*( s3 -2.0*s4 +s5 )**2 + 0.25*( 3.0*s3 -4.0*s4     +s5 )**2
 
        w1=0.10/(epsilon+b1)**2
        w2=0.60/(epsilon+b2)**2
        w3=0.30/(epsilon+b3)**2
 
        dum(i,j,k)=rrv(i,j,k)*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo
      enddo

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
! Advection in z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k,s1,s2,s3,s4,s5,f1,f2,f3,b1,b2,b3,w1,w2,w3)
      do k=4,nk-2
      do j=1,nj
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          s1=s(i,j,k-3)
          s2=s(i,j,k-2)
          s3=s(i,j,k-1)
          s4=s(i,j,k  )
          s5=s(i,j,k+1)
        else
          s1=s(i,j,k+2)
          s2=s(i,j,k+1)
          s3=s(i,j,k  )
          s4=s(i,j,k-1)
          s5=s(i,j,k-2)
        endif

        f1=( 2.0*s1 -7.0*s2 +11.0*s3 )*onedsix
        f2=(    -s2 +5.0*s3  +2.0*s4 )*onedsix
        f3=( 2.0*s3 +5.0*s4      -s5 )*onedsix

        b1=thdtw*( s1 -2.0*s2 +s3 )**2 + 0.25*(     s1 -4.0*s2 +3.0*s3 )**2
        b2=thdtw*( s2 -2.0*s3 +s4 )**2 + 0.25*(     s2             -s4 )**2
        b3=thdtw*( s3 -2.0*s4 +s5 )**2 + 0.25*( 3.0*s3 -4.0*s4     +s5 )**2
 
        w1=0.10/(epsilon+b1)**2
        w2=0.60/(epsilon+b2)**2
        w3=0.30/(epsilon+b3)**2
 
        dum(i,j,k)=rrw(i,j,k)*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
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


