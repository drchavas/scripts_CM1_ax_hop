!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bc2d(s)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie,jb:je) :: s

      integer i,j

!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
      if(wbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j)
        do j=jb,je
          s(   0,j)=s(ni  ,j)
          s(  -1,j)=s(ni-1,j)
          s(  -2,j)=s(ni-2,j)
        enddo
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j)
        do j=jb,je
          s(   0,j)=s( 1,j)
          s(  -1,j)=s( 1,j)
          s(  -2,j)=s( 1,j)
        enddo
#ifdef MPI
      elseif(ibw.eq.1.and.wbc.eq.3)then
#else
      elseif(wbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j)
        do j=jb,je
          s(   0,j)=s(   1,j)
          s(  -1,j)=s(   2,j)
          s(  -2,j)=s(   3,j)
        enddo
      endif

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
      if(ebc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j)
        do j=jb,je
          s(ni+1,j)=s(   1,j)
          s(ni+2,j)=s(   2,j)
          s(ni+3,j)=s(   3,j)
        enddo
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j)
        do j=jb,je
          s(ni+1,j)=s(ni,j)
          s(ni+2,j)=s(ni,j)
          s(ni+3,j)=s(ni,j)
        enddo
#ifdef MPI
      elseif(ibe.eq.1.and.ebc.eq.3)then
#else
      elseif(ebc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j)
        do j=jb,je
          s(ni+1,j)=s(  ni,j)
          s(ni+2,j)=s(ni-1,j)
          s(ni+3,j)=s(ni-2,j)
        enddo
      endif

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
      if(sbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i)
        do i=ib,ie
          s(i,   0)=s(i,nj  )
          s(i,  -1)=s(i,nj-1)
          s(i,  -2)=s(i,nj-2)
        enddo
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i)
        do i=ib,ie
          s(i,   0)=s(i, 1)
          s(i,  -1)=s(i, 1)
          s(i,  -2)=s(i, 1)
        enddo
#ifdef MPI
      elseif(ibs.eq.1.and.sbc.eq.3)then
#else
      elseif(sbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i)
        do i=ib,ie
          s(i,   0)=s(i, 1)
          s(i,  -1)=s(i, 2)
          s(i,  -2)=s(i, 3)
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

#ifndef MPI
      if(nbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i)
        do i=ib,ie
          s(i,nj+1)=s(i,   1)
          s(i,nj+2)=s(i,   2)
          s(i,nj+3)=s(i,   3)
        enddo
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i)
        do i=ib,ie
          s(i,nj+1)=s(i,nj)
          s(i,nj+2)=s(i,nj)
          s(i,nj+3)=s(i,nj)
        enddo
#ifdef MPI
      elseif(ibn.eq.1.and.nbc.eq.3)then
#else
      elseif(nbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i)
        do i=ib,ie
          s(i,nj+1)=s(i,nj  )
          s(i,nj+2)=s(i,nj-1)
          s(i,nj+3)=s(i,nj-2)
        enddo
      endif

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcs(s)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: s

      integer i,j,k

!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
      if(wbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=jb,je
          s(   0,j,k)=s(ni  ,j,k)
          s(  -1,j,k)=s(ni-1,j,k)
          s(  -2,j,k)=s(ni-2,j,k)
        enddo
        enddo
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          s(   0,j,k)=s( 1,j,k)
          s(  -1,j,k)=s( 1,j,k)
          s(  -2,j,k)=s( 1,j,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibw.eq.1.and.wbc.eq.3)then
#else
      elseif(wbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          s(   0,j,k)=s(   1,j,k)
          s(  -1,j,k)=s(   2,j,k)
          s(  -2,j,k)=s(   3,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
      if(ebc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=jb,je
          s(ni+1,j,k)=s(   1,j,k)
          s(ni+2,j,k)=s(   2,j,k)
          s(ni+3,j,k)=s(   3,j,k)
        enddo
        enddo
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          s(ni+1,j,k)=s(ni,j,k)
          s(ni+2,j,k)=s(ni,j,k)
          s(ni+3,j,k)=s(ni,j,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibe.eq.1.and.ebc.eq.3)then
#else
      elseif(ebc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          s(ni+1,j,k)=s(  ni,j,k)
          s(ni+2,j,k)=s(ni-1,j,k)
          s(ni+3,j,k)=s(ni-2,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
      if(sbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=ib,ie
          s(i,   0,k)=s(i,nj  ,k)
          s(i,  -1,k)=s(i,nj-1,k)
          s(i,  -2,k)=s(i,nj-2,k)
        enddo
        enddo
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+1
          s(i,   0,k)=s(i, 1,k)
          s(i,  -1,k)=s(i, 1,k)
          s(i,  -2,k)=s(i, 1,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibs.eq.1.and.sbc.eq.3)then
#else
      elseif(sbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+1
          s(i,   0,k)=s(i, 1,k)
          s(i,  -1,k)=s(i, 2,k)
          s(i,  -2,k)=s(i, 3,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

#ifndef MPI
      if(nbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=ib,ie
          s(i,nj+1,k)=s(i,   1,k)
          s(i,nj+2,k)=s(i,   2,k)
          s(i,nj+3,k)=s(i,   3,k)
        enddo
        enddo
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+1
          s(i,nj+1,k)=s(i,nj,k)
          s(i,nj+2,k)=s(i,nj,k)
          s(i,nj+3,k)=s(i,nj,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibn.eq.1.and.nbc.eq.3)then
#else
      elseif(nbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+1
          s(i,nj+1,k)=s(i,nj  ,k)
          s(i,nj+2,k)=s(i,nj-1,k)
          s(i,nj+3,k)=s(i,nj-2,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  top/bottom boundary condition

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=0,nj+1
      do i=0,ni+1
        s(i,j, 0)=s(i,j,1)
        s(i,j,nk+1)=s(i,j,nk  )
      enddo
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcp(p)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: p

      integer i,j,k

!-----------------------------------------------------------------------
!  west boundary condition

      if(wbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          p(   0,j,k)=p(ni  ,j,k)
          p(  -1,j,k)=p(ni-1,j,k)
          p(  -2,j,k)=p(ni-2,j,k)
        enddo
        enddo
      elseif(wbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          p(  -1,j,k)=p( 0,j,k)
          p(  -2,j,k)=p( 0,j,k)
        enddo
        enddo
      elseif(wbc.eq.3)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          p(   0,j,k)=p(   1,j,k)
          p(  -1,j,k)=p(   2,j,k)
          p(  -2,j,k)=p(   3,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  east boundary condition

      if(ebc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          p(ni+1,j,k)=p(   1,j,k)
          p(ni+2,j,k)=p(   2,j,k)
          p(ni+3,j,k)=p(   3,j,k)
        enddo
        enddo
      elseif(ebc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          p(ni+2,j,k)=p(ni+1,j,k)
          p(ni+3,j,k)=p(ni+1,j,k)
        enddo
        enddo
      elseif(ebc.eq.3)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          p(ni+1,j,k)=p(  ni,j,k)
          p(ni+2,j,k)=p(ni-1,j,k)
          p(ni+3,j,k)=p(ni-2,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  south boundary condition

      if(sbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=ib,ie
          p(i,   0,k)=p(i,nj  ,k)
          p(i,  -1,k)=p(i,nj-1,k)
          p(i,  -2,k)=p(i,nj-2,k)
        enddo
        enddo
      elseif(sbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+1
          p(i,  -1,k)=p(i, 0,k)
          p(i,  -2,k)=p(i, 0,k)
        enddo
        enddo
      elseif(sbc.eq.3)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+1
          p(i,   0,k)=p(i, 1,k)
          p(i,  -1,k)=p(i, 2,k)
          p(i,  -2,k)=p(i, 3,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

      if(nbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=ib,ie
          p(i,nj+1,k)=p(i,   1,k)
          p(i,nj+2,k)=p(i,   2,k)
          p(i,nj+3,k)=p(i,   3,k)
        enddo
        enddo
      elseif(nbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+1
          p(i,nj+2,k)=p(i,nj+1,k)
          p(i,nj+3,k)=p(i,nj+1,k)
        enddo
        enddo
      elseif(nbc.eq.3)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+1
          p(i,nj+1,k)=p(i,nj  ,k)
          p(i,nj+2,k)=p(i,nj-1,k)
          p(i,nj+3,k)=p(i,nj-2,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  top/bottom boundary condition

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=jb,je
      do i=ib,ie
        p(i,j,   0)=p(i,j, 1)
        p(i,j,nk+1)=p(i,j,nk)
      enddo
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcu(u)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'  

      real, dimension(ib:ie+1,jb:je,kb:ke) :: u

      integer i,j,k  

!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
      if(wbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=jb,je
          u(   0,j,k)=u(ni  ,j,k)
          u(  -1,j,k)=u(ni-1,j,k)
          u(  -2,j,k)=u(ni-2,j,k)
        enddo
        enddo
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          u(   0,j,k)=u(1,j,k)
          u(  -1,j,k)=u(1,j,k)
          u(  -2,j,k)=u(1,j,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibw.eq.1.and.wbc.eq.3)then
#else
      elseif(wbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          u(   1,j,k)=0.
          u(   0,j,k)=-u(   2,j,k)
          u(  -1,j,k)=-u(   3,j,k)
          u(  -2,j,k)=-u(   4,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
      if(ebc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=jb,je
          u(ni+2,j,k)=u(2,j,k)
          u(ni+3,j,k)=u(3,j,k)
          u(ni+4,j,k)=u(4,j,k)
        enddo
        enddo
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          u(ni+2,j,k)=u(ni+1,j,k)
          u(ni+3,j,k)=u(ni+1,j,k)
          u(ni+4,j,k)=u(ni+1,j,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibe.eq.1.and.ebc.eq.3)then
#else
      elseif(ebc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+1
          u(ni+1,j,k)=0.
          u(ni+2,j,k)=-u(ni  ,j,k)
          u(ni+3,j,k)=-u(ni-1,j,k)
          u(ni+4,j,k)=-u(ni-2,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
      if(sbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=ib,ie+1
          u(i,   0,k)=u(i,nj  ,k)
          u(i,  -1,k)=u(i,nj-1,k)
          u(i,  -2,k)=u(i,nj-2,k)
        enddo
        enddo
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+2
          u(i,   0,k)=u(i, 1,k)
          u(i,  -1,k)=u(i, 1,k)
          u(i,  -2,k)=u(i, 1,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibs.eq.1.and.sbc.eq.3)then
#else
      elseif(sbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+2
          u(i,   0,k)=u(i, 1,k)
          u(i,  -1,k)=u(i, 2,k)
          u(i,  -2,k)=u(i, 3,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

#ifndef MPI
      if(nbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=ib,ie+1
          u(i,nj+1,k)=u(i,   1,k)
          u(i,nj+2,k)=u(i,   2,k)
          u(i,nj+3,k)=u(i,   3,k)
        enddo
        enddo
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+2
          u(i,nj+1,k)=u(i,nj,k)
          u(i,nj+2,k)=u(i,nj,k)
          u(i,nj+3,k)=u(i,nj,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibn.eq.1.and.nbc.eq.3)then
#else
      elseif(nbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+2
          u(i,nj+1,k)=u(i,nj  ,k)
          u(i,nj+2,k)=u(i,nj-1,k)
          u(i,nj+3,k)=u(i,nj-2,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  top/bottom boundary condition

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=0,nj+1
      do i=0,ni+2
        u(i,j, 0)=u(i,j,1)
        u(i,j,nk+1)=u(i,j,nk  )
      enddo
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcv(v)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v
 
      integer i,j,k

!-----------------------------------------------------------------------
!  south boundary condition
 
#ifndef MPI
      if(sbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=ib,ie
          v(i,   0,k)=v(i,nj  ,k)
          v(i,  -1,k)=v(i,nj-1,k)
          v(i,  -2,k)=v(i,nj-2,k)
        enddo
        enddo
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+1
          v(i,   0,k)=v(i,1,k)
          v(i,  -1,k)=v(i,1,k)
          v(i,  -2,k)=v(i,1,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibs.eq.1.and.sbc.eq.3)then
#else
      elseif(sbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+1
          v(i,   1,k)=0.
          v(i,   0,k)=-v(i,2,k)
          v(i,  -1,k)=-v(i,3,k)
          v(i,  -2,k)=-v(i,4,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

      IF(axisymm.eq.1)THEN
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=ib,ie
          v(i,2,k)=v(i,1,k)
        enddo
        enddo
      ENDIF

#ifndef MPI
      if(nbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=ib,ie
          v(i,nj+2,k)=v(i,   2,k)
          v(i,nj+3,k)=v(i,   3,k)
          v(i,nj+4,k)=v(i,   4,k)
        enddo
        enddo
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+1
          v(i,nj+2,k)=v(i,nj+1,k)
          v(i,nj+3,k)=v(i,nj+1,k)
          v(i,nj+4,k)=v(i,nj+1,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibn.eq.1.and.nbc.eq.3)then
#else
      elseif(nbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=0,ni+1
          v(i,nj+1,k)=0.
          v(i,nj+2,k)=-v(i,nj  ,k)
          v(i,nj+3,k)=-v(i,nj-1,k)
          v(i,nj+4,k)=-v(i,nj-2,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
      if(wbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=jb,je+1
          v(   0,j,k)=v(ni  ,j,k)
          v(  -1,j,k)=v(ni-1,j,k)
          v(  -2,j,k)=v(ni-2,j,k)
        enddo
        enddo
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+2
          v(   0,j,k)=v( 1,j,k)
          v(  -1,j,k)=v( 1,j,k)
          v(  -2,j,k)=v( 1,j,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibw.eq.1.and.wbc.eq.3)then
#else
      elseif(wbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+2
          v(   0,j,k)=v( 1,j,k)
          v(  -1,j,k)=v( 2,j,k)
          v(  -2,j,k)=v( 3,j,k)
        enddo
        enddo
      endif

      IF(axisymm.eq.1)THEN
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+2
          v( 0,j,k) = -v(1,j,k)
          v(-1,j,k) = -v(2,j,k)
          v(-2,j,k) = -v(3,j,k)
        enddo
        enddo
      ENDIF

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
      if(ebc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=jb,je+1
          v(ni+1,j,k)=v(   1,j,k)
          v(ni+2,j,k)=v(   2,j,k)
          v(ni+3,j,k)=v(   3,j,k)
        enddo
        enddo
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+2
          v(ni+1,j,k)=v(ni,j,k)
          v(ni+2,j,k)=v(ni,j,k)
          v(ni+3,j,k)=v(ni,j,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibe.eq.1.and.ebc.eq.3)then
#else
      elseif(ebc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=0,nj+2
          v(ni+1,j,k)=v(ni  ,j,k)
          v(ni+2,j,k)=v(ni-1,j,k)
          v(ni+3,j,k)=v(ni-2,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  top/bottom boundary condition

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=0,nj+2
      do i=0,ni+1
        v(i,j, 0)=v(i,j,1)
        v(i,j,nk+1)=v(i,j,nk  )
      enddo
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcw(w,flag)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w
      integer flag
 
      integer i,j,k

!-----------------------------------------------------------------------
!  west boundary condition
 
#ifndef MPI
      if(wbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=2,nk
        do j=jb,je
          w(   0,j,k)=w(ni  ,j,k)
          w(  -1,j,k)=w(ni-1,j,k)
          w(  -2,j,k)=w(ni-2,j,k)
        enddo
        enddo
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=2,nk
        do j=0,nj+1
          w(   0,j,k)=w( 1,j,k)
          w(  -1,j,k)=w( 1,j,k)
          w(  -2,j,k)=w( 1,j,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibw.eq.1.and.wbc.eq.3)then
#else
      elseif(wbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=2,nk
        do j=0,nj+1
          w(   0,j,k)=w( 1,j,k)
          w(  -1,j,k)=w( 2,j,k)
          w(  -2,j,k)=w( 3,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
      if(ebc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=2,nk
        do j=jb,je
          w(ni+1,j,k)=w(   1,j,k)
          w(ni+2,j,k)=w(   2,j,k)
          w(ni+3,j,k)=w(   3,j,k)
        enddo
        enddo
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=2,nk
        do j=0,nj+1
          w(ni+1,j,k)=w(ni,j,k)
          w(ni+2,j,k)=w(ni,j,k)
          w(ni+3,j,k)=w(ni,j,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibe.eq.1.and.ebc.eq.3)then
#else
      elseif(ebc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=2,nk
        do j=0,nj+1
          w(ni+1,j,k)=w(ni  ,j,k)
          w(ni+2,j,k)=w(ni-1,j,k)
          w(ni+3,j,k)=w(ni-2,j,k)
        enddo 
        enddo 
      endif

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
      if(sbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=2,nk
        do i=ib,ie
          w(i,   0,k)=w(i,nj  ,k)
          w(i,  -1,k)=w(i,nj-1,k)
          w(i,  -2,k)=w(i,nj-2,k)
        enddo
        enddo
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=2,nk
        do i=0,ni+1
          w(i,   0,k)=w(i, 1,k)
          w(i,  -1,k)=w(i, 1,k)
          w(i,  -2,k)=w(i, 1,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibs.eq.1.and.sbc.eq.3)then
#else
      elseif(sbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=2,nk
        do i=0,ni+1
          w(i,   0,k)=w(i, 1,k)
          w(i,  -1,k)=w(i, 2,k)
          w(i,  -2,k)=w(i, 3,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

#ifndef MPI
      if(nbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=2,nk
        do i=ib,ie
          w(i,nj+1,k)=w(i,   1,k)
          w(i,nj+2,k)=w(i,   2,k)
          w(i,nj+3,k)=w(i,   3,k)
        enddo
        enddo
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=2,nk
        do i=0,ni+1
          w(i,nj+1,k)=w(i,nj,k)
          w(i,nj+2,k)=w(i,nj,k)
          w(i,nj+3,k)=w(i,nj,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibn.eq.1.and.nbc.eq.3)then
#else
      elseif(nbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=2,nk
        do i=0,ni+1
          w(i,nj+1,k)=w(i,nj  ,k)
          w(i,nj+2,k)=w(i,nj-1,k)
          w(i,nj+3,k)=w(i,nj-2,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  top/bottom boundary condition

    IF(flag.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=0,nj+1
      do i=0,ni+1
        w(i,j, 0)=-w(i,j,2)
        w(i,j, 1)=0.0
        w(i,j,nk+1)=0.0
        w(i,j,nk+2)=-w(i,j,nk  )
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcwsfc(dzdx,dzdy,u,v,w)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(itb:ite,jtb:jte) :: dzdx,dzdy
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w

      integer i,j

!-----------------------------------------------------------------------

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=0,nj+1
      do i=0,ni+1
        w(i,j,1)=0.5*( (u(i,j,1)+u(i+1,j,1))*dzdx(i,j)    &
                      +(v(i,j,1)+v(i,j+1,1))*dzdy(i,j) )
      enddo
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bct(t)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie,jb:je,kb:ke+1) :: t

      integer i,j,k

!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
      if(wbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk+1
        do j=jb,je
          t(   0,j,k)=t(ni  ,j,k)
          t(  -1,j,k)=t(ni-1,j,k)
          t(  -2,j,k)=t(ni-2,j,k)
        enddo
        enddo
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk+1
        do j=0,nj+1
          t(   0,j,k)=t( 1,j,k)
          t(  -1,j,k)=t( 1,j,k)
          t(  -2,j,k)=t( 1,j,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibw.eq.1.and.wbc.eq.3)then
#else
      elseif(wbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk+1
        do j=0,nj+1
          t(   0,j,k)=t(   1,j,k)
          t(  -1,j,k)=t(   2,j,k)
          t(  -2,j,k)=t(   3,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
      if(ebc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk+1
        do j=jb,je
          t(ni+1,j,k)=t(   1,j,k)
          t(ni+2,j,k)=t(   2,j,k)
          t(ni+3,j,k)=t(   3,j,k)
        enddo
        enddo
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk+1
        do j=0,nj+1
          t(ni+1,j,k)=t(ni,j,k)
          t(ni+2,j,k)=t(ni,j,k)
          t(ni+3,j,k)=t(ni,j,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibe.eq.1.and.ebc.eq.3)then
#else
      elseif(ebc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk+1
        do j=0,nj+1
          t(ni+1,j,k)=t(  ni,j,k)
          t(ni+2,j,k)=t(ni-1,j,k)
          t(ni+3,j,k)=t(ni-2,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
      if(sbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk+1
        do i=ib,ie
          t(i,   0,k)=t(i,nj  ,k)
          t(i,  -1,k)=t(i,nj-1,k)
          t(i,  -2,k)=t(i,nj-2,k)
        enddo
        enddo
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk+1
        do i=0,ni+1
          t(i,   0,k)=t(i, 1,k)
          t(i,  -1,k)=t(i, 1,k)
          t(i,  -2,k)=t(i, 1,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibs.eq.1.and.sbc.eq.3)then
#else
      elseif(sbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk+1
        do i=0,ni+1
          t(i,   0,k)=t(i, 1,k)
          t(i,  -1,k)=t(i, 2,k)
          t(i,  -2,k)=t(i, 3,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

#ifndef MPI
      if(nbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk+1
        do i=ib,ie
          t(i,nj+1,k)=t(i,   1,k)
          t(i,nj+2,k)=t(i,   2,k)
          t(i,nj+3,k)=t(i,   3,k)
        enddo
        enddo
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk+1
        do i=0,ni+1
          t(i,nj+1,k)=t(i,nj,k)
          t(i,nj+2,k)=t(i,nj,k)
          t(i,nj+3,k)=t(i,nj,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibn.eq.1.and.nbc.eq.3)then
#else
      elseif(nbc.eq.3)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk+1
        do i=0,ni+1
          t(i,nj+1,k)=t(i,nj  ,k)
          t(i,nj+2,k)=t(i,nj-1,k)
          t(i,nj+3,k)=t(i,nj-2,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  top/bottom boundary condition

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=0,nj+1
      do i=0,ni+1
        t(i,j, 0)=t(i,j,1)
        t(i,j,nk+2)=t(i,j,nk+1)
      enddo
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


#ifdef MPI
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcs2(s)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real s(ib:ie,jb:je,kb:ke)

      integer i,j,k

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            s(  -2,j,k)=s( 1,j,k)
            s(  -1,j,k)=s( 1,j,k)
            s(   0,j,k)=s( 1,j,k)
          enddo
        endif

        if(patchnww)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            s(  -2,j,k)=s( 1,j,k)
            s(  -1,j,k)=s( 1,j,k)
            s(   0,j,k)=s( 1,j,k)
          enddo
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            s(ni+1,j,k)=s(ni,j,k)
            s(ni+2,j,k)=s(ni,j,k)
            s(ni+3,j,k)=s(ni,j,k)
          enddo
        endif

        if(patchnee)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            s(ni+1,j,k)=s(ni,j,k)
            s(ni+2,j,k)=s(ni,j,k)
            s(ni+3,j,k)=s(ni,j,k)
          enddo
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            s(i,  -2,k)=s(i, 1,k)
            s(i,  -1,k)=s(i, 1,k)
            s(i,   0,k)=s(i, 1,k)
          enddo
        endif

        if(patchses)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            s(i,  -2,k)=s(i, 1,k)
            s(i,  -1,k)=s(i, 1,k)
            s(i,   0,k)=s(i, 1,k)
          enddo
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            s(i,nj+1,k)=s(i,nj,k)
            s(i,nj+2,k)=s(i,nj,k)
            s(i,nj+3,k)=s(i,nj,k)
          enddo
        endif

        if(patchnen)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            s(i,nj+1,k)=s(i,nj,k)
            s(i,nj+2,k)=s(i,nj,k)
            s(i,nj+3,k)=s(i,nj,k)
          enddo
        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bct2(t)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real t(ib:ie,jb:je,kb:ke+1)

      integer i,j,k

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            t(  -2,j,k)=t( 1,j,k)
            t(  -1,j,k)=t( 1,j,k)
            t(   0,j,k)=t( 1,j,k)
          enddo
        endif

        if(patchnww)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            t(  -2,j,k)=t( 1,j,k)
            t(  -1,j,k)=t( 1,j,k)
            t(   0,j,k)=t( 1,j,k)
          enddo
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            t(ni+1,j,k)=t(ni,j,k)
            t(ni+2,j,k)=t(ni,j,k)
            t(ni+3,j,k)=t(ni,j,k)
          enddo
        endif

        if(patchnee)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            t(ni+1,j,k)=t(ni,j,k)
            t(ni+2,j,k)=t(ni,j,k)
            t(ni+3,j,k)=t(ni,j,k)
          enddo
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            t(i,  -2,k)=t(i, 1,k)
            t(i,  -1,k)=t(i, 1,k)
            t(i,   0,k)=t(i, 1,k)
          enddo
        endif

        if(patchses)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            t(i,  -2,k)=t(i, 1,k)
            t(i,  -1,k)=t(i, 1,k)
            t(i,   0,k)=t(i, 1,k)
          enddo
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            t(i,nj+1,k)=t(i,nj,k)
            t(i,nj+2,k)=t(i,nj,k)
            t(i,nj+3,k)=t(i,nj,k)
          enddo
        endif

        if(patchnen)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            t(i,nj+1,k)=t(i,nj,k)
            t(i,nj+2,k)=t(i,nj,k)
            t(i,nj+3,k)=t(i,nj,k)
          enddo
        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcu2(u)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real u(ib:ie+1,jb:je,kb:ke)

      integer i,j,k

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            u(  -2,j,k)=u( 1,j,k)
            u(  -1,j,k)=u( 1,j,k)
            u(   0,j,k)=u( 1,j,k)
          enddo
        endif

        if(patchnww)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            u(  -2,j,k)=u( 1,j,k)
            u(  -1,j,k)=u( 1,j,k)
            u(   0,j,k)=u( 1,j,k)
          enddo
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            u(ni+2,j,k)=u(ni+1,j,k)
            u(ni+3,j,k)=u(ni+1,j,k)
            u(ni+4,j,k)=u(ni+1,j,k)
          enddo
        endif

        if(patchnee)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            u(ni+2,j,k)=u(ni+1,j,k)
            u(ni+3,j,k)=u(ni+1,j,k)
            u(ni+4,j,k)=u(ni+1,j,k)
          enddo
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            u(i,  -2,k)=u(i, 1,k)
            u(i,  -1,k)=u(i, 1,k)
            u(i,   0,k)=u(i, 1,k)
          enddo
        endif

        if(patchses)then
          i=ni+2
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            u(i,  -2,k)=u(i, 1,k)
            u(i,  -1,k)=u(i, 1,k)
            u(i,   0,k)=u(i, 1,k)
          enddo
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            u(i,nj+1,k)=u(i,nj,k)
            u(i,nj+2,k)=u(i,nj,k)
            u(i,nj+3,k)=u(i,nj,k)
          enddo
        endif

        if(patchnen)then
          i=ni+2
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            u(i,nj+1,k)=u(i,nj,k)
            u(i,nj+2,k)=u(i,nj,k)
            u(i,nj+3,k)=u(i,nj,k)
          enddo
        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcv2(v)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real v(ib:ie,jb:je+1,kb:ke)

      integer i,j,k

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            v(  -2,j,k)=v( 1,j,k)
            v(  -1,j,k)=v( 1,j,k)
            v(   0,j,k)=v( 1,j,k)
          enddo
        endif

        if(patchnww)then
          j=nj+2
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            v(  -2,j,k)=v( 1,j,k)
            v(  -1,j,k)=v( 1,j,k)
            v(   0,j,k)=v( 1,j,k)
          enddo
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            v(ni+1,j,k)=v(ni,j,k)
            v(ni+2,j,k)=v(ni,j,k)
            v(ni+3,j,k)=v(ni,j,k)
          enddo
        endif

        if(patchnee)then
          j=nj+2
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            v(ni+1,j,k)=v(ni,j,k)
            v(ni+2,j,k)=v(ni,j,k)
            v(ni+3,j,k)=v(ni,j,k)
          enddo
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            v(i,  -2,k)=v(i, 1,k)
            v(i,  -1,k)=v(i, 1,k)
            v(i,   0,k)=v(i, 1,k)
          enddo
        endif

        if(patchses)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            v(i,  -2,k)=v(i, 1,k)
            v(i,  -1,k)=v(i, 1,k)
            v(i,   0,k)=v(i, 1,k)
          enddo
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            v(i,nj+2,k)=v(i,nj+1,k)
            v(i,nj+3,k)=v(i,nj+1,k)
            v(i,nj+4,k)=v(i,nj+1,k)
          enddo
        endif

        if(patchnen)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk
            v(i,nj+2,k)=v(i,nj+1,k)
            v(i,nj+3,k)=v(i,nj+1,k)
            v(i,nj+4,k)=v(i,nj+1,k)
          enddo
        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcw2(w)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real w(ib:ie,jb:je,kb:ke+1)

      integer i,j,k

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            w(  -2,j,k)=w( 1,j,k)
            w(  -1,j,k)=w( 1,j,k)
            w(   0,j,k)=w( 1,j,k)
          enddo
        endif

        if(patchnww)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            w(  -2,j,k)=w( 1,j,k)
            w(  -1,j,k)=w( 1,j,k)
            w(   0,j,k)=w( 1,j,k)
          enddo
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            w(ni+1,j,k)=w(ni,j,k)
            w(ni+2,j,k)=w(ni,j,k)
            w(ni+3,j,k)=w(ni,j,k)
          enddo
        endif

        if(patchnee)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            w(ni+1,j,k)=w(ni,j,k)
            w(ni+2,j,k)=w(ni,j,k)
            w(ni+3,j,k)=w(ni,j,k)
          enddo
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            w(i,  -2,k)=w(i, 1,k)
            w(i,  -1,k)=w(i, 1,k)
            w(i,   0,k)=w(i, 1,k)
          enddo
        endif

        if(patchses)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            w(i,  -2,k)=w(i, 1,k)
            w(i,  -1,k)=w(i, 1,k)
            w(i,   0,k)=w(i, 1,k)
          enddo
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            w(i,nj+1,k)=w(i,nj,k)
            w(i,nj+2,k)=w(i,nj,k)
            w(i,nj+3,k)=w(i,nj,k)
          enddo
        endif

        if(patchnen)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(k)
          do k=1,nk+1
            w(i,nj+1,k)=w(i,nj,k)
            w(i,nj+2,k)=w(i,nj,k)
            w(i,nj+3,k)=w(i,nj,k)
          enddo
        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcs2_2d(s)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real s(ib:ie,jb:je)

      integer i,j

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
          s(  -2,j)=s( 1,j)
          s(  -1,j)=s( 1,j)
          s(   0,j)=s( 1,j)
        endif

        if(patchnww)then
          j=nj+1
          s(  -2,j)=s( 1,j)
          s(  -1,j)=s( 1,j)
          s(   0,j)=s( 1,j)
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          j=0
          s(ni+1,j)=s(ni,j)
          s(ni+2,j)=s(ni,j)
          s(ni+3,j)=s(ni,j)
        endif

        if(patchnee)then
          j=nj+1
          s(ni+1,j)=s(ni,j)
          s(ni+2,j)=s(ni,j)
          s(ni+3,j)=s(ni,j)
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          i=0
          s(i,  -2)=s(i, 1)
          s(i,  -1)=s(i, 1)
          s(i,   0)=s(i, 1)
        endif

        if(patchses)then
          i=ni+1
          s(i,  -2)=s(i, 1)
          s(i,  -1)=s(i, 1)
          s(i,   0)=s(i, 1)
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          i=0
          s(i,nj+1)=s(i,nj)
          s(i,nj+2)=s(i,nj)
          s(i,nj+3)=s(i,nj)
        endif

        if(patchnen)then
          i=ni+1
          s(i,nj+1)=s(i,nj)
          s(i,nj+2)=s(i,nj)
          s(i,nj+3)=s(i,nj)
        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine radbcew(radbcw,radbce,ua)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, dimension(jb:je,kb:ke) :: radbcw,radbce
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
 
      integer j,k
      real cbcw,cbce
 
      if(ibw.eq.1.and.wbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(j,k,cbcw)
        do k=1,nk
        do j=1,nj
          cbcw=ua(1,j,k)-30.0
          if(cbcw.lt.0.0)then
            radbcw(j,k)=cbcw
          else
            radbcw(j,k)=0.
          endif
        enddo
        enddo
      endif
 
      if(ibe.eq.1.and.ebc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(j,k,cbce)
        do k=1,nk
        do j=1,nj
          cbce=ua(ni+1,j,k)+30.0
          if(cbce.gt.0.0)then
            radbce(j,k)=cbce
          else
            radbce(j,k)=0.
          endif
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine radbcns(radbcs,radbcn,va)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
 
      integer i,k
      real cbcs,cbcn

      if(ibs.eq.1.and.sbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(i,k,cbcs)
        do k=1,nk
        do i=1,ni
          cbcs=va(i,1,k)-30.0
          if(cbcs.lt.0.0)then
            radbcs(i,k)=cbcs
          else
            radbcs(i,k)=0.
          endif
        enddo
        enddo
      endif
 
      if(ibn.eq.1.and.nbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(i,k,cbcn)
        do k=1,nk
        do i=1,ni
          cbcn=va(i,nj+1,k)+30.0
          if(cbcn.gt.0.0)then
            radbcn(i,k)=cbcn
          else
            radbcn(i,k)=0.
          endif
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine radbcew4(ruf,radbcw,radbce,u1,u2,dt)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie+1) :: ruf
      real, dimension(jb:je,kb:ke) :: radbcw,radbce
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u1,u2
      real dt

      integer j,k
      real umax,avgw,avge,foo,cbcw,cbce

      umax=350.0

      if(ibw.eq.1.and.wbc.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(j,k,foo,avgw,cbcw)
        do j=1,nj
          avgw=0.
          do k=1,nk
            foo=(u1(3,j,k)-u1(2,j,k))
            cbcw=dx*ruf(2)*(u1(2,j,k)-u2(2,j,k))   &
                   /(dt*(sign(1.e-10,foo)+foo))
            cbcw=max(min(cbcw,0.0),-umax)
            avgw=avgw+cbcw
          enddo
          avgw=avgw/float(nk)
          do k=1,nk
            radbcw(j,k)=avgw
          enddo
        enddo
      endif

      if(ibe.eq.1.and.ebc.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(j,k,foo,avge,cbce)
        do j=1,nj
          avge=0.
          do k=1,nk
            foo=(u1(ni+1-1,j,k)-u1(ni+1-2,j,k))
            cbce=dx*ruf(ni+1-1)*(u1(ni+1-1,j,k)-u2(ni+1-1,j,k))   &
                   /(dt*(sign(1.e-10,foo)+foo))
            cbce=min(max(cbce,0.0),umax)
            avge=avge+cbce
          enddo
          avge=avge/float(nk)
          do k=1,nk
            radbce(j,k)=avge
          enddo
        enddo
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine radbcns4(rvf,radbcs,radbcn,v1,v2,dt)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(jb:je+1) :: rvf
      real, dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v1,v2
      real dt

      integer i,k
      real umax,avgs,avgn,foo,cbcs,cbcn

      umax=350.0

      if(ibs.eq.1.and.sbc.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(i,k,avgs,foo,cbcs)
        do i=1,ni
          avgs=0.
          do k=1,nk
            foo=(v1(i,3,k)-v1(i,2,k))
            cbcs=dy*rvf(2)*(v1(i,2,k)-v2(i,2,k))   &
                   /(dt*(sign(1.e-10,foo)+foo))
            cbcs=max(min(cbcs,0.0),-umax)
            avgs=avgs+cbcs
          enddo
          avgs=avgs/float(nk)
          do k=1,nk
            radbcs(i,k)=avgs
          enddo
        enddo
      endif

      if(ibn.eq.1.and.nbc.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(i,k,avgn,foo,cbcn)
        do i=1,ni
          avgn=0.
          do k=1,nk
            foo=(v1(i,nj+1-1,k)-v1(i,nj+1-2,k))
            cbcn=dy*rvf(nj+1-1)*(v1(i,nj+1-1,k)-v2(i,nj+1-1,k))   &
                   /(dt*(sign(1.e-10,foo)+foo))
            cbcn=min(max(cbcn,0.0),umax)
            avgn=avgn+cbcn
          enddo
          avgn=avgn/float(nk)
          do k=1,nk
            radbcn(i,k)=avgn
          enddo
        enddo
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine restrict_openbc_we(rvh,rmh,rho0,u3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u3d

      integer i,j,k
      real*8 :: fluxout,fluxin,tem,u1,t3
      real*8, dimension(nk) :: temout,temin

!$omp parallel do default(shared)   &
!$omp private(k)
      do k=1,nk
        temout(k) = 0.0d0
        temin(k)  = 0.0d0
      enddo

      if(wbc.eq.2.and.ibw.eq.1)then
        i=1
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          temout(k)=temout(k)-min(0.0,rho0(1,j,k)*u3d(i,j,k)*rvh(j)*rmh(1,j,k))
          temin(k) =temin(k) +max(0.0,rho0(1,j,k)*u3d(i,j,k)*rvh(j)*rmh(1,j,k))
        enddo
        enddo
      endif

      if(ebc.eq.2.and.ibe.eq.1)then
        i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          temout(k)=temout(k)+max(0.0,rho0(ni,j,k)*u3d(i,j,k)*rvh(j)*rmh(ni,j,k))
          temin(k) =temin(k) -min(0.0,rho0(ni,j,k)*u3d(i,j,k)*rvh(j)*rmh(ni,j,k))
        enddo
        enddo
      endif

      fluxout = 0.0d0
      fluxin  = 0.0d0

      do k=1,nk
        fluxout = fluxout + temout(k)
        fluxin  = fluxin  + temin(k)
      enddo

#ifdef MPI
      tem=0.0d0
      call MPI_ALLREDUCE(fluxout,tem,1,MPI_DOUBLE_PRECISION,MPI_SUM,   &
                         MPI_COMM_WORLD,ierr)
      fluxout=tem
      tem=0.0d0
      call MPI_ALLREDUCE(fluxin,tem,1,MPI_DOUBLE_PRECISION,MPI_SUM,   &
                         MPI_COMM_WORLD,ierr)
      fluxin=tem
#endif

      t3=(fluxin+1.0d-20)/(fluxout+1.0d-20)

      if(wbc.eq.2.and.ibw.eq.1)then
        i=1
!$omp parallel do default(shared)   &
!$omp private(j,k,u1)
        do k=1,nk
        do j=1,nj
          u1=rho0(1,j,k)*u3d(i,j,k)
          if(u1.lt.0.0)then
            u3d(i,j,k)=u1*t3/rho0(1,j,k)
          endif
        enddo
        enddo
      endif

      if(ebc.eq.2.and.ibe.eq.1)then
        i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k,u1)
        do k=1,nk
        do j=1,nj
          u1=rho0(ni,j,k)*u3d(i,j,k)
          if(u1.gt.0.0)then
            u3d(i,j,k)=u1*t3/rho0(ni,j,k)
          endif
        enddo
        enddo
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine restrict_openbc_sn(ruh,rmh,rho0,v3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      real, dimension(ib:ie) :: ruh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v3d

      integer i,j,k
      real*8 :: fluxout,fluxin,tem,u1,t3
      real*8, dimension(nk) :: temout,temin

!$omp parallel do default(shared)   &
!$omp private(k)
      do k=1,nk
        temout(k) = 0.0d0
        temin(k)  = 0.0d0
      enddo

      if(sbc.eq.2.and.ibs.eq.1)then
        j=1
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          temout(k)=temout(k)-min(0.0,rho0(i,1,k)*v3d(i,j,k)*ruh(i)*rmh(i,1,k))
          temin(k) =temin(k) +max(0.0,rho0(i,1,k)*v3d(i,j,k)*ruh(i)*rmh(i,1,k))
        enddo
        enddo
      endif

      if(nbc.eq.2.and.ibn.eq.1)then
        j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          temout(k)=temout(k)+max(0.0,rho0(i,nj,k)*v3d(i,j,k)*ruh(i)*rmh(i,nj,k))
          temin(k) =temin(k) -min(0.0,rho0(i,nj,k)*v3d(i,j,k)*ruh(i)*rmh(i,nj,k))
        enddo
        enddo
      endif

      fluxout = 0.0d0
      fluxin  = 0.0d0

      do k=1,nk
        fluxout = fluxout + temout(k)
        fluxin  = fluxin  + temin(k)
      enddo

#ifdef MPI
      tem=0.0d0
      call MPI_ALLREDUCE(fluxout,tem,1,MPI_DOUBLE_PRECISION,MPI_SUM,   &
                         MPI_COMM_WORLD,ierr)
      fluxout=tem
      tem=0.0d0
      call MPI_ALLREDUCE(fluxin,tem,1,MPI_DOUBLE_PRECISION,MPI_SUM,   &
                         MPI_COMM_WORLD,ierr)
      fluxin=tem
#endif

      t3=(fluxin+1.0d-20)/(fluxout+1.0d-20)

      if(sbc.eq.2.and.ibs.eq.1)then
        j=1
!$omp parallel do default(shared)   &
!$omp private(i,k,u1)
        do k=1,nk
        do i=1,ni
          u1=rho0(i,1,k)*v3d(i,j,k)
          if(u1.lt.0.0)then
            v3d(i,j,k)=u1*t3/rho0(i,1,k)
          endif
        enddo
        enddo
      endif

      if(nbc.eq.2.and.ibn.eq.1)then
        j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k,u1)
        do k=1,nk
        do i=1,ni
          u1=rho0(i,nj,k)*v3d(i,j,k)
          if(u1.gt.0.0)then
            v3d(i,j,k)=u1*t3/rho0(i,nj,k)
          endif
        enddo
        enddo
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine movesfc(rmax,dt,uh,vh,sfc,s,dum1,dum2)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real :: rmax,dt
      real, dimension(ib:ie) :: uh
      real, dimension(jb:je) :: vh
      real, dimension(ib:ie,jb:je) :: sfc,s,dum1,dum2

      integer i,j,nrk
      real s1,s2,s3,s4,s5,f1,f2,f3,b1,b2,b3,w1,w2,w3,tem

      real epsilon,onedsix,thdtw
      parameter(epsilon=1.0e-8)
      parameter(onedsix=1.0/6.0)
      parameter(thdtw=13.0/12.0)

#ifdef MPI
      real, dimension(3,nj) :: west,newwest,east,neweast
      real, dimension(ni,3) :: south,newsouth,north,newnorth
      integer reqs(8)
#endif

!------------------------------------------------------------

!$omp parallel do default(shared)  &
!$omp private(i,j)
    do j=1,nj
    do i=1,ni
      s(i,j)=sfc(i,j)
    enddo
    enddo

!$omp parallel do default(shared)  &
!$omp private(i,j)
    do j=0,nj+1
    do i=0,ni+1
      dum1(i,j)=0.0
      dum2(i,j)=0.0
    enddo
    enddo

    DO nrk=1,3

!-------------------------
!  set boundary conditions

#ifdef MPI
      call comm_2d_start(s,west,newwest,east,neweast,   &
                           south,newsouth,north,newnorth,reqs)
#else
      if(wbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j)
        do j=1,nj
          s( 0,j)=s(ni  ,j)
          s(-1,j)=s(ni-1,j)
          s(-2,j)=s(ni-2,j)
        enddo
      endif
      if(ebc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j)
        do j=1,nj
          s(ni+1,j)=s(1,j)
          s(ni+2,j)=s(2,j)
          s(ni+3,j)=s(3,j)
        enddo
      endif
      if(sbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i)
        do i=1,ni
          s(i, 0)=s(i,nj  )
          s(i,-1)=s(i,nj-1)
          s(i,-2)=s(i,nj-2)
        enddo
      endif
      if(nbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i)
        do i=1,ni
          s(i,nj+1)=s(i,1)
          s(i,nj+2)=s(i,2)
          s(i,nj+3)=s(i,3)
        enddo
      endif
#endif

      if(ibw.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j)
        do j=1,nj
          s( 0,j)=s(1,j)
          s(-1,j)=s(1,j)
          s(-2,j)=s(1,j)
        enddo
      endif

      if(ibe.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j)
        do j=1,nj
          s(ni+1,j)=s(ni,j)
          s(ni+2,j)=s(ni,j)
          s(ni+3,j)=s(ni,j)
        enddo
      endif

      if(ibs.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i)
        do i=1,ni
          s(i, 0)=s(i,1)
          s(i,-1)=s(i,1)
          s(i,-2)=s(i,1)
        enddo
      endif

      if(ibn.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i)
        do i=1,ni
          s(i,nj+1)=s(i,nj)
          s(i,nj+2)=s(i,nj)
          s(i,nj+3)=s(i,nj)
        enddo
      endif

!-------------------------

#ifdef MPI
      call comm_2dew_end(s,west,newwest,east,neweast,reqs)
#endif

!-------------------------

    if(abs(umove).gt.0.01)then
!$omp parallel do default(shared)  &
!$omp private(i,j,s1,s2,s3,s4,s5,f1,f2,f3,b1,b2,b3,w1,w2,w3)
      do j=1,nj
      do i=1,ni+1
        if(umove.ge.0.0)then
          s1=s(i+2,j)
          s2=s(i+1,j)
          s3=s(i  ,j)
          s4=s(i-1,j)
          s5=s(i-2,j)
        else
          s1=s(i-3,j)
          s2=s(i-2,j)
          s3=s(i-1,j)
          s4=s(i  ,j)
          s5=s(i+1,j)
        endif

        f1=( 2.0*s1 -7.0*s2 +11.0*s3 )*onedsix
        f2=(    -s2 +5.0*s3  +2.0*s4 )*onedsix
        f3=( 2.0*s3 +5.0*s4      -s5 )*onedsix

        b1=thdtw*( s1 -2.0*s2 +s3 )**2 + 0.25*(     s1 -4.0*s2 +3.0*s3 )**2
        b2=thdtw*( s2 -2.0*s3 +s4 )**2 + 0.25*(     s2             -s4 )**2
        b3=thdtw*( s3 -2.0*s4 +s5 )**2 + 0.25*( 3.0*s3 -4.0*s4     +s5 )**2

        w1=0.10/((epsilon+b1)**2)
        w2=0.60/((epsilon+b2)**2)
        w3=0.30/((epsilon+b3)**2)

        dum1(i,j)=((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo
    endif

!-------------------------

#ifdef MPI
      call comm_2dns_end(s,south,newsouth,north,newnorth,reqs)
#endif

!-------------------------

    if(abs(vmove).gt.0.01)then
!$omp parallel do default(shared)  &
!$omp private(i,j,s1,s2,s3,s4,s5,f1,f2,f3,b1,b2,b3,w1,w2,w3)
      do j=1,nj+1
      do i=1,ni
        if(vmove.ge.0.0)then
          s1=s(i,j+2)
          s2=s(i,j+1)
          s3=s(i,j  )
          s4=s(i,j-1)
          s5=s(i,j-2)
        else
          s1=s(i,j-3)
          s2=s(i,j-2)
          s3=s(i,j-1)
          s4=s(i,j  )
          s5=s(i,j+1)
        endif

        f1=( 2.0*s1 -7.0*s2 +11.0*s3 )*onedsix
        f2=(    -s2 +5.0*s3  +2.0*s4 )*onedsix
        f3=( 2.0*s3 +5.0*s4      -s5 )*onedsix

        b1=thdtw*( s1 -2.0*s2 +s3 )**2 + 0.25*(     s1 -4.0*s2 +3.0*s3 )**2
        b2=thdtw*( s2 -2.0*s3 +s4 )**2 + 0.25*(     s2             -s4 )**2
        b3=thdtw*( s3 -2.0*s4 +s5 )**2 + 0.25*( 3.0*s3 -4.0*s4     +s5 )**2

        w1=0.10/((epsilon+b1)**2)
        w2=0.60/((epsilon+b2)**2)
        w3=0.30/((epsilon+b3)**2)

        dum2(i,j)=((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo
    endif

      tem=dt/(4-nrk)

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        s(i,j)=max( rmax , sfc(i,j)+tem*(                          &
                          umove*(dum1(i+1,j)-dum1(i,j))*rdx*uh(i)  &
                         +vmove*(dum2(i,j+1)-dum2(i,j))*rdy*vh(j) ) )
      enddo
      enddo

!-------------------------

    ENDDO

!$omp parallel do default(shared)  &
!$omp private(i,j)
    do j=1,nj
    do i=1,ni
      sfc(i,j)=s(i,j)
    enddo
    enddo

!----------------------------------------------------------------

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bsx(dt,bsq,rvh,rmh,dum)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real :: dt
      real*8 bsq
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: dum

      integer j,k
      real*8 tem
      real*8 bud(nk)

!-----------------------------------------------------------------------

!$omp parallel do default(shared)  &
!$omp private(k)
      do k=1,nk
        bud(k)=0.0d0
      enddo

      if(ibw.eq.1 .and. wbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          bud(k)=bud(k)+dum(1,j,k)*rvh(j)*rmh(1,j,k)
        enddo
        enddo
      endif

      if(ibe.eq.1 .and. ebc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          bud(k)=bud(k)-dum(ni+1,j,k)*rvh(j)*rmh(ni+1,j,k)
        enddo
        enddo
      endif

      tem=dt*dy*dz
      do k=1,nk
        bsq=bsq+bud(k)*tem
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_stat=time_stat+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bsy(dt,bsq,ruh,rmh,dum)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real :: dt
      real*8 bsq
      real, dimension(ib:ie) :: ruh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: dum

      integer i,k
      real*8 tem
      real*8 bud(nk)

!-----------------------------------------------------------------------

!$omp parallel do default(shared)  &
!$omp private(k)
      do k=1,nk
        bud(k)=0.0d0
      enddo

      if(ibs.eq.1 .and. sbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          bud(k)=bud(k)+dum(i,1,k)*ruh(i)*rmh(i,1,k)
        enddo
        enddo
      endif

      if(ibn.eq.1 .and. nbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          bud(k)=bud(k)-dum(i,nj+1,k)*ruh(i)*rmh(i,nj+1,k)
        enddo
        enddo
      endif

      tem=dt*dx*dz
      do k=1,nk
        bsq=bsq+bud(k)*tem
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_stat=time_stat+mytime()

      return
      end


