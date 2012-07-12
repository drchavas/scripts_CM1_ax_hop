

!-----------------------------------------------------------------------
!  message passing routines
!-----------------------------------------------------------------------


      function nabor(i,j,nx,ny)
      implicit none
      integer i,j,nx,ny,nabor
      integer newi,newj

      newi=i
      newj=j

      if ( newi .lt.  1 ) newi = nx
      if ( newi .gt.  nx) newi = 1

      if ( newj .lt.  1 ) newj = ny
      if ( newj .gt.  ny) newj = 1

      nabor = (newi-1) + (newj-1)*nx

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine prepcorners(s)
      implicit none

      include 'input.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: s

      integer :: i,j

      do j=0,nj+1
      do i=0,ni+1
        s(i,j,0)    = s(i,j,1)
        s(i,j,nk+1) = s(i,j,nk)
      enddo
      enddo

      return
      end


