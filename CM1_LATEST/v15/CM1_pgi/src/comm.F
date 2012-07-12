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

#ifdef MPI
      subroutine getcorner(s)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real s(ib:ie,jb:je,kb:ke)
 
      integer k
      integer reqs(8),status(MPI_STATUS_SIZE)
      real, dimension(nk) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer tag,count

!-----

      tag=5001
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(2),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          se1(k)=s(ni,1,k)
        enddo
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(1),ierr)
      endif
 
!-----

      tag=5002
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(4),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          ne1(k)=s(ni,nj,k)
        enddo
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(3),ierr)
      endif
 
!-----

      tag=5003
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          sw1(k)=s(1,1,k)
        enddo
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif
 
!-----

      tag=5004
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          nw1(k)=s(1,nj,k)
        enddo
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif
 
!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          s(0,nj+1,k)=nw2(k)
        enddo
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          s(0,0,k)=sw2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          s(ni+1,nj+1,k)=ne2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          s(ni+1,0,k)=se2(k)
        enddo
      endif

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcornert(t)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real t(ib:ie,jb:je,kb:ke+1)

      integer k
      integer reqs(8),status(MPI_STATUS_SIZE)
      real, dimension(nk+1) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer tag,count

!-----

      tag=5001
      count=nkp1

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(2),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          se1(k)=t(ni,1,k)
        enddo
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(1),ierr)
      endif

!-----

      tag=5002
      count=nkp1

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(4),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          ne1(k)=t(ni,nj,k)
        enddo
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(3),ierr)
      endif

!-----

      tag=5003
      count=nkp1

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          sw1(k)=t(1,1,k)
        enddo
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif

!-----

      tag=5004
      count=nkp1

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          nw1(k)=t(1,nj,k)
        enddo
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif

!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          t(0,nj+1,k)=nw2(k)
        enddo
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          t(0,0,k)=sw2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          t(ni+1,nj+1,k)=ne2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          t(ni+1,0,k)=se2(k)
        enddo
      endif

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcorneru(u)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real u(ib:ie+1,jb:je,kb:ke)

      integer k
      integer reqs(8),status(MPI_STATUS_SIZE)
      real, dimension(nk) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer tag,count

!-----

      tag=5001
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(2),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          se1(k)=u(ni,1,k)
        enddo
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(1),ierr)
      endif

!-----

      tag=5002
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(4),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          ne1(k)=u(ni,nj,k)
        enddo
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(3),ierr)
      endif

!-----

      tag=5003
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          sw1(k)=u(2,1,k)
        enddo
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif

!-----

      tag=5004
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          nw1(k)=u(2,nj,k)
        enddo
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif

!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          u(0,nj+1,k)=nw2(k)
        enddo
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          u(0,0,k)=sw2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          u(ni+2,nj+1,k)=ne2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          u(ni+2,0,k)=se2(k)
        enddo
      endif

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcornerv(v)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real v(ib:ie,jb:je+1,kb:ke)

      integer k
      integer reqs(8),status(MPI_STATUS_SIZE)
      real, dimension(nk) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer tag,count

!-----

      tag=5001
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(2),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          se1(k)=v(ni,2,k)
        enddo
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(1),ierr)
      endif

!-----

      tag=5002
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(4),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          ne1(k)=v(ni,nj,k)
        enddo
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(3),ierr)
      endif

!-----

      tag=5003
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          sw1(k)=v(1,2,k)
        enddo
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif

!-----

      tag=5004
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          nw1(k)=v(1,nj,k)
        enddo
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif

!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          v(0,nj+2,k)=nw2(k)
        enddo
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          v(0,0,k)=sw2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          v(ni+1,nj+2,k)=ne2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          v(ni+1,0,k)=se2(k)
        enddo
      endif

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcornerw(w)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real w(ib:ie,jb:je,kb:ke+1)

      integer k
      integer reqs(8),status(MPI_STATUS_SIZE)
      real, dimension(nk+1) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer tag,count

!-----

      tag=5001
      count=nkp1

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(2),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          se1(k)=w(ni,1,k)
        enddo
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(1),ierr)
      endif

!-----

      tag=5002
      count=nkp1

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(4),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          ne1(k)=w(ni,nj,k)
        enddo
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(3),ierr)
      endif

!-----

      tag=5003
      count=nkp1

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          sw1(k)=w(1,1,k)
        enddo
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif

!-----

      tag=5004
      count=nkp1

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          nw1(k)=w(1,nj,k)
        enddo
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif

!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          w(0,nj+1,k)=nw2(k)
        enddo
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          w(0,0,k)=sw2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          w(ni+1,nj+1,k)=ne2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          w(ni+1,0,k)=se2(k)
        enddo
      endif

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2we_start(s,west,newwest,east,neweast,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real s(ib:ie,jb:je,kb:ke)
      real west(2,nj,nk),newwest(2,nj,nk)
      real east(2,nj,nk),neweast(2,nj,nk)
      integer reqs(4)

      integer i,j,k
      integer tag,count

!------------------------------------------------

      count=cs2we
      nf=nf+1
      tag=nf

      ! receive east
      if(ibe.eq.0)then
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(2),ierr)
      endif

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,2
          west(i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(1),ierr)
      endif

!----------

      count=cs2we
      nf=nf+1
      tag=nf

      ! receive west
      if(ibw.eq.0)then
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(4),ierr)
      endif

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,2
          east(i,j,k)=s(ni-2+i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(3),ierr)
      endif

!----------

      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2we_end(s,west,newwest,east,neweast,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real s(ib:ie,jb:je,kb:ke)
      real west(2,nj,nk),newwest(2,nj,nk)
      real east(2,nj,nk),neweast(2,nj,nk)
      integer reqs(4)

      integer i,j,k
      integer status(MPI_STATUS_SIZE)

!-------------------------------------------------------------------

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,2
          s(ni+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      endif

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,2
          s(i-2,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      endif

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2sn_start(s,south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real s(ib:ie,jb:je,kb:ke)
      real south(ni,2,nk),newsouth(ni,2,nk)
      real north(ni,2,nk),newnorth(ni,2,nk)
      integer reqs(4)

      integer i,j,k
      integer tag,count

!----------

      count=cs2sn
      nf=nf+1
      tag=nf

      ! receive south
      if(ibs.eq.0)then
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(2),ierr)
      endif

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,2
        do i=1,ni
          north(i,j,k)=s(i,nj-2+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(1),ierr)
      endif

!----------

      count=cs2sn
      nf=nf+1
      tag=nf

      ! receive north
      if(ibn.eq.0)then
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(4),ierr)
      endif

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,2
        do i=1,ni
          south(i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(3),ierr)
      endif

!----------

      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2sn_end(s,south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real s(ib:ie,jb:je,kb:ke)
      real south(ni,2,nk),newsouth(ni,2,nk)
      real north(ni,2,nk),newnorth(ni,2,nk)
      integer reqs(4)

      integer i,j,k
      integer status(MPI_STATUS_SIZE)

!----------

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,2
        do i=1,ni
          s(i,j-2,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,2
        do i=1,ni
          s(i,nj+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_3s_start(s,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real s(ib:ie,jb:je,kb:ke)
      real west(3,nj,nk),newwest(3,nj,nk)
      real east(3,nj,nk),neweast(3,nj,nk)
      real south(ni,3,nk),newsouth(ni,3,nk)
      real north(ni,3,nk),newnorth(ni,3,nk)
      integer reqs(8)
 
      integer i,j,k
      integer tag,count
 
!------------------------------------------------

      count=cs3we
      nf=nf+1
      tag=nf

      ! receive east
      if(ibe.eq.0)then
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(2),ierr)
      endif

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          west(i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(1),ierr)
      endif
 
!----------
 
      count=cs3we
      nf=nf+1
      tag=nf

      ! receive west
      if(ibw.eq.0)then
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(4),ierr)
      endif
 
      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          east(i,j,k)=s(ni-3+i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(3),ierr)
      endif
 
!----------
 
      count=cs3sn
      nf=nf+1
      tag=nf

      ! receive south
      if(ibs.eq.0)then
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(6),ierr)
      endif
 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          north(i,j,k)=s(i,nj-3+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif
 
!----------
 
      count=cs3sn
      nf=nf+1
      tag=nf

      ! receive north
      if(ibn.eq.0)then
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(8),ierr)
      endif
 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          south(i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif
 
!----------
 
      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3s_end(s,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real s(ib:ie,jb:je,kb:ke)
      real west(3,nj,nk),newwest(3,nj,nk)
      real east(3,nj,nk),neweast(3,nj,nk)
      real south(ni,3,nk),newsouth(ni,3,nk)
      real north(ni,3,nk),newnorth(ni,3,nk)
      integer reqs(8)
 
      integer i,j,k
      integer status(MPI_STATUS_SIZE)

!-------------------------------------------------------------------
 
      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          s(ni+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      endif
 
      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          s(i-3,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      endif
 
!----------

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          s(i,j-3,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      endif
 
      if(ibn.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          s(i,nj+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      endif

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()
 
!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then

          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(1,0,k)
            enddo
          endif

          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(1,nj+1,k)
            enddo
          endif

        endif

        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni+1,nj,k)
            enddo
          endif
 
        endif
 
      endif

!-----------------------------------------------------------
!  Mirror b.c. patch
!
!      if(ibn.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=1,nk
!        do i=0,ni+1
!          s(i,nj+1,k)=s(i,nj  ,k)
!          s(i,nj+2,k)=s(i,nj-1,k)
!          s(i,nj+3,k)=s(i,nj-2,k)
!        enddo
!        enddo
!      endif
!
!      if(ibw.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=1,nk
!        do j=0,nj+1
!          s(-2,j,k)=s(3,j,k)
!          s(-1,j,k)=s(2,j,k)
!          s( 0,j,k)=s(1,j,k)
!        enddo
!        enddo
!      endif
!
!-----------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_3t_start(t,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real t(ib:ie,jb:je,kb:ke+1)
      real west(3,nj,nk+1),newwest(3,nj,nk+1)
      real east(3,nj,nk+1),neweast(3,nj,nk+1)
      real south(ni,3,nk+1),newsouth(ni,3,nk+1)
      real north(ni,3,nk+1),newnorth(ni,3,nk+1)
      integer reqs(8)

      integer i,j,k
      integer tag,count

!------------------------------------------------

      count=ct3we
      nf=nf+1
      tag=nf

      ! receive east
      if(ibe.eq.0)then
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(2),ierr)
      endif

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,nj
        do i=1,3
          west(i,j,k)=t(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(1),ierr)
      endif

!----------

      count=ct3we
      nf=nf+1
      tag=nf

      ! receive west
      if(ibw.eq.0)then
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(4),ierr)
      endif

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,nj
        do i=1,3
          east(i,j,k)=t(ni-3+i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(3),ierr)
      endif

!----------

      count=ct3sn
      nf=nf+1
      tag=nf

      ! receive south
      if(ibs.eq.0)then
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(6),ierr)
      endif

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,3
        do i=1,ni
          north(i,j,k)=t(i,nj-3+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif

!----------

      count=ct3sn
      nf=nf+1
      tag=nf

      ! receive north
      if(ibn.eq.0)then
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(8),ierr)
      endif

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,3
        do i=1,ni
          south(i,j,k)=t(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif

!----------

      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_3t_end(t,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real t(ib:ie,jb:je,kb:ke+1)
      real west(3,nj,nk+1),newwest(3,nj,nk+1)
      real east(3,nj,nk+1),neweast(3,nj,nk+1)
      real south(ni,3,nk+1),newsouth(ni,3,nk+1)
      real north(ni,3,nk+1),newnorth(ni,3,nk+1)
      integer reqs(8)

      integer i,j,k
      integer status(MPI_STATUS_SIZE)

!-------------------------------------------------------------------

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,nj
        do i=1,3
          t(ni+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      endif

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,nj
        do i=1,3
          t(i-3,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      endif

!----------

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,3
        do i=1,ni
          t(i,j-3,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,3
        do i=1,ni
          t(i,nj+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      endif

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------
!  patch for corner

     if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then

        if(ibw.eq.1)then

          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,0,k)=t(1,0,k)
            enddo
          endif

          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,nj+1,k)=t(1,nj+1,k)
            enddo
          endif

        endif

        if(ibe.eq.1)then

          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,0,k)=t(ni,0,k)
            enddo
          endif

          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,nj+1,k)=t(ni,nj+1,k)
            enddo
          endif

        endif

      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then

        if(ibs.eq.1)then

          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,0,k)=t(0,1,k)
            enddo
          endif

          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,0,k)=t(ni+1,1,k)
            enddo
          endif

        endif

        if(ibn.eq.1)then

          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,nj+1,k)=t(0,nj,k)
            enddo
          endif

          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,nj+1,k)=t(ni+1,nj,k)
            enddo
          endif

        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3u_start(u,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real u(ib:ie+1,jb:je,kb:ke)
      real west(3,nj,nk),newwest(3,nj,nk)
      real east(3,nj,nk),neweast(3,nj,nk)
      real south(ni+1,3,nk),newsouth(ni+1,3,nk)
      real north(ni+1,3,nk),newnorth(ni+1,3,nk)
      integer reqs(8)
 
      integer i,j,k
      integer tag,count
 
!------------------------------------------------

      nu=nu+1
      tag=1000+nu
      count=cs3we

      ! receive east
      if(ibe.eq.0)then
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(2),ierr)
      endif
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          west(i,j,k)=u(i+1,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(1),ierr)
      endif
 
!----------
 
      nu=nu+1
      tag=1000+nu
      count=cs3we
 
      ! receive west
      if(ibw.eq.0)then
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(4),ierr)
      endif

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          east(i,j,k)=u(ni-3+i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(3),ierr)
      endif

!----------
 
      nu=nu+1
      tag=1000+nu
      count=cu3sn

      ! receive north
      if(ibn.eq.0)then
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(6),ierr)
      endif
 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni+1
          south(i,j,k)=u(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif
 
!----------
 
      nu=nu+1
      tag=1000+nu
      count=cu3sn

      ! receive south
      if(ibs.eq.0)then
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(8),ierr)
      endif
 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni+1
          north(i,j,k)=u(i,nj-3+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif
 
!----------

      if(timestats.ge.1) time_mpu1=time_mpu1+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3u_end(u,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real u(ib:ie+1,jb:je,kb:ke)
      real west(3,nj,nk),newwest(3,nj,nk)
      real east(3,nj,nk),neweast(3,nj,nk)
      real south(ni+1,3,nk),newsouth(ni+1,3,nk)
      real north(ni+1,3,nk),newnorth(ni+1,3,nk)
      integer reqs(8)
 
      integer i,j,k
      integer status(MPI_STATUS_SIZE)

!----------
 
      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          u(ni+1+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      endif
 
      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          u(i-3,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      endif
 
      if(ibn.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni+1
          u(i,nj+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      endif
 
      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni+1
          u(i,j-3,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      endif

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif
 
      if(timestats.ge.1) time_mpu2=time_mpu2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,0,k)=u(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,nj+1,k)=u(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,0,k)=u(ni+1,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,nj+1,k)=u(ni+1,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,0,k)=u(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,0,k)=u(ni+2,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,nj+1,k)=u(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,nj+1,k)=u(ni+2,nj,k)
            enddo
          endif
 
        endif
 
      endif

!-----------------------------------------------------------
!  Mirror b.c. patch
!
!      if(ibn.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=1,nk
!        do i=0,ni+2
!          u(i,nj+1,k)=u(i,nj  ,k)
!          u(i,nj+2,k)=u(i,nj-1,k)
!          u(i,nj+3,k)=u(i,nj-2,k)
!        enddo
!        enddo
!      endif
!
!      if(ibw.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=1,nk
!        do j=0,nj+1
!          u(-2,j,k)=-u(4,j,k)
!          u(-1,j,k)=-u(3,j,k)
!          u( 0,j,k)=-u(2,j,k)
!          u( 1,j,k)=0.0
!        enddo
!        enddo
!      endif
!
!-----------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
!----------
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3v_start(v,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real v(ib:ie,jb:je+1,kb:ke)
      real west(3,nj+1,nk),newwest(3,nj+1,nk)
      real east(3,nj+1,nk),neweast(3,nj+1,nk)
      real south(ni,3,nk),newsouth(ni,3,nk)
      real north(ni,3,nk),newnorth(ni,3,nk)
      integer reqs(8)
 
      integer i,j,k
      integer tag,count
 
!------------------------------------------------

      nv=nv+1
      tag=2000+nv
      count=cv3we

      ! receive east
      if(ibe.eq.0)then
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(2),ierr)
      endif
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,3
          west(i,j,k)=v(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(1),ierr)
      endif
 
!----------
 
      nv=nv+1
      tag=2000+nv
      count=cv3we
 
      ! receive west
      if(ibw.eq.0)then
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(4),ierr)
      endif

      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,3
          east(i,j,k)=v(ni-3+i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(3),ierr)
      endif

!----------
 
      nv=nv+1
      tag=2000+nv
      count=cs3sn
 
      ! receive north
      if(ibn.eq.0)then
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(6),ierr)
      endif

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          south(i,j,k)=v(i,j+1,k)
        enddo
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif
 
!----------
 
      nv=nv+1
      tag=2000+nv
      count=cs3sn
 
      ! receive south
      if(ibs.eq.0)then
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(8),ierr)
      endif

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          north(i,j,k)=v(i,nj-3+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif
 
!----------

      if(timestats.ge.1) time_mpv1=time_mpv1+mytime()
 
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3v_end(v,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real v(ib:ie,jb:je+1,kb:ke)
      real west(3,nj+1,nk),newwest(3,nj+1,nk)
      real east(3,nj+1,nk),neweast(3,nj+1,nk)
      real south(ni,3,nk),newsouth(ni,3,nk)
      real north(ni,3,nk),newnorth(ni,3,nk)
      integer reqs(8)
 
      integer i,j,k
      integer status(MPI_STATUS_SIZE)
 
!--------
 
      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,3
          v(ni+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      endif

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,3
          v(i-3,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      endif
 
      if(ibn.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          v(i,nj+1+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      endif
 
      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          v(i,j-3,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      endif

!--------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif
 
      if(timestats.ge.1) time_mpv2=time_mpv2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,0,k)=v(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,nj+2,k)=v(1,nj+2,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,0,k)=v(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,nj+2,k)=v(ni,nj+2,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,0,k)=v(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,0,k)=v(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,nj+2,k)=v(0,nj+1,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,nj+2,k)=v(ni+1,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

!-----------------------------------------------------------
!  Mirror b.c. patch
!
!      if(ibn.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=1,nk
!        do i=0,ni+1
!          v(i,nj+1,k)=0.0
!          v(i,nj+2,k)=-v(i,nj  ,k)
!          v(i,nj+3,k)=-v(i,nj-1,k)
!          v(i,nj+4,k)=-v(i,nj-2,k)
!        enddo
!        enddo
!      endif
!
!      if(ibw.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=1,nk
!        do j=0,nj+2
!          v(-2,j,k)=v(3,j,k)
!          v(-1,j,k)=v(2,j,k)
!          v( 0,j,k)=v(1,j,k)
!        enddo
!        enddo
!      endif
!
!-----------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
!----------
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3w_start(w,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real w(ib:ie,jb:je,kb:ke+1)
      real west(3,nj,nk-1),newwest(3,nj,nk-1)
      real east(3,nj,nk-1),neweast(3,nj,nk-1)
      real south(ni,3,nk-1),newsouth(ni,3,nk-1)
      real north(ni,3,nk-1),newnorth(ni,3,nk-1)
      integer reqs(8)
 
      integer i,j,k
      integer tag,count
 
!------------------------------------------------

      nw=nw+1
      tag=3000+nw
      count=cw3we

      ! receive east
      if(ibe.eq.0)then
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(2),ierr)
      endif
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,3
          west(i,j,k-1)=w(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,    &
                       MPI_COMM_WORLD,reqs(1),ierr)
      endif
 
!----------
 
      nw=nw+1
      tag=3000+nw
      count=cw3we
 
      ! receive west
      if(ibw.eq.0)then
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(4),ierr)
      endif

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,3
          east(i,j,k-1)=w(ni-3+i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(3),ierr)
      endif

!----------

      nw=nw+1
      tag=3000+nw
      count=cw3sn
 
      ! receive north
      if(ibn.eq.0)then
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(6),ierr)
      endif

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,3
        do i=1,ni
          south(i,j,k-1)=w(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif
 
!----------
 
      nw=nw+1
      tag=3000+nw
      count=cw3sn
 
      ! receive south
      if(ibs.eq.0)then
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(8),ierr)
      endif

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,3
        do i=1,ni
          north(i,j,k-1)=w(i,nj-3+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif
 
!----------

      if(timestats.ge.1) time_mpw1=time_mpw1+mytime()
 
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3w_end(w,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real w(ib:ie,jb:je,kb:ke+1)
      real west(3,nj,nk-1),newwest(3,nj,nk-1)
      real east(3,nj,nk-1),neweast(3,nj,nk-1)
      real south(ni,3,nk-1),newsouth(ni,3,nk-1)
      real north(ni,3,nk-1),newnorth(ni,3,nk-1)
      integer reqs(8)
 
      integer i,j,k
      integer status(MPI_STATUS_SIZE)

!--------
 
      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,3
          w(ni+i,j,k)=neweast(i,j,k-1)
        enddo
        enddo
        enddo
      endif
 
      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,3
          w(i-3,j,k)=newwest(i,j,k-1)
        enddo
        enddo
        enddo
      endif
 
      if(ibn.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,3
        do i=1,ni
          w(i,nj+j,k)=newnorth(i,j,k-1)
        enddo
        enddo
        enddo
      endif
 
      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,3
        do i=1,ni
          w(i,j-3,k)=newsouth(i,j,k-1)
        enddo
        enddo
        enddo
      endif
 
!--------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(timestats.ge.1) time_mpw2=time_mpw2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,0,k)=w(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,nj+1,k)=w(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,0,k)=w(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,nj+1,k)=w(ni,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,0,k)=w(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,0,k)=w(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,nj+1,k)=w(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,nj+1,k)=w(ni+1,nj,k)
            enddo
          endif
 
        endif
 
      endif

!-----------------------------------------------------------
!  Mirror b.c. patch
!
!      if(ibn.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=2,nk
!        do i=0,ni+1
!          w(i,nj+1,k)=w(i,nj  ,k)
!          w(i,nj+2,k)=w(i,nj-1,k)
!          w(i,nj+3,k)=w(i,nj-2,k)
!        enddo
!        enddo
!      endif
!
!      if(ibw.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=2,nk
!        do j=0,nj+1
!          w(-2,j,k)=w(3,j,k)
!          w(-1,j,k)=w(2,j,k)
!          w( 0,j,k)=w(1,j,k)
!        enddo
!        enddo
!      endif
!
!-----------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
!----------
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1s_start(s,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real s(ib:ie,jb:je,kb:ke)
      real west(nj,nk),newwest(nj,nk)
      real east(nj,nk),neweast(nj,nk)
      real south(ni,nk),newsouth(ni,nk)
      real north(ni,nk),newnorth(ni,nk)
      integer reqs(8)
 
      integer i,j,k
      integer tag,count
 
!------------------------------------------------

      nf=nf+1
      tag=nf
      count=cs1we

      ! receive east
      if(ibe.eq.0)then
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(2),ierr)
      endif
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          west(j,k)=s(1,j,k)
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(1),ierr)
      endif
 
!----------

      nf=nf+1
      tag=nf
      count=cs1we

      ! receive west
      if(ibw.eq.0)then
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(4),ierr)
      endif
 
      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          east(j,k)=s(ni,j,k)
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(3),ierr)
      endif
 
!----------

      nf=nf+1
      tag=nf
      count=cs1sn

      ! receive south
      if(ibs.eq.0)then
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(6),ierr)
      endif
 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          north(i,k)=s(i,nj,k)
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif
 
!----------

      nf=nf+1
      tag=nf
      count=cs1sn

      ! receive north
      if(ibn.eq.0)then
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(8),ierr)
      endif
 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          south(i,k)=s(i,1,k)
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif

!----------
 
      if(timestats.ge.1) time_mps1=time_mps1+mytime()
 
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1s_end(s,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real s(ib:ie,jb:je,kb:ke)
      real west(nj,nk),newwest(nj,nk)
      real east(nj,nk),neweast(nj,nk)
      real south(ni,nk),newsouth(ni,nk)
      real north(ni,nk),newnorth(ni,nk)
      integer reqs(8)
 
      integer i,j,k
      integer status(MPI_STATUS_SIZE)

!---------------------------------------------------------------------
 
      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          s(ni+1,j,k)=neweast(j,k)
        enddo
        enddo
      endif
 
      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          s(0,j,k)=newwest(j,k)
        enddo
        enddo
      endif
 
      if(ibs.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          s(i,0,k)=newsouth(i,k)
        enddo
        enddo
      endif
 
      if(ibn.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          s(i,nj+1,k)=newnorth(i,k)
        enddo
        enddo
      endif

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni+1,nj,k)
            enddo
          endif
 
        endif
 
      endif

!-----------------------------------------------------------
!  Mirror b.c. patch
!
!      if(ibn.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=1,nk
!        do i=0,ni+1
!          s(i,nj+1,k)=s(i,nj  ,k)
!          s(i,nj+2,k)=s(i,nj-1,k)
!          s(i,nj+3,k)=s(i,nj-2,k)
!        enddo
!        enddo
!      endif
!
!      if(ibw.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=1,nk
!        do j=0,nj+1
!          s(-2,j,k)=s(3,j,k)
!          s(-1,j,k)=s(2,j,k)
!          s( 0,j,k)=s(1,j,k)
!        enddo
!        enddo
!      endif
!
!-----------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
!----------
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1s_end_we(s,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real s(ib:ie,jb:je,kb:ke)
      real west(nj,nk),newwest(nj,nk)
      real east(nj,nk),neweast(nj,nk)
      real south(ni,nk),newsouth(ni,nk)
      real north(ni,nk),newnorth(ni,nk)
      integer reqs(8)

      integer i,j,k
      integer status(MPI_STATUS_SIZE)

!---------------------------------------------------------------------

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          s(ni+1,j,k)=neweast(j,k)
        enddo
        enddo
      endif

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          s(0,j,k)=newwest(j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1s_end_sn(s,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real s(ib:ie,jb:je,kb:ke)
      real west(nj,nk),newwest(nj,nk)
      real east(nj,nk),neweast(nj,nk)
      real south(ni,nk),newsouth(ni,nk)
      real north(ni,nk),newnorth(ni,nk)
      integer reqs(8)

      integer i,j,k
      integer status(MPI_STATUS_SIZE)

!-----------------------------------------------------------

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          s(i,0,k)=newsouth(i,k)
        enddo
        enddo
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          s(i,nj+1,k)=newnorth(i,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 

      subroutine comm_1s_end_misc(s,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real s(ib:ie,jb:je,kb:ke)
      real west(nj,nk),newwest(nj,nk)
      real east(nj,nk),neweast(nj,nk)
      real south(ni,nk),newsouth(ni,nk)
      real north(ni,nk),newnorth(ni,nk)
      integer reqs(8)

      integer i,j,k
      integer status(MPI_STATUS_SIZE)

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------
!  patch for corner

      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then

        if(ibw.eq.1)then

          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(1,0,k)
            enddo
          endif

          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(1,nj+1,k)
            enddo
          endif

        endif

        if(ibe.eq.1)then

          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni,0,k)
            enddo
          endif

          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni,nj+1,k)
            enddo
          endif

        endif

      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then

        if(ibs.eq.1)then

          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(0,1,k)
            enddo
          endif

          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni+1,1,k)
            enddo
          endif

        endif

        if(ibn.eq.1)then

          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(0,nj,k)
            enddo
          endif

          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni+1,nj,k)
            enddo
          endif

        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1t_start(t,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real t(ib:ie,jb:je,kb:ke+1)
      real west(nj,nk+1),newwest(nj,nk+1)
      real east(nj,nk+1),neweast(nj,nk+1)
      real south(ni,nk+1),newsouth(ni,nk+1)
      real north(ni,nk+1),newnorth(ni,nk+1)
      integer reqs(8)

      integer i,j,k
      integer tag,count

!------------------------------------------------

      nf=nf+1
      tag=nf
      count=ct1we

      ! receive east
      if(ibe.eq.0)then
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(2),ierr)
      endif

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nkp1
        do j=1,nj
          west(j,k)=t(1,j,k)
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(1),ierr)
      endif

!----------

      nf=nf+1
      tag=nf
      count=ct1we

      ! receive west
      if(ibw.eq.0)then
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(4),ierr)
      endif

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nkp1
        do j=1,nj
          east(j,k)=t(ni,j,k)
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(3),ierr)
      endif

!----------

      nf=nf+1
      tag=nf
      count=ct1sn

      ! receive south
      if(ibs.eq.0)then
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(6),ierr)
      endif

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nkp1
        do i=1,ni
          north(i,k)=t(i,nj,k)
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif

!----------

      nf=nf+1
      tag=nf
      count=ct1sn

      ! receive north
      if(ibn.eq.0)then
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(8),ierr)
      endif

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nkp1
        do i=1,ni
          south(i,k)=t(i,1,k)
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif

!----------

      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1t_end(t,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real t(ib:ie,jb:je,kb:ke+1)
      real west(nj,nk+1),newwest(nj,nk+1)
      real east(nj,nk+1),neweast(nj,nk+1)
      real south(ni,nk+1),newsouth(ni,nk+1)
      real north(ni,nk+1),newnorth(ni,nk+1)
      integer reqs(8)

      integer i,j,k
      integer status(MPI_STATUS_SIZE)

!---------------------------------------------------------------------

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nkp1
        do j=1,nj
          t(ni+1,j,k)=neweast(j,k)
        enddo
        enddo
      endif

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nkp1
        do j=1,nj
          t(0,j,k)=newwest(j,k)
        enddo
        enddo
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nkp1
        do i=1,ni
          t(i,0,k)=newsouth(i,k)
        enddo
        enddo
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nkp1
        do i=1,ni
          t(i,nj+1,k)=newnorth(i,k)
        enddo
        enddo
      endif

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------
!  patch for corner

      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then

        if(ibw.eq.1)then

          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,0,k)=t(1,0,k)
            enddo
          endif

          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,nj+1,k)=t(1,nj+1,k)
            enddo
          endif

        endif

        if(ibe.eq.1)then

          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,0,k)=t(ni,0,k)
            enddo
          endif

          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,nj+1,k)=t(ni,nj+1,k)
            enddo
          endif

        endif

      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then

        if(ibs.eq.1)then

          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,0,k)=t(0,1,k)
            enddo
          endif

          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,0,k)=t(ni+1,1,k)
            enddo
          endif

        endif

        if(ibn.eq.1)then

          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,nj+1,k)=t(0,nj,k)
            enddo
          endif

          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,nj+1,k)=t(ni+1,nj,k)
            enddo
          endif

        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1u_start(u,uw1,uw2,ue1,ue2,   &
                                 us1,us2,un1,un2,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real u(ib:ie+1,jb:je,kb:ke)
      real uw1(nj,nk),uw2(nj,nk)
      real ue1(nj,nk),ue2(nj,nk)
      real us1(ni+1,nk),us2(ni+1,nk)
      real un1(ni+1,nk),un2(ni+1,nk)
      integer reqs(8)
 
      integer i,j,k
      integer tag,count
 
!----

      nu=nu+1
      tag=1000+nu
      count=cs1we
 
        ! receive east
        if(ibe.eq.0)then
          call mpi_irecv(ue2,count,MPI_REAL,myeast,tag,   &
                        MPI_COMM_WORLD,reqs(2),ierr)
        endif

        ! send west
        if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            uw1(j,k)=u(2,j,k)
          enddo
          enddo
          call mpi_isend(uw1,count,MPI_REAL,mywest,tag,   &
                         MPI_COMM_WORLD,reqs(1),ierr)
        endif

!----

      nu=nu+1
      tag=1000+nu
      count=cs1we
 
        ! receive west
        if(ibw.eq.0)then
          call mpi_irecv(uw2,count,MPI_REAL,mywest,tag,   &
                        MPI_COMM_WORLD,reqs(4),ierr)
        endif

        ! send east
        if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            ue1(j,k)=u(ni,j,k)
          enddo
          enddo
          call mpi_isend(ue1,count,MPI_REAL,myeast,tag,   &
                         MPI_COMM_WORLD,reqs(3),ierr)
        endif

!----

      nu=nu+1
      tag=1000+nu
      count=cu1sn
 
        ! receive north
        if(ibn.eq.0)then
          call mpi_irecv(un2,count,MPI_REAL,mynorth,tag,   &
                        MPI_COMM_WORLD,reqs(6),ierr)
        endif

        ! send south
        if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni+1
            us1(i,k)=u(i,1,k)
          enddo
          enddo
          call mpi_isend(us1,count,MPI_REAL,mysouth,tag,   &
                         MPI_COMM_WORLD,reqs(5),ierr)
        endif

!----

      nu=nu+1
      tag=1000+nu
      count=cu1sn
 
        ! receive south
        if(ibs.eq.0)then
          call mpi_irecv(us2,count,MPI_REAL,mysouth,tag,   &
                        MPI_COMM_WORLD,reqs(8),ierr)
        endif

        ! send north
        if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni+1
            un1(i,k)=u(i,nj,k)
          enddo
          enddo
          call mpi_isend(un1,count,MPI_REAL,mynorth,tag,   &
                         MPI_COMM_WORLD,reqs(7),ierr)
        endif

!----

      if(timestats.ge.1) time_mpu1=time_mpu1+mytime()
 
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1u_end(u,uw1,uw2,ue1,ue2,   &
                               us1,us2,un1,un2,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real u(ib:ie+1,jb:je,kb:ke)
      real uw1(nj,nk),uw2(nj,nk)
      real ue1(nj,nk),ue2(nj,nk)
      real us1(ni+1,nk),us2(ni+1,nk)
      real un1(ni+1,nk),un2(ni+1,nk)
      integer reqs(8)
 
      integer i,j,k
      integer status(MPI_STATUS_SIZE)

!-----
 
      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          u(ni+2,j,k)=ue2(j,k)
        enddo
        enddo
      endif
 
      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          u(0,j,k)=uw2(j,k)
        enddo
        enddo
      endif

!-----
 
      if(ibn.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni+1
          u(i,nj+1,k)=un2(i,k)
        enddo
        enddo
      endif
 
      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni+1
          u(i,0,k)=us2(i,k)
        enddo
        enddo
      endif

!-----

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(timestats.ge.1) time_mpu2=time_mpu2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,0,k)=u(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,nj+1,k)=u(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,0,k)=u(ni+1,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,nj+1,k)=u(ni+1,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,0,k)=u(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,0,k)=u(ni+2,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,nj+1,k)=u(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,nj+1,k)=u(ni+2,nj,k)
            enddo
          endif
 
        endif
 
      endif

!-----------------------------------------------------------
!  Mirror b.c. patch
!
!      if(ibn.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=1,nk
!        do i=0,ni+2
!          u(i,nj+1,k)=u(i,nj  ,k)
!          u(i,nj+2,k)=u(i,nj-1,k)
!          u(i,nj+3,k)=u(i,nj-2,k)
!        enddo
!        enddo
!      endif
!
!      if(ibw.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=1,nk
!        do j=0,nj+1
!          u(-2,j,k)=-u(4,j,k)
!          u(-1,j,k)=-u(3,j,k)
!          u( 0,j,k)=-u(2,j,k)
!          u( 1,j,k)=0.0
!        enddo
!        enddo
!      endif
!
!-----------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
!----------
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1v_start(v,vw1,vw2,ve1,ve2,   &
                                 vs1,vs2,vn1,vn2,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real v(ib:ie,jb:je+1,kb:ke)
      real vw1(nj+1,nk),vw2(nj+1,nk)
      real ve1(nj+1,nk),ve2(nj+1,nk)
      real vs1(ni,nk),vs2(ni,nk)
      real vn1(ni,nk),vn2(ni,nk)
      integer reqs(8)
 
      integer i,j,k
      integer tag,count

!-----
 
      nv=nv+1
      tag=2000+nv
      count=cv1we
 
        ! receive east
        if(ibe.eq.0)then
          call mpi_irecv(ve2,count,MPI_REAL,myeast,tag,   &
                        MPI_COMM_WORLD,reqs(2),ierr)
        endif

        ! send west
        if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj+1
            vw1(j,k)=v(1,j,k)
          enddo
          enddo
          call mpi_isend(vw1,count,MPI_REAL,mywest,tag,   &
                         MPI_COMM_WORLD,reqs(1),ierr)
        endif

!-----

 
      nv=nv+1
      tag=2000+nv
      count=cv1we
 
        ! receive west
        if(ibw.eq.0)then
          call mpi_irecv(vw2,count,MPI_REAL,mywest,tag,   &
                        MPI_COMM_WORLD,reqs(4),ierr)
        endif

        ! send east
        if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj+1
            ve1(j,k)=v(ni,j,k)
          enddo
          enddo
          call mpi_isend(ve1,count,MPI_REAL,myeast,tag,   &
                         MPI_COMM_WORLD,reqs(3),ierr)
        endif

!-----

 
      nv=nv+1
      tag=2000+nv
      count=cs1sn

        ! receive north
        if(ibn.eq.0)then
          call mpi_irecv(vn2,count,MPI_REAL,mynorth,tag,   &
                        MPI_COMM_WORLD,reqs(6),ierr)
        endif
 
        ! send south
        if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            vs1(i,k)=v(i,2,k)
          enddo
          enddo
          call mpi_isend(vs1,count,MPI_REAL,mysouth,tag,   &
                         MPI_COMM_WORLD,reqs(5),ierr)
        endif
 
!-----

 
      nv=nv+1
      tag=2000+nv
      count=cs1sn
 
        ! receive south
        if(ibs.eq.0)then
          call mpi_irecv(vs2,count,MPI_REAL,mysouth,tag,   &
                        MPI_COMM_WORLD,reqs(8),ierr)
        endif

        ! send north
        if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            vn1(i,k)=v(i,nj,k)
          enddo
          enddo
          call mpi_isend(vn1,count,MPI_REAL,mynorth,tag,   &
                         MPI_COMM_WORLD,reqs(7),ierr)
        endif
 
!-----

      if(timestats.ge.1) time_mpv1=time_mpv1+mytime()
 
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1v_end(v,vw1,vw2,ve1,ve2,   &
                               vs1,vs2,vn1,vn2,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real v(ib:ie,jb:je+1,kb:ke)
      real vw1(nj+1,nk),vw2(nj+1,nk)
      real ve1(nj+1,nk),ve2(nj+1,nk)
      real vs1(ni,nk),vs2(ni,nk)
      real vn1(ni,nk),vn2(ni,nk)
      integer reqs(8)
 
      integer i,j,k
      integer status(MPI_STATUS_SIZE)

!-----

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj+1
          v(ni+1,j,k)=ve2(j,k)
        enddo
        enddo
      endif
 
      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj+1
          v(0,j,k)=vw2(j,k)
        enddo
        enddo
      endif
 
!-----
 
      if(ibn.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          v(i,nj+2,k)=vn2(i,k)
        enddo
        enddo
      endif
 
      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          v(i,0,k)=vs2(i,k)
        enddo
        enddo
      endif

!-----

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(timestats.ge.1) time_mpv2=time_mpv2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,0,k)=v(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,nj+2,k)=v(1,nj+2,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,0,k)=v(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,nj+2,k)=v(ni,nj+2,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,0,k)=v(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,0,k)=v(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,nj+2,k)=v(0,nj+1,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,nj+2,k)=v(ni+1,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

!-----------------------------------------------------------
!  Mirror b.c. patch
!
!      if(ibn.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=1,nk
!        do i=0,ni+1
!          v(i,nj+1,k)=0.0
!          v(i,nj+2,k)=-v(i,nj  ,k)
!          v(i,nj+3,k)=-v(i,nj-1,k)
!          v(i,nj+4,k)=-v(i,nj-2,k)
!        enddo
!        enddo
!      endif
!
!      if(ibw.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=1,nk
!        do j=0,nj+2
!          v(-2,j,k)=v(3,j,k)
!          v(-1,j,k)=v(2,j,k)
!          v( 0,j,k)=v(1,j,k)
!        enddo
!        enddo
!      endif
!
!-----------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
!----------

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1w_start(w,ww1,ww2,we1,we2,   &
                                 ws1,ws2,wn1,wn2,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real w(ib:ie,jb:je,kb:ke+1)
      real ww1(nj,nk-1),ww2(nj,nk-1)
      real we1(nj,nk-1),we2(nj,nk-1)
      real ws1(ni,nk-1),ws2(ni,nk-1)
      real wn1(ni,nk-1),wn2(ni,nk-1)
      integer reqs(8)
 
      integer i,j,k
      integer tag,count

!-----

      nw=nw+1
      tag=3000+nw
      count=cw1we
 
        ! receive east
        if(ibe.eq.0)then
          call mpi_irecv(we2,count,MPI_REAL,myeast,tag,   &
                        MPI_COMM_WORLD,reqs(2),ierr)
        endif

        ! send west
        if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=2,nk
          do j=1,nj
            ww1(j,k-1)=w(1,j,k)
          enddo
          enddo
          call mpi_isend(ww1,count,MPI_REAL,mywest,tag,   &
                         MPI_COMM_WORLD,reqs(1),ierr)
        endif

!-----

      nw=nw+1
      tag=3000+nw
      count=cw1we
 
        ! receive west
        if(ibw.eq.0)then
          call mpi_irecv(ww2,count,MPI_REAL,mywest,tag,   &
                        MPI_COMM_WORLD,reqs(4),ierr)
        endif

        ! send east
        if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=2,nk
          do j=1,nj
            we1(j,k-1)=w(ni,j,k)
          enddo
          enddo
          call mpi_isend(we1,count,MPI_REAL,myeast,tag,   &
                         MPI_COMM_WORLD,reqs(3),ierr)
        endif

!-----

      nw=nw+1
      tag=3000+nw
      count=cw1sn
 
        ! receive north
        if(ibn.eq.0)then
          call mpi_irecv(wn2,count,MPI_REAL,mynorth,tag,   &
                        MPI_COMM_WORLD,reqs(6),ierr)
        endif

        ! send south
        if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=2,nk
          do i=1,ni
            ws1(i,k-1)=w(i,1,k)
          enddo
          enddo
          call mpi_isend(ws1,count,MPI_REAL,mysouth,tag,   &
                         MPI_COMM_WORLD,reqs(5),ierr)
        endif

!-----

      nw=nw+1
      tag=3000+nw
      count=cw1sn
 
        ! receive south
        if(ibs.eq.0)then
          call mpi_irecv(ws2,count,MPI_REAL,mysouth,tag,   &
                        MPI_COMM_WORLD,reqs(8),ierr)
        endif

        ! send north
        if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=2,nk
          do i=1,ni
            wn1(i,k-1)=w(i,nj,k)
          enddo
          enddo
          call mpi_isend(wn1,count,MPI_REAL,mynorth,tag,   &
                         MPI_COMM_WORLD,reqs(7),ierr)
        endif

!-----

      if(timestats.ge.1) time_mpw1=time_mpw1+mytime()
 
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1w_end(w,ww1,ww2,we1,we2,   &
                               ws1,ws2,wn1,wn2,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real w(ib:ie,jb:je,kb:ke+1)
      real ww1(nj,nk-1),ww2(nj,nk-1)
      real we1(nj,nk-1),we2(nj,nk-1)
      real ws1(ni,nk-1),ws2(ni,nk-1)
      real wn1(ni,nk-1),wn2(ni,nk-1)
      integer reqs(8)
 
      integer i,j,k
      integer status(MPI_STATUS_SIZE)
 
!-----

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=2,nk
        do j=1,nj
          w(ni+1,j,k)=we2(j,k-1)
        enddo
        enddo
      endif
 
      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=2,nk
        do j=1,nj
          w(0,j,k)=ww2(j,k-1)
        enddo
        enddo
      endif
 
!-----
 
      if(ibn.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=2,nk
        do i=1,ni
          w(i,nj+1,k)=wn2(i,k-1)
        enddo
        enddo
      endif
 
      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=2,nk
        do i=1,ni
          w(i,0,k)=ws2(i,k-1)
        enddo
        enddo
      endif

!-----

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(timestats.ge.1) time_mpw2=time_mpw2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,0,k)=w(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,nj+1,k)=w(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,0,k)=w(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,nj+1,k)=w(ni,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,0,k)=w(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,0,k)=w(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,nj+1,k)=w(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,nj+1,k)=w(ni+1,nj,k)
            enddo
          endif
 
        endif
 
      endif

!-----------------------------------------------------------
!  Mirror b.c. patch
!
!      if(ibn.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=2,nk
!        do i=0,ni+1
!          w(i,nj+1,k)=w(i,nj  ,k)
!          w(i,nj+2,k)=w(i,nj-1,k)
!          w(i,nj+3,k)=w(i,nj-2,k)
!        enddo
!        enddo
!      endif
!
!      if(ibw.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=2,nk
!        do j=0,nj+1
!          w(-2,j,k)=w(3,j,k)
!          w(-1,j,k)=w(2,j,k)
!          w( 0,j,k)=w(1,j,k)
!        enddo
!        enddo
!      endif
!
!-----------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
!----------

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2d_start(s,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real, dimension(ib:ie,jb:je) :: s
      real, dimension(3,nj) :: west,newwest,east,neweast
      real, dimension(ni,3) :: south,newsouth,north,newnorth
      integer reqs(8)

      integer i,j
      integer tag,count

!------------------------------------------------

      count=3*nj
      nf=nf+1
      tag=nf

      ! receive east
      if(ibe.eq.0)then
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(2),ierr)
      endif

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,3
          west(i,j)=s(i,j)
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(1),ierr)
      endif

!----------

      count=3*nj
      nf=nf+1
      tag=nf

      ! receive west
      if(ibw.eq.0)then
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(4),ierr)
      endif

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,3
          east(i,j)=s(ni-3+i,j)
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(3),ierr)
      endif

!----------

      count=ni*3
      nf=nf+1
      tag=nf

      ! receive south
      if(ibs.eq.0)then
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(6),ierr)
      endif

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,3
        do i=1,ni
          north(i,j)=s(i,nj-3+j)
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif

!----------

      count=ni*3
      nf=nf+1
      tag=nf

      ! receive north
      if(ibn.eq.0)then
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(8),ierr)
      endif

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,3
        do i=1,ni
          south(i,j)=s(i,j)
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif

!----------

      if(timestats.ge.1) time_mpq1=time_mpq1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2dew_end(s,west,newwest,east,neweast,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real, dimension(ib:ie,jb:je) :: s
      real, dimension(3,nj) :: west,newwest,east,neweast
      integer reqs(8)

      integer i,j
      integer status(MPI_STATUS_SIZE)

!-------------------------------------------------------------------

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,3
          s(ni+i,j)=neweast(i,j)
        enddo
        enddo
      endif

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,3
          s(i-3,j)=newwest(i,j)
        enddo
        enddo
      endif

!----------

      if(timestats.ge.1) time_mpq2=time_mpq2+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2dns_end(s,south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real, dimension(ib:ie,jb:je) :: s
      real, dimension(ni,3) :: south,newsouth,north,newnorth
      integer reqs(8)

      integer i,j
      integer status(MPI_STATUS_SIZE)

!-------------------------------------------------------------------

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,3
        do i=1,ni
          s(i,j-3)=newsouth(i,j)
        enddo
        enddo
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,3
        do i=1,ni
          s(i,nj+j)=newnorth(i,j)
        enddo
        enddo
      endif

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

!----------

      if(timestats.ge.1) time_mpq2=time_mpq2+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2d_corner(s)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real, dimension(ib:ie,jb:je) :: s

      integer reqs(8),status(MPI_STATUS_SIZE)
      real :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer tag,count

!-----

      tag=5001
      count=1

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(2),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        se1=s(ni,1)
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(1),ierr)
      endif

!-----

      tag=5002
      count=1

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(4),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        ne1=s(ni,nj)
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(3),ierr)
      endif

!-----

      tag=5003
      count=1

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        sw1=s(1,1)
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif

!-----

      tag=5004
      count=1

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nw1=s(1,nj)
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif

!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(2),status,ierr)
        s(0,nj+1)=nw2
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(4),status,ierr)
        s(0,0)=sw2
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(6),status,ierr)
        s(ni+1,nj+1)=ne2
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(8),status,ierr)
        s(ni+1,0)=se2
      endif

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine fincomm(qa,q3d,qw1,qw2,qe1,qe2,   &
                                qs1,qs2,qn1,qn2,reqs_q)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: qa,q3d
      real, dimension(3,nj,nk) :: qw1,qw2,qe1,qe2
      real, dimension(ni,3,nk) :: qs1,qs2,qn1,qn2
      integer, dimension(8) :: reqs_q

      integer i,j,k

!-----------------------------------------------------------------

      call comm_3s_end(q3d,qw1,qw2,qe1,qe2,   &
                           qs1,qs2,qn1,qn2,reqs_q)

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        qa(i,j,k)=q3d(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()

!-----------------------------------------------------------------

      return
      end


#endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine prepcorners(s)
      implicit none

      include 'input.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: s

      integer :: i,j

#ifdef MPI
      call getcorner(s)
      call bcs2(s)
#endif
      do j=0,nj+1
      do i=0,ni+1
        s(i,j,0)    = s(i,j,1)
        s(i,j,nk+1) = s(i,j,nk)
      enddo
      enddo

      return
      end


