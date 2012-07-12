

      subroutine maxmin(izz,jzz,kzz,f,nstat,rstat,amax,amin)
      implicit none

      include 'input.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      integer :: izz,jzz,kzz,nstat
      real, dimension(stat_out) :: rstat
      real, dimension(1-ngxy:izz+ngxy,1-ngxy:jzz+ngxy,1-ngz:kzz+ngz) :: f
      character*6 :: amax,amin

!-----------------------------------------------------------------------

      integer :: i,j,k
      integer :: imax,jmax,kmax,imin,jmin,kmin
      integer, dimension(nk+1) :: imaxt,jmaxt,kmaxt,imint,jmint,kmint
      real, dimension(nk+1) :: tmax,tmin
      real :: fmax,fmin,rmax,rmin
      integer :: loc
      real, dimension(2) :: mmax,nmax,mmin,nmin

!-----------------------------------------------------------------------

!$omp parallel do default(shared)    &
!$omp private(i,j,k)
      do k=1,kzz
        tmax(k)=-99999999.
        tmin(k)= 99999999.
        do j=1,jzz
        do i=1,izz
          if(f(i,j,k).gt.tmax(k))then
            tmax(k)=f(i,j,k)
            imaxt(k)=i
            jmaxt(k)=j
            kmaxt(k)=k
          endif
          if(f(i,j,k).lt.tmin(k))then
            tmin(k)=f(i,j,k)
            imint(k)=i
            jmint(k)=j
            kmint(k)=k
          endif
        enddo
        enddo
      enddo

      fmax=-99999999.
      fmin= 99999999.
      do k=1,kzz
        if(tmax(k).gt.fmax)then
          fmax=tmax(k)
          imax=imaxt(k)
          jmax=jmaxt(k)
          kmax=kmaxt(k)
        endif
        if(tmin(k).lt.fmin)then
          fmin=tmin(k)
          imin=imint(k)
          jmin=jmint(k)
          kmin=kmint(k)
        endif 
      enddo

#ifdef MPI
      mmax(1)=fmax
      mmax(2)=myid
      call MPI_ALLREDUCE(mmax,nmax,1,MPI_2REAL,MPI_MAXLOC,   &
                         MPI_COMM_WORLD,ierr)
      loc=nint(nmax(2))
      imax=imax+(myi-1)*ni
      jmax=jmax+(myj-1)*nj
      call MPI_BCAST(imax,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jmax,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kmax,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)

      mmin(1)=fmin
      mmin(2)=myid
      call MPI_ALLREDUCE(mmin,nmin,1,MPI_2REAL,MPI_MINLOC,   &
                         MPI_COMM_WORLD,ierr)
      loc=nint(nmin(2))
      imin=imin+(myi-1)*ni
      jmin=jmin+(myj-1)*nj
      call MPI_BCAST(imin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jmin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kmin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)

      fmax=nmax(1)
      fmin=nmin(1)

    if(myid.eq.0)then
#endif
      write(6,100) amax,fmax,imax,jmax,kmax,    &
                   amin,fmin,imin,jmin,kmin
100   format(2x,a6,':',1x,f13.6,i5,i5,i5,    &
             4x,a6,':',1x,f13.6,i5,i5,i5)

      nstat = nstat + 1
      rstat(nstat) = fmax
      nstat = nstat + 1
      rstat(nstat) = fmin
#ifdef MPI
    endif
#endif

      if(timestats.ge.1) time_stat=time_stat+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine maxmin2d(izz,jzz,f,nstat,rstat,amax,amin)
      implicit none
        
      include 'input.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      integer :: izz,jzz,nstat
      real, dimension(stat_out) :: rstat
      real, dimension(1-ngxy:izz+ngxy,1-ngxy:jzz+ngxy) :: f
      character*6 :: amax,amin
        
!-----------------------------------------------------------------------
          
      integer :: i,j
      integer :: imax,jmax,imin,jmin
      integer, dimension(jzz) :: imaxt,jmaxt,imint,jmint
      real, dimension(jzz) :: tmax,tmin
      real :: fmax,fmin,rmax,rmin
      integer :: loc
      real, dimension(2) :: mmax,nmax,mmin,nmin
          
!-----------------------------------------------------------------------

!$omp parallel do default(shared)    &
!$omp private(i,j)
      do j=1,jzz
        tmax(j)=-99999999.
        tmin(j)= 99999999.
        do i=1,izz
          if(f(i,j).gt.tmax(j))then
            tmax(j)=f(i,j)
            imaxt(j)=i
            jmaxt(j)=j
          endif
          if(f(i,j).lt.tmin(j))then
            tmin(j)=f(i,j)
            imint(j)=i
            jmint(j)=j
          endif
        enddo
      enddo

      fmax=-99999999.
      fmin= 99999999.
      do j=1,jzz
        if(tmax(j).gt.fmax)then
          fmax=tmax(j)
          imax=imaxt(j)
          jmax=jmaxt(j)
        endif
        if(tmin(j).lt.fmin)then
          fmin=tmin(j)
          imin=imint(j)
          jmin=jmint(j)
        endif
      enddo

#ifdef MPI
      mmax(1)=fmax
      mmax(2)=myid
      call MPI_ALLREDUCE(mmax,nmax,1,MPI_2REAL,MPI_MAXLOC,   &
                         MPI_COMM_WORLD,ierr)
      loc=nint(nmax(2))
      imax=imax+(myi-1)*ni
      jmax=jmax+(myj-1)*nj
      call MPI_BCAST(imax,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jmax,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)

      mmin(1)=fmin
      mmin(2)=myid
      call MPI_ALLREDUCE(mmin,nmin,1,MPI_2REAL,MPI_MINLOC,   &
                         MPI_COMM_WORLD,ierr)
      loc=nint(nmin(2))
      imin=imin+(myi-1)*ni
      jmin=jmin+(myj-1)*nj
      call MPI_BCAST(imin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jmin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)

      fmax=nmax(1)
      fmin=nmin(1)

    if(myid.eq.0)then
#endif
      write(6,100) amax,fmax,imax,jmax,1,    &
                   amin,fmin,imin,jmin,1
100   format(2x,a6,':',1x,f13.6,i5,i5,i5,    &
             4x,a6,':',1x,f13.6,i5,i5,i5)

      nstat = nstat + 1
      rstat(nstat) = fmax
      nstat = nstat + 1
      rstat(nstat) = fmin
#ifdef MPI
    endif
#endif

      if(timestats.ge.1) time_stat=time_stat+mytime()

      return
      end


