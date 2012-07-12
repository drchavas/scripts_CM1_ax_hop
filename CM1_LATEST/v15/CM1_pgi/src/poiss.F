

      subroutine poiss(uh,vh,mh,rmh,mf,rmf,pi0,thv0,rho0,rf0,   &
                       def,divx,ppi,uten,vten,wten,             &
                       cfb,cfa,cfc,d1,d2,pdt,deft,rhs,trans,dttmp)

      use singleton
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, dimension(ib:ie) :: uh
      real, dimension(jb:je) :: vh
      real, dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,rho0,rf0
      real, dimension(ib:ie,jb:je,kb:ke) :: def,divx,ppi
      real, dimension(ib:ie+1,jb:je,kb:ke) :: uten
      real, dimension(ib:ie,jb:je+1,kb:ke) :: vten
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wten
      real, dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: cfb
      real, dimension(kpb:kpe) :: cfa,cfc,d1,d2
      complex, dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: pdt,deft
      complex, dimension(ipb:ipe,jpb:jpe) :: rhs,trans
      real dttmp

      integer i,j,k
      real :: tem
      real, dimension(0:nk+1) :: r1
      complex, dimension(nk) :: lgbth,lgbph

!!!
!!!  Get the forcing from the exact numerical representation
!!!  of the right-hand side of the momentum equation
!!!
!!!  The divx term reduces accumulation of numerical errors over time
!!!

      IF(axisymm.eq.1)THEN
        if(myid.eq.0)then
          print *
          print *,'  The anelastic/incompressible solver cannot be '
          print *,'  used with the axisymmetric model(yet)'
          print *
          print *,'  Stopping model ...'
          print *
        endif
        call stopcm1
      ENDIF

      tem = 1.0/dttmp

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        def(i,j,k)=rho0(1,1,k)*(                                  &
                   (uten(i+1,j,k)-uten(i,j,k))*rdx*uh(i)          &
                  +(vten(i,j+1,k)-vten(i,j,k))*rdy*vh(j) )        &
                  +( rf0(i,j,k+1)*wten(i,j,k+1)                   &
                    -rf0(i,j,k  )*wten(i,j,k  ) )*rdz*mh(i,j,k)   &
                  +divx(i,j,k)*tem
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------
!!!
!!! fourier transform the total forcing
!!!

      DO k=1,nk

!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          rhs(i,j)=cmplx(def(i,j,k)*d1(k),0.0)
        enddo
        enddo

        if(imirror.eq.1)then

!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            rhs(ipe+1-i,j)=rhs(i,j)
          enddo
          enddo

        endif

        if(jmirror.eq.1)then

!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            rhs(i,jpe+1-j)=rhs(i,j)
          enddo
          enddo

        endif

        trans=fft(rhs)

!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=jpb,jpe
        do i=ipb,ipe
          deft(i,j,k)=trans(i,j)
        enddo
        enddo

      ENDDO

!-----------------------------------------------------------------------

!!!
!!! solve the tri-diagonal matrix
!!!

!$omp parallel do default(shared)  &
!$omp private(i,j,k,r1,lgbth,lgbph)
      DO j=jpb,jpe
      DO i=ipb,ipe

        if(i.eq.1.and.j.eq.1)then
          r1(nk+1)=0.0
          r1(nk)=0.0
          do k=nk,2,-1
            r1(k-1)=(deft(i,j,k)-cfc(k)*r1(k+1)-cfb(i,j,k)*r1(k))/cfa(k)
          enddo
          do k=1,nk
            pdt(i,j,k)=cmplx( r1(k) , 0.0 )
          enddo
        else
          lgbth(1)=-cfc(1)/cfb(i,j,1)
          lgbph(1)= deft(i,j,1)/cfb(i,j,1)
          do k=2,nk
            lgbth(k)=-cfc(k)/(cfa(k)*lgbth(k-1)+cfb(i,j,k))
            lgbph(k)=(deft(i,j,k)-cfa(k)*lgbph(k-1))/(cfa(k)*lgbth(k-1)+cfb(i,j,k))
          enddo
          pdt(i,j,nk)=lgbph(nk)
          do k=nk-1,1,-1
            pdt(i,j,k)=lgbth(k)*pdt(i,j,k+1)+lgbph(k)
          enddo
        endif

      ENDDO
      ENDDO

!---------------------------------------------

!!!
!!! reverse fourier transform and we're done, not so bad after all!
!!!

      DO k=1,nk

!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=jpb,jpe
        do i=ipb,ipe
          rhs(i,j)=pdt(i,j,k)
        enddo
        enddo

        trans=fft(rhs,inv=.true.)

!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          ppi(i,j,k)=real(trans(i,j))*d2(k)
        enddo
        enddo

      ENDDO

      if(timestats.ge.1) time_poiss=time_poiss+mytime()

!-----------------------------------------------------------------------

      call bcp(ppi)

!-----------------------------------------------------------------------
!  All done.

      return
      end


