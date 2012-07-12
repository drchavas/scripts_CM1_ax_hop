

      subroutine getqli(ip,q,ql,qi)
      implicit none

      include 'input.incl'
      include 'timestat.incl'

      integer :: ip
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q
      real, dimension(ib:ie,jb:je,kb:ke) :: ql,qi

      integer :: i,j,k,n

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1-ip,nj+ip
      do i=1-ip,ni+ip
        ql(i,j,k)=0.0
        qi(i,j,k)=0.0
      enddo
      enddo
      enddo

      do n=nql1,nql2
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1-ip,nj+ip
        do i=1-ip,ni+ip
          ql(i,j,k)=ql(i,j,k)+q(i,j,k,n)
        enddo
        enddo
        enddo
      enddo

      IF(iice.eq.1)THEN
        do n=nqs1,nqs2
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1-ip,nj+ip
          do i=1-ip,ni+ip
            qi(i,j,k)=qi(i,j,k)+q(i,j,k,n)
          enddo
          enddo
          enddo
        enddo
      ENDIF

      if(timestats.ge.1) time_misc=time_misc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefq(rmax,asq,ruh,rvh,rmh,rho,q3d)
      implicit none

      include 'input.incl'
      include 'timestat.incl'

      real rmax
      real*8 asq
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho,q3d

      integer i,j,k
      real*8 t1,t2,t3
      real*8 a1,a2,tem
      real*8, dimension(nj) :: budj
      real*8, dimension(nk) :: budk

!----------------------------------------------------------------------

      tem = dx*dy*dz

      IF(pdscheme.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(j)
        do j=1,nj
          budj(j)=0.0d0
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k,t1,t2,t3,a1,a2)
        do j=1,nj
        do i=1,ni
          t1=0.0d0
          t2=0.0d0
          a1=0.0d0
          a2=0.0d0
          do k=1,nk
            t1=t1+rho(i,j,k)*q3d(i,j,k)
            a1=a1+rho(i,j,k)*q3d(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)
!!!            q3d(i,j,k)=max(0.0,q3d(i,j,k))
            if(q3d(i,j,k).lt.rmax) q3d(i,j,k)=0.0
            t2=t2+rho(i,j,k)*q3d(i,j,k)
          enddo
          t3=(t1+1.0d-20)/(t2+1.0d-20)
          if(t3.lt.0.0) t3=1.0d0
          do k=1,nk
            q3d(i,j,k)=t3*q3d(i,j,k)
            a2=a2+rho(i,j,k)*q3d(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)
          enddo
          budj(j)=budj(j)+a2-a1
        enddo
        enddo

        do j=1,nj
          asq=asq+budj(j)*tem
        enddo

      ELSE

!$omp parallel do default(shared)  &
!$omp private(k)
        do k=1,nk
          budk(k)=0.0d0
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k,a1,a2)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          a1=rho(i,j,k)*q3d(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)
!!!          q3d(i,j,k)=max(0.0,q3d(i,j,k))
          if(q3d(i,j,k).lt.rmax) q3d(i,j,k)=0.0
          a2=rho(i,j,k)*q3d(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)
          budk(k)=budk(k)+a2-a1
        enddo
        enddo
        enddo

        do k=1,nk
          asq=asq+budk(k)*tem
        enddo

      ENDIF

!----------------------------------------------------------------------

      if(timestats.ge.1) time_misc=time_misc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine calcprs(pi0,prs,pp3d)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0
      real, dimension(ib:ie,jb:je,kb:ke) :: prs,pp3d
 
      integer i,j,k
 
!----------------------------------------------------------------------
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
      enddo
      enddo
      enddo
 
!----------------------------------------------------------------------
 
      if(timestats.ge.1) time_misc=time_misc+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calcrho(pi0,th0,rho,prs,pp3d,th3d,q3d)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,th0
      real, dimension(ib:ie,jb:je,kb:ke) :: rho,prs,pp3d,th3d
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d

      integer i,j,k

!----------------------------------------------------------------------

      IF(imoist.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          rho(i,j,k)=prs(i,j,k)                         &
             /( rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))     &
                  *(1.0+max(0.0,q3d(i,j,k,nqv))*reps) )
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          rho(i,j,k)=prs(i,j,k)   &
             /(rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k)))
        enddo
        enddo
        enddo

      ENDIF

!----------------------------------------------------------------------

      if(timestats.ge.1) time_misc=time_misc+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calcdbz(rho,qr,qs,qg,dbz)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'goddard.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: rho,qr,qs,qg,dbz

      integer :: i,j,k
      real :: n0r,n0g,n0s,rhor,rhog,rhos,gamma,zer,zeg,zes
      real, parameter :: epp = 1.0e-8

  IF(ptype.eq.2)THEN

    rhor = 1000.0 * ROQR
    rhog = 1000.0 * ROQG
    rhos = 1000.0 * ROQS

    n0r = 1.0e8 * TNW
    n0g = 1.0e8 * TNG
    n0s = 1.0e8 * TNSS

!!!    print *,'  rhor,rhog,rhos = ',rhor,rhog,rhos
!!!    print *,'  n0r,n0g,n0s    = ',n0r,n0g,n0s

!$omp parallel do default(shared)  &
!$omp private(i,j,k,gamma,zer,zeg,zes)
    do k=1,nk
    do j=1,nj
    do i=1,ni

    if(qr(i,j,k).ge.epp)then
      !--- rain ---
      gamma=(3.14159*n0r*rhor/(rho(i,j,k)*qr(i,j,k)))**0.25
      zer=720.0*n0r*(gamma**(-7))
    else
      zer=0.0
    endif

    if(qg(i,j,k).ge.epp)then
      !--- graupel/hail ---
      gamma=(3.14159*n0g*rhog/(rho(i,j,k)*qg(i,j,k)))**0.25
      zeg=720.0*n0g*(gamma**(-7))*((rhog/rhor)**2)*0.224
    else
      zeg=0.0
    endif

    if(qs(i,j,k).ge.epp)then
      !--- snow ---
      gamma=(3.14159*n0s*rhos/(rho(i,j,k)*qs(i,j,k)))**0.25
      zes=720.0*n0s*(gamma**(-7))*((rhos/rhor)**2)*0.224
    else
      zes=0.0
    endif

      !--- dbz ---

    if( (zer+zeg+zes).gt.1.0e-18 )then
      dbz(i,j,k)=10.0*log10((zer+zeg+zes)*1.0e18)
    else
      dbz(i,j,k)=0.0
    endif

    enddo
    enddo
    enddo

  ELSE

    write(outfile,*)
    write(outfile,*) ' ptype = ',ptype
    write(outfile,*)
    write(outfile,*) ' calcdbz is not valid for this value of ptype'
    write(outfile,*)
    call stopcm1

  ENDIF

      if(timestats.ge.1) time_write=time_write+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calcuh(uf,vf,zh,zf,ua,va,wa,uh,zeta)
      implicit none

      ! Subroutine to calculate vertically integrated updraft helicity
      ! Reference:  Kain et al, 2008, WAF, p 931

      include 'input.incl'
      include 'timestat.incl'

      ! note:  need zh,zf Above Ground Level

      real, intent(in), dimension(ib:ie+1) :: uf
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je) :: uh
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: zeta

      real, parameter :: zz0 = 1000.0     ! bottom of integration layer (m AGL)
      real, parameter :: zzt = 6000.0     ! top of integration layer (m AGL)

      integer :: i,j,k
      real :: wbar,zbar

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
    DO k=1,nk
    DO j=1,nj+1
    DO i=1,ni+1
      zeta(i,j,k) = (va(i,j,k)-va(i-1,j,k))*rdx*uf(i)   &
                   -(ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)
    ENDDO
    ENDDO
    ENDDO

!$omp parallel do default(shared)  &
!$omp private(i,j,k,wbar,zbar)
    DO j=1,nj
    DO i=1,ni
      uh(i,j) = 0.0
      DO k=1,nk
        IF( zh(i,j,k).ge.zz0 .and. zh(i,j,k).le.zzt )THEN
          ! note:  only consider cyclonically rotating updrafts
          !        (so, w and zeta must both be positive)
          wbar = max( 0.0 , 0.5*(wa(i,j,k)+wa(i,j,k+1)) )
          zbar = max( 0.0 , 0.25*(zeta(i,j,k)+zeta(i+1,j,k)   &
                                 +zeta(i,j+1,k)+zeta(i+1,j+1,k)) )
          uh(i,j) = uh(i,j) + (min(zf(i,j,k+1),zzt)-max(zf(i,j,k),zz0))*wbar*zbar
        ENDIF
      ENDDO
    ENDDO
    ENDDO

      return
      end subroutine calcuh


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calcvort(xh,xf,uf,vf,zh,mf,zf,ua,va,wa,xvort,yvort,zvort,tem)
      implicit none

      ! Subroutine to calculate 3 components of vorticity
      ! at scalar points.

      include 'input.incl'
      include 'timestat.incl'

      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,zf
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: xvort,yvort,zvort,tem

      integer :: i,j,k

        !---
        ! x-vort:
        tem=0.0
        if(axisymm.eq.0)then
          ! Cartesian grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=2,nk
          do j=1,nj+1
          do i=1,ni
            tem(i,j,k) = (wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j)   &
                        -(va(i,j,k)-va(i,j,k-1))*rdz*0.5*(mf(i,j-1,k)+mf(i,j,k))
          enddo
          enddo
          enddo
          IF( (iturb.ge.1.and.bcturbu.eq.2).or.(dns.eq.1.and.bc_wind.eq.1) )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do j=1,nj+1
            do i=1,ni
              tem(i,j,1)=tem(i,j,2)
              tem(i,j,nk+1)=tem(i,j,nk)
            enddo
            enddo
          ENDIF
          IF( (iturb.ge.1.and.bcturbu.eq.3).or.(dns.eq.1.and.bc_wind.eq.2).or.(ipbl.eq.1) )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do j=1,nj+1
            do i=1,ni
              tem(i,j,1)=-2.0*va(i,j,1)*rdz*0.5*(mf(i,j-1,1)+mf(i,j,1))
              tem(i,j,nk+1)=2.0*va(i,j,nk)*rdz*0.5*(mf(i,j-1,nk+1)+mf(i,j,nk+1))
            enddo
            enddo
          ENDIF
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            xvort(i,j,k) = 0.25*(tem(i,j,k)+tem(i,j+1,k)+tem(i,j,k+1)+tem(i,j+1,k+1))
          enddo
          enddo
          enddo
        else
          ! Axisymmetric grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            tem(i,j,k) = -(va(i,j,k)-va(i,j,k-1))*rdz*0.5*(mf(i,j-1,k)+mf(i,j,k))
          enddo
          enddo
          enddo
          IF( (iturb.ge.1.and.bcturbu.eq.2).or.(dns.eq.1.and.bc_wind.eq.1) )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do j=1,nj
            do i=1,ni
              tem(i,j,1)=tem(i,j,2)
              tem(i,j,nk+1)=tem(i,j,nk)
            enddo
            enddo
          ENDIF
          IF( (iturb.ge.1.and.bcturbu.eq.3).or.(dns.eq.1.and.bc_wind.eq.2).or.(ipbl.eq.1) )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do j=1,nj
            do i=1,ni
              tem(i,j,1)=-2.0*va(i,j,1)*rdz*0.5*(mf(i,j-1,1)+mf(i,j,1))
              tem(i,j,nk+1)=2.0*va(i,j,nk)*rdz*0.5*(mf(i,j-1,nk+1)+mf(i,j,nk+1))
            enddo
            enddo
          ENDIF
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            xvort(i,j,k) = 0.5*(tem(i,j,k)+tem(i,j,k+1))
          enddo
          enddo
          enddo
        endif


        !---
        ! y-vort:
        tem=0.0
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni+1
            tem(i,j,k) = (ua(i,j,k)-ua(i,j,k-1))*rdz*0.5*(mf(i-1,j,k)+mf(i,j,k))   &
                        -(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)
          enddo
          enddo
          enddo
          IF( (iturb.ge.1.and.bcturbu.eq.2).or.(dns.eq.1.and.bc_wind.eq.1) )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do j=1,nj
            do i=1,ni+1
              tem(i,j,1)=tem(i,j,2)
              tem(i,j,nk+1)=tem(i,j,nk)
            enddo
            enddo
          ENDIF
          IF( (iturb.ge.1.and.bcturbu.eq.3).or.(dns.eq.1.and.bc_wind.eq.2).or.(ipbl.eq.1) )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do j=1,nj
            do i=1,ni+1
              tem(i,j,1)=2.0*ua(i,j,1)*rdz*0.5*(mf(i-1,j,1)+mf(i,j,1))
              tem(i,j,nk+1)=-2.0*ua(i,j,nk)*rdz*0.5*(mf(i-1,j,nk+1)+mf(i,j,nk+1))
            enddo
            enddo
          ENDIF
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            yvort(i,j,k) = 0.25*(tem(i,j,k)+tem(i+1,j,k)+tem(i,j,k+1)+tem(i+1,j,k+1))
          enddo
          enddo
          enddo


        !---
        ! z-vort:
        tem=0.0
        if(axisymm.eq.0)then
          ! Cartesian grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
            do j=1,nj+1
            do i=1,ni+1
              tem(i,j,k) = (va(i,j,k)-va(i-1,j,k))*rdx*uf(i)   &
                          -(ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)
            enddo
            enddo
            do j=1,nj
            do i=1,ni
              zvort(i,j,k) = 0.25*(tem(i,j,k)+tem(i+1,j,k)+tem(i,j+1,k)+tem(i+1,j+1,k))
            enddo
            enddo
          enddo
        else
          ! Axisymmetric grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
            do j=1,nj
            tem(1,j,k) = 0.0
            do i=2,ni+1
              tem(i,j,k) = (va(i,j,k)*xh(i)-va(i-1,j,k)*xh(i-1))*rdx*uf(i)/xf(i)
            enddo
            enddo
            do j=1,nj
            do i=1,ni
              zvort(i,j,k) = 0.5*(tem(i,j,k)+tem(i+1,j,k))
            enddo
            enddo
          enddo
        endif

      return
      end subroutine calcvort


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calccpch(zf,th0,qv0,cpc,cph,tha,qa)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      real, intent(in),    dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in),    dimension(ib:ie,jb:je,kb:ke) :: th0,qv0
      real, intent(inout), dimension(ib:ie,jb:je) :: cpc,cph
      real, intent(in),    dimension(ib:ie,jb:je,kb:ke) :: tha
      real, intent(in),    dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa

      integer :: i,j,k,n
      real :: ql
      real, dimension(nk) :: bb

      ! defines top of cold pool / location to stop calculation of C
      real, parameter :: bcrit = 0.0

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,ql,bb)
    DO j=1,nj
    DO i=1,ni
      cpc(i,j) = 0.0
      cph(i,j) = 0.0
    ! only calculate cpc/cph if surface theta-pert is at least -1 K
    IF( tha(i,j,1).le.-1.0 )THEN
      bb = 0.0
      if(imoist.eq.1)then
      do k=1,nk
        ql = 0.0
        do n=nql1,nql2
          ql=ql+qa(i,j,k,n)
        enddo
        if(iice.eq.1)then
          do n=nqs1,nqs2
            ql=ql+qa(i,j,k,n)
          enddo
        endif
        bb(k) = g*( tha(i,j,k)/th0(i,j,k)      &
                   + repsm1*(qa(i,j,k,nqv)-qv0(i,j,k)) - ql )
      enddo
      endif
      k = 1
      do while( bb(k).lt.bcrit .and. k.lt.nk )
        if( cpc(i,j).lt.0.0 ) cpc(i,j) = 0.0
        cpc(i,j) = cpc(i,j) - 2.0*bb(k)*(zf(i,j,k+1)-zf(i,j,k))
        k = k + 1
      enddo
      if( cpc(i,j).gt.0.0 )then
        cpc(i,j) = sqrt(cpc(i,j))
        cph(i,j) = zf(i,j,k-1) + (zf(i,j,k)-zf(i,j,k-1))*(bcrit-bb(k-1))   &
                                                        /(bb(k)-bb(k-1))
        ! account for terrain:
        cph(i,j) = cph(i,j) - zf(i,j,1)
      endif
    ENDIF
    ENDDO
    ENDDO

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calccref(cref,dbz)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(inout), dimension(ib:ie,jb:je) :: cref
      real, intent(in),    dimension(ib:ie,jb:je,kb:ke) :: dbz

      integer :: i,j,k

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do j=1,nj
      do i=1,ni
        cref(i,j)=0.0
        do k=1,nk
          cref(i,j)=max(cref(i,j),dbz(i,j,k))
        enddo
      enddo
      enddo


      if(timestats.ge.1) time_write=time_write+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calcthe(zh,pi0,th0,the,rh,prs,ppi,tha,qa)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: zh,pi0,th0
      real, dimension(ib:ie,jb:je,kb:ke) :: the,rh,prs,ppi,tha
      real, dimension(ib:ie,jb:je,kb:ke,numq) :: qa

      integer i,j,k,n
      real tx,cpm
      real rslf
      real, parameter :: l0 = 2.555e6

! Reference:  Bryan, 2008, MWR, p. 5239

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,tx,cpm)
      do j=1,nj
      do k=1,nk
      do i=1,ni
        if(zh(i,j,k).le.10000.)then
          tx=(th0(i,j,k)+tha(i,j,k))*(pi0(i,j,k)+ppi(i,j,k))
          cpm=cp
          the(i,j,k)=tx                                              &
            *((p00*(1.0+qa(i,j,k,nqv)*reps)/prs(i,j,k))**(rd/cpm))   &
            *(rh(i,j,k)**(-qa(i,j,k,nqv)*rv/cpm))                    &
            *exp(l0*qa(i,j,k,nqv)/(cpm*tx))
        else
          the(i,j,k)=the(i,j,k-1)
        endif
      enddo
      enddo
      enddo

      if(timestats.ge.1) time_stat=time_stat+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine cloud(nstat,rstat,zh,qci)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif
 
      integer nstat
      real, dimension(stat_out) :: rstat
      real, dimension(ib:ie,jb:je,kb:ke) :: zh
      real, dimension(ib:ie,jb:je,kb:ke) :: qci

      integer i,j,k
      real qcbot(nk),qctop(nk),bot,top,var

!$omp parallel do default(shared)  &
!$omp private(k)
      do k=1,nk
        qcbot(k)=maxz
        qctop(k)=0.0
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          if(qci(i,j,k).ge.clwsat)then
            qctop(k)=max(qctop(k),zh(i,j,k))
            qcbot(k)=min(qcbot(k),zh(i,j,k))
          endif
        enddo
        enddo
      enddo

      top=0.0
      do k=1,nk
        top=max(top,qctop(k))
      enddo

      bot=maxz
      do k=1,nk
        bot=min(bot,qcbot(k))
      enddo

#ifdef MPI
      call MPI_REDUCE(bot,var,1,MPI_REAL,MPI_MIN,0,MPI_COMM_WORLD,ierr)
      bot=var
      call MPI_REDUCE(top,var,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      top=var
      if(myid.eq.0)then
#endif

      if(bot.eq.maxz) bot=0.0

      write(6,100) 'QCTOP ',top,1,1,1,   &
                   'QCBOT ',bot,1,1,1
100   format(2x,a6,':',1x,f13.6,i5,i5,i5,   &
             4x,a6,':',1x,f13.6,i5,i5,i5)
 
      nstat = nstat + 1
      rstat(nstat) = top
      nstat = nstat + 1
      rstat(nstat) = bot

#ifdef MPI
      endif
#endif

      if(timestats.ge.1) time_stat=time_stat+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine vertvort(nstat,rstat,xh,xf,uf,vf,zh,ua,va)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      integer nstat
      real, dimension(stat_out) :: rstat
      real, dimension(ib:ie) :: xh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: zh
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va

      integer i,j,k,n,n1km,n2km,n3km,n4km,n5km
      real vort,vmax,var
      character*6 text

!-----
!  note:  does not account for terrain

      n1km=0
      n2km=0
      n3km=0
      n4km=0
      n5km=0

      do k=nk,1,-1
        if(zh(1,1,k).ge.1000.0) n1km=k
        if(zh(1,1,k).ge.2000.0) n2km=k
        if(zh(1,1,k).ge.3000.0) n3km=k
        if(zh(1,1,k).ge.4000.0) n4km=k
        if(zh(1,1,k).ge.5000.0) n5km=k
      enddo

      do n=1,6
        if(n.eq.1)then
          k=1
          text='VORSFC'
        elseif(n.eq.2)then
          k=n1km
          text='VOR1KM'
        elseif(n.eq.3)then
          k=n2km
          text='VOR2KM'
        elseif(n.eq.4)then
          k=n3km
          text='VOR3KM'
        elseif(n.eq.5)then
          k=n4km
          text='VOR4KM'
        elseif(n.eq.6)then
          k=n5km
          text='VOR5KM'
        endif
        vmax=-9999999.
      IF( axisymm.ne.1 )THEN
        do j=1+ibs,nj+1-ibn
        do i=1+ibw,ni+1-ibe
          vort=(va(i,j,k)-va(i-1,j,k))*rdx*uf(i)   &
              -(ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)
          vmax=max(vmax,vort)
        enddo
        enddo
      ELSE
        do j=1,nj+1
        do i=2,ni+1
          vort=(xh(i)*va(i,j,k)-xh(i-1)*va(i-1,j,k))*rdx*uf(i)/xf(i)
          vmax=max(vmax,vort)
        enddo
        enddo
      ENDIF
#ifdef MPI
        call MPI_REDUCE(vmax,var,1,MPI_REAL,MPI_MAX,0,   &
                        MPI_COMM_WORLD,ierr)
        vmax=var
        if(myid.eq.0)then
#endif
        write(6,100) text,vmax
        nstat = nstat + 1
        rstat(nstat) = vmax
#ifdef MPI
        endif
#endif
      ENDDO

100   format(2x,a6,':',1x,e13.6)

      if(timestats.ge.1) time_stat=time_stat+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 
      subroutine calccfl(nstat,rstat,dt,acfl,uf,vf,mf,ua,va,wa,writeit)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif
 
      integer nstat
      real, dimension(stat_out) :: rstat
      real :: dt
      real*8 :: acfl
      real, dimension(ib:ie+1) :: uf
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa
      integer :: writeit
 
      integer i,j,k
      integer imax,jmax,kmax
      integer imaxt(nk),jmaxt(nk),kmaxt(nk)
      real dtdx,dtdy,dtdz,cfl(nk),fmax
      integer :: loc
      real, dimension(2) :: mmax,nmax
 
      dtdx=dt*rdx
      dtdy=dt*rdy
      dtdz=dt*rdz

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
        cfl(k)=-99999.0
        do j=1,nj
        do i=1,ni+1
          if( abs(ua(i,j,k)*dtdx*uf(i)).gt.cfl(k) )then
            cfl(k)=abs(ua(i,j,k)*dtdx*uf(i))
            imaxt(k)=i
            jmaxt(k)=j
            kmaxt(k)=k
          endif
        enddo
        enddo
      if(axisymm.eq.0)then
        do j=1,nj+1
        do i=1,ni
          if( abs(va(i,j,k)*dtdy*vf(j)).gt.cfl(k) )then
            cfl(k)=abs(va(i,j,k)*dtdy*vf(j))
            imaxt(k)=i
            jmaxt(k)=j
            kmaxt(k)=k
          endif
        enddo
        enddo
      endif
        do j=1,nj
        do i=1,ni
          if( abs(wa(i,j,k)*dtdz*mf(i,j,k)).gt.cfl(k) )then
            cfl(k)=abs(wa(i,j,k)*dtdz*mf(i,j,k))
            imaxt(k)=i
            jmaxt(k)=j
            kmaxt(k)=k
          endif
        enddo
        enddo
      enddo

      fmax=-99999999.
      imax=1
      jmax=1
      kmax=1
      do k=1,nk
        if(cfl(k).gt.fmax)then
          fmax=cfl(k)
          imax=imaxt(k)
          jmax=jmaxt(k)
          kmax=kmaxt(k)
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
      fmax=nmax(1)

      if(myid.eq.0)then
#endif

    IF(writeit.eq.1)THEN
      nstat = nstat + 1
      IF( adapt_dt.eq.1 )THEN
        write(6,100) 'CFLMAX',sngl(acfl),imax,jmax,kmax
        rstat(nstat) = sngl(acfl)
        acfl = 0.0
      ELSE
        write(6,100) 'CFLMAX',fmax,imax,jmax,kmax
        rstat(nstat) = fmax
      ENDIF
100   format(2x,a6,':',1x,f13.6,i5,i5,i5)
    ENDIF

#ifdef MPI
      endif
#endif

      cflmax = fmax

      if(fmax.ge.1.50) stopit=.true.
 
      if(timestats.ge.1) time_stat=time_stat+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calcksmax(nstat,rstat,dt,uh,vh,mf,kmh,kmv,khh,khv)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      integer nstat
      real, dimension(stat_out) :: rstat
      real :: dt
      real, dimension(ib:ie) :: uh
      real, dimension(jb:je) :: vh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv

      integer i,j,k
      integer imaxh,jmaxh,kmaxh
      integer imaxv,jmaxv,kmaxv
      integer imaxth(nk),jmaxth(nk),kmaxth(nk)
      integer imaxtv(nk),jmaxtv(nk),kmaxtv(nk)
      real dtdx,dtdy,dtdz,tem,ksh(nk),ksv(nk),fhmax,fvmax
      integer :: loc
      real, dimension(2) :: mmax,nmax

    IF(iturb.ge.1)THEN

      dtdx=dt*rdx*rdx
      dtdy=dt*rdy*rdy
      dtdz=dt*rdz*rdz

!$omp parallel do default(shared)  &
!$omp private(i,j,k,tem)
      do k=2,nk
        ksh(k)=-99999.0
        ksv(k)=-99999.0
        do j=1,nj
        do i=1,ni
!!!          tem = max( abs(kmh(i,j,k))*dtdx*uh(i)*uh(i) ,   &
!!!                     abs(khh(i,j,k))*dtdx*uh(i)*uh(i) )
          tem = khh(i,j,k)*dtdx*uh(i)*uh(i)
          if( tem.gt.ksh(k) )then
            ksh(k)=tem
            imaxth(k)=i
            jmaxth(k)=j
            kmaxth(k)=k
          endif
!!!          tem = max( abs(kmh(i,j,k))*dtdy*vh(j)*vh(j) ,   &
!!!                     abs(khh(i,j,k))*dtdy*vh(j)*vh(j) )
          tem = khh(i,j,k)*dtdy*vh(j)*vh(j)
          if( tem.gt.ksh(k) )then
            ksh(k)=tem
            imaxth(k)=i
            jmaxth(k)=j
            kmaxth(k)=k
          endif
!!!          tem = max( abs(kmv(i,j,k))*dtdz*mf(i,j,k)*mf(i,j,k) ,   &
!!!                     abs(khv(i,j,k))*dtdz*mf(i,j,k)*mf(i,j,k) )
          tem = khv(i,j,k)*dtdz*mf(i,j,k)*mf(i,j,k)
          if( tem.gt.ksv(k) )then
            ksv(k)=tem
            imaxtv(k)=i
            jmaxtv(k)=j
            kmaxtv(k)=k
          endif
        enddo
        enddo
      enddo

      fhmax=-99999999.
      fvmax=-99999999.
      imaxh=1
      jmaxh=1
      kmaxh=1
      imaxv=1
      jmaxv=1
      kmaxv=1
      do k=2,nk
        if(ksh(k).gt.fhmax)then
          fhmax=ksh(k)
          imaxh=imaxth(k)
          jmaxh=jmaxth(k)
          kmaxh=kmaxth(k)
        endif
        if(ksv(k).gt.fvmax)then
          fvmax=ksv(k)
          imaxv=imaxtv(k)
          jmaxv=jmaxtv(k)
          kmaxv=kmaxtv(k)
        endif
      enddo

    ELSE

      fhmax = 0.0
      fvmax = 0.0
      imaxh = 0
      jmaxh = 0
      kmaxh = 0
      imaxv = 0
      jmaxv = 0
      kmaxv = 0

    ENDIF

#ifdef MPI
      mmax(1)=fhmax
      mmax(2)=myid
      call MPI_ALLREDUCE(mmax,nmax,1,MPI_2REAL,MPI_MAXLOC,   &
                         MPI_COMM_WORLD,ierr)
      loc=nint(nmax(2))
      imaxh=imaxh+(myi-1)*ni
      jmaxh=jmaxh+(myj-1)*nj
      call MPI_BCAST(imaxh,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jmaxh,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kmaxh,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      fhmax=nmax(1)

      mmax(1)=fvmax
      mmax(2)=myid
      call MPI_ALLREDUCE(mmax,nmax,1,MPI_2REAL,MPI_MAXLOC,   &
                         MPI_COMM_WORLD,ierr)
      loc=nint(nmax(2))
      imaxv=imaxv+(myi-1)*ni
      jmaxv=jmaxv+(myj-1)*nj
      call MPI_BCAST(imaxv,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jmaxv,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kmaxv,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      fvmax=nmax(1)

      if(myid.eq.0)then
#endif

      write(6,100) 'KSHMAX',fhmax,imaxh,jmaxh,kmaxh

      nstat = nstat + 1
      rstat(nstat) = fhmax

      write(6,100) 'KSVMAX',fvmax,imaxv,jmaxv,kmaxv

      nstat = nstat + 1
      rstat(nstat) = fvmax

100   format(2x,a6,':',1x,f13.6,i5,i5,i5)

#ifdef MPI
      endif
#endif

      if(timestats.ge.1) time_stat=time_stat+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getrmw(nstat,rstat,xh,va)
      implicit none

      include 'input.incl'
      include 'timestat.incl'

      integer :: nstat
      real, dimension(stat_out) :: rstat
      real, dimension(ib:ie) :: xh
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va

      integer :: i,k,imax,jmax,kmax
      real :: rmax,vmax

      integer, dimension(nk) :: imaxt,kmaxt
      real, dimension(nk) :: rmaxt,vmaxt

!$omp parallel do default(shared)  &
!$omp private(i,k)
      do k=1,nk
        vmaxt(k) = 0.0
        do i=1,ni
          IF( va(i,1,k).gt.vmaxt(k) )THEN
            vmaxt(k) = va(i,1,k)
            rmaxt(k) = xh(i)
            imaxt(k) = i
            kmaxt(k) = k
          ENDIF
        enddo
      enddo

      vmax = 0.0
      do k=1,nk
        IF( vmaxt(k).gt.vmax )THEN
          vmax = vmaxt(k)
          rmax = rmaxt(k)
          imax = imaxt(k)
          kmax = kmaxt(k)
        ENDIF
      enddo

      jmax = 1

      write(6,131) 'RMW   ',rmax,imax,jmax,kmax
131   format(2x,a6,':',1x,f13.6,i5,i5,i5)
      nstat = nstat + 1
      rstat(nstat) = rmax

      if(timestats.ge.1) time_stat=time_stat+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calcmass(nstat,rstat,ruh,rvh,rmh,rho)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif
 
      integer nstat
      real, dimension(stat_out) :: rstat
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho
 
      integer i,j,k
      real*8 foo(nk),tmass,var
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
        foo(k)=0.0d0
        do j=1,nj
        do i=1,ni
          foo(k)=foo(k)+rho(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)
        enddo
        enddo
      enddo
 
      tmass=0.0d0
      do k=1,nk
        tmass=tmass+foo(k)
      enddo

#ifdef MPI
      var=0.0d0
      call MPI_REDUCE(tmass,var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,   &
                      MPI_COMM_WORLD,ierr)
      tmass=var
      if(myid.eq.0)then
#endif

      tmass=tmass*(dx*dy*dz)
 
      write(6,100) 'TMASS ',tmass
100   format(2x,a6,':',1x,e13.6)
 
      nstat = nstat + 1
      rstat(nstat) = tmass

#ifdef MPI
      endif
#endif
 
      if(timestats.ge.1) time_stat=time_stat+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine totmois(nstat,rstat,train,ruh,rvh,rmh,qv,ql,qi,rho)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      integer nstat
      real, dimension(stat_out) :: rstat
      real*8 :: train
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: qv,ql,qi,rho
 
      integer i,j,k
      real*8 foo(nk),tmass,var,totrain

!$omp parallel do default(shared)  &
!$omp private(k)
      do k=1,nk
        foo(k)=0.0d0
      enddo
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          foo(k)=foo(k)+rho(i,j,k)*(qv(i,j,k)+ql(i,j,k)+qi(i,j,k))*ruh(i)*rvh(j)*rmh(i,j,k)
        enddo
        enddo
      enddo
 
      tmass=0.0d0
      do k=1,nk
        tmass=tmass+foo(k)
      enddo

      totrain=train
#ifdef MPI
      totrain=0.0d0
      call MPI_REDUCE(train,totrain,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,   &
                      MPI_COMM_WORLD,ierr)
#endif

#ifdef MPI
      var=0.0d0
      call MPI_REDUCE(tmass,var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,   &
                      MPI_COMM_WORLD,ierr)
      tmass=var
      if(myid.eq.0)then
#endif

      tmass=tmass*(dx*dy*dz)+totrain

      write(6,100) 'TMOIS ',tmass
100   format(2x,a6,':',1x,e13.6)
 
      nstat = nstat + 1
      rstat(nstat) = tmass

#ifdef MPI
      endif
#endif
 
      if(timestats.ge.1) time_stat=time_stat+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine totq(nstat,rstat,ruh,rvh,rmh,q,rho,aname)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      integer nstat
      real, dimension(stat_out) :: rstat
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: q,rho
      character*6 aname

      integer i,j,k
      real*8 foo(nk),tmass,var

!$omp parallel do default(shared)  &
!$omp private(k)
      do k=1,nk
        foo(k)=0.0d0
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        foo(k)=foo(k)+rho(i,j,k)*q(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)
      enddo
      enddo
      enddo

      tmass=0.0d0
      do k=1,nk
        tmass=tmass+foo(k)
      enddo

#ifdef MPI
      var=0.0d0
      call MPI_REDUCE(tmass,var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,   &
                      MPI_COMM_WORLD,ierr)
      tmass=var
      if(myid.eq.0)then
#endif

      tmass=tmass*(dx*dy*dz)

      write(6,100) aname,tmass
100   format(2x,a6,':',1x,e13.6)

      nstat = nstat + 1
      rstat(nstat) = tmass

#ifdef MPI
      endif
#endif

      if(timestats.ge.1) time_stat=time_stat+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine calcener(nstat,rstat,ruh,rvh,zh,rmh,pi0,th0,rho,ua,va,wa,ppi,tha,   &
                          qv,ql,qi,vr)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif
 
      integer nstat
      real, dimension(stat_out) :: rstat
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,rmh,pi0,th0
      real, dimension(ib:ie,jb:je,kb:ke) :: rho
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,tha,qv,ql,qi,vr
 
      integer i,j,k
      real*8 :: u,v,w,tmp,qtot,ek,ei,ep,et,le,var,tem
      real*8, dimension(nk) :: foo1,foo2,foo3,foo4

!$omp parallel do default(shared)  &
!$omp private(i,j,k,u,v,w,tmp,qtot,tem)
      do k=1,nk
        foo1(k)=0.0d0      ! = ek
        foo2(k)=0.0d0      ! = ei
        foo3(k)=0.0d0      ! = ep
        foo4(k)=0.0d0      ! = le
        do j=1,nj
        do i=1,ni
          tem=ruh(i)*rvh(j)*rmh(i,j,k)
          u=umove+0.5*(ua(i,j,k)+ua(i+1,j,k))
          v=vmove+0.5*(va(i,j,k)+va(i,j+1,k))
          w=0.5*(wa(i,j,k)+wa(i,j,k+1))
          qtot=qv(i,j,k)+ql(i,j,k)+qi(i,j,k)
          foo1(k)=foo1(k)+rho(i,j,k)*tem*(1.0+qtot)*0.5*(u*u+v*v+w*w)   &
               +ql(i,j,k)*rho(i,j,k)*tem*0.5*(vr(i,j,k)**2-2.0*w*vr(i,j,k))
          tmp=(th0(i,j,k)+tha(i,j,k))*(pi0(i,j,k)+ppi(i,j,k))
          foo2(k)=foo2(k)+rho(i,j,k)*tem*(cv+cvv*qv(i,j,k))*tmp
          foo3(k)=foo3(k)+rho(i,j,k)*tem*(1.0+qtot)*g*zh(i,j,k)
          foo4(k)=foo4(k)+rho(i,j,k)*tem*ql(i,j,k)*(cpl*tmp-lv1)   &
                         +rho(i,j,k)*tem*qi(i,j,k)*(cpi*tmp-ls1)
        enddo
        enddo
      enddo

      ek=0.0d0
      ei=0.0d0
      ep=0.0d0
      le=0.0d0
 
      do k=1,nk
        ek=ek+foo1(k)
        ei=ei+foo2(k)
        ep=ep+foo3(k)
        le=le+foo4(k)
      enddo

      ek=ek*(dx*dy*dz)
      ei=ei*(dx*dy*dz)
      ep=ep*(dx*dy*dz)
      le=le*(dx*dy*dz)

#ifdef MPI
      var=0.0d0
      call MPI_REDUCE(ek,var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,   &
                      MPI_COMM_WORLD,ierr)
      ek=var
      var=0.0d0
      call MPI_REDUCE(ei,var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,   &
                      MPI_COMM_WORLD,ierr)
      ei=var
      var=0.0d0
      call MPI_REDUCE(ep,var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,   &
                      MPI_COMM_WORLD,ierr)
      ep=var
      var=0.0d0
      call MPI_REDUCE(le,var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,   &
                      MPI_COMM_WORLD,ierr)
      le=var
      if(myid.eq.0)then
#endif

      et=ek+ei+ep+le
 
      write(6,100) 'TENERG',et
100   format(2x,a6,':',1x,e13.6)

      nstat = nstat + 1
      rstat(nstat) = ek
      nstat = nstat + 1
      rstat(nstat) = ei
      nstat = nstat + 1
      rstat(nstat) = ep
      nstat = nstat + 1
      rstat(nstat) = le
      nstat = nstat + 1
      rstat(nstat) = et

#ifdef MPI
      endif
#endif
 
      if(timestats.ge.1) time_stat=time_stat+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine calcmoe(nstat,rstat,ruh,rvh,rmh,rho,ua,va,wa,qv,ql,qi,vr)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif
 
      integer nstat
      real, dimension(stat_out) :: rstat
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, dimension(ib:ie,jb:je,kb:ke) :: qv,ql,qi,vr
 
      integer i,j,k
      real*8 :: tmu,tmv,tmw,qtot,var,tem
      real*8, dimension(nk) :: foo1,foo2,foo3

!$omp parallel do default(shared)  &
!$omp private(k)
      do k=1,nk
        foo1(k)=0.0d0
        foo2(k)=0.0d0
        foo3(k)=0.0d0
      enddo
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k,qtot,tem)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          qtot=qv(i,j,k)+ql(i,j,k)+qi(i,j,k)
          tem=ruh(i)*rvh(j)*rmh(i,j,k)
          foo1(k)=foo1(k)   &
                +rho(i,j,k)*tem*(1.0+qtot)*( umove+0.5*(ua(i,j,k)+ua(i+1,j,k)) )
          foo2(k)=foo2(k)   &
                +rho(i,j,k)*tem*(1.0+qtot)*( vmove+0.5*(va(i,j,k)+va(i,j+1,k)) )
          foo3(k)=foo3(k)                                                &
                +rho(i,j,k)*tem*(1.0+qtot)*( 0.5*(wa(i,j,k)+wa(i,j,k+1)) )   &
                -rho(i,j,k)*tem*ql(i,j,k)*vr(i,j,k)
        enddo
        enddo
      enddo

      tmu=0.0d0
      tmv=0.0d0
      tmw=0.0d0
      do k=1,nk
        tmu=tmu+foo1(k)
        tmv=tmv+foo2(k)
        tmw=tmw+foo3(k)
      enddo

#ifdef MPI
      var=0.0d0
      call MPI_REDUCE(tmu,var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                      MPI_COMM_WORLD,ierr)
      tmu=var
      var=0.0d0
      call MPI_REDUCE(tmv,var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                      MPI_COMM_WORLD,ierr)
      tmv=var
      var=0.0d0
      call MPI_REDUCE(tmw,var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                      MPI_COMM_WORLD,ierr)
      tmw=var
      if(myid.eq.0)then
#endif

      tmu=tmu*(dx*dy*dz)
      tmv=tmv*(dx*dy*dz)
      tmw=tmw*(dx*dy*dz)
 
      write(6,100) 'TMU   ',tmu
      write(6,100) 'TMV   ',tmv
      write(6,100) 'TMW   ',tmw
100   format(2x,a6,':',1x,e13.6)
 
      nstat = nstat + 1
      rstat(nstat) = tmu
      nstat = nstat + 1
      rstat(nstat) = tmv
      nstat = nstat + 1
      rstat(nstat) = tmw

#ifdef MPI
      endif
#endif
 
      if(timestats.ge.1) time_stat=time_stat+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine tmf(nstat,rstat,ruh,rvh,rho,wa)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      integer nstat
      real, dimension(stat_out) :: rstat
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa

      integer i,j,k
      real*8 :: tmfu,tmfd,mf,var
      real*8, dimension(nk) :: foo1,foo2

!$omp parallel do default(shared)  &
!$omp private(i,j,k,mf)
      do k=1,nk
        foo1(k)=0.0d0
        foo2(k)=0.0d0
        do j=1,nj
        do i=1,ni
          mf=rho(i,j,k)*0.5*(wa(i,j,k)+wa(i,j,k+1))*ruh(i)*rvh(j)
          foo1(k)=foo1(k)+max(mf,0.0d0)
          foo2(k)=foo2(k)+min(mf,0.0d0)
        enddo
        enddo
      enddo

      tmfu=0.0d0
      tmfd=0.0d0
      do k=1,nk
        tmfu=tmfu+foo1(k)
        tmfd=tmfd+foo2(k)
      enddo

#ifdef MPI
      var=0.0d0
      call MPI_REDUCE(tmfu,var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                      MPI_COMM_WORLD,ierr)
      tmfu=var
      var=0.0d0
      call MPI_REDUCE(tmfd,var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                      MPI_COMM_WORLD,ierr)
      tmfd=var
      if(myid.eq.0)then
#endif

      tmfu=tmfu*dx*dy
      tmfd=tmfd*dx*dy

      write(6,100) 'TMFU  ',tmfu
      write(6,100) 'TMFD  ',tmfd
100   format(2x,a6,':',1x,e13.6)

      nstat = nstat + 1
      rstat(nstat) = tmfu
      nstat = nstat + 1
      rstat(nstat) = tmfd

#ifdef MPI
      endif
#endif

      if(timestats.ge.1) time_stat=time_stat+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine zinterp(sigma,zs,zh,dum1,dum2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(kb:ke) :: sigma
      real, dimension(ib:ie,jb:je) :: zs
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,dum1,dum2

      integer i,j,k,kk,kup,kdn
      real, dimension(nk) :: zref

      do k=1,nk
!!!        zref(k)=(k*dz-0.5*dz)
        zref(k)=sigma(k)
      enddo

      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum2(i,j,k)=dum1(i,j,k)
      enddo
      enddo
      enddo

      do k=1,nk
      do j=1,nj
      do i=1,ni
        if( (zref(k).lt.zh(i,j,1)).or.(zref(k).gt.zh(i,j,nk)) )then
          dum1(i,j,k)=-99999999.
        elseif(zs(i,j).lt.0.1 .or. zref(k).eq.zh(i,j,1))then
          dum1(i,j,k)=dum2(i,j,k)
        else
          kup=0
          kdn=0
          do kk=1,nk
            if(zref(k).gt.zh(i,j,kk)) kdn=kk
          enddo
          kup=kdn+1
          if(kup.le.0.or.kdn.le.0.or.kup.ge.nk+1.or.kdn.ge.nk+1)then
            print *,kdn,kup
            print *,zs(i,j),zh(i,j,kdn),zref(k),zh(i,j,kup)
            print *,i,j,k
            call stopcm1
          endif
          dum1(i,j,k)=dum2(i,j,kdn)+(dum2(i,j,kup)-dum2(i,j,kdn))   &
                                   *(  zref(k  )  -zh(i,j,kdn))     &
                                   /(  zh(i,j,kup)-zh(i,j,kdn))
        endif
      enddo
      enddo
      enddo

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!
! ---------------------------------------------------------------------
! THIS FUNCTION CALCULATES THE LIQUID SATURATION VAPOR MIXING RATIO AS
! A FUNCTION OF TEMPERATURE AND PRESSURE
!
      FUNCTION RSLF(P,T)

      IMPLICIT NONE
      include 'input.incl'
      include 'constants.incl'
      REAL ESL,RSLF,X,T,P,C0,C1,C2,C3,C4,C5,C6,C7,C8
      PARAMETER (C0= .611583699E03)
      PARAMETER (C1= .444606896E02)
      PARAMETER (C2= .143177157E01)
      PARAMETER (C3= .264224321E-1)
      PARAMETER (C4= .299291081E-3)
      PARAMETER (C5= .203154182E-5)
      PARAMETER (C6= .702620698E-8)
      PARAMETER (C7= .379534310E-11)
      PARAMETER (C8=-.321582393E-13)

!  Note to self ... this should be changed, somehow, in the future.
!  GHB 060806

    if(ptype.eq.1.or.ptype.eq.2.or.ptype.eq.3.or.ptype.eq.5.or.ptype.eq.6)then

      ! from Bolton (1980, MWR)
      esl=611.2 * EXP( 17.67 * ( T  - 273.15 ) / ( T  - 29.65 ) )
      rslf= eps * ESL /(P-ESL)

    elseif(ptype.eq.4)then

      rslf=380.00*exp(17.2693882-4097.8531/(t-35.86))/p

!    elseif(ptype.eq.3)then
!
!      X=MAX(-80.,T-273.16)
!      ESL=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
!      RSLF=eps*ESL/(P-ESL)

    else

      esl=611.2 * EXP( 17.67 * ( T  - 273.15 ) / ( T  - 29.65 ) )
      rslf= eps * ESL /(P-ESL)

    endif

      RETURN
      END

!
! ---------------------------------------------------------------------
! THIS FUNCTION CALCULATES THE ICE SATURATION VAPOR MIXING RATIO AS A
! FUNCTION OF TEMPERATURE AND PRESSURE
!
      FUNCTION RSIF(P,T)

      IMPLICIT NONE
      include 'input.incl'
      include 'constants.incl'
      REAL ESI,RSIF,X,T,P,C0,C1,C2,C3,C4,C5,C6,C7,C8
      PARAMETER (C0= .609868993E03)
      PARAMETER (C1= .499320233E02)
      PARAMETER (C2= .184672631E01)
      PARAMETER (C3= .402737184E-1)
      PARAMETER (C4= .565392987E-3)
      PARAMETER (C5= .521693933E-5)
      PARAMETER (C6= .307839583E-7)
      PARAMETER (C7= .105785160E-9)
      PARAMETER (C8= .161444444E-12)

!  Note to self ... this should be changed, somehow, in the future.
!  GHB 060806

    if(ptype.eq.1.or.ptype.eq.2.or.ptype.eq.3.or.ptype.eq.5.or.ptype.eq.6)then

      ! from Tao et al (1989, MWR)
      esi=611.2 * EXP( 21.8745584 * ( T  - 273.15 ) / ( T  - 7.66 ) )
      rsif= eps * ESI /(P-ESI)

    elseif(ptype.eq.4)then

      rsif=380.00*exp(21.87455-5807.4743/(t-7.66))/p

!    elseif(ptype.eq.3)then
!
!      X=MAX(-80.,T-273.16)
!      ESI=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
!      RSIF=eps*ESI/(P-ESI)

    else

      esi=611.2 * EXP( 21.8745584 * ( T  - 273.15 ) / ( T  - 7.66 ) )
      rsif= eps * ESI /(P-ESI)

    endif

      RETURN
      END


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


! REV PROCESSED 213 LINES OF CODE. PROGRAM DONE.
      REAL FUNCTION GAMMA(X)
      implicit none
!D    DOUBLE PRECISION FUNCTION DGAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY
      REAL               &
!D    DOUBLE PRECISION   &
          C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE,   &
          TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/,   &
           SQRTPI/0.9189385332046727417803297E0/,                      &
           PI/3.1415926535897932384626434E0/
!D    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
!D   1     SQRTPI/0.9189385332046727417803297D0/,
!D   2     PI/3.1415926535897932384626434D0/
!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,   &
           XINF/3.4E38/
!D    DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
!D   1     XINF/1.79D308/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,   &
             -3.79804256470945635097577E+2,6.29331155312818442661052E+2,   &
             8.66966202790413211295064E+2,-3.14512729688483675254357E+4,   &
             -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,   &
            -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,   &
              2.25381184209801510330112E+4,4.75584627752788110767815E+3,   &
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!D    DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
!D   1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
!D   2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
!D   3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
!D    DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
!D   1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
!D   2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
!D   3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      DATA C/-1.910444077728E-03,8.4171387781295E-04,                     &
           -5.952379913043012E-04,7.93650793500350248E-04,                &
           -2.777777777777681622553E-03,8.333333333333333331554247E-02,   &
            5.7083835261E-03/
!D    DATA C/-1.910444077728D-03,8.4171387781295D-04,
!D   1     -5.952379913043012D-04,7.93650793500350248D-04,
!D   2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
!D   3      5.7083835261D-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I)
!D    CONV(I) = DBLE(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO 260 I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
  260   CONTINUE
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          RES=RES/Y1
        ELSEIF(Y1.GT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          DO 290 I=1,N
            RES=RES*Y
            Y=Y+ONE
  290     CONTINUE
        ENDIF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO 350 I=1,6
            SUM=SUM/YSQ+C(I)
  350     CONTINUE
          SUM=SUM/Y-Y+SQRTPI
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE)RES=FACT/RES
  900 GAMMA=RES
!D900 DGAMMA = RES
      RETURN
! ---------- LAST LINE OF GAMMA ----------
      END


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function mytime()
      implicit none

      include 'timestat.incl'

      integer count,rate,max
      real time_current,rcount

!----------------------------------------------------------
!  Platform-independent timer

      call system_clock(count,rate,max)
      if( count.lt.count_last )then
        ! simple kludge ... do nothing
        ! fix some other day   (GHB, 101018)
!!!        rcount = float(count+max)
!!!        time_current=rcount*clock_rate
!!!        mytime=time_current-time_last
!!!        rcount = float(count)
!!!        time_current=rcount*clock_rate
        rcount = float(count)
        time_current=rcount*clock_rate
        mytime=0.0
      else
        rcount = float(count)
        time_current=rcount*clock_rate
        mytime=time_current-time_last
      endif
      time_last=time_current
      count_last=count

!----------------------------------------------------------

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine stopcm1()
      implicit none

#ifdef MPI
      include 'mpif.h'

      call mpi_abort( mpi_comm_world, 1 )
#endif

      stop

      return
      end


