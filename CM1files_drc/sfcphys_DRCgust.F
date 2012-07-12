!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcecd(cdu,cdv,ce,u0,v0,rf0,u1,v1,w1,u,v)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(inout), dimension(ib:ie,jb:je) :: cdu,cdv,ce
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf0
      real, intent(inout), dimension(ib:ie,jb:je) :: u1,v1,w1
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0,u
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0,v

      integer i,j
      real usfc,vsfc,wsp

      real, parameter :: dcd1  =  0.5e-3
      real, parameter :: dcd2  =  2.4e-3
      real, parameter :: dwsp1 =  5.0
      real, parameter :: dwsp2 = 35.0

!-----------------------------------------------------------------------
! for pertflx=1, the base-state is not included in calculation of wind speed

      IF(pertflx.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          u1(i,j) = u(i,j,1)-u0(i,j,1)
          v1(i,j) = v(i,j,1)-v0(i,j,1)
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          u1(i,j) = u(i,j,1)
          v1(i,j) = v(i,j,1)
        enddo
        enddo

      ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,usfc,vsfc)
      do j=1,nj
      do i=1,ni
        usfc=0.5*( u1(i,j)+u1(i+1,j) )
        vsfc=0.5*( v1(i,j)+v1(i,j+1) )
        w1(i,j)= 3.0 + sqrt(usfc**2+vsfc**2)   !DRC 05-22-11 gustiness
      enddo
      enddo

!-----------------------------------------------------------------------
!  Get Cd on u pts:

  IF(idrag.eq.1)THEN

    IF(cecd.eq.1)THEN
      ! constant value:

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1
        cdu(i,j) = cnstcd
      enddo
      enddo

    ELSEIF(cecd.eq.2)THEN
      ! Deacon's formula:  see Rotunno and Emanuel (1987, JAS, p. 547)

!$omp parallel do default(shared)   &
!$omp private(i,j,usfc,vsfc,wsp)
      do j=1,nj
      do i=1,ni+1
        usfc=u1(i,j)
        vsfc=0.25*( ( v1(i  ,j)+v1(i  ,j+1) )   &
                   +( v1(i-1,j)+v1(i-1,j+1) ) )
        wsp=sqrt(usfc**2+vsfc**2)
        cdu(i,j) = 1.1e-3+(4.0e-5*wsp)
      enddo
      enddo

    ELSEIF(cecd.eq.3)THEN
      ! Based on Donelan (2004, GRL)

!$omp parallel do default(shared)   &
!$omp private(i,j,usfc,vsfc,wsp)
      do j=1,nj
      do i=1,ni+1
        usfc=u1(i,j)
        vsfc=0.25*( ( v1(i  ,j)+v1(i  ,j+1) )   &
                   +( v1(i-1,j)+v1(i-1,j+1) ) )
        wsp=sqrt(usfc**2+vsfc**2)
        cdu(i,j) = dcd1+( wsp - dwsp1 )*( dcd2 - dcd1 )/( dwsp2 - dwsp1 )
        cdu(i,j) = min(cdu(i,j),dcd2)
        cdu(i,j) = max(cdu(i,j),dcd1)
      enddo
      enddo

    ELSE

      write(outfile,*)
      write(outfile,*) '  Undefined value for cecd '
      write(outfile,*)
      call stopcm1

    ENDIF

!-----------------------------------------------------------------------
!  Get Cd on v pts:

    IF(cecd.eq.1)THEN
      ! constant value:

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        cdv(i,j) = cnstcd
      enddo
      enddo

    ELSEIF(cecd.eq.2)THEN
      ! Deacon's formula:  see Rotunno and Emanuel (1987, JAS, p. 547)

!$omp parallel do default(shared)   &
!$omp private(i,j,usfc,vsfc,wsp)
      do j=1,nj+1
      do i=1,ni
        usfc=0.25*( ( u1(i,j  )+u1(i+1,j  ) )   &
                   +( u1(i,j-1)+u1(i+1,j-1) ) )
        vsfc=v1(i,j)
        wsp=sqrt(usfc**2+vsfc**2)
        cdv(i,j) = 1.1e-3+(4.0e-5*wsp)
      enddo
      enddo

    ELSEIF(cecd.eq.3)THEN
      ! Based on Donelan (2004, GRL)

!$omp parallel do default(shared)   &
!$omp private(i,j,usfc,vsfc,wsp)
      do j=1,nj+1
      do i=1,ni
        usfc=0.25*( ( u1(i,j  )+u1(i+1,j  ) )   &
                   +( u1(i,j-1)+u1(i+1,j-1) ) )
        vsfc=v1(i,j)
        wsp=sqrt(usfc**2+vsfc**2)
        cdv(i,j) = dcd1+( wsp - dwsp1 )*( dcd2 - dcd1 )/( dwsp2 - dwsp1 )
        cdv(i,j) = min(cdv(i,j),dcd2)
        cdv(i,j) = max(cdv(i,j),dcd1)
      enddo
      enddo

    ELSE

      write(outfile,*)
      write(outfile,*) '  Undefined value for cecd '
      write(outfile,*)
      call stopcm1

    ENDIF

  ENDIF   ! endif for idrag=1 check

!-----------------------------------------------------------------------
!  Get Ce:

  IF(isfcflx.eq.1)THEN

    IF(cecd.eq.1)THEN
      ! constant value:

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        ce(i,j) = cnstce
      enddo
      enddo

    ELSEIF(cecd.eq.2)THEN
      ! Deacon's formula:  see Rotunno and Emanuel (1987, JAS, p. 547)

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        ce(i,j) = 1.1e-3+(4.0e-5*w1(i,j))
      enddo
      enddo

    ELSEIF(cecd.eq.3)THEN
      ! Constant, based on Drennan et al. (2007, JAS, p. 1103)

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        ce(i,j) = 1.18e-3
      enddo
      enddo

    ELSE

      write(outfile,*)
      write(outfile,*) '  Undefined value for cecd '
      write(outfile,*)
      call stopcm1

    ENDIF

  ENDIF

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine sfcdrag(cdu,cdv,u0,v0,rf0,u1,v1,t13,t23,u,v)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in), dimension(ib:ie,jb:je) :: cdu,cdv
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf0
      real, intent(in), dimension(ib:ie,jb:je) :: u1,v1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t13,t23
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0,u
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0,v

      integer i,j
      real usfc,vsfc,wspd

!-----------------------------------------------------------------------
!  x component:

!$omp parallel do default(shared)   &
!$omp private(i,j,usfc,vsfc,wspd)
      do j=1,nj
      do i=1,ni+1
        usfc=u1(i,j)
        vsfc=0.25*( ( v1(i  ,j)+v1(i  ,j+1) )   &
                   +( v1(i-1,j)+v1(i-1,j+1) ) )
        wspd=sqrt(usfc**2+vsfc**2)
        t13(i,j,1) = cdu(i,j)*usfc*wspd
      enddo
      enddo

!-----------------------------------------------------------------------
!  y component:

!$omp parallel do default(shared)   &
!$omp private(i,j,usfc,vsfc,wspd)
      do j=1,nj+1
      do i=1,ni
        usfc=0.25*( ( u1(i,j  )+u1(i+1,j  ) )   &
                   +( u1(i,j-1)+u1(i+1,j-1) ) )
        vsfc=v1(i,j)
        wspd=sqrt(usfc**2+vsfc**2)
        t23(i,j,1) = cdv(i,j)*vsfc*wspd
      enddo
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine sfcflux(dt,ruh,xf,rvh,pi0s,ce,zh,pi0,thv0,th0,u0,v0,tsk,thflux,qvflux,mavail,   &
                         rho,u1,v1,w1,u,v,ppi,tha,qva,qsfc)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: rvh
      real, intent(in), dimension(ib:ie,jb:je) :: pi0s,ce
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,pi0,thv0,th0
      real, intent(in), dimension(ib:ie,jb:je) :: tsk
      real, intent(inout), dimension(ib:ie,jb:je) :: thflux,qvflux
      real, intent(in), dimension(ibl:iel,jbl:jel) :: mavail
      real, intent(in), dimension(ib:ie,jb:je) :: u1,v1,w1
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0,u
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0,v
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,ppi,tha
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem) :: qva
      real*8, intent(inout) :: qsfc

      integer i,j
      real pisfc,usfc,vsfc,qvsat

      real :: rhosfc,tem
      real rslf
      real*8, dimension(nj) :: bud1

!-----------------------------------------------------------------------
!  sensible heat flux:

!$omp parallel do default(shared)   &
!$omp private(i,j,pisfc)
      do j=1,nj
      do i=1,ni
        pisfc = pi0s(i,j) + ppi(i,j,1)
        thflux(i,j)=ce(i,j)*w1(i,j)*(tsk(i,j)/pisfc-th0(i,j,1)-tha(i,j,1))
      enddo
      enddo

!-----------------------------------------------------------------------
!  latent heat flux:

    IF(imoist.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,pisfc,qvsat)
      do j=1,nj
      do i=1,ni
        pisfc = pi0s(i,j) + ppi(i,j,1)
        qvsat=rslf(p00*(pisfc**cpdrd),tsk(i,j))
        qvflux(i,j)=ce(i,j)*w1(i,j)*(qvsat-qva(i,j,1))*mavail(i,j)
      enddo
      enddo

!-----------------------------------------------------------------------
!  some budget calculations (only calculated if imoist=1):

!$omp parallel do default(shared)  &
!$omp private(j)
      do j=1,nj
        bud1(j)=0.0d0
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,rhosfc)
      do j=1,nj
      do i=1,ni
        rhosfc=rho(i,j,1)-zh(i,j,1)*(rho(i,j,2)-rho(i,j,1))   &
                                   /( zh(i,j,2)- zh(i,j,1))
        if(axisymm.eq.1) rhosfc=rhosfc*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
        bud1(j)=bud1(j)+qvflux(i,j)*ruh(i)*rvh(j)*rhosfc
      enddo
      enddo

      tem = dt*dx*dy

      do j=1,nj
        qsfc=qsfc+bud1(j)*tem
      enddo

    ENDIF

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine sfcdiags(tsk,thflux,qvflux,cdu,cdv,ce,u1,v1,w1,   &
                          xland,psfc,qsfc,u10,v10,hfx,qfx,cda,znt,ust,gz1oz0,   &
                          psim,psih,br,zol,mol,hpbl,wspd,dsxy,th2,t2,q2, &
                          zs,zh,pi0s,pi0,th0,ppi,tha,rho,qa,ua,va)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in), dimension(ib:ie,jb:je) :: tsk,thflux,qvflux,   &
                                                  cdu,cdv,ce,u1,v1,w1
      real, intent(in), dimension(ibl:iel,jbl:jel) :: xland
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: psfc,qsfc,u10,v10,hfx,qfx, &
                                    cda,znt,ust,gz1oz0,psim,psih,br,zol,mol,hpbl,wspd,dsxy,th2,t2,q2
      real, intent(in), dimension(itb:ite,jtb:jte) :: zs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je) :: pi0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,th0,ppi,tha,rho
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va

      integer :: i,j
      real :: pisfc,thgb,thx,thvx,tskv,za,govrth,dthvdz,vconv,vsgd,dthvm,   &
              rhosfc,val,fluxc
      real :: rslf

      REAL    , PARAMETER ::  VCONVC=1.
      REAL    , PARAMETER ::  CZO=0.0185
      REAL    , PARAMETER ::  OZO=1.59E-5
      REAL    , PARAMETER ::  KARMAN = 0.4
      REAL    , PARAMETER ::  EP1 = rv/rd - 1.0
      REAL,     PARAMETER ::  PRT=1.

      ! surface layer diagnostics:

!$omp parallel do default(shared)   &
!$omp private(i,j,pisfc,thgb,thx,thvx,tskv,za,govrth,dthvdz,vconv,vsgd,   &
!$omp dthvm,rhosfc,val,fluxc)
      do j=1,nj
      do i=1,ni
        pisfc = pi0s(i,j) + ppi(i,j,1)
        psfc(i,j) = p00*(pisfc**cpdrd)
        thgb = tsk(i,j)/pisfc
        thx = th0(i,j,1)+tha(i,j,1)
        thvx = thx*(1.+EP1*qa(i,j,1,nqv))
        qsfc(i,j) = rslf(psfc(i,j),tsk(i,j))
        tskv = thgb*(1.0+ep1*qsfc(i,j))
        za = zh(i,j,1)
        if( terrain_flag ) za = za - zs(i,j)
        govrth = g/thx
        ! use old znt (as in sfclay):
        gz1oz0(i,j) = alog(za/znt(i,j))
        DTHVDZ = THVX-TSKV
        if (xland(i,j).lt.1.5) then
          ! land:
          fluxc = max(thflux(i,j) + ep1*tskv*qvflux(i,j),0.)
          VCONV = vconvc*(g/tsk(i,j)*hpbl(i,j)*fluxc)**.33
        else
          ! ocean:
          IF(-DTHVDZ.GE.0)THEN
            DTHVM=-DTHVDZ
          ELSE
            DTHVM=0.
          ENDIF
          VCONV = 2.*SQRT(DTHVM)
        endif
! Mahrt and Sun low-res correction
        VSGD = 0.32 * (max(dsxy(i,j)/5000.-1.,0.))**.33
        wspd(i,j) = sqrt( w1(i,j)*w1(i,j) + vconv*vconv + vsgd*vsgd )
        wspd(i,j) = max(0.1,wspd(i,j))
        br(i,j) = govrth*za*DTHVDZ/(wspd(i,j)**2)
        rhosfc=rho(i,j,1)-zh(i,j,1)*(rho(i,j,2)-rho(i,j,1))   &
                                   /( zh(i,j,2)- zh(i,j,1))
        hfx(i,j) = thflux(i,j)*cp*rhosfc
        qfx(i,j) = qvflux(i,j)*rhosfc
        cda(i,j) = 0.25*( (cdu(i,j)+cdu(i+1,j))+(cdv(i,j)+cdv(i,j+1)) )
        ! new ust/znt:
        ! TO PREVENT OSCILLATIONS AVERAGE WITH OLD VALUE
        ust(i,j) = 0.5*( ust(i,j) + sqrt(cda(i,j)*wspd(i,j)*wspd(i,j)) )
        if( (xland(i,j)-1.5).ge.0.0 )then
          znt(i,j) = CZO*ust(i,j)*ust(i,j)/G+OZO
        endif
        ! impose neutral sfc layer:
        psim(i,j) = 0.0
        psih(i,j) = 0.0
        zol(i,j) = 0.0
        mol(i,j) = 0.0
        val = alog(10.0/znt(i,j))/alog(za/znt(i,j))
        u10(i,j) = 0.5*(ua(i,j,1)+ua(i+1,j,1))*val
        v10(i,j) = 0.5*(va(i,j,1)+va(i,j+1,1))*val
        val = alog(2.0/znt(i,j))/alog(za/znt(i,j))
        th2(i,j) = thgb+(thx-thgb)*val
        q2(i,j) = qsfc(i,j)+(qa(i,j,1,nqv)-qsfc(i,j))*val
        t2(i,j) = th2(i,j)*pisfc
      enddo
      enddo

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()
      end subroutine sfcdiags


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine gethpbl(zh,th0,tha,qa,hpbl)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,th0,tha
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: hpbl

      integer :: i,j,kk
      real :: thx,thvx,thv,thvlast,thcrit

      REAL    , PARAMETER ::  EP1 = rv/rd - 1.0

      ! (NEEDED BY SFCLAY ... THIS IS A ROUGH ESTIMATE ONLY)
      ! (ONLY NEEDED WHEN IPBL=0)
      ! (USE WITH CAUTION)
      ! extraordinarily simple calculation:  define pbl depth as 
      ! level where thv is first greater than thv at lowest model level
      ! 110104:  add 0.5 K, for the sake of slightly stable PBLs

!$omp parallel do default(shared)   &
!$omp private(i,j,kk,thx,thvx,thv,thvlast,thcrit)
      do j=1,nj
      do i=1,ni
        hpbl(i,j) = 0.0
        kk = 1
        thx = th0(i,j,1)+tha(i,j,1)
        thvx = thx*(1.+EP1*qa(i,j,1,nqv))
        thvlast = thvx
        thcrit = thvx+0.5
        do while( hpbl(i,j).lt.1.0e-12 .and. kk.lt.nk )
          kk = kk + 1
          thv = (th0(i,j,kk)+tha(i,j,kk))*(1.0+EP1*qa(i,j,kk,nqv))
          if( thv.ge.thcrit )then
            hpbl(i,j) = zh(i,j,kk-1)+(zh(i,j,kk)-zh(i,j,kk-1))   &
                                    *(thcrit-thvlast)/(thv-thvlast)
          endif
          thvlast = thv
        enddo
        if( kk.gt.(nk-1) .or. hpbl(i,j).lt.1.0e-12 ) hpbl(i,j) = 0.0
      enddo
      enddo

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()
      end subroutine gethpbl


