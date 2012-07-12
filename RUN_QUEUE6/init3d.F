

      subroutine init3d(num_soil_layers,qbudget,asq,bsq,                  &
                        xh,rxh,uh,ruh,xf,rxf,uf,ruf,yh,vh,rvh,yf,vf,rvf,  &
                        xfref,yfref,                                      &
                        zh,mh,rmh,zf,mf,rmf,rho0s,pi0s,prs0s,pi0,prs0,rho0,thv0,th0,t0,qv0,   &
                        u0,v0,rh0,qc0,ql0,rr0,rf0,rrf0,                   &
                        zs,gz,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,cdu,cdv,ce,    &
                        radbcw,radbce,radbcs,radbcn,                      &
                        dum1,dum2,dum3,dum4,divx,rho,prs,                 &
                        t11,t12,t13,t22,t23,t33,                          &
                        rru,ua,u3d,uten,uten1,rrv,va,v3d,vten,vten1,      &
                        rrw,wa,w3d,wten,wten1,ppi,pp3d,ppten,sten,        &
                        tha,th3d,thten,thten1,thterm,tk,qa,q3d,qten,      &
                        kmh,kmv,khh,khv,tkea,tke3d,tketen,                &
                        pta,pt3d,ptten,                                   &
                        pdata,cfb,cfa,cfc,ad1,ad2,pdt,deft,rhs,trans)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'

      integer :: num_soil_layers
      real*8, dimension(nbudget) :: qbudget
      real*8, dimension(numq) :: asq,bsq
      real, dimension(ib:ie) :: xh,rxh,uh,ruh
      real, dimension(ib:ie+1) :: xf,rxf,uf,ruf
      real, dimension(jb:je) :: yh,vh,rvh
      real, dimension(jb:je+1) :: yf,vf,rvf
      real, dimension(-2:nx+4) :: xfref
      real, dimension(-2:ny+4) :: yfref
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf,rmf
      real, dimension(ib:ie,jb:je) :: rho0s,pi0s,prs0s
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,prs0,rho0,thv0,th0,t0,qv0,rh0
      real, dimension(ib:ie,jb:je,kb:ke) :: qc0,ql0,rr0,rf0,rrf0
      real, dimension(itb:ite,jtb:jte) :: zs,gz
      real, dimension(ib:ie,jb:je,nrain) :: rain,sws,svs,sps,srs,sgs,sus,shs
      real, dimension(ib:ie,jb:je) :: thflux,qvflux,cdu,cdv,ce
      real, dimension(jb:je,kb:ke) :: radbcw,radbce
      real, dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4
      real, dimension(ib:ie,jb:je,kb:ke) :: divx,rho,prs
      real, dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u0,rru,ua,u3d,uten,uten1
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v0,rrv,va,v3d,vten,vten1
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw,wa,w3d,wten,wten1
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppten,sten
      real, dimension(ib:ie,jb:je,kb:ke) :: tha,th3d,thten,thten1,thterm,tk
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa,q3d,qten
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d,tketen
      real, dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta,pt3d,ptten
      real, dimension(npvals,nparcels) :: pdata
      real, dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: cfb
      real, dimension(kpb:kpe) :: cfa,cfc,ad1,ad2
      complex, dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: pdt,deft
      complex, dimension(ipb:ipe,jpb:jpe) :: rhs,trans
 
!-----------------------------------------------------------------------

      integer i,j,k,l,n,nn,nbub,nloop
      integer ic,jc,ifoo,jfoo
      integer :: shape_qv
      real ric,rjc
      real xc,yc,zc,bhrad,bvrad,bptpert,beta,omega,tmp,zdep
      real thvnew(nk),pinew(nk)
      real thl,ql,qt,th1,t1,ql2,rm,cpm,v1,v2,th2

      real, dimension(:), allocatable :: rref
      real, dimension(:,:), allocatable :: vref,piref,thref,thvref,qvref,qvsat,qv_bub
      real*8 :: rmax,vmax,frac,angle
      real :: r0,zdd,dd2,dd1,vr,rr,diff,xref,yref,xmax,ymax,zmin,zmax
      real :: r0drmax,dz0,dz0_qv,z_bltop
      real :: rc_qv,r0_qv,rhpert_max,zc_qv,zmax_qv,zmin_qv,Lh_qv,Lz_qv,relh
      real :: c1,c2,mult,nominal_dx
      integer :: ival,ni1,ni2,ni3
      integer :: i1,i2,ii,nref

      real rmin,foo1,foo2,umax,umin,vmin
      real :: rand,amplitude
      real*8 :: dpi

      logical :: setppi,maintain_rh

!--------------------------

      real rslf

      write(outfile,*) 'Inside INIT3D'
      write(outfile,*)

      convinit = 0
      setppi = .true.
      maintain_rh = .false.

!------------------------------------------------------------------
!  Initialize arrays with base state values.

      do n=1,nbudget
        qbudget(n)=0.0d0
      enddo

      do n=1,numq
        asq(n)=0.0d0
        bsq(n)=0.0d0
      enddo

      do n=1,nrain
      do j=jb,je
      do i=ib,ie
        ! these are all positive-definite, so set initial value to zero:
        rain(i,j,n)=0.0
        sws(i,j,n)=0.0
        srs(i,j,n)=0.0
        sgs(i,j,n)=0.0
        shs(i,j,n)=0.0
        ! for sps, we want to get a MINIMUM value at the surface, so...
        ! set sps to an absurdly large number:
        sps(i,j,n)=200000.0
        ! svs and sus can be negative or positive, 
        ! but we want to get a MAXIMUM value, so...
        ! set svs and sus to an absurdly low (negative) number:
        svs(i,j,n)=-1000.0
        sus(i,j,n)=-1000.0
      enddo
      enddo
      enddo

      do j=jb,je
      do i=ib,ie
        thflux(i,j)=0.0
        qvflux(i,j)=0.0
           cdu(i,j)=0.0
           cdv(i,j)=0.0
            ce(i,j)=0.0
      enddo
      enddo

      do k=kb,ke
      do j=jb,je
      do i=ib,ie+1
        ua(i,j,k)=u0(i,j,k)
        u3d(i,j,k)=u0(i,j,k)
        rru(i,j,k)=0.0
        uten(i,j,k)=0.0
        uten1(i,j,k)=0.0
      enddo
      enddo
      enddo

      do k=kb,ke
      do j=jb,je+1
      do i=ib,ie
        va(i,j,k)=v0(i,j,k)
        v3d(i,j,k)=v0(i,j,k)
        rrv(i,j,k)=0.0
        vten(i,j,k)=0.0
        vten1(i,j,k)=0.0
      enddo
      enddo
      enddo

      do k=kb,ke+1
      do j=jb,je
      do i=ib,ie
        rrw(i,j,k)=0.0
        wa(i,j,k)=0.0
        w3d(i,j,k)=0.0
        wten(i,j,k)=0.0
        wten1(i,j,k)=0.0
      enddo
      enddo
      enddo

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        dum1(i,j,k)=0.0
        dum2(i,j,k)=0.0
        dum3(i,j,k)=0.0
        dum4(i,j,k)=0.0
        divx(i,j,k)=0.0
        rho(i,j,k)=rho0(i,j,k)
        prs(i,j,k)=prs0(i,j,k)
        ppi(i,j,k)=0.0
        pp3d(i,j,k)=0.0
        ppten(i,j,k)=0.0
        sten(i,j,k)=0.0
        tha(i,j,k)=0.0
        th3d(i,j,k)=0.0
        thten(i,j,k)=0.0
        thten1(i,j,k)=0.0
        thterm(i,j,k)=0.0
        tk(i,j,k)=0.0
        t11(i,j,k)=0.0
        t12(i,j,k)=0.0
        t13(i,j,k)=0.0
        t22(i,j,k)=0.0
        t23(i,j,k)=0.0
        t33(i,j,k)=0.0
      enddo
      enddo
      enddo

      IF(thsmall.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=jb,je
        do i=ib,ie
          tk(i,j,k)=-0.5*rdz*mf(i,j,k)*rf0(i,j,k)*(th0(i,j,k)-th0(i,j,k-1))
        enddo
        enddo
        enddo
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=jb,je
        do i=ib,ie
          tk(i,j,1)=-0.5*rf0(i,j,1)*(th0(i,j,2)-th0(i,j,1))   &
                                   /( zh(i,j,2)- zh(i,j,1))
          tk(i,j,nk+1)=-0.5*rf0(i,j,nk+1)*(th0(i,j,nk)-th0(i,j,nk-1))   &
                                         /( zh(i,j,nk)- zh(i,j,nk-1))
        enddo
        enddo
      ENDIF

!-----

      ! Start by defining all q arrays to be zero

      do n=1,numq
      do k=kbm,kem
      do j=jbm,jem
      do i=ibm,iem
          qa(i,j,k,n)=0.0
         q3d(i,j,k,n)=0.0
        qten(i,j,k,n)=0.0
      enddo
      enddo
      enddo
      enddo

      ! Now, define the qv arrays

    IF(imoist.eq.1)THEN

      do k=kbm,kem
      do j=jbm,jem
      do i=ibm,iem
        qa(i,j,k,nqv)=qv0(i,j,k)
      enddo
      enddo
      enddo

!---- This is here to ensure that certain idealized cases work ----

      IF( (isnd.eq.4 .or. isnd.eq.9 .or. isnd.eq.10 .or. isnd.eq.11) )THEN

        do k=kbm,kem
        do j=jbm,jem
        do i=ibm,iem
          qa(i,j,k,2)=qc0(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

      IF( (isnd.eq.4 .or. isnd.eq.9 .or. isnd.eq.10) .and. iice.eq.1 )THEN

        do k=kbm,kem
        do j=jbm,jem
        do i=ibm,iem
          qa(i,j,k,4)=ql0(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

    ENDIF

      ! this array must be zero for all cases

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        ql0(i,j,k)=0.0
      enddo
      enddo
      enddo

!-----

      do k=kbc,kec
      do j=jbc,jec
      do i=ibc,iec
        kmh(i,j,k)=0.0
        kmv(i,j,k)=0.0
        khh(i,j,k)=0.0
        khv(i,j,k)=0.0
      enddo
      enddo
      enddo

      do k=kbt,ket
      do j=jbt,jet
      do i=ibt,iet
        tkea(i,j,k)=0.0
        tke3d(i,j,k)=0.0
        tketen(i,j,k)=0.0
      enddo
      enddo
      enddo

!-----

      do n=1,npt
      do k=kbp,kep
      do j=jbp,jep
      do i=ibp,iep
        pta(i,j,k,n)=0.0
        pt3d(i,j,k,n)=0.0
        ptten(i,j,k,n)=0.0
      enddo
      enddo
      enddo
      enddo

    IF(iptra.eq.1)THEN
      ! define concentrations for passive fluid tracers here:
      do n=1,npt
      do k=kbp,kep
      do j=jbp,jep
      do i=ibp,iep
        if(n.eq.1)then
          pta(i,j,k,n)=0.0
          if(zh(i,j,k).lt.3000.0) pta(i,j,k,n)=0.001
        endif
        if(n.eq.2)then
          pta(i,j,k,n)=0.0
          if(zh(i,j,k).gt.3000.0.and.zh(i,j,k).lt.6000.0) pta(i,j,k,n)=0.001
        endif
        if(n.eq.3)then
          pta(i,j,k,n)=0.0
          if(zh(i,j,k).gt.6000.0.and.zh(i,j,k).lt.9000.0) pta(i,j,k,n)=0.001
        endif
      enddo
      enddo
      enddo
      enddo
    ENDIF

!-----

      do n=1,nparcels
      do i=1,npvals
        pdata(i,n)=0.0
      enddo
      enddo

      IF(iprcl.eq.1)THEN
        ! define initial locations of parcels here:

        write(outfile,*)
        write(outfile,*) '  Parcels ! '
        write(outfile,*) '  npvals,nparcels = ',npvals,nparcels
        write(outfile,*) '  Initial parcel locations (x,y,z):'
        n = 0
        do k=1,10
        do j=1,60
        do i=1,60
          n = n + 1
          if(n.gt.nparcels)then
            write(outfile,*)
            write(outfile,*) ' You are trying to define too many parcels'
            write(outfile,*)
            write(outfile,*) ' Increase the value of nparcels in namelist.input'
            write(outfile,*)
            call stopcm1
          endif
          pdata(1,n) = minx + 2000.0*(i-1)
          pdata(2,n) = miny + 2000.0*(j-1)
          pdata(3,n) = zh(1,1,1) + 1000.0*(k-1)
          write(outfile,*) n,pdata(1,n),pdata(2,n),pdata(3,n)
        enddo
        enddo
        enddo
        write(outfile,*)

        do i=1,nparcels
          pdata(ifx,i)=pdata(1,i)
          pdata(ify,i)=pdata(2,i)
          pdata(ifz,i)=pdata(3,i)
        enddo

      ENDIF

!-----

      do k=kb,ke
      do j=jb,je
        radbcw(j,k)=0.0
        radbce(j,k)=0.0
      enddo
      enddo

      do k=kb,ke
      do i=ib,ie
        radbcs(i,k)=0.0
        radbcn(i,k)=0.0
      enddo
      enddo

!-----------------------------------------------------------------------
!  iinit = 1
!  Warm bubble
!  reference:

      IF(iinit.eq.1)THEN

        write(outfile,*)
        write(outfile,*) '  Warm bubble'
        write(outfile,*)

        ric     =      0.0  ! center of bubble in x-direction (m)
        rjc     =      0.0  ! center of bubble in y-direction (m)
        zc = 4375.0 ! height of center of bubble above ground (m)
        bhrad   =  10000.0  ! horizontal radius of bubble (m)
        bvrad   =   1400.0  ! vertical radius of bubble (m)
        bptpert =      1.0  ! max potential temp perturbation (K)

        ! By default, CM1 sets qv=constant at a constant height level for 
        ! this value of iinit.  If you would rather have rh=constant at 
        ! a constant height level, then set this to .true.
        maintain_rh = .false.

        do k=1,nk
        do j=1,nj
        do i=1,ni
          beta=sqrt(                             &
                    ((xh(i)-ric)/bhrad)**2       &
                   +((yh(j)-rjc)/bhrad)**2       &
                   +((zh(i,j,k)-zc)/bvrad)**2)
          if(beta.lt.1.0)then
            tha(i,j,k)=bptpert*(cos(0.5*pi*beta)**2)
          else
            tha(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iinit = 2
!  Cold pool (dam break style)
!  reference:  

      ELSEIF(iinit.eq.2)THEN

        write(outfile,*)
        write(outfile,*) '  Cold pool .... periodic in N-S'
        write(outfile,*)

        ric      =  200000.0   ! eastern edge of cold pool
        zdep     =    2500.0   ! depth of cold pool (m)
        bptpert  =      -6.0   ! max temp perturbation at sfc (K)

        ! By default, CM1 sets qv=constant at a constant height level for 
        ! this value of iinit.  If you would rather have rh=constant at 
        ! a constant height level, then set this to .true.
        maintain_rh = .true.

        do k=1,nk
        do j=1,nj
        do i=1,ni
          if( (xh(i).le.ric).and.(zh(i,j,k).lt.zdep) )then
            tha(i,j,k)=bptpert*(zdep-zh(i,j,k))/zdep
          else
            tha(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo


!-----------------------------------------------------------------------
!  iinit = 3
!  Line of warm bubbles
!  reference:  

      ELSEIF(iinit.eq.3)THEN

        write(outfile,*)
        write(outfile,*) '  Line of warm bubbles'
        write(outfile,*)
 
        nbub    =      3     ! number of warm bubbles
        ric     =  30000.0   ! center of bubble in x-direction (m)
        zc = 4375.0 ! height of center of bubble above ground (m)
        bhrad   =  10000.0   ! horizontal radius of bubble (m)
        bvrad   =   1400.0   ! vertical radius of bubble (m)
        bptpert =      2.0   ! max potential temp perturbation (K)

        ! By default, CM1 sets qv=constant at a constant height level for 
        ! this value of iinit.  If you would rather have rh=constant at 
        ! a constant height level, then set this to .true.
        maintain_rh = .false.

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          tha(i,j,k)=0.0
        enddo
        enddo
        enddo

        do n=1,nbub

          if(n.eq.1) rjc=  3000.0
          if(n.eq.2) rjc= 33000.0
          if(n.eq.3) rjc= 63000.0

          write(outfile,*) '  ric,rjc=',n,ric,rjc
 
          do k=kb,ke
          do j=jb,je
          do i=ib,ie
            beta=sqrt(                        &
                    ((xh(i)-ric)/bhrad)**2    &
                   +((yh(j)-rjc)/bhrad)**2    &
                   +((zh(i,j,k)-zc)/bvrad)**2)
            if(beta.lt.1.0)then
              tha(i,j,k)=bptpert*(cos(0.5*pi*beta)**2)
            else
              tha(i,j,k)=max(0.0,tha(i,j,k))
            endif
          enddo
          enddo
          enddo

        enddo


!-----------------------------------------------------------------------
!  iinit = 4
!  moist bubble for moist benchmark
!  reference:  Bryan and Fritsch, 2002, MWR, 130, 2917-2928.

      ELSEIF(iinit.eq.4)THEN

        ! parameters for dry counterpart bubble

        ric      =      0.0       ! x-location of bubble center (m)
        zc = 4375.0 ! z-location of bubble center (m)
        bhrad    =   2000.0       ! horizontal radius of bubble (m)
        bvrad    =   2000.0       ! vertical radius of bubble (m)
        bptpert  =      2.0       ! maximum potential temp. pert. (K)

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          beta=sqrt( ((xh(i)-ric)/bhrad)**2    &
                    +((zh(i,j,k)-zc)/bvrad)**2)
          if(beta.lt.1.0)then
            dum1(i,j,k)=bptpert*(cos(0.5*pi*beta)**2)
          else
            dum1(i,j,k)=0.
          endif
        enddo
        enddo
        enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          tha(i,j,k)=0.
          ppi(i,j,k)=0.
        enddo
        enddo
        enddo

        do nn=1,30
          do k=kb,ke
          do j=jb,je
          do i=ib,ie
            qa(i,j,k,nqv)=rh0(i,j,k)*rslf(prs0(i,j,k),(th0(i,j,k)+tha(i,j,k))*pi0(i,j,k))
          enddo
          enddo
          enddo

          do k=kb,ke
          do j=jb,je
          do i=ib,ie
            qa(i,j,k,2)=max(qt_mb-qa(i,j,k,nqv),0.0)
          enddo
          enddo
          enddo

          do k=kb,ke
          do j=jb,je
          do i=ib,ie
            tha(i,j,k)=( (dum1(i,j,k)/300.)+(1.0+qt_mb)/(1.0+qa(i,j,k,nqv)) )  &
               *thv0(i,j,k)*(1.0+qa(i,j,k,nqv))/(1.0+reps*qa(i,j,k,nqv)) - th0(i,j,k)
            if(abs(tha(i,j,k)).lt.1.e-4) tha(i,j,k)=0.
          enddo
          enddo
          enddo
        enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          qa(i,j,k,nqv)=rslf(prs0(i,j,k),(th0(i,j,k)+tha(i,j,k))*pi0(i,j,k))
          qa(i,j,k,2  )=max(qt_mb-qa(i,j,k,nqv),0.0)
        enddo
        enddo
        enddo

!-----------------------------------------------------------------
!  iinit = 5
!  density current sim

      ELSEIF(iinit.eq.5)THEN

        write(outfile,*)
        write(outfile,*) '  Cold pool (elipse, following Straka)'
        write(outfile,*)

        ric     =     0.0
        rjc     =     0.0
        zc      =  3000.0
        bhrad   =  4000.0
        bvrad   =  2000.0
        bptpert =   -15.0

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          beta=sqrt(                           &
                     ((xh(i)-ric)/bhrad)**2    &
!!!                    +((yh(j)-rjc)/bhrad)**2    &
                    +((zh(i,j,k)-zc)/bvrad)**2)
          if(beta.lt.1.0)then
            dum1(i,j,k)=bptpert*(cos(pi*beta)+1.0)*0.5
          else
            dum1(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          tmp=t0(i,j,k)+dum1(i,j,k)
          tha(i,j,k)=tmp/pi0(i,j,k)-th0(i,j,k)
          if(abs(tha(i,j,k)).lt.1.e-4) tha(i,j,k)=0.0
          ppi(i,j,k)=0.0
        enddo
        enddo
        enddo

!------------------------------------------------------------------
!  Rotunno-Emanuel tropical cyclone vortex
!  (see Rotunno and Emanuel, 1987, JAS, for more information)

      ELSEIF(iinit.eq.7)THEN

        r0     =   412500.0
        rmax   =    82500.0
        vmax   =       15.0
        zdd    =    20000.0

        dd2 = 2.0 * rmax / ( r0 + rmax )

        allocate(  rref(nx)       )
        allocate(  vref(nx,0:nk+1))
        allocate( piref(nx,0:nk+1))
        allocate( thref(nx,0:nk+1))
        allocate(thvref(nx,0:nk+1))
        allocate( qvref(nx,0:nk+1))

          rref=0.0
          vref=0.0
         piref=0.0
         thref=0.0
        thvref=0.0
         qvref=0.0

        IF(ibalance.ne.0)THEN
          write(outfile,*)
          write(outfile,*) ' Please use ibalance = 0 with iinit=7'
          write(outfile,*)
          write(outfile,*) ' ... stopping inside init3d ... '
          write(outfile,*)
          call stopcm1
        ENDIF
        IF(terrain_flag)THEN
          write(outfile,*)
          write(outfile,*) ' iinit=7 is not setup for use with terrain'
          write(outfile,*)
          write(outfile,*) ' ... stopping inside init3d ... '
          write(outfile,*)
          call stopcm1
        ENDIF

        IF(axisymm.eq.1)THEN
          nref = nx
          xref = 0.0
          do i=1,nref
            rref(i) = 0.5*(xfref(i)+xfref(i+1))
          enddo
        ELSE
          nref = nx/2+1
          xref = xfref(nx/2+1)
          yref = yfref(ny/2+1)
          xmax = 0.5*(xfref(nx)+xfref(nx+1))
          ymax = 0.5*(yfref(ny)+yfref(ny+1))
          do i=1,nref
            rref(i) = 0.5*(xfref(nx/2+i)+xfref(nx/2+i+1))-xref
          enddo
        ENDIF

!!!        print *
!!!        print *,'  v:'
        do k=1,nk
        do i=1,nref
          if(rref(i).lt.r0)then
            dd1 = 2.0 * rmax / ( rref(i) + rmax )
            vr = sqrt( vmax**2 * (rref(i)/rmax)**2     &
            * ( dd1 ** 3 - dd2 ** 3 ) + 0.25*fcor*fcor*rref(i)*rref(i) )   &
                    - 0.5 * fcor * rref(i)
          else
            vr = 0.0
          endif
          if(zh(1,1,k).lt.zdd)then
            vref(i,k) = vr * (zdd-zh(1,1,k))/(zdd-0.0)
          else
            vref(i,k) = 0.0
          endif
!!!          if(k.eq.1) print *,i,xh(ni/2+i),rref(i),vref(i,k)
        enddo
        enddo
!!!        print *

      ! need to iterate for qv to converge:
      DO nloop=1,20

        do k=1,nk
        do i=1,nref
          if(imoist.eq.1)   &
          qvref(i,k) = rh0(1,1,k)*rslf(p00*((pi0(1,1,k)+piref(i,k))**cpdrd),   &
                             (pi0(1,1,k)+piref(i,k))*(th0(1,1,k)+thref(i,k)) )
          thvref(i,k)=(th0(1,1,k)+thref(i,k))*(1.0+reps*qvref(i,k))   &
                                             /(1.0+qvref(i,k))
        enddo
        enddo

!!!        print *,'  pi:'
        do k=1,nk
          piref(nref,k)=0.0
          do i=nref,2,-1
            piref(i-1,k) = piref(i,k)                                       &
         + (rref(i-1)-rref(i))/(cp*0.5*(thvref(i-1,k)+thvref(i,k))) * 0.5 * &
             ( vref(i  ,k)*vref(i  ,k)/rref(i)                              &
              +vref(i-1,k)*vref(i-1,k)/rref(i-1)                            &
               + fcor * ( vref(i,k) + vref(i-1,k) ) )
!!!            if(k.eq.1) print *,i-1,rref(i-1),piref(i-1,k)
          enddo
        enddo
!!!        print *

        do i=1,nref
          piref(i,   0) = piref(i, 1)
          piref(i,nk+1) = piref(i,nk)
        enddo

        do k=2,nk
        do i=1,nref
          thref(i,k) = 0.5*( cp*0.5*(thvref(i,k)+thvref(i,k+1))*(piref(i,k+1)-piref(i,k))*rdz*mf(1,1,k+1)     &
                            +cp*0.5*(thvref(i,k)+thvref(i,k-1))*(piref(i,k)-piref(i,k-1))*rdz*mf(1,1,k) )   &
                          *thv0(1,1,k)/g
          thref(i,k)=(thv0(1,1,k)+thref(i,k))*(1.0+qvref(i,k))/(1.0+reps*qvref(i,k))-th0(1,1,k)
        enddo
        enddo

        k=1
        do i=1,nref
          thref(i,k) = ( cp*0.5*(thvref(i,k)+thvref(i,k+1))*(piref(i,k+1)-piref(i,k))*rdz*mf(1,1,k+1) )   &
                          *thv0(1,1,k)/g
          thref(i,k)=(thv0(1,1,k)+thref(i,k))*(1.0+qvref(i,k))/(1.0+reps*qvref(i,k))-th0(1,1,k)
        enddo

        write(outfile,*) nloop,thref(1,1),qvref(1,1),piref(1,1)

      ENDDO   ! enddo for iteration

        IF(axisymm.eq.1)THEN

          do k=1,nk
          do i=1,ni
             va(i,1,k) =  vref(i,k)
            ppi(i,1,k) = piref(i,k)
            tha(i,1,k) = thref(i,k)
            if(imoist.eq.1) qa(i,1,k,nqv) = qvref(i,k)
          enddo
          enddo

        ELSE

          do j=1,nj+1
          do i=1,ni+1
            ! scalar points:
            rr = sqrt( (xh(i)-xref)**2 + (yh(j)-yref)**2 )
            rr = min( rr , xmax-xref )
            ! need to account for grid stretching.  Do simple search:
            diff = -1.0e20
            ii = 0
            do while( diff.lt.0.0 )
              ii = ii + 1
              if( ii.gt.nref )then
                write(6,*)
                write(6,*) ' ii,nref = ',ii,nref
                write(6,*) ' rr      = ',rr,xmax,xref
                write(6,*) ' rref    = ',rref(ii-1),rref(ii-1)-rr
                write(6,*)
                call stopcm1
              endif
              diff = rref(ii)-rr
            enddo
            i2 = ii
            i1 = i2-1
            frac = (      rr-rref(i1))   &
                  /(rref(i2)-rref(i1))
            do k=1,nk
              ppi(i,j,k) = piref(i1,k)+(piref(i2,k)-piref(i1,k))*frac
              tha(i,j,k) = thref(i1,k)+(thref(i2,k)-thref(i1,k))*frac
              if(imoist.eq.1) qa(i,j,k,nqv) = qvref(i1,k)+(qvref(i2,k)-qvref(i1,k))*frac
            enddo

            ! u:
            rr = sqrt( (xf(i)-xref)**2 + (yh(j)-yref)**2 )
            rr = min( rr , xmax-xref )
            ! need to account for grid stretching.  Do simple search:
            diff = -1.0e20
            ii = 0
            do while( diff.lt.0.0 )
              ii = ii + 1
              if( ii.gt.nref )then
                write(6,*)
                write(6,*) ' ii,nref = ',ii,nref
                write(6,*) ' rr      = ',rr,xmax,xref
                write(6,*)
                call stopcm1
              endif
              diff = rref(ii)-rr
            enddo
            if( abs(rr-rref(ii)).lt.tsmall .and. ii.eq.1 ) ii = 2
            i2 = ii
            i1 = i2-1
            frac = (      rr-rref(i1))   &
                  /(rref(i2)-rref(i1))
            do k=1,nk
              angle = datan2(dble(yh(j)-yref),dble(xf(i)-xref))
              ua(i,j,k) = -( vref(i1,k)+( vref(i2,k)- vref(i1,k))*frac )*sin(angle)
            enddo

            ! v:
            rr = sqrt( (yf(j)-yref)**2 + (xh(i)-xref)**2 )
            rr = min( rr , xmax-xref )
            ! need to account for grid stretching.  Do simple search:
            diff = -1.0e20
            ii = 0
            do while( diff.lt.0.0 )
              ii = ii + 1
              if( ii.gt.nref )then
                write(6,*)
                write(6,*) ' ii,nref = ',ii,nref
                write(6,*) ' rr      = ',rr,xmax,xref
                write(6,*)
                call stopcm1
              endif
              diff = rref(ii)-rr
            enddo
            if( abs(rr-rref(ii)).lt.tsmall .and. ii.eq.1 ) ii = 2
            i2 = ii
            i1 = i2-1
            frac = (      rr-rref(i1))   &
                  /(rref(i2)-rref(i1))
            do k=1,nk
              angle = datan2(dble(yf(j)-yref),dble(xh(i)-xref))
              va(i,j,k) = (vref(i1,k)+( vref(i2,k)- vref(i1,k))*frac )*cos(angle)
            enddo
          enddo
          enddo

!!!          print *
!!!          print *,'  symmtest:'
!!!          j = nj/2 + 5
!!!          k = 1
!!!          do i=1,nref
!!!            print *,i,j,ua(i,j,k),va(j,i,k),ua(i,j,k)+va(j,i,k)
!!!          enddo
!!!          print *

        ENDIF

        call bcu(ua)
        call bcv(va)
        call bcs(ppi)
        call bcs(tha)

        call calcprs(pi0,prs,ppi)

        deallocate(  rref)
        deallocate(  vref)
        deallocate( piref)
        deallocate( thref)
        deallocate(thvref)
        deallocate( qvref)

        setppi = .false.

!-----------------------------------------------------------------------
!  iinit = 8
!  Line thermal with random small-amplitude perturbations

      ELSEIF(iinit.eq.8)THEN

        write(outfile,*)
        write(outfile,*) '  Warm bubble'
        write(outfile,*)

        ric     = 150000.0  ! center of bubble in x-direction (m)
        zc = 4375.0 ! height of center of bubble above ground (m)
        bhrad   =  10000.0  ! horizontal radius of bubble (m)
        bvrad   =   1500.0  ! vertical radius of bubble (m)
        bptpert =      2.0  ! max potential temp perturbation (K)

        ! By default, CM1 sets qv=constant at a constant height level for 
        ! this value of iinit.  If you would rather have rh=constant at 
        ! a constant height level, then set this to .true.
        maintain_rh = .false.

        call random_seed

        do n=1,myid
        do k=1,nk
        do j=1,nj
        do i=1,ni
          call random_number(rand)
        enddo
        enddo
        enddo
        enddo

        do k=1,nk
        do j=1,nj
        do i=1,ni
          call random_number(rand)
          beta=sqrt(                             &
                    ((xh(i)-ric)/bhrad)**2       &
                   +((zh(i,j,k)-zc)/bvrad)**2)
          if(beta.lt.1.0)then
            tha(i,j,k)=bptpert*(cos(0.5*pi*beta)**2)   &
                      +0.2*(2.0*rand-1.0)
          else
            tha(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo

!------------------------------------------------------------------
!  iinit = 9
!  Forced convergence
!  Reference:  Loftus et al, 2008: MWR, v. 136, pp. 2408--2421.

      ELSEIF(iinit.eq.9)THEN

        ! User-defined settings:
        Dmax     =  -1.0e-3     ! maximum divergence (s^{-1})
        zdeep    =  2000.0      ! depth (m) of forced convergence
        lamx     = 10000.0      ! Loftus et al lambda_x parameter
        lamy     = 10000.0      ! Loftus at al lambda_y parameter
        xcent    =     0.0      ! x-location (m)
        ycent    =     0.0      ! y-location (m)
        convtime =   900.0      ! time (s) at beginning of simulation over
                                ! which convergence is applied

        ! Don't change anything below here:
        convinit = 1
        IF( ny.eq.1 )THEN
          ! 2D (x-z):
          Aconv = (-0.5*Dmax)/( (1.0/(lamx**2)) )
          lamy = 1.0e20
        ELSEIF( nx.eq.1 )THEN
          ! 2D (y-z):
          Aconv = (-0.5*Dmax)/( (1.0/(lamy**2)) )
          lamx = 1.0e20
        ELSE
          ! 3D:
          Aconv = (-0.5*Dmax)/( (1.0/(lamx**2))+(1.0/(lamy**2)) )
        ENDIF

!------------------------------------------------------------------
!  DRC test init 1
!  Desc: playing with the initial condition moist vortex

      ELSEIF(iinit.eq.10)THEN

        allocate(  rref(nx)       )
        allocate(  vref(nx,0:nk+1))
        allocate( piref(nx,0:nk+1))
        allocate( thref(nx,0:nk+1))
        allocate(thvref(nx,0:nk+1))
        allocate( qvref(nx,0:nk+1))
        allocate( qvsat(nx,0:nk+1))
        allocate( qv_bub(nx,0:nk+1))

		!values for initial condition (i.e. does NOT include background!)
          rref=0.0	!?
          vref=0.0	!winds?
         piref=0.0	!exner: T=theta*pi
         thref=0.0	!pot temp
        thvref=0.0	!wet-bulb pot temp
         qvref=0.0	!water vapor mixing ratio
         
        IF(ibalance.ne.0)THEN
          write(outfile,*)
          write(outfile,*) ' Please use ibalance = 0 with iinit=7'
          write(outfile,*)
          write(outfile,*) ' ... stopping inside init3d ... '
          write(outfile,*)
          call stopcm1
        ENDIF
        IF(terrain_flag)THEN
          write(outfile,*)
          write(outfile,*) ' iinit=7 is not setup for use with terrain'
          write(outfile,*)
          write(outfile,*) ' ... stopping inside init3d ... '
          write(outfile,*)
          call stopcm1
        ENDIF

        IF(axisymm.eq.1)THEN
          nref = nx	!number of points from center to edge
          xref = 0.0	!TC center is at start of x-array
          do i=1,nref
            rref(i) = 0.5*(xfref(i)+xfref(i+1))	!r is same as x but staggered
          enddo
        ELSE	!3D
          nref = nx/2+1	!number of points from center to edge
          xref = xfref(nx/2+1)	!TC center is in middle of x- and y-array
          yref = yfref(ny/2+1)
          xmax = 0.5*(xfref(nx)+xfref(nx+1))	!edge is half of domain width
          ymax = 0.5*(yfref(ny)+yfref(ny+1))
          do i=1,nref
            rref(i) = 0.5*(xfref(nx/2+i)+xfref(nx/2+i+1))-xref
          enddo
        ENDIF

!!! START INITIAL MID-LEVEL VORTEX
!!!        print *
!!!        print *,'  v:'
		! calculate v(r,z) at each point using eqn (37) of RE87

        z_bltop = 1500.0 !height of boundary layer top [m]
		
        !Initial vortex
        r0 = 400000.0 ! outer radius of vortex, where V=0 [m]
        rmax = 54600.00 ! rmax [m]
        vmax = 1.6336994E-04 ! max wind speed in vortex [m/s]
        zc = 4375.0 ! height of center of vortex above ground [m]
        dz0 = 2875.0 ! vertical scale of vortex edge from center [m]

		!convert units
!		r0 = 400000.0 ![km]-->[m]
!		rmax = 54600.00 ![km]-->[m]

		!additional length scales
        zmin = max(zc-dz0,z_bltop) ! height of bottom edge of vortex above ground [m]; must stay above boundary layer
        zmax = zc+dz0 ! height of top edge of vortex above ground [m]

!        Lhv   = .5*(r0_qv-rc_qv)  ! vertical decay scale of vortex winds [m]
!        Lz   = .5*(.5*((zmax-zc)+(zc-zmin)))   ! horizontal decay scale of vortex winds [m]

        dd2 = 2.0 * rmax / ( r0 + rmax )	!term in radial part of v(r,z) equation
		
        do k=1,nk
        do i=1,nref
          if(rref(i).lt.r0)then
            dd1 = 2.0 * rmax / ( rref(i) + rmax )	!1st term in parentheses of (37)
            vr = sqrt( vmax**2 * (rref(i)/rmax)**2     &
            * ( dd1 ** 3 - dd2 ** 3 ) + 0.25*fcor*fcor*rref(i)*rref(i) )   &
                    - 0.5 * fcor * rref(i)			!v(r), no z-dependence yet

          else	!r>=r0
            vr = 0.0
          endif
          if(zh(1,1,k).le.zmax.and.zh(1,1,k).ge.zc)then
            vref(i,k) = vr * (zmax-zh(1,1,k))/(zmax-zc)	!v(r) decays linearly upwards from zc with height
		  elseif(zh(1,1,k).le.zc.and.zh(1,1,k).ge.zmin)then
            vref(i,k) = vr * (zh(1,1,k)-zmin)/(zc-zmin)	!v(r) decays linearly downwards from zc with height
          else						
            vref(i,k) = 0.0	!v(r) = 0 at and above z_sponge
          endif
!!!          if(k.eq.1) print *,i,xh(ni/2+i),rref(i),vref(i,k)
        enddo
        enddo
!!!        print *
!!! END INITIAL MID-LEVEL VORTEX

!!! START RELATIVE HUMIDITY BUBBLE DRC 03-08-11
!!!        print *,'  qv:'

		!!!! calculate qv(r,z) at each point

      	!basic state from sounding
        do k=1,nk
        do i=1,nref        
			qvsat(i,k) = rslf(p00*((pi0(1,1,k)+piref(i,k))**cpdrd),   &
                             (pi0(1,1,k)+piref(i,k))*(th0(1,1,k)+thref(i,k)) )
	        qvref(i,k) = rh0(1,1,k)*qvsat(i,k)	!this is the basic state based on input_sounding alone	    
        enddo
        enddo
	       
        !Initial relative humidity anomaly
        rc_qv = 0.0 ! radius of center of anomaly [m]
        r0_qv = 200000.0 ! outer radius of anomaly [m]
        rhpert_max = 1.9 ! max relative humidity perturbation (%); 0.01455 kg/kg --> RH=80% at lowest level
        zc_qv = 4375.0 ! height of center of bubble above ground [m]
        dz0_qv = 5000.0 ! vertical scale of bubble edge from center [m]
        shape_qv = 2 ! 1=gaussian, 2=constant

		!convert units
        zmax_qv = zc_qv+dz0_qv ! height of top edge of bubble above ground [m]
        zmin_qv = max(zc_qv-dz0_qv,z_bltop) ! height of bottom edge of bubble above ground; must stay above boundary layer
!		rc_qv = 0.0 ![km]-->[m]
!		r0_qv = 200000.0 ![km]-->[m]
!		rhpert_max = rhpert_max/100

		!additional length scales		
        Lh_qv   = .5*(r0_qv-rc_qv)  ! vertical decay scale of bubble [m]
        Lz_qv   = .5*(.5*((zmax_qv-zc_qv)+(zc_qv-zmin_qv)))   ! horizontal decay scale of bubble [m]
        
        
        do k=1,nk
        do i=1,nref
          if(rref(i).lt.(rc_qv+r0_qv).and.rref(i).ge.(max(0.0,rc_qv-r0_qv)).and.zh(1,1,k).lt.zmax_qv.and.zh(1,1,k).gt.zmin_qv)then
			if(rhpert_max.lt.(1-(qvref(i,k)/qvsat(i,k))))then	!perturbation will not saturate gridbox
				if(shape_qv.eq.1)then
		            qv_bub(i,k)=rhpert_max*qvsat(i,k)*(exp(-((rref(i)-rc_qv)/Lh_qv)**2))*(exp(-((zh(1,1,k)-zc_qv)/Lz_qv)**2))
		        elseif(shape_qv.eq.2)then
		            qv_bub(i,k)=rhpert_max*qvsat(i,k)
		    	endif
			else	!perturbation will saturate gridbox, thus set qv_bub that will cause RH=100%
				if(shape_qv.eq.1)then
					qv_bub(i,k)=(1-(qvref(i,k)/qvsat(i,k)))*qvsat(i,k)*(exp(-((rref(i)-rc_qv)/Lh_qv)**2))*(exp(-((zh(1,1,k)-zc_qv)/Lz_qv)**2))
		        elseif(shape_qv.eq.2)then
					qv_bub(i,k)=(1-(qvref(i,k)/qvsat(i,k)))*qvsat(i,k)
				endif
			endif
          else
            qv_bub(i,k)=0.0
          endif

          qvref(i,k)=qvref(i,k)+qv_bub(i,k)
!          qvref(i,k)=qvref(i,k)/qvsat(i,k)
!          qvref(i,k)=qv_bub(i,k)

        enddo
        enddo
		
		write(outfile,*) qv_bub(1,1:10), zh(1,1,1:10)
		
!          if(rref(i).lt.r0)then
!        	qv_bub = rhpert_max* (rref(i)/rmax)**2     &
!            * ( dd1 ** 3 - dd2 ** 3 ) + 0.25*fcor*fcor*rref(i)*rref(i) )   &
!                    - 0.5 * fcor * rref(i)			!v(r), no z-dependence yet                    
!		          else	!r>=r0
!            qvr = 0.0
!          endif
!!!          if(k.eq.1) print *,i,xh(ni/2+i),rref(i),vref(i,k)

!!!        print *

!!! END MOISTURE BUBBLE DRC 03-08-11

      !!! iterate to calculate piref corresponding to v(r,z) in GWB and thvref corresponding to piref given qvref
      DO nloop=1,20

        do k=1,nk
        do i=1,nref
          if(imoist.eq.1)   &	!rh0=initial rel humidity (calculated in isnd=7 code)
          thvref(i,k)=(th0(1,1,k)+thref(i,k))*(1.0+reps*qvref(i,k))   &
                                             /(1.0+qvref(i,k))
                                             
                    !reps = 1/eps, eps = rd/rv = .622
                    !rslf = saturation mixing ratio wrt liquid water function (i.e. qv_sat(P,T))--located in misclibs.F
                    !p00=100000 Pa
                    !cpdrd=cpd/rd (ratio of specific heat of dry air at constant pressure to dry gas constant)
        enddo
        enddo

!!!        print *,'  pi:'
        do k=1,nk
          piref(nref,k)=0.0
          do i=nref,2,-1
            piref(i-1,k) = piref(i,k)                                       &
         + (rref(i-1)-rref(i))/(cp*0.5*(thvref(i-1,k)+thvref(i,k))) * 0.5 * &
             ( vref(i  ,k)*vref(i  ,k)/rref(i)                              &
              +vref(i-1,k)*vref(i-1,k)/rref(i-1)                            &
               + fcor * ( vref(i,k) + vref(i-1,k) ) )
!!!            if(k.eq.1) print *,i-1,rref(i-1),piref(i-1,k)
          enddo
        enddo
!!!        print *

        do i=1,nref
          piref(i,   0) = piref(i, 1)
          piref(i,nk+1) = piref(i,nk)
        enddo

        do k=2,nk
        do i=1,nref
          thref(i,k) = 0.5*( cp*0.5*(thvref(i,k)+thvref(i,k+1))*(piref(i,k+1)-piref(i,k))*rdz*mf(1,1,k+1)     &
                            +cp*0.5*(thvref(i,k)+thvref(i,k-1))*(piref(i,k)-piref(i,k-1))*rdz*mf(1,1,k) )   &
                          *thv0(1,1,k)/g
          thref(i,k)=(thv0(1,1,k)+thref(i,k))*(1.0+qvref(i,k))/(1.0+reps*qvref(i,k))-th0(1,1,k)
        enddo
        enddo

        k=1
        do i=1,nref
          thref(i,k) = ( cp*0.5*(thvref(i,k)+thvref(i,k+1))*(piref(i,k+1)-piref(i,k))*rdz*mf(1,1,k+1) )   &
                          *thv0(1,1,k)/g
          thref(i,k)=(thv0(1,1,k)+thref(i,k))*(1.0+qvref(i,k))/(1.0+reps*qvref(i,k))-th0(1,1,k)
        enddo
		
        write(outfile,*) nloop,thref(1,1),qvref(1,1),piref(1,1),vref(7,1),relh

      ENDDO   ! enddo for iteration

!		do i=1,nk
!        write(outfile,*) rh0(1,1,i)
!        enddo
        
        IF(axisymm.eq.1)THEN

          do k=1,nk
          do i=1,ni
             va(i,1,k) =  vref(i,k)
            ppi(i,1,k) = piref(i,k)
            tha(i,1,k) = thref(i,k)
            if(imoist.eq.1) qa(i,1,k,nqv) = qvref(i,k)
          enddo
          enddo

        ELSE

          do j=1,nj+1
          do i=1,ni+1
            ! scalar points:
            rr = sqrt( (xh(i)-xref)**2 + (yh(j)-yref)**2 )
            rr = min( rr , xmax-xref )
            ! need to account for grid stretching.  Do simple search:
            diff = -1.0e20
            ii = 0
            do while( diff.lt.0.0 )
              ii = ii + 1
              if( ii.gt.nref )then
                write(6,*)
                write(6,*) ' ii,nref = ',ii,nref
                write(6,*) ' rr      = ',rr,xmax,xref
                write(6,*) ' rref    = ',rref(ii-1),rref(ii-1)-rr
                write(6,*)
                call stopcm1
              endif
              diff = rref(ii)-rr
            enddo
            i2 = ii
            i1 = i2-1
            frac = (      rr-rref(i1))   &
                  /(rref(i2)-rref(i1))
            do k=1,nk
              ppi(i,j,k) = piref(i1,k)+(piref(i2,k)-piref(i1,k))*frac
              tha(i,j,k) = thref(i1,k)+(thref(i2,k)-thref(i1,k))*frac
              if(imoist.eq.1) qa(i,j,k,nqv) = qvref(i1,k)+(qvref(i2,k)-qvref(i1,k))*frac
            enddo

            ! u:
            rr = sqrt( (xf(i)-xref)**2 + (yh(j)-yref)**2 )
            rr = min( rr , xmax-xref )
            ! need to account for grid stretching.  Do simple search:
            diff = -1.0e20
            ii = 0
            do while( diff.lt.0.0 )
              ii = ii + 1
              if( ii.gt.nref )then
                write(6,*)
                write(6,*) ' ii,nref = ',ii,nref
                write(6,*) ' rr      = ',rr,xmax,xref
                write(6,*)
                call stopcm1
              endif
              diff = rref(ii)-rr
            enddo
            if( abs(rr-rref(ii)).lt.tsmall .and. ii.eq.1 ) ii = 2
            i2 = ii
            i1 = i2-1
            frac = (      rr-rref(i1))   &
                  /(rref(i2)-rref(i1))
            do k=1,nk
              angle = datan2(dble(yh(j)-yref),dble(xf(i)-xref))
              ua(i,j,k) = -( vref(i1,k)+( vref(i2,k)- vref(i1,k))*frac )*sin(angle)
            enddo

            ! v:
            rr = sqrt( (yf(j)-yref)**2 + (xh(i)-xref)**2 )
            rr = min( rr , xmax-xref )
            ! need to account for grid stretching.  Do simple search:
            diff = -1.0e20
            ii = 0
            do while( diff.lt.0.0 )
              ii = ii + 1
              if( ii.gt.nref )then
                write(6,*)
                write(6,*) ' ii,nref = ',ii,nref
                write(6,*) ' rr      = ',rr,xmax,xref
                write(6,*)
                call stopcm1
              endif
              diff = rref(ii)-rr
            enddo
            if( abs(rr-rref(ii)).lt.tsmall .and. ii.eq.1 ) ii = 2
            i2 = ii
            i1 = i2-1
            frac = (      rr-rref(i1))   &
                  /(rref(i2)-rref(i1))
            do k=1,nk
              angle = datan2(dble(yf(j)-yref),dble(xh(i)-xref))
              va(i,j,k) = (vref(i1,k)+( vref(i2,k)- vref(i1,k))*frac )*cos(angle)
            enddo
          enddo
          enddo

!!!          print *
!!!          print *,'  symmtest:'
!!!          j = nj/2 + 5
!!!          k = 1
!!!          do i=1,nref
!!!            print *,i,j,ua(i,j,k),va(j,i,k),ua(i,j,k)+va(j,i,k)
!!!          enddo
!!!          print *

        ENDIF

        call bcu(ua)
        call bcv(va)
        call bcs(ppi)
        call bcs(tha)

        call calcprs(pi0,prs,ppi)

        deallocate(  rref)
        deallocate(  vref)
        deallocate( piref)
        deallocate( thref)
        deallocate(thvref)
        deallocate( qvref)

        setppi = .false.

!------------------------------------------------------------------

      ENDIF    ! end of iinit options

!------------------------------------------------------------------

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Random perturbations:

      IF( irandp.eq.1 )THEN

        ! this is the amplitude of the theta perturbations
        ! (plus or minus this value in K)
        amplitude = 0.5

        ! initialize the random number generator
        call random_seed

        ! this makes sure that each processor has a different set
        ! of random numbers for MPI runs
        do n=1,myid
        do k=1,nk
        do j=1,nj
        do i=1,ni
          call random_number(rand)
        enddo
        enddo
        enddo
        enddo

        ! random numbers added here
        ! (can be modified to only place perturbations in certain
        !  locations, but this default code simply puts them
        !  everywhere)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          call random_number(rand)
          tha(i,j,k)=tha(i,j,k)+amplitude*(2.0*rand-1.0)
        enddo
        enddo
        enddo

      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!----------------------------------------------
!  arrays for elliptic solver

      do k=kpb,kpe
      do j=jpb,jpe
      do i=ipb,ipe
        cfb(i,j,k) = 0.0
      enddo
      enddo
      enddo


      do k=kpb,kpe
        cfa(k) = 0.0
        cfc(k) = 0.0
        ad1(k) = 0.0
        ad2(k) = 0.0
      enddo

      IF( (ibalance.eq.2).or.(psolver.eq.4).or.(psolver.eq.5) )THEN

        dpi = 4.0d0*datan(1.0d0)
        write(outfile,*) '  dpi = ',dpi

        IF(psolver.le.3)THEN
          do k=1,nk
            cfa(k)=mh(1,1,k)*mf(1,1,k  )*rf0(1,1,k  )*0.5*(thv0(1,1,k-1)+thv0(1,1,k))/(dz*dz*rho0(1,1,k)*thv0(1,1,k))
            cfc(k)=mh(1,1,k)*mf(1,1,k+1)*rf0(1,1,k+1)*0.5*(thv0(1,1,k)+thv0(1,1,k+1))/(dz*dz*rho0(1,1,k)*thv0(1,1,k))
            ad1(k) = 1.0/(cp*rho0(1,1,k)*thv0(1,1,k))
            ad2(k) = 1.0
          enddo
          cfa( 1) = 0.0
          cfc(nk) = 0.0
          do j=jpb,jpe
          do i=ipb,ipe
            do k=1,nk
              cfb(i,j,k)=2.0d0*( dcos(2.0d0*dpi*dble(i-1)/dble(ipe))          &
                                +dcos(2.0d0*dpi*dble(j-1)/dble(jpe))          &
                                -2.0d0)/(dx*dx) - cfa(k) - cfc(k)
            enddo
          enddo
          enddo
        ELSE
          do k=1,nk
            cfa(k)=mh(1,1,k)*mf(1,1,k  )*rf0(1,1,k  )/(dz*dz*rho0(1,1,k-1))
            cfc(k)=mh(1,1,k)*mf(1,1,k+1)*rf0(1,1,k+1)/(dz*dz*rho0(1,1,k+1))
            ad1(k) = 1.0
            ad2(k) = 1.0/rho0(1,1,k)
          enddo
          cfa( 1) = 0.0
          cfc(nk) = 0.0
          do j=jpb,jpe
          do i=ipb,ipe
            do k=2,nk-1
              cfb(i,j,k)=2.0d0*( dcos(2.0d0*dpi*dble(i-1)/dble(ipe))          &
                                +dcos(2.0d0*dpi*dble(j-1)/dble(jpe))          &
                                -2.0d0)/(dx*dx)                               &
                    -mh(1,1,k)*mf(1,1,k+1)*rf0(1,1,k+1)/(dz*dz*rho0(1,1,k))   &
                    -mh(1,1,k)*mf(1,1,k  )*rf0(1,1,k  )/(dz*dz*rho0(1,1,k))
            enddo
            cfb(i,j,1)=2.0d0*( dcos(2.0d0*dpi*dble(i-1)/dble(ipe))          &
                              +dcos(2.0d0*dpi*dble(j-1)/dble(jpe))          &
                              -2.0d0)/(dx*dx)                               &
                  -mh(1,1,1)*mf(1,1,2  )*rf0(1,1,2  )/(dz*dz*rho0(1,1,1))
            cfb(i,j,nk)=2.0d0*( dcos(2.0d0*dpi*dble(i-1)/dble(ipe))          &
                              +dcos(2.0d0*dpi*dble(j-1)/dble(jpe))          &
                              -2.0d0)/(dx*dx)                               &
                  -mh(1,1,nk)*mf(1,1,nk  )*rf0(1,1,nk  )/(dz*dz*rho0(1,1,nk))
          enddo
          enddo
        ENDIF

      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!-----------------------------------------------------------------
!  Get 3d pressure
        
      if(imoist.eq.1 .and. maintain_rh)then

        !! maintain rh
        write(outfile,*)
        write(outfile,*) '  Constant rh across domain:'
        write(outfile,*)

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          qa(i,j,k,nqv)=rh0(i,j,k)*rslf(prs0(i,j,k),(th0(i,j,k)+tha(i,j,k))*pi0(i,j,k))
        enddo
        enddo
        enddo

      endif


    IF(setppi)THEN

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        ppi(i,j,k)=0.0
      enddo
      enddo
      enddo

      IF(ibalance.eq.1)THEN

        ! hydrostatic balance ... integrate top-down

        do j=1,nj
        do i=1,ni
          ! virtual potential temperature

          if(imoist.eq.1)then
            do k=1,nk
              qt=0.0
              do n=nql1,nql2
                qt=qt+qa(i,j,k,n)
              enddo
              if(iice.eq.1)then
                do n=nqs1,nqs2
                  qt=qt+qa(i,j,k,n)
                enddo
              endif
              thvnew(k)=(th0(i,j,k)+tha(i,j,k))*(1.0+reps*qa(i,j,k,nqv))   &
                                               /(1.0+qa(i,j,k,nqv)+qt)
            enddo
          else
            do k=1,nk
              thvnew(k)=th0(i,j,k)+tha(i,j,k)
            enddo
          endif

          ! non-dimensional pressure
          pinew(nk)=pi0(i,j,nk)
          do k=nk-1,1,-1
            pinew(k)=pinew(k+1)+g*(zh(i,j,k+1)-zh(i,j,k))   &
                    /(cp*0.5*(thvnew(k+1)+thvnew(k)))
          enddo

          ! new pressure
          do k=1,nk
            ppi(i,j,k)=pinew(k)-pi0(i,j,k)
            if(abs(ppi(i,j,k)).lt.1.0e-6) ppi(i,j,k)=0.0
          enddo

        enddo
        enddo

      ELSEIF(ibalance.eq.2)THEN

        write(outfile,*)
        write(outfile,*) '  ibalance = 2'
        write(outfile,*)

        if(stretch_x.ge.1.or.stretch_y.ge.1)then
          print *,'  this option not supported with horizontal grid stretching'
          print *,'  (yet)'
          call stopcm1
        endif

#ifdef MPI
        print *,'  This option is not (yet) supported in MPI mode'
        print *,'  (sorry)'
        print *
        call stopcm1
#endif

        ! buoyancy pressure

        ! th3d stores theta-v

        if(imoist.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qt=0.0
            do n=nql1,nql2
              qt=qt+qa(i,j,k,n)
            enddo
            if(iice.eq.1)then
              do n=nqs1,nqs2
                qt=qt+qa(i,j,k,n)
              enddo
            endif
            th3d(i,j,k)=(th0(i,j,k)+tha(i,j,k))*(1.0+reps*qa(i,j,k,nqv))   &
                       /(1.0+qa(i,j,k,nqv)+qt)
          enddo
          enddo
          enddo
        else
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            th3d(i,j,k)=th0(i,j,k)+tha(i,j,k)
          enddo
          enddo
          enddo
        endif

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          th3d(i,j,0   ) = th3d(i,j,1)
          th3d(i,j,nk+1) = th3d(i,j,nk)
        enddo
        enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum4(i,j,k)=g*( th3d(i,j,k)/thv0(i,j,k)-1.0 )
        enddo
        enddo
        enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          dum4(i,j,0   ) = -dum4(i,j,1)
          dum4(i,j,nk+1) = -dum4(i,j,nk)
        enddo
        enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          wten(i,j,k)=0.5*( dum4(i,j,k-1)+dum4(i,j,k) )
        enddo
        enddo
        enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          ppi(i,j,k)=0.0
          dum3(i,j,k)=0.0
          divx(i,j,k)=0.0
          uten(i,j,k)=0.0
          vten(i,j,k)=0.0
        enddo
        enddo
        enddo

        call poiss(uh,vh,mh,rmh,mf,rmf,pi0,thv0,rho0,rf0,    &
                   dum3,divx,ppi,uten,vten,wten,             &
                   cfb,cfa,cfc,ad1,ad2,pdt,deft,rhs,trans,dtl)

        IF(psolver.eq.4.or.psolver.eq.5)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=kb,ke
          do j=jb,je
          do i=ib,ie
            ppi(i,j,k)=((prs0(1,1,k)+ppi(i,j,k)*rho0(1,1,k))*rp00)**rovcp   &
                      -pi0(1,1,k)
            pp3d(i,j,k)=ppi(i,j,k)
          enddo
          enddo
          enddo

        ENDIF

        call bcs(ppi)

      ENDIF

    ENDIF

!------------------------------------------------------------------

      write(outfile,*)
      write(outfile,*) 'Leaving INIT3D'

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getset(dzdx,dzdy,pi0,th0,rho0,prs0,rho,prs,               &
                        ua,u3d,va,v3d,wa,w3d,ppi,pp3d,                     &
                        tha,th3d,qa,q3d,tkea,tke3d,pta,pt3d,               &
                        reqs_u,reqs_v,reqs_w,reqs_s,reqs_tk,               &
                        uw31,uw32,ue31,ue32,us31,us32,un31,un32,           &
                        vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,           &
                        ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,           &
                        sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,           &
                        tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      real, dimension(itb:ite,jtb:jte) :: dzdx,dzdy
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,th0,rho0,prs0
      real, dimension(ib:ie,jb:je,kb:ke) :: rho,prs
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua,u3d
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va,v3d
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa,w3d
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d
      real, dimension(ib:ie,jb:je,kb:ke) :: tha,th3d
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa,q3d
      real, dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d
      real, dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta,pt3d
      integer, dimension(rmp) :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_tk
      real, dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, dimension(imp,cmp,kmp)   :: vs31,vs32,vn31,vn32
      real, dimension(cmp,jmp,kmp-1) :: ww31,ww32,we31,we32
      real, dimension(imp,cmp,kmp-1) :: ws31,ws32,wn31,wn32
      real, dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, dimension(cmp,jmp,kmt)   :: tkw1,tkw2,tke1,tke2
      real, dimension(imp,cmp,kmt)   :: tks1,tks2,tkn1,tkn2

!----------
 
      integer i,j,k,n

!------------------------------------------------------------------
!  Make sure boundary values are set properly

      write(outfile,*) 'Inside GETSET'
      write(outfile,*)

      call bcu(ua)
      call bcv(va)
      call bcw(wa,1)
      call bcs(ppi)
      call bcs(tha)
      if(imoist.eq.1)then
        do n=1,numq
          call bcs(qa(ibm,jbm,kbm,n))
        enddo
      endif
      if(iturb.eq.1)then
        call bct(tkea)
      endif
      if(iptra.eq.1)then
        do n=1,npt
          call bcs(pta(ib,jb,kb,n))
        enddo
      endif

#ifdef MPI
!------------------------------------------------------------------

      nf=0
      nu=0
      nv=0
      nw=0

      call comm_3u_start(ua,uw31,uw32,ue31,ue32,   &
                            us31,us32,un31,un32,reqs_u)
      call comm_3u_end(ua,uw31,uw32,ue31,ue32,   &
                          us31,us32,un31,un32,reqs_u)

      call comm_3v_start(va,vw31,vw32,ve31,ve32,   &
                            vs31,vs32,vn31,vn32,reqs_v)
      call comm_3v_end(va,vw31,vw32,ve31,ve32,   &
                          vs31,vs32,vn31,vn32,reqs_v)

      call comm_3w_start(wa,ww31,ww32,we31,we32,   &
                            ws31,ws32,wn31,wn32,reqs_w)
      call comm_3w_end(wa,ww31,ww32,we31,we32,   &
                          ws31,ws32,wn31,wn32,reqs_w)

      call comm_3s_start(ppi,sw31,sw32,se31,se32,   &
                             ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(ppi,sw31,sw32,se31,se32,   &
                           ss31,ss32,sn31,sn32,reqs_s)

      call comm_3s_start(tha,sw31,sw32,se31,se32,   &
                             ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(tha,sw31,sw32,se31,se32,   &
                           ss31,ss32,sn31,sn32,reqs_s)

      IF(imoist.eq.1)THEN
        do n=1,numq
          call comm_3s_start(qa(ibm,jbm,kbm,n),sw31,sw32,se31,se32,   &
                                               ss31,ss32,sn31,sn32,reqs_s)
          call comm_3s_end(qa(ibm,jbm,kbm,n),sw31,sw32,se31,se32,   &
                                             ss31,ss32,sn31,sn32,reqs_s)
        enddo
      ENDIF

      IF(iturb.eq.1)THEN
        call comm_3t_start(tkea,tkw1,tkw2,tke1,tke2,   &
                                tks1,tks2,tkn1,tkn2,reqs_tk)
        call comm_3t_end(tkea,tkw1,tkw2,tke1,tke2,   &
                              tks1,tks2,tkn1,tkn2,reqs_tk)
      ENDIF

      IF(iptra.eq.1)THEN
        do n=1,npt
          call comm_3s_start(pta(ib,jb,kb,n),sw31,sw32,se31,se32,   &
                                             ss31,ss32,sn31,sn32,reqs_s)
          call comm_3s_end(pta(ib,jb,kb,n),sw31,sw32,se31,se32,   &
                                           ss31,ss32,sn31,sn32,reqs_s)
        enddo
      ENDIF

      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

#endif
      if(terrain_flag)then
        call bcwsfc(dzdx,dzdy,ua,va,wa)
        call bc2d(wa(ib,jb,1))
      endif
!------------------------------------------------------------------
!  Get stuff

    IF(psolver.eq.4.or.psolver.eq.5)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        rho(i,j,k)=rho0(i,j,k)
        prs(i,j,k)=prs0(i,j,k)
      enddo
      enddo
      enddo

    ELSE

      call calcprs(pi0,prs,ppi)
 
      call calcrho(pi0,th0,rho,prs,ppi,tha,qa)

    ENDIF

!------------------------------------------------------------------

      do k=kb,ke
      do j=jb,je
      do i=ib,ie+1
        u3d(i,j,k)=ua(i,j,k)
      enddo
      enddo
      enddo
 
      do k=kb,ke
      do j=jb,je+1
      do i=ib,ie
        v3d(i,j,k)=va(i,j,k)
      enddo
      enddo
      enddo
 
      do k=kb,ke+1
      do j=jb,je
      do i=ib,ie
        w3d(i,j,k)=wa(i,j,k)
      enddo
      enddo
      enddo

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        pp3d(i,j,k)=ppi(i,j,k)
        th3d(i,j,k)=tha(i,j,k)
      enddo
      enddo
      enddo

      if(imoist.eq.1)then
        do n=1,numq
        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          q3d(i,j,k,n)=qa(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
      endif

      if(iturb.eq.1)then
        do k=kbt,ket
        do j=jbt,jet
        do i=ibt,iet
          tke3d(i,j,k)=tkea(i,j,k)
        enddo
        enddo
        enddo
      endif

      if(iptra.eq.1)then
        do n=1,npt
        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          pt3d(i,j,k,n)=pta(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
      endif

      write(outfile,*)
      write(outfile,*) 'Leaving GETSET'
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine convinitu(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibw,ibe,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xf,yh,zh,u0,u3d)
      implicit none

      integer, intent(in) :: myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibw,ibe
      real, intent(in) :: zdeep,lamx,lamy,xcent,ycent,aconv
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in),    dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u3d

      integer :: i,j,k
      real :: term1,term2,term3,term4,umo

!!!      if(myid.eq.0) print *,'    convinitu '
!$omp parallel do default(shared)   &
!$omp private(i,j,k,term1,term2,term3,term4,umo)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        term4 = (zdeep-0.5*(zh(i-1,j,k)+zh(i,j,k)))/zdeep
        if (term4 .gt. 0.0) then
          term1 = -(2.0*Aconv*(xf(i)-xcent))/(lamx**2)
          term2 = -((xf(i)-xcent)/lamx)**2
          term3 = -((yh(j)-ycent)/lamy)**2
          umo = term1*(exp(term2)*exp(term3))*term4
          if( abs(umo).gt.0.01 ) u3d(i,j,k) = u0(i,j,k)+umo
        endif
      enddo
      enddo
      enddo

      return
      end subroutine convinitu


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine convinitv(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibs,ibn,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xh,yf,zh,v0,v3d)
      implicit none

      integer, intent(in) :: myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibs,ibn
      real, intent(in) :: zdeep,lamx,lamy,xcent,ycent,aconv
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in),    dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v3d

      integer :: i,j,k
      real :: term1,term2,term3,term4,vmo

!!!      if(myid.eq.0) print *,'    convinitv '
!$omp parallel do default(shared)   &
!$omp private(i,j,k,term1,term2,term3,term4,vmo)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        term4 = (zdeep-0.5*(zh(i,j-1,k)+zh(i,j,k)))/zdeep
        if (term4 .gt. 0.0) then
          term1 = -(2.0*Aconv*(yf(j)-ycent))/(lamy**2)
          term2 = -((xh(i)-xcent)/lamx)**2
          term3 = -((yf(j)-ycent)/lamy)**2
          vmo = term1*(exp(term2)*exp(term3))*term4
          if( abs(vmo).gt.0.01 ) v3d(i,j,k) = v0(i,j,k)+vmo
        endif
      enddo
      enddo
      enddo

      return
      end subroutine convinitv


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
