
!-----------------------------------------------------------------------
!  CM1 Numerical Model, Release 15  (cm1r15)
!  13 January 2011
!  http://www.mmm.ucar.edu/people/bryan/cm1/
!-----------------------------------------------------------------------
!  Quick Index:
!    ua/u3d     = velocity in x-direction (m/s)
!    va/v3d     = velocity in y-direction (m/s)
!    wa/w3d     = velocity in z-direction (m/s)
!    tha/th3d   = perturbation potential temperature (K)
!    ppi/pp3d   = perturbation nondimensional pressure ("Exner function")
!    qa/q3d     = mixing ratios of moisture (kg/kg)
!    tkea/tke3d = SUBGRID turbulence kinetic energy (m^2/s^2)
!    kmh/kmv    = turbulent diffusion coefficients for momentum (m^2/s)
!    khh/khv    = turbulent diffusion coefficients for scalars (m^2/s)
!                 (h = horizontal, v = vertical)
!    prs        = pressure (Pa)
!    rho        = density (kg/m^3)
!
!    th0,pi0,prs0,etc = base-state arrays
!
!    xh         = x (m) at scalar points
!    xf         = x (m) at u points
!    yh         = y (m) at scalar points
!    yf         = y (m) at v points
!    zh         = z (m above sea level) of scalar points (aka, "half levels")
!    zf         = z (m above sea level) of w points (aka, "full levels")
!
!    For the axisymmetric model (axisymm=1), xh and xf are radius (m).
!
!  See "The governing equations for CM1" for more details:
!        http://www.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf
!-----------------------------------------------------------------------
!  Some notes:
!
!  - Upon entering solve, the arrays ending in "a" (eg, ua,wa,tha,qa,etc)
!    are equivalent to the arrays ending in "3d" (eg, u3d,w3d,th3d,q3d,etc).
!  - The purpose of solve is to update the variables from time "t" to time
!    "t+dt".  Values at time "t+dt" are stored in the "3d" arrays.
!  - The "ghost zones" (boundaries beyond the computational subdomain) are
!    filled out completely (3 rows/columns) for the "3d" arrays.  To save 
!    unnecessary computations, starting with cm1r15 the "ghost zones" of 
!    the "a" arrays are only filled out to 1 row/column.  Hence, if you 
!    need to do calculations that use a large stencil, you must use the 
!    "3d" arrays (not the "a") arrays.
!  - Arrays named "ten" store tendencies.  Those ending "ten1" store
!    pre-RK tendencies that are calculated once and then held fixed during
!    the RK (Runge-Kutta) sub-steps. 
!  - CM1 uses a low-storage three-step Runge-Kutta scheme.  See Wicker
!    and Skamarock (2002, MWR, p 2088) for more information.
!  - CM1 uses a staggered C grid.  Hence, u arrays have one more grid point
!    in the i direction, v arrays have one more grid point in the j 
!    direction, and w arrays have one more grid point in the k direction
!    (compared to scalar arrays).
!  - CM1 assumes the subgrid turbulence parameters (tke,km,kh) are located
!    at the w points. 
!-----------------------------------------------------------------------

      subroutine solve(nstep,nrec,prec,nwrite,nrst,rbufsz,num_soil_layers,ndt, &
                   dt,dtlast,mtime,stattim,taptim,rsttim,radtim,adt,acfl, &
                   dodrag,dosfcflx,cloudvar,rhovar,qname,budname,bud,bud2,qbudget,asq,bsq, &
                   xh,rxh,uh,ruh,xf,rxf,uf,ruf,yh,vh,rvh,yf,vf,rvf,  &
                   sigma,sigmaf,tauh,taus,zh,mh,rmh,tauf,zf,mf,rmf,  &
                   rstat,rho0s,pi0s,prs0s,rth0s,pi0,rho0,prs0,thv0,th0,qv0,qc0, &
                   ql0,rr0,rf0,rrf0,rru0,rrv0,                       &
                   zs,gz,dzdx,dzdy,rain,sws,svs,sps,srs,sgs,sus,shs, &
                   doimpl,tsk,thflux,qvflux,cdu,cdv,ce,u1,v1,w1,     &
                   radbcw,radbce,radbcs,radbcn,                      &
                   dum1,dum2,dum3,dum4,divx,rho,prs,                 &
                   t11,t12,t13,t22,t23,t33,                          &
                   gx,u0,rru,ua,u3d,uten,uten1,                      &
                   gy,v0,rrv,va,v3d,vten,vten1,                      &
                   rrw,wa,w3d,wten,wten1,ppi,pp3d,ppten,sten,        &
                   tha,th3d,thten,thten1,thterm,tk,qa,q3d,qten,zvdarray, &
                   kmh,kmv,khh,khv,tkea,tke3d,tketen,                &
                   dissten,thpten,qvpten,qcpten,qipten,upten,vpten,  &
                   swten,lwten,o30,radsw,rnflx,radswnet,radlwin,rad2d,   &
                   x,y,z,za,zp,                                      &
                   lu_index,kpbl2d,psfc,u10,v10,hfx,qfx,xland,znt,ust,  &
                   hpbl,wspd,psim,psih,gz1oz0,br,                    &
                   CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                    &
                   MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                 &
                   CK,CKA,CD,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD, &
                   f2d,gsw,glw,chklowq,capg,snowc,dsxy,              &
                   slab_zs,slab_dzs,tslb,tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,       &
                   pta,pt3d,ptten,                                   &
                   pdata,cfb,cfa,cfc,ad1,ad2,pdt,deft,rhs,trans,     &
                   reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_tk,reqs_q,reqs_t, &
                   ww1,ww2,we1,we2,ws1,ws2,wn1,wn2,                  &
                   pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                  &
                   vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,                  &
                   uw31,uw32,ue31,ue32,us31,us32,un31,un32,          &
                   vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,          &
                   ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,          &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,          &
                   pw31,pw32,pe31,pe32,ps31,ps32,pn31,pn32,          &
                   tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,          &
                   qw1,qw2,qe1,qe2,qs1,qs2,qn1,qn2,                  &
                   tw1,tw2,te1,te2,ts1,ts2,tn1,tn2,ploc,packet)
      use module_mp_thompson
      use module_mp_morr_two_moment
      use module_sf_sfclay
      use module_bl_ysu
      use module_sf_slab
      use module_sf_oml
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'radcst.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif

!-----------------------------------------------------------------------
! Arrays and variables passed into solve

      integer, intent(in) :: nstep
      integer, intent(inout) :: nrec,prec,nwrite,nrst
      integer, intent(in) :: rbufsz,num_soil_layers
      integer, intent(inout) :: ndt
      real, intent(inout) :: dt,dtlast
      real*8, intent(inout) :: mtime
      real*8, intent(inout) :: stattim,taptim,rsttim,radtim,adt,acfl
      logical, intent(in) :: dodrag,dosfcflx
      logical, intent(in), dimension(maxq) :: cloudvar,rhovar
      character*3, intent(in), dimension(maxq) :: qname
      character*6, intent(in), dimension(maxq) :: budname
      real*8, intent(inout), dimension(nk) :: bud
      real*8, intent(inout), dimension(nj) :: bud2
      real*8, intent(inout), dimension(nbudget) :: qbudget
      real*8, intent(inout), dimension(numq) :: asq,bsq
      real, intent(in), dimension(ib:ie) :: xh,rxh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,uf,ruf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf,rvf
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: tauh,taus,zh,mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: tauf,zf,mf,rmf
      real, intent(inout), dimension(stat_out) :: rstat
      real, intent(in), dimension(ib:ie,jb:je) :: rho0s,pi0s,prs0s,rth0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,prs0,thv0,th0,qv0,qc0
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: ql0,rr0,rf0,rrf0,rru0,rrv0
      real, intent(in), dimension(itb:ite,jtb:jte) :: zs,gz,dzdx,dzdy
      real, intent(in), dimension(itb:ite+1,jtb:jte,ktb:kte) :: gx
      real, intent(in), dimension(itb:ite,jtb:jte+1,ktb:kte) :: gy
      real, intent(inout), dimension(ib:ie,jb:je,nrain) :: rain,sws,svs,sps,srs,sgs,sus,shs
      logical, intent(inout), dimension(ib:ie,jb:je) :: doimpl
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,thflux,qvflux,cdu,cdv,ce,u1,v1,w1
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,divx,rho,prs
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: rru,ua,u3d,uten,uten1
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: rrv,va,v3d,vten,vten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,wa,w3d,wten,wten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppten,sten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tha,th3d,thten,thten1,thterm
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: tk
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa,q3d,qten
      real, intent(inout), dimension(ibzvd:iezvd,jbzvd:jezvd,kbzvd:kezvd,nqzvd) :: zvdarray
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d,tketen
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dissten
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten
      real, intent(in), dimension(ibr:ier,jbr:jer,kbr:ker) :: o30
      real, intent(inout), dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin
      real, intent(inout), dimension(ni,nj,nrad2d) :: rad2d
      real, intent(inout), dimension(ni+1) :: x
      real, intent(inout), dimension(nj+1) :: y
      real, intent(inout), dimension(nk+3) :: z,za
      real, intent(inout), dimension(ni,nj,nk+3) :: zp
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: lu_index
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: kpbl2d
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: psfc,u10,v10,hfx,qfx,xland,znt,ust, &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,       &
                                      CK,CKA,CD,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      f2d,gsw,glw,chklowq,capg,snowc,dsxy
      real, intent(in), dimension(num_soil_layers) :: slab_zs,slab_dzs
      real, intent(inout), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta,pt3d,ptten
      real, intent(inout), dimension(npvals,nparcels) :: pdata
      real, intent(in), dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: cfb
      real, intent(in), dimension(kpb:kpe) :: cfa,cfc,ad1,ad2
      complex, intent(inout), dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: pdt,deft
      complex, intent(inout), dimension(ipb:ipe,jpb:jpe) :: rhs,trans
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_tk
      integer, intent(inout), dimension(rmp,numq) :: reqs_q
      integer, intent(inout), dimension(rmp,npt) :: reqs_t
      real, intent(inout), dimension(jmp,kmp-1) :: ww1,ww2,we1,we2
      real, intent(inout), dimension(imp,kmp-1) :: ws1,ws2,wn1,wn2
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(jmp,kmp) :: vw1,vw2,ve1,ve2
      real, intent(inout), dimension(imp,kmp) :: vs1,vs2,vn1,vn2
      real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, intent(inout), dimension(imp,cmp,kmp)   :: vs31,vs32,vn31,vn32
      real, intent(inout), dimension(cmp,jmp,kmp-1) :: ww31,ww32,we31,we32
      real, intent(inout), dimension(imp,cmp,kmp-1) :: ws31,ws32,wn31,wn32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: pw31,pw32,pe31,pe32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ps31,ps32,pn31,pn32
      real, intent(inout), dimension(cmp,jmp,kmt)   :: tkw1,tkw2,tke1,tke2
      real, intent(inout), dimension(imp,cmp,kmt)   :: tks1,tks2,tkn1,tkn2
      real, intent(inout), dimension(cmp,jmp,kmp,numq) :: qw1,qw2,qe1,qe2
      real, intent(inout), dimension(imp,cmp,kmp,numq) :: qs1,qs2,qn1,qn2
      real, intent(inout), dimension(cmp,jmp,kmp,npt) :: tw1,tw2,te1,te2
      real, intent(inout), dimension(imp,cmp,kmp,npt) :: ts1,ts2,tn1,tn2
      real, intent(inout), dimension(3,nparcels) :: ploc
      real, intent(inout), dimension(npvals+1,nparcels) :: packet

!-----------------------------------------------------------------------
! Arrays and variables defined inside solve

      integer i,j,k,n,nrk,bflag,pdef,nn,fnum
      real :: delqv,delpi,delth,delt,fac
      real :: foo1,foo2

      logical :: simple_comm,qcom,tcom,dorad

      real :: tout,cfl_limit,max_change,dtsm

!      real dttmp,rtime,rdt,tem,thrad,prad,ql
      real dttmp,rtime,rdt,tem,thrad,prad,ql,thrad_drc,T_drc  !DRC 06-02-11
      real T_tpp
      real :: cpm,cvm
      real*8 afoo,bfoo
      logical :: getdbz

      logical :: doirrad,dosorad
      real :: saltitude,sazimuth,zen
      real, dimension(2) :: x1
      real, dimension(2) :: y1
      real, dimension(rbufsz) :: radbuf
      real, dimension(nkr) :: swtmp,lwtmp
      real, dimension(nkr) :: tem1,tem2,tem3,tem4,tem5,   &
                              tem6,tem7,tem8,tem9,tem10,   &
                              tem11,tem12,tem13,tem14,tem15,   &
                              tem16,tem17
      real, dimension(nkr) :: ptprt,pprt,qv,qc,qr,qi,qs,qh,   &
                              ptbar,pbar,appi,rhostr,zpp,o31

      real :: rad2deg,albedo,albedoz,tema,temb,frac_snowcover
      real :: dtsfc0,dtsfc

      logical :: flag_qi
      integer :: isfflx
      real :: ep1,ep2,karman,rovg
      real :: SVP1,SVP2,SVP3,SVPT0,p1000mb,eomeg,stbolt
      integer :: ifsnow
      real :: dtmin

#ifdef MPI
      ! for mpi comm of cda:
      integer, dimension(8) :: reqs
      real, dimension(3,nj) :: west,newwest,east,neweast
      real, dimension(ni,3) :: south,newsouth,north,newnorth
#endif

!--------------------------------------------------------------------

#ifdef MPI
      nf=0
      nu=0
      nv=0
      nw=0
#endif

      afoo=0.0d0
      bfoo=0.0d0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Adaptive timestepping:
!   (assumes cflmax has already been calculated)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF( (adapt_dt.eq.1) .and. (myid.eq.0) )THEN
        ! only processor 0 does this:

        cfl_limit  = 1.00    ! maximum CFL allowed  (actually a "target" value)
        max_change = 0.25    ! maximum (percentage) change in timestep

        IF( cflmax.gt.cfl_limit )THEN
          ! decrease timestep timestep:
          dttmp = max( 1.0-max_change , cfl_limit/cflmax )*dt
        ELSE
          ! increase timestep:
          dttmp = min( 1.0+max_change , cfl_limit/cflmax )*dt
        ENDIF

        ! don't allow dt to exceed twice initial timestep
        dttmp = min( dttmp , 2.0*dtl )

        ! ramp-down timestep when approaching output time
      IF( taptim.gt.0.0 )THEN
        tout = sngl( taptim - mtime )
        if( tout.gt.(2.0*dttmp) .and. tout.le.(3.0*dttmp)  )then
          dttmp = 0.33333333*tout
        elseif( tout.gt.dttmp .and. tout.le.(2.0*dttmp)  )then
          dttmp = 0.5*tout
        elseif( tout.le.dttmp )then
          dttmp = tout
        endif
      ENDIF

      IF( rsttim.gt.0.0 )THEN
        ! ramp-down timestep when approaching restart time
        tout = sngl( rsttim - mtime )
        if( tout.gt.(2.0*dttmp) .and. tout.le.(3.0*dttmp)  )then
          dttmp = 0.33333333*tout
        elseif( tout.gt.dttmp .and. tout.le.(2.0*dttmp)  )then
          dttmp = 0.5*tout
        elseif( tout.le.dttmp )then
          dttmp = tout
        endif
      ENDIF

      IF( stattim.gt.0.0 )THEN
        ! ramp-down timestep when approaching stat time
        tout = sngl( stattim - mtime )
        if( tout.gt.(2.0*dttmp) .and. tout.le.(3.0*dttmp)  )then
          dttmp = 0.33333333*tout
        elseif( tout.gt.dttmp .and. tout.le.(2.0*dttmp)  )then
          dttmp = 0.5*tout
        elseif( tout.le.dttmp )then
          dttmp = tout
        endif
      ENDIF

        dt = dttmp

        ! Algorithm to determine number of small steps:
        IF( psolver.eq.2 )THEN
          ! check dx,dy,dz:
          IF( ny.eq.1 )THEN
            ! 2D sims (x-z):
            dtsm = 0.50*min( min_dx , min_dz )/350.0
          ELSEIF( nx.eq.1 )THEN
            ! 2D sims (y-z):
            dtsm = 0.50*min( min_dy , min_dz )/350.0
          ELSE
            ! 3D sims:
            dtsm = 0.50*min( min_dx , min_dy , min_dz )/350.0
          ENDIF
        ELSEIF( psolver.eq.3 )THEN
          ! check dx,dy:
          IF( ny.eq.1 )THEN
            ! 2D sims (x-z):
            dtsm = 0.60*min_dx/350.0
          ELSEIF( nx.eq.1 )THEN
            ! 2D sims (y-z):
            dtsm = 0.60*min_dy/350.0
          ELSE
            ! 3D sims:
            dtsm = 0.60*min( min_dx , min_dy )/350.0
          ENDIF
        ENDIF
        nsound = max( nint( dt/dtsm ) , 4 )
        if( mod(nsound,2).ne.0 ) nsound = nsound + 1
        if( dt/float(nsound).gt.dtsm ) nsound = nsound + 2

        if( nsound.gt.24 )then
          nsound = 24
          dt = 24*dtsm
        endif

        print *,'cflmax,dt,nsound:',cflmax,dt,nsound

        ! end of processor 0 stuff
      ENDIF

      IF( adapt_dt.eq.1 )THEN
        ! all processors:
#ifdef MPI
        call MPI_BCAST(dt    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nsound,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
        ndt = ndt + 1
        adt = adt + dt
        acfl = acfl + cflmax
        if(timestats.ge.1) time_misc=time_misc+mytime()
        IF( dt.ne.dtlast )THEN
          IF( (imoist.eq.1).and.(ptype.eq.2) )then
            call consat2(dt)
            if(timestats.ge.1) time_microphy=time_microphy+mytime()
          ENDIF
          IF( (imoist.eq.1).and.(ptype.eq.4) )then
            call lfoice_init(dt)
            if(timestats.ge.1) time_microphy=time_microphy+mytime()
          ENDIF
          dtlast = dt
        ENDIF
      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   radiation  ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF( radopt.eq.1 )THEN

        ! time at end of timestep:
        rtime=sngl(mtime+dt)
        dorad = .false.
        IF( rtime.ge.sngl(radtim) ) dorad = .true.
        dtrad = max( dtrad , dt )

        IF( dorad )THEN
!$omp parallel do default(shared)  &
!$omp private(i)
          do i=1,ni+1
            x(i)=xf(i)
          enddo
!$omp parallel do default(shared)  &
!$omp private(j)
          do j=1,nj+1
            y(j)=yf(j)
          enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk+3
          do j=1,nj
          do i=1,ni
            zp(i,j,k) =   zf(i,j,k-1)
          enddo
          enddo
          enddo
          i = 1
          j = 1
!$omp parallel do default(shared)  &
!$omp private(k)
          do k=1,nk+3
            z(k)=zf(1,1,k-1)
            za(k)=zh(1,1,min(k-1,ke))
          enddo
          rtime=sngl(mtime+dt)
          write(outfile,*) '  Calculating radiation tendency:'
          if(timestats.ge.1) time_rad=time_rad+mytime()
          call bcs(prs)
          CALL zenangl( ni,nj, x,y, zp(1,1,2),    &
                rad2d(1,1,ncosz), rad2d(1,1,ncosss), radsw,              &
                dum1(1,1,1),dum1(1,1,2),dum1(1,1,3),dum1(1,1,4),        &
                dum2(1,1,1),dum2(1,1,2),dum2(1,1,3),dum2(1,1,4),        &
                saltitude,sazimuth,dx,dy,dt,rtime,                     &
                ctrlat,ctrlon,year,month,day,hour,minute,second,jday )
          if(myid.eq.0)then
            print *,'    solar zenith angle  (degrees) = ',   &
                                   acos(rad2d(ni,nj,ncosz))*degdpi
            print *,'    solar azimuth angle (degrees) = ',sazimuth*degdpi
          endif
!-----------------------------------------------------------------------
!
!  Calculate surface albedo which is dependent on solar zenith angle
!  and soil moisture. Set the albedo for different types of solar
!  flux to be same.
!
!    rsirbm   Solar IR surface albedo for beam radiation
!    rsirdf   Solar IR surface albedo for diffuse radiation
!    rsuvbm   Solar UV surface albedo for beam radiation
!    rsuvdf   Solar UV surface albedo for diffuse radiation
!
!-----------------------------------------------------------------------
!
  rad2deg = 180.0/3.141592654

!$omp parallel do default(shared)  &
!$omp private(i,j,albedo,albedoz,frac_snowcover,tema)
  DO j=1,nj
    DO i=1,ni

      ! let's just use MM5/WRF value, instead:
      albedo = albd(i,j)

      ! arps code for albedo:
      ! (not sure I trust this.....)

!      albedoz = 0.01 * ( EXP( 0.003286         & ! zenith dependent albedo
!          * SQRT( ( ACOS(rad2d(i,j,ncosz))*rad2deg ) ** 3 ) ) - 1.0 )
!
!      IF ( soilmodel == 0 ) THEN             ! soil type not defined
!!!!        stop 12321
!        tema = 0
!      ELSE
!        tema = qsoil(i,j,1)/wsat(soiltyp(i,j))
!      END IF
!
!      frac_snowcover = MIN(snowdpth(i,j)/snowdepth_crit, 1.0)
!
!      IF ( tema > 0.5 ) THEN
!        albedo = albedoz + (1.-frac_snowcover)*0.14                     &
!                         + frac_snowcover*snow_albedo
!      ELSE
!        albedo = albedoz + (1.-frac_snowcover)*(0.31 - 0.34 * tema)     &
!                         + frac_snowcover*snow_albedo
!      END IF
!        albedo = albedoz

      rad2d(i,j,nrsirbm) = albedo
      rad2d(i,j,nrsirdf) = albedo
      rad2d(i,j,nrsuvbm) = albedo
      rad2d(i,j,nrsuvdf) = albedo

    END DO
  END DO
          ! big OpenMP parallelization loop:
!$omp parallel do default(shared)  &
!$omp private(i,j,k,ptprt,pprt,qv,qc,qr,qi,qs,qh,appi,o31,                 &
!$omp tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,tem10,        &
!$omp tem11,tem12,tem13,tem14,tem15,tem16,tem17,radbuf,swtmp,lwtmp,   &
!$omp doirrad,dosorad,z,za,zpp,ptbar,pbar,rhostr,x1,y1)
        do j=1,nj
        do i=1,ni
          swtmp = 0.0
          lwtmp = 0.0
          do k=1,nk+2
            ptprt(k) =  tha(i,j,k-1)
             pprt(k) =  prs(i,j,k-1) - prs0(i,j,k-1)
               qv(k) =   qa(i,j,k-1,nqv)
               qc(k) =   qa(i,j,k-1,nqc)
               qr(k) =   qa(i,j,k-1,nqr)
               qi(k) =   qa(i,j,k-1,nqi)
               qs(k) =   qa(i,j,k-1,nqs)
               qh(k) =   qa(i,j,k-1,nqg)
             appi(k) =  pi0(i,j,k-1) + ppi(i,j,k-1)
              o31(k) =  o30(i,j,k-1)
          enddo
          ptprt(1) = ptprt(2)
           pprt(1) =  pprt(2)
          ptprt(nk+2) = ptprt(nk+1)
           pprt(nk+2) =  pprt(nk+1)
          x1(1) = xf(i)
          x1(2) = xf(i+1)
          y1(1) = yf(j)
          y1(2) = yf(j+1)
          do k=1,nk+3
            z(k)=zf(i,j,k-1)
            za(k)=zh(i,j,min(k-1,ke))
            zpp(k) =   zp(i,j,k)
          enddo
          do k=2,nk+2
            ptbar(k) =  th0(i,j,k-1)
             pbar(k) = prs0(i,j,k-1)
           rhostr(k) = rho0(i,j,k-1)
          enddo
            ptbar(1) = rth0s(i,j)**(-1)
             pbar(1) = prs0s(i,j)
           rhostr(1) = rho0s(i,j)
            doirrad = .true.
            dosorad = .true.
          CALL radtrns(nir,njr,nkr, rbufsz, 0,myid,dx,dy,            &
                 ib,ie,jb,je,kb,ke,xh,yh,prs0s(i,j),                  &
                 ptprt,pprt,qv,qc,qr,qi,qs,qh,                          &
                 ptbar,pbar,appi,o31,rhostr, tsk(i,j), zpp ,                                 &
                 radsw(i,j),rnflx(i,j),radswnet(i,j),radlwin(i,j), rad2d(i,j,ncosss),            &
                 rad2d(i,j,nrsirbm),rad2d(i,j,nrsirdf),rad2d(i,j,nrsuvbm),                       &
                 rad2d(i,j,nrsuvdf), rad2d(i,j,ncosz),sazimuth,                                  &
                 rad2d(i,j,nfdirir),rad2d(i,j,nfdifir),rad2d(i,j,nfdirpar),rad2d(i,j,nfdifpar),  &
                 tem1, tem2, tem3, tem4, tem5,                &
                 tem6, tem7, tem8, tem9, tem10,               &
                 tem11,tem12,tem13,tem14,tem15,tem16,  &
                 radbuf(1), tem17,swtmp,lwtmp,doirrad,dosorad)
          do k=1,nk
            swten(i,j,k) = swtmp(k+1)
            lwten(i,j,k) = lwtmp(k+1)
          enddo
        enddo
        enddo
          radtim=radtim+dtrad
        ENDIF
        if(timestats.ge.1) time_rad=time_rad+mytime()

      ENDIF


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   surface  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!-------------------------------------------------------------------
!  prepare some arrays for WRF surface/pbl physics:

      ! between here and call to ysu:
      ! DO NOT CHANGE:  dum1,dum2,dum3,dum4,sten,t11,t23

      IF((oceanmodel.eq.2).or.(ipbl.eq.1).or.(sfcmodel.eq.2))THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=0.5*(ua(i,j,k)+ua(i+1,j,k))
          dum2(i,j,k)=0.5*(va(i,j,k)+va(i,j+1,k))
          dum3(i,j,k)=th0(i,j,k)+tha(i,j,k)
          sten(i,j,k)=pi0(i,j,k)+ppi(i,j,k)
          dum4(i,j,k)=dum3(i,j,k)*sten(i,j,k)
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          t11(i,j,k) = dz*rmf(i,j,k)
          t23(i,j,k) = prs(i,j,k-1)+(prs(i,j,k)-prs(i,j,k-1))   &
                                   *( zf(i,j,k)- zh(i,j,k-1))   &
                                   /( zh(i,j,k)- zh(i,j,k-1))
        enddo
        enddo
        enddo

        ! values at surface, top of model:
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          t23(i,j,1) =  prs(i,j,1)-zh(i,j,1)*( prs(i,j,2)- prs(i,j,1))   &
                                            /(  zh(i,j,2)-  zh(i,j,1))
          t23(i,j,nk+1)= prs(i,j,nk)+(zf(i,j,nk+1)-zh(i,j,nk))       &
                                    *( prs(i,j,nk)- prs(i,j,nk-1))   &
                                    /(  zh(i,j,nk)-  zh(i,j,nk-1))
          psfc(i,j) = t23(i,j,1)
        enddo
        enddo

        ep1 = rv/rd - 1.0
        ep2 = rd/rv
        karman = 0.4
        rovg = rd/g

        ! dum1 = u at scalars
        ! dum2 = v at scalars
        ! dum3 = th
        ! dum4 = t
        ! sten = pi
        ! t11 = dz8w
        ! t12 = qvten
        ! t13 = qcten
        ! t22 = qiten
        ! t23 = p3di
        ! t33 = exch_h
        ! divx = uten
        ! thterm = vten

        isfflx = 1
        SVP1=0.6112
        SVP2=17.67
        SVP3=29.65
        SVPT0=273.15
        p1000mb      = 100000.
        EOMEG=7.2921E-5
        STBOLT=5.67051E-8

        IF(radopt.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            gsw(i,j)=radsw(i,j)
            glw(i,j)=radlwin(i,j)
          enddo
          enddo
        ELSE
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            gsw(i,j)=0.0
            glw(i,j)=0.0
          enddo
          enddo
        ENDIF

        if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()
      ENDIF

      IF( sfcmodel.ge.1 .and. ipbl.eq.0 )THEN

        call gethpbl(zh,th0,tha,qa,hpbl)

      ENDIF

!-------------------------------------------------------------------
! surface schemes:

!---------------------------------------------------------------------------
! original CM1 formulation:

    IF(sfcmodel.eq.1)THEN

      if(isfcflx.eq.1.or.idrag.eq.1)then
        call getcecd(cdu,cdv,ce,u0,v0,rf0,u1,v1,w1,ua,va)
      endif

      ! get surface flux
      if(isfcflx.eq.1)then
        call sfcflux(dt,ruh,xf,rvh,pi0s,ce,zh,pi0,thv0,th0,u0,v0,tsk,thflux,qvflux,mavail, &
                     rho,u1,v1,w1,ua,va,ppi,tha,qa(ibm,jbm,kbm,nqv), &
                     qbudget(8))
      endif

      call sfcdiags(tsk,thflux,qvflux,cdu,cdv,ce,u1,v1,w1,   &
                    xland,psfc,qsfc,u10,v10,hfx,qfx,cda,znt,ust,gz1oz0,   &
                    psim,psih,br,zol,mol,hpbl,wspd,dsxy,th2,t2,q2,      &
                    zs,zh,pi0s,pi0,th0,ppi,tha,rho,qa,ua,va)

!-------------------------------------------------------------------

    ELSEIF(sfcmodel.eq.2)THEN

      ! surface layer:
      ! (needed by sfcmodel=2 and ipbl=1)
      call SFCLAY(dum1,dum2,dum4,qa(ib,jb,kb,nqv),prs,t11,       &
                   CP,G,ROVCP,RD,XLV,PSFC,CHS,CHS2,CQS2,CPMM,    &
                   ZNT,UST,hpbl,MAVAIL,ZOL,MOL,REGIME,PSIM,PSIH, &
                   XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,QGH,QSFC,RMOL, &
                   U10,V10,TH2,T2,Q2,                            &
                   GZ1OZ0,WSPD,BR,ISFFLX,dsxy,                   &
                   SVP1,SVP2,SVP3,SVPT0,EP1,EP2,                 &
                   KARMAN,EOMEG,STBOLT,                          &
                   P1000mb,                                      &
                   1  ,ni+1 , 1  ,nj+1 , 1  ,nk+1 ,                    &
                   ib ,ie , jb ,je , kb ,ke ,                    &
                   1  ,ni , 1  ,nj , 1  ,nk ,                    &
                   ustm,ck,cka,cd,cda,isftcflx,iz0tlnd           )

      ifsnow = 0
      dtmin = dt/60.0

      ! slab scheme (MM5/WRF):
      call SLAB(dum4,qa(ib,jb,kb,nqv),prs,FLHC,FLQC,                      &
                   PSFC,XLAND,TMN,HFX,QFX,LH,TSK,QSFC,CHKLOWQ,  &
                   GSW,GLW,CAPG,THC,SNOWC,EMISS,MAVAIL,         &
                   DT,ROVCP,XLV,DTMIN,IFSNOW,               &
                   SVP1,SVP2,SVP3,SVPT0,EP2,                    &
                   KARMAN,EOMEG,STBOLT,                         &
                   TSLB,slab_ZS,slab_DZS,num_soil_layers, .true. ,       &
                   P1000mb,                                     &
                     1, ni+1,   1, nj+1,   1, nk+1,             &
                    ib, ie,  jb, je,  kb, ke,                   &
                     1, ni,   1, nj,   1, nk                    )

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()
    ENDIF

!-------------------------------------------------------------------
! simple ocean mixed layer model based Pollard, Rhines and Thompson (1973)
!   (from WRF)

    IF(oceanmodel.eq.2)THEN

        CALL oceanml(tml,t0ml,hml,h0ml,huml,hvml,ust,dum1,dum2, &
                     tmoml,f2d,g,oml_gamma,                     &
                     xland,hfx,lh,tsk,gsw,glw,emiss,            &
                     dt,STBOLT,                                 &
                       1, ni+1,   1, nj+1,   1, nk+1,           &
                      ib, ie,  jb, je,  kb, ke,                 &
                       1, ni,   1, nj,   1, nk                  )

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()
    ENDIF

!-------------------------------------------------------------------
!  PBL scheme:

      IF(ipbl.eq.1)THEN

        divx = 0.0
        thterm = 0.0
        thten = 0.0
        t12 = 0.0
        t13 = 0.0
        t22 = 0.0

        if( iice.eq.1 )then
          flag_qi = .true.
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,nqv) = qa(i,j,k,nqi)
          enddo
          enddo
          enddo
        else
          flag_qi = .false.
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,nqv) = 0.0
          enddo
          enddo
          enddo
        endif

        IF(output_km.eq.1.or.output_kh.eq.1)THEN
          ! ppten = exch_m
          t33=0.0
          ppten=0.0
        ENDIF

        ! PBL:
        call ysu(dum1,dum2,dum3,dum4,qa(ib,jb,kb,nqv),         &
                  qa(ib,jb,kb,nqc),qten(ib,jb,kb,nqv),prs,t23,sten,  &
                  divx,thterm,thten,                           &
                  t12,t13,t22,flag_qi,                         &
                  cp,g,rovcp,rd,rovg,ep1,ep2,karman,xlv,rv,    &
                  t11 ,psfc,                                   &
!!!                  znu,znw,mut,p_top,                        &
                  znt,ust,hpbl,psim,psih,                      &
                  xland,hfx,qfx,gz1oz0,wspd,br,                &
                  dt,kpbl2d,                                   &
                  t33,ppten,                                   &
                  u10,v10,                                     &
                  1  ,ni+1 , 1  ,nj+1 , 1  ,nk+1 ,             &
                  ib ,ie , jb ,je , kb ,ke ,                   &
                  1  ,ni , 1  ,nj , 1  ,nk ,                   &
                  regime                                       )
        if(timestats.ge.1) time_pbl=time_pbl+mytime()

        call bcs(divx)
#ifdef MPI
        call comm_1s_start(divx,pw1,pw2,pe1,pe2,   &
                                ps1,ps2,pn1,pn2,reqs_s)
#endif
        call bcs(thterm)
#ifdef MPI
        call comm_1s_start(thterm,vw1,vw2,ve1,ve2,   &
                                  vs1,vs2,vn1,vn2,reqs_p)
#endif
        IF(output_km.eq.1.or.output_kh.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            khv(i,j,k) = t33(i,j,k)
            kmv(i,j,k) = ppten(i,j,k)
          enddo
          enddo
          enddo
        ENDIF
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          thpten(i,j,k) = thten(i,j,k)
          qvpten(i,j,k) =   t12(i,j,k)
          qcpten(i,j,k) =   t13(i,j,k)
          qipten(i,j,k) =   t22(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
#ifdef MPI
        call comm_1s_end(divx,pw1,pw2,pe1,pe2,   &
                              ps1,ps2,pn1,pn2,reqs_s)
#endif
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=0,ni+1
           upten(i,j,k) =  divx(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
#ifdef MPI
        call comm_1s_end(thterm,vw1,vw2,ve1,ve2,   &
                                vs1,vs2,vn1,vn2,reqs_p)
#endif
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=1,ni
           vpten(i,j,k) =thterm(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()

      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   subgrid turbulence schemes  cccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Misc prep:
!  Also, set surface stresses:

      IF( sfcmodel.ge.2 )THEN
        ! put WRF parameters into CM1 arrays:

        IF( dosfcflx )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            thflux(i,j) = hfx(i,j)/(cp*rho(i,j,1))
            qvflux(i,j) = qfx(i,j)/rho(i,j,1)
          enddo
          enddo
        ENDIF
        IF( dodrag )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=0,nj+1
          do i=0,ni+1
            u1(i,j) = 0.5*(ua(i,j,1)+ua(i+1,j,1))
            v1(i,j) = 0.5*(va(i,j,1)+va(i,j+1,1))
            w1(i,j) = sqrt( u1(i,j)**2 + v1(i,j)**2 )
            ce(i,j) = cka(i,j)
          enddo
          enddo
          call bc2d(cda)
          call bc2d(ust)
#ifdef MPI
          call comm_2d_start(cda,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
          call comm_2dew_end(cda,west,newwest,east,neweast,reqs)
          call comm_2dns_end(cda,south,newsouth,north,newnorth,reqs)
          call comm_2d_start(ust,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
          call comm_2dew_end(ust,west,newwest,east,neweast,reqs)
          call comm_2dns_end(ust,south,newsouth,north,newnorth,reqs)
#endif
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj+1
          do i=1,ni+1
            cdu(i,j) = 0.5*(cda(i-1,j)+cda(i,j))
            cdv(i,j) = 0.5*(cda(i,j-1)+cda(i,j))
            t13(i,j,1) = ((0.5*(ust(i-1,j)+ust(i,j)))**2)*ua(i,j,1)/max(0.5*(w1(i-1,j)+w1(i,j)),0.1)
            t23(i,j,1) = ((0.5*(ust(i,j-1)+ust(i,j)))**2)*va(i,j,1)/max(0.5*(w1(i,j-1)+w1(i,j)),0.1)
          enddo
          enddo
        ENDIF

      ENDIF

      IF(sfcmodel.eq.1)THEN
        ! get surface drag
        if(idrag.eq.1)then
          call sfcdrag(cdu,cdv,u0,v0,rf0,u1,v1,t13,t23,ua,va)
        endif
      ENDIF

!--------------------------------------------------------------------
!                 For turbulence section only:
!  dum1 = squared Brunt-Vaisala frequency (N_m^2) (nm)
!  dum2 = Vertical deformation terms (S_v^2) (defsq)
!  dum3 = Horizontal deformation terms (S_h^2) (defh)
!
!  Arrays available for temporary storage:
!  dum4,divx,ppten,sten,thterm

      IF(iturb.ge.1)THEN

        ! squared Brunt-Vaisala frequency:
        call calcnm(mf,pi0,thv0,th0,cloudvar,dum1,dum2,dum3,dum4,divx,   &
                    prs,ppi,tha,qa)

        ! deformation:
        call calcdef(dodrag,xh,rxh,uh,xf,rxf,uf,vh,vf,mh,mf,dum2,dum3,   &
                     divx,ppten,ua,va,wa,t11,t12,t13,t22,t23,t33,gx,gy)

      ENDIF

!--------------------------------------------------------------------
!  iturb=1:  tke scheme  (for large eddy simulation)
!    Reference:  Deardorff, 1980, Bound Layer Meteor, p. 495
!                see also Stevens, Moeng, Sullivan, 1999, JAS, p. 3963

      IF(iturb.eq.1)THEN
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nkt
        do j=1,nj
        do i=1,ni
          tketen(i,j,k)=0.0
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()

        call turbtke(dt,dodrag,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,   &
                     dum1,dum2,dum3,dum4,divx,ppten,sten,thterm,   &
                     kmh,kmv,khh,khv,tkea,tketen,t13,t23,ua,va)

!-------------------------------------------------
!  iturb=2:  Smagorinsky scheme  (for large eddy simulation)
!    Reference:  see, e.g., Stevens, Moeng, Sullivan, 1999, JAS, p. 3963

      ELSEIF(iturb.eq.2)THEN

        call turbsmag(dt,dodrag,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,  &
                      dum1,dum2,dum3,dum4,divx,sten,               &
                      kmh,kmv,khh,khv,t13,t23,ua,va)

!-------------------------------------------------
!  iturb=3:  parameterized turbulence  (no explicit turbulence)
!    Reference:  Rotunno and Emanuel, 1987, JAS, p. 542
!                Bryan and Rotunno, 2009, MWR, p. 1770

      ELSEIF(iturb.eq.3)THEN

        IF( l_v.le.tsmall .and. ipbl.eq.1 )THEN
          ! save kmv,khv:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            sten(i,j,k)=kmv(i,j,k)
            ppten(i,j,k)=khv(i,j,k)
          enddo
          enddo
          enddo
        ENDIF

        call turbparam(nstep,zf,dt,dodrag,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,  &
                         dum1,dum2,dum3,kmh,kmv,khh,khv,t13,t23,ua,va)

        IF( l_v.le.tsmall .and. ipbl.eq.1 )THEN
          ! restore kmv,khv:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            kmv(i,j,k)=sten(i,j,k)
            khv(i,j,k)=ppten(i,j,k)
          enddo
          enddo
          enddo
!!!          ! diagnostic:  effective khv
!!!          !   dum1 is theta flux:
!!!          do j=1,nj
!!!          do i=1,ni
!!!            dum1(i,j,nk+1) = 0.0
!!!            do k=nk,1,-1
!!!              dum1(i,j,k) = dum1(i,j,k+1)+thpten(i,j,k)*rho0(i,j,k)*dz*rmh(i,j,k)
!!!            enddo
!!!            do k=1,nk+1
!!!              khv(i,j,k) = -dum1(i,j,k)/((th0(i,j,k)+tha(i,j,k))-(th0(i,j,k-1)+tha(i,j,k-1)))*rdz*mf(i,j,k)*rf0(i,j,k)
!!!            enddo
!!!          enddo
!!!          enddo
        ENDIF

!-------------------------------------------------

      ELSEIF(iturb.ne.0)THEN

        print *,'  unknown turbulence setting ... '
        call stopcm1

      ENDIF

!-------------------------------------------------
!  check for columns that need vertically implicit diffusion:

      IF(iturb.ge.1)THEN

        tem = 0.125*dz*dz/dt

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=0,nj+1
        do i=0,ni+1
          doimpl(i,j) = .false.
          k = 2
          do while( ( .not. doimpl(i,j) ) .and. (k.le.nk) )
            if( khv(i,j,k) .gt. tem*rmf(i,j,k)*rmf(i,j,k) )then
              doimpl(i,j) = .true.
            endif
            k = k + 1
          enddo
        enddo
        enddo
        if(timestats.ge.1) time_turb=time_turb+mytime()

      ENDIF

!-------------------------------------------------
!  some more calculations for TKE scheme:

      IF(iturb.eq.1)THEN

        call turbt(dt,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho0,rr0,rf0,rrf0,   &
                   dum1,dum2,dum3,dum4,divx,tkea,tketen,kmh,kmv,gx,gy,doimpl)

        if(idiff.eq.1)then
          if(difforder.eq.2)then
            call diff2w(1,0,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,thterm,tkea,tketen)
          elseif(difforder.eq.6)then
            ! for diff6, use '3d' array
            call diff6w(dt,dum1,dum2,dum3,tke3d,tketen)
          endif
        endif

      ENDIF

!-------------------------------------------------
!  Get turbulent stresses:

      IF(iturb.ge.1)THEN

        call gettau(dodrag,xf,rxf,rho0,rf0,kmh,kmv,t11,t12,t13,t22,t23,t33,ua)

      ENDIF

!--------------------------------------------------------------------
!  Dissipative heating term:

      IF(idiss.eq.1)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dissten(i,j,k) = 0.0
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_turb=time_turb+mytime()
      ENDIF

      IF(iturb.ge.1.and.idiss.eq.1)THEN

        call getepsilon(rxh,uh,xf,rxf,uf,yh,vh,yf,vf,mh,mf,rr0,rrf0,   &
                        dum1,dum2,dum3,dum4,divx,ppten,sten,dissten,    &
                        t11,t12,t13,t22,t23,t33,ua,va,wa,gx,gy)

      ENDIF

      IF(ipbl.eq.1.and.idiss.eq.1)THEN
        ! Dissipative heating from ysu scheme:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          ! assume t13,t23 are zero at top of domain:
          t13(i,j,nk+1) = 0.0
          t23(i,j,nk+1) = 0.0
          do k=nk,1,-1
            t13(i,j,k) = t13(i,j,k+1)-upten(i,j,k)*rho0(i,j,k)*dz*rmh(i,j,k)
            t23(i,j,k) = t23(i,j,k+1)-vpten(i,j,k)*rho0(i,j,k)*dz*rmh(i,j,k)
          enddo
          do k=2,nk
            dum2(i,j,k)=0.5*((ua(i,j,k  )+ua(i+1,j,k  ))  &
                            -(ua(i,j,k-1)+ua(i+1,j,k-1)))*rdz*mf(i,j,k)
            dum3(i,j,k)=0.5*((va(i,j,k  )+va(i,j+1,k  ))  &
                            -(va(i,j,k-1)+va(i,j+1,k-1)))*rdz*mf(i,j,k)
          enddo
          dum2(i,j,1)=2.0*0.5*(ua(i,j,1)+ua(i+1,j,1))*rdz*mf(i,j,1)
          dum3(i,j,1)=2.0*0.5*(va(i,j,1)+va(i,j+1,1))*rdz*mf(i,j,1)
          dum2(i,j,nk+1)=0.0
          dum3(i,j,nk+1)=0.0
          do k=1,nk
            dissten(i,j,k)=dissten(i,j,k)+rr0(i,j,k)*0.5*(  &
                             ( t13(i,j,k  )*dum2(i,j,k  )   &
                              +t13(i,j,k+1)*dum2(i,j,k+1) ) &
                            +( t23(i,j,k  )*dum3(i,j,k  )   &
                              +t23(i,j,k+1)*dum3(i,j,k+1) ) )
          enddo
        enddo
        enddo
        t13 = 0.0
        t23 = 0.0

      ENDIF

!--------------------------------------------------------------------


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   Pre-RK calculations   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
!-------------------------------------------------------------------
!  Parcel update

      if(iprcl.eq.1)then
        ! rtime valid at beginning of time step
        rtime=sngl(mtime)
        call parcel_driver(prec,dt,xh,uh,ruh,yh,vh,rvh,zh,mh,rmh,mf,        &
                           pi0,thv0,th0,dum1,dum2,dum3,dum4,divx,prs,    &
                           ua,va,wa,ppi,thten,tha,qa,khv,pdata,rtime,    &
                           ploc,packet,reqs_p,                           &
                           pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2)
      endif

!--------------------------------------------------------------------
!  radbc
 
      if(irbc.eq.1)then

        if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,ua)
 
        if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,va)

      endif

!--------------------------------------------------------------------
!  U-equation
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
!!!        uten1(i,j,k)=0.
        uten1(i,j,k)=-rdalpha*0.5*(tauh(i-1,j,k)+tauh(i,j,k))*(ua(i,j,k)-u0(i,j,k))
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      if(idiff.ge.1)then
        if(difforder.eq.2)then
          call diff2u(1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,dissten,ua,uten1)
        elseif(difforder.eq.6)then
          ! for diff6, use '3d' array
          call diff6u(dt,u0,dum1,dum2,dum3,u3d,uten1)
        endif
      endif

      if(dns.eq.1)then
        call diff2u(2,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,dissten,ua,uten1)
      endif
 
      if(iturb.ge.1)then
        call turbu(dt,dodrag,xh,ruh,xf,rxf,uf,vh,mh,mf,rmf,rho0,rf0,rru0,   &
           dum1,dum2,dum3,dum4,ua,uten1,wa,t11,t12,t13,t22,kmv,gx,gy,doimpl)
      endif

      if(ipbl.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          uten1(i,j,k) = uten1(i,j,k) + 0.5*( upten(i-1,j,k)+ upten(i,j,k))
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
      endif

!--------------------------------------------------------------------
!  V-equation
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
!!!        vten1(i,j,k)=0.
        vten1(i,j,k)=-rdalpha*0.5*(tauh(i,j-1,k)+tauh(i,j,k))*(va(i,j,k)-v0(i,j,k))
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      if(idiff.ge.1)then
        if(difforder.eq.2)then
          call diff2v(1,xh,uh,rxf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,dissten,va,vten1)
        elseif(difforder.eq.6)then
          ! for diff6, use '3d' array
          call diff6v(dt,v0,dum1,dum2,dum3,v3d,vten1)
        endif
      endif

      if(dns.eq.1)then
        call diff2v(2,xh,uh,rxf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,dissten,va,vten1)
      endif
 
      if(iturb.ge.1)then
        call turbv(dt,dodrag,xh,rxh,uh,xf,rvh,vf,mh,mf,rho0,rf0,rrv0,   &
            dum1,dum2,dum3,dum4,va,vten1,wa,t12,t22,t23,kmv,gx,gy,doimpl)
      endif

      if(ipbl.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vten1(i,j,k) = vten1(i,j,k) + 0.5*( vpten(i,j-1,k)+ vpten(i,j,k))
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
      endif
 
!--------------------------------------------------------------------
!  W-equation
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk+1
      do j=1,nj
      do i=1,ni
!!!        wten1(i,j,k)=0.0
        wten1(i,j,k)=-rdalpha*tauf(i,j,k)*wa(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      if(idiff.ge.1)then
        if(difforder.eq.2)then
          call diff2w(1,1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,dissten,wa,wten1)
        elseif(difforder.eq.6)then
          ! for diff6, use '3d' array
          call diff6w(dt,dum1,dum2,dum3,w3d,wten1)
        endif
      endif

      if(dns.eq.1)then
        call diff2w(2,1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,dissten,wa,wten1)
      endif
 
      if(iturb.ge.1)then
        call turbw(dt,xh,rxh,uh,xf,vh,mh,mf,rho0,rf0,rrf0,   &
               dum1,dum2,dum3,dum4,wa,wten1,t13,t23,t33,kmv,gx,gy,doimpl)
      endif

!--------------------------------------------------------------------
!  THETA-equation
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
!!!        thten1(i,j,k)=0.0
        thten1(i,j,k)=-rdalpha*taus(i,j,k)*tha(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      if(idiff.eq.1)then
        if(difforder.eq.2)then
          call diff2s(1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,tha,thten1)
        elseif(difforder.eq.6)then
          ! for diff6, use '3d' array
          call diff6s(dt,ql0,dum1,dum2,dum3,th3d,thten1)
        endif
      endif

!----- cvm (if needed) -----!

      IF( neweqts.ge.1 .and. (idiss.eq.1.or.rterm.eq.1) )THEN
        ! store cvm in dum1:
        ! store ql  in dum2:
        ! store qi  in dum3:
        call getqli(0,qa,dum2,dum3)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=cv+cvv*qa(i,j,k,nqv)+cpl*dum2(i,j,k)+cpi*dum3(i,j,k)
        enddo
        enddo
        enddo
      ENDIF

!----- store appropriate rho for budget calculations in dum2 -----!

      IF(axisymm.eq.1)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = rho(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
        enddo
        enddo
        enddo
      ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = rho(i,j,k)
        enddo
        enddo
        enddo
      ENDIF

!---- Dissipative heating term:

      IF(idiss.eq.1)THEN
        ! use dissten array to store epsilon
        if(imoist.eq.1.and.neweqts.ge.1)then
          ! moist, new equations:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            thten1(i,j,k)=thten1(i,j,k)   &
                        +dissten(i,j,k)/( cpdcv*dum1(i,j,k)*(pi0(i,j,k)+ppi(i,j,k)) )
          enddo
          enddo
          enddo
        else
          ! traditional equations:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            thten1(i,j,k)=thten1(i,j,k)   &
                        +dissten(i,j,k)/( cp*(pi0(i,j,k)+ppi(i,j,k)) )
          enddo
          enddo
          enddo
        endif
      ENDIF

!---- Rotunno-Emanuel "radiation" term
!---- (currently capped at 2 K/day ... see RE87 p 546)

      IF(rterm.eq.1)THEN
!$omp parallel do default(shared)  &
!$omp private(k)
        do k=1,nk
          bud(k)=0.0d0
        enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k,thrad,prad)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ! NOTE:  thrad is a POTENTIAL TEMPERATURE tendency

        !DRC 05-23-11
!          thrad = -tha(i,j,k)/(12.0*3600.0)
!          if( tha(i,j,k).gt. 1.0 ) thrad = -1.0/(12.0*3600.0)
!          if( tha(i,j,k).lt.-1.0 ) thrad =  1.0/(12.0*3600.0)
        T_drc = (th0(i,j,k)+tha(i,j,k))*(prs(i,j,k)/p00)**(rd/cp) 
        T_tpp = 200
        if(T_drc.gt.T_tpp)then  !assume cooling occurs only where warm enough
                                !(strat should be isothermal at 200K!)
          thrad_drc = -0.5 / (12.0*3600.0)
          thrad = thrad_drc
        else
          thrad = -tha(i,j,k)/(80*12.0*3600.0)
        endif
        !DRC 05-23-11; 1 K/day in troposphere only; 40-day period relaxation to
        !reference in stratosphere
        !END DRC 05-23-11
          thten1(i,j,k)=thten1(i,j,k)+thrad
          ! associated pressure tendency:
          prad = (pi0(i,j,k)+ppi(i,j,k))*rddcv*thrad/(th0(i,j,k)+tha(i,j,k))
          ! budget:
          bud(k) = bud(k) + dum1(i,j,k)*dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*( &
                            thrad*(pi0(i,j,k)+ppi(i,j,k))    &
                           + prad*(th0(i,j,k)+tha(i,j,k)) )
        enddo
        enddo
        enddo
        tem = dt*dx*dy*dz
        do k=1,nk
          qbudget(10) = qbudget(10) + tem*bud(k)
        enddo
      ENDIF
      if(timestats.ge.1) time_misc=time_misc+mytime()

      if( (iturb.ge.1).or.(dns.eq.1) )then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
          sten(i,j,k)=th0(i,j,k)+tha(i,j,k)
        enddo
        enddo
        enddo
      endif

      if(dns.eq.1)then
        call diff2s(2,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,sten,thten1)
      endif

      IF( radopt.eq.1 )THEN
        ! tendency from radiation scheme:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          thten1(i,j,k)=thten1(i,j,k)+(swten(i,j,k)+lwten(i,j,k))
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_rad=time_rad+mytime()
      ENDIF

      IF( ipbl.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          thten1(i,j,k) = thten1(i,j,k) + thpten(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
      ENDIF

      if(iturb.ge.1)then
        call turbs(1,dt,dosfcflx,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho0,rr0,rf0,thflux,   &
                   dum1,dum2,dum3,dum4,divx,sten,thten1,khh,khv,gx,gy,doimpl)
      endif

!-------------------------------------------------------------------
!  contribution to pressure tendency from potential temperature:

      IF(neweqts.ge.1.or.imoist.eq.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ppten(i,j,k) = thten1(i,j,k)*rddcv   &
                        *(pi0(i,j,k)+ppi(i,j,k))/(th0(i,j,k)+tha(i,j,k))
        enddo
        enddo
        enddo
      ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ppten(i,j,k) = 0.0
        enddo
        enddo
        enddo
      ENDIF
      if(timestats.ge.1) time_misc=time_misc+mytime()

!-------------------------------------------------------------------
!  budget calculations:

      if(dosfcflx.and.imoist.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(j)
        do j=1,nj
          bud2(j) = 0.0d0
        enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k,delpi,delth,delqv,delt,n)
        do j=1,nj
        do i=1,ni
          k = 1
          delth = rf0(i,j,1)*rr0(i,j,1)*rdz*mh(i,j,1)*thflux(i,j)
          delqv = rf0(i,j,1)*rr0(i,j,1)*rdz*mh(i,j,1)*qvflux(i,j)
          delpi = rddcv*(pi0(i,j,1)+ppi(i,j,1))*(           &
                                delqv/(eps+qa(i,j,1,nqv))   &
                               +delth/(th0(i,j,1)+tha(i,j,1))  )
          delt = (pi0(i,j,k)+ppi(i,j,k))*delth   &
                +(th0(i,j,k)+tha(i,j,k))*delpi
          bud2(j) = bud2(j) + dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*(        &
                  cv*delt                                                   &
                + cvv*qa(i,j,k,nqv)*delt                                    &
                + cvv*(pi0(i,j,k)+ppi(i,j,k))*(th0(i,j,k)+tha(i,j,k))*delqv &
                + g*zh(i,j,k)*delqv   )
          do n=nql1,nql2
            bud2(j) = bud2(j) + dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*cpl*qa(i,j,k,n)*delt
          enddo
          if(iice.eq.1)then
            do n=nqs1,nqs2
              bud2(j) = bud2(j) + dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*cpi*qa(i,j,k,n)*delt
            enddo
          endif
        enddo
        enddo
        tem = dt*dx*dy*dz
        do j=1,nj
          qbudget(9) = qbudget(9) + tem*bud2(j)
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()
      endif

!-------------------------------------------------------------------
!  Passive Tracers

      if(iptra.eq.1)then
        do n=1,npt
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ptten(i,j,k,n)=0.0
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_misc=time_misc+mytime()
          if(idiff.eq.1)then
            if(difforder.eq.2)then
              call diff2s(1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,pta(ib,jb,kb,n),ptten(ib,jb,kb,n))
            elseif(difforder.eq.6)then
              ! for diff6, use '3d' array
              call diff6s(dt,ql0,dum1,dum2,dum3,pt3d(ib,jb,kb,n),ptten(ib,jb,kb,n))
            endif
          endif
          if(iturb.ge.1)then
            call turbs(0,dt,dosfcflx,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho0,rr0,rf0,qvflux,   &
                       dum1,dum2,dum3,dum4,divx,pta(ib,jb,kb,n),ptten(ib,jb,kb,n),khh,khv,gx,gy,doimpl)
          endif
        enddo
      endif

!-------------------------------------------------------------------
!  Moisture

      if(imoist.eq.1)then
        DO n=1,numq
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,n)=0.0
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_misc=time_misc+mytime()
!---------------------------
          ! qv:
          if(n.eq.nqv)then
            if(idiff.eq.1)then
              if(difforder.eq.2)then
                call diff2s(1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,qa(ib,jb,kb,n),qten(ib,jb,kb,n))
              elseif(difforder.eq.6)then
                ! for diff6, use '3d' array
                call diff6s(dt,qv0,dum1,dum2,dum3,q3d(ib,jb,kb,n),qten(ib,jb,kb,n))
              endif
            endif
            if(iturb.ge.1)then
              call turbs(1,dt,dosfcflx,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho0,rr0,rf0,qvflux,   &
                         dum1,dum2,dum3,dum4,divx,qa(ib,jb,kb,n),qten(ib,jb,kb,n),khh,khv,gx,gy,doimpl)
            endif
            if(ipbl.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
              do k=1,nk
              do j=1,nj
              do i=1,ni
                qten(i,j,k,nqv) = qten(i,j,k,nqv) + qvpten(i,j,k)
              enddo
              enddo
              enddo
              if(timestats.ge.1) time_pbl=time_pbl+mytime()
            endif
            IF(neweqts.ge.1)THEN
              ! contribution to pressure tendency from water vapor:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
              do k=1,nk
              do j=1,nj
              do i=1,ni
                ppten(i,j,k)=ppten(i,j,k)+qten(i,j,k,n)   &
                            *rddcv*(pi0(i,j,k)+ppi(i,j,k))/(eps+qa(i,j,k,n))
              enddo
              enddo
              enddo
              if(timestats.ge.1) time_misc=time_misc+mytime()
            ENDIF
!---------------------------
          ! not qv:
          else
            if(idiff.eq.1)then
              if(difforder.eq.2)then
                call diff2s(1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,qa(ib,jb,kb,n),qten(ib,jb,kb,n))
              elseif(difforder.eq.6)then
                ! for diff6, use '3d' array
                call diff6s(dt,ql0,dum1,dum2,dum3,q3d(ib,jb,kb,n),qten(ib,jb,kb,n))
              endif
            endif
            if(iturb.ge.1)then
              call turbs(0,dt,dosfcflx,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho0,rr0,rf0,qvflux,   &
                         dum1,dum2,dum3,dum4,divx,qa(ib,jb,kb,n),qten(ib,jb,kb,n),khh,khv,gx,gy,doimpl)
            endif
          endif
!---------------------------
        ENDDO
        IF(ipbl.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if(nqc.ne.0)   &
            qten(i,j,k,nqc) = qten(i,j,k,nqc) + qcpten(i,j,k)
            if(nqi.ne.0)   &
            qten(i,j,k,nqi) = qten(i,j,k,nqi) + qipten(i,j,k)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_pbl=time_pbl+mytime()
        ENDIF
      endif

!-------------------------------------------------------------------


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   Begin RK section   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ! time at end of full timestep:
      rtime=sngl(mtime+dt)

!--------------------------------------------------------------------
! RK3 begin

      DO NRK=1,3

        dttmp=dt/float(4-nrk)

!--------------------------------------------------------------------

      ! terms in theta and pi equations for proper mass/energy conservation
      ! Reference:  Bryan and Fritsch (2002, MWR)
      IF(neweqts.ge.1 .and. imoist.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum4(i,j,k)=cpl*q3d(i,j,k,nql1)
        enddo
        enddo
        enddo

        IF(nql2.ge.2)THEN
        DO n=(nql1+1),nql2
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum4(i,j,k)=dum4(i,j,k)+cpl*q3d(i,j,k,n)
          enddo
          enddo
          enddo
        ENDDO
        ENDIF

        if(iice.eq.1)then
          DO n=nqs1,nqs2
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum4(i,j,k)=dum4(i,j,k)+cpi*q3d(i,j,k,n)
            enddo
            enddo
            enddo
          ENDDO
        endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k,cpm,cvm)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          cpm=cp+cpv*q3d(i,j,k,nqv)+dum4(i,j,k)
          cvm=cv+cvv*q3d(i,j,k,nqv)+dum4(i,j,k)
          thterm(i,j,k)=( rd+rv*q3d(i,j,k,nqv)-rovcp*cpm )/cvm
          t22(i,j,k)=rovcp*cpm/cvm
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          t22(i,j,k)=rddcv
        enddo
        enddo
        enddo

      ENDIF   ! for neweqts=1
      if(timestats.ge.1) time_misc=time_misc+mytime()

!--------------------------------------------------------------------
        IF(nrk.ge.2)THEN
#ifdef MPI
          call comm_3u_end(u3d,uw31,uw32,ue31,ue32,   &
                               us31,us32,un31,un32,reqs_u)
          call comm_3v_end(v3d,vw31,vw32,ve31,ve32,   &
                               vs31,vs32,vn31,vn32,reqs_v)
          call comm_3w_end(w3d,ww31,ww32,we31,we32,   &
                               ws31,ws32,wn31,wn32,reqs_w)
#endif
          if(terrain_flag)then
            call bcwsfc(dzdx,dzdy,u3d,v3d,w3d)
            call bc2d(w3d(ib,jb,1))
          endif
        ENDIF
!--------------------------------------------------------------------

        IF(.not.terrain_flag)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=2,nk
          do j=0,nj+1
          do i=0,ni+1
            rrw(i,j,k)=rf0(i,j,k)*w3d(i,j,k)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_advs=time_advs+mytime()

        ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            rrw(i,j,k)=rf0(i,j,k)*( w3d(i,j,k)*gz(i,j)                 &
                 +0.25*(sigma(k)-zt)/(zt-zs(i,j))*(                    &
                       ( (u3d(i+1,j,k-1)+u3d(i  ,j,k-1))               &
                        +(u3d(i+1,j,k  )+u3d(i  ,j,k  )) )*dzdx(i,j)   &
                      +( (v3d(i,j+1,k-1)+v3d(i,j  ,k-1))               &
                        +(v3d(i,j+1,k  )+v3d(i,j  ,k  )) )*dzdy(i,j) ) )
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_advs=time_advs+mytime()

          call bcw(rrw,0)
#ifdef MPI
          call comm_1w_start(rrw,ww1,ww2,we1,we2,   &
                                 ws1,ws2,wn1,wn2,reqs_w)
#endif

        ENDIF

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=0,ni+2
          rru(i,j,k)=0.5*(rho0(i-1,j,k)+rho0(i,j,k))*u3d(i,j,k)
        enddo
        enddo
        enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+2
        do i=0,ni+1
          rrv(i,j,k)=0.5*(rho0(i,j-1,k)+rho0(i,j,k))*v3d(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_advs=time_advs+mytime()

#ifdef MPI
        if(terrain_flag)then
          call comm_1w_end(rrw,ww1,ww2,we1,we2,   &
                               ws1,ws2,wn1,wn2,reqs_w)
        endif
#endif

        IF(.not.terrain_flag)THEN

!------------
        IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            divx(i,j,k)=(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)        &
                       +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)        &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(i,j,k)
          enddo
          enddo
          enddo

        ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            divx(i,j,k)=(xf(i+1)*rru(i+1,j,k)-xf(i)*rru(i,j,k))*rdx*uh(i)*rxh(i)   &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(i,j,k)
          enddo
          enddo
          enddo

        ENDIF
!------------

        ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            divx(i,j,k)=(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)        &
                       +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)        &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(i,j,k)/gz(i,j)
          enddo
          enddo
          enddo

        ENDIF
        if(timestats.ge.1) time_divx=time_divx+mytime()


!--------------------------------------------------------------------
!  TKE advection
 
        IF(iturb.eq.1)THEN

          ! use wten for tke tendency, step tke forward:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nkt
          do j=1,nj
          do i=1,ni
            wten(i,j,k)=tketen(i,j,k)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_misc=time_misc+mytime()

#ifdef MPI
        IF(nrk.ge.2)THEN
          call comm_3t_end(tke3d,tkw1,tkw2,tke1,tke2,   &
                                 tks1,tks2,tkn1,tkn2,reqs_tk)
        ENDIF
#endif

          if(hadvorder.eq.5)then
            call adv5w(xh,rxh,uh,xf,vh,gz,mf,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,  &
                       rru,rrv,rrw,tke3d,wten)
          elseif(hadvorder.eq.6)then
            call adv6w(xh,rxh,uh,xf,vh,gz,mf,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,  &
                       rru,rrv,rrw,tke3d,wten)
          endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nkt
          do j=1,nj
          do i=1,ni
            tke3d(i,j,k)=tkea(i,j,k)+dttmp*wten(i,j,k)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_misc=time_misc+mytime()

          if(nrk.eq.3)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
            do k=1,nkt
            do j=1,nj
            do i=1,ni
              if(tke3d(i,j,k).lt.1.0e-6) tke3d(i,j,k)=0.0
            enddo
            enddo
            enddo
          endif
          if(timestats.ge.1) time_integ=time_integ+mytime()


          call bct(tke3d)
#ifdef MPI
          call comm_3t_start(tke3d,tkw1,tkw2,tke1,tke2,   &
                                   tks1,tks2,tkn1,tkn2,reqs_tk)
#endif

        ENDIF

!--------------------------------------------------------------------
!  Passive Tracers

        if(iptra.eq.1)then
          DO n=1,npt

      ! t33 = dummy

      bflag=0
      if(stat_qsrc.eq.1 .and. nrk.eq.3) bflag=1

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=ptten(i,j,k,n)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_misc=time_misc+mytime()


      if(nrk.eq.3)then
        pdef = 1
      else
        pdef = 0
      endif

#ifdef MPI
          IF(nrk.ge.2)THEN
            call comm_3s_end(pt3d(ib,jb,kb,n),                           &
                  tw1(1,1,1,n),tw2(1,1,1,n),te1(1,1,1,n),te2(1,1,1,n),   &
                  ts1(1,1,1,n),ts2(1,1,1,n),tn1(1,1,1,n),tn2(1,1,1,n),   &
                  reqs_t(1,n))
          ENDIF
#endif


    IF( (advweno.eq.1) .or. (advweno.eq.2.and.nrk.eq.3) )THEN
        call wenos(bflag,bfoo,xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,       &
                   rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,   &
                   divx,rru,rrv,rrw,pt3d(ib,jb,kb,n),sten,dt)
    ELSE
      if(hadvorder.eq.5.or.advweno.ge.1)then
        call adv5s(bflag,bfoo,xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,       &
                   rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,   &
                   rru,rrv,rrw,pta(ib,jb,kb,n),pt3d(ib,jb,kb,n),sten,pdef,dttmp)
      elseif(hadvorder.eq.6)then
        call adv6s(bflag,bfoo,xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,       &
                   rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,   &
                   rru,rrv,rrw,pta(ib,jb,kb,n),pt3d(ib,jb,kb,n),sten,pdef,dttmp)
      endif
    ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        pt3d(i,j,k,n)=pta(i,j,k,n)+dttmp*sten(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()

        if(nrk.eq.3) call pdefq(0.0,afoo,ruh,rvh,rmh,rho,pt3d(ib,jb,kb,n))

        call bcs(pt3d(ib,jb,kb,n))
#ifdef MPI
        call comm_3s_start(pt3d(ib,jb,kb,n)   &
                     ,tw1(1,1,1,n),tw2(1,1,1,n),te1(1,1,1,n),te2(1,1,1,n)     &
                     ,ts1(1,1,1,n),ts2(1,1,1,n),tn1(1,1,1,n),tn2(1,1,1,n)     &
                     ,reqs_t(1,n) )
#endif

          ENDDO
        endif

!--------------------------------------------------------------------
!  finish comms for q/theta:
#ifdef MPI
        IF(nrk.ge.2)THEN
        if(imoist.eq.1)then
          DO n=numq,1,-1
            call comm_3s_end(q3d(ib,jb,kb,n),                            &
                  qw1(1,1,1,n),qw2(1,1,1,n),qe1(1,1,1,n),qe2(1,1,1,n),   &
                  qs1(1,1,1,n),qs2(1,1,1,n),qn1(1,1,1,n),qn2(1,1,1,n),   &
                  reqs_q(1,n))
          ENDDO
        endif
          call comm_3s_end(th3d,sw31,sw32,se31,se32,   &
                                ss31,ss32,sn31,sn32,reqs_s)
        ENDIF
#endif
!--------------------------------------------------------------------
!  Calculate misc. variables
!
!    These arrays store variables that are used later in the
!    SOUND subroutine.  Do not modify t11 or t22 until after sound!
!    dum1 = vapor
!    dum2 = all liquid
!    dum3 = all solid

        IF(imoist.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            dum1(i,j,k)=max(q3d(i,j,k,nqv),0.0)
          enddo
          enddo
          enddo

          call getqli(1,q3d,dum2,dum3)

        ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            dum1(i,j,k)=0.0
            dum2(i,j,k)=0.0
            dum3(i,j,k)=0.0
          enddo
          enddo
          enddo

        ENDIF

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          t11(i,j,k)=(th0(i,j,k)+th3d(i,j,k))*(1.0+reps*dum1(i,j,k))     &
                     /(1.0+dum1(i,j,k)+max(0.0,dum2(i,j,k))+max(0.0,dum3(i,j,k)))
        enddo
        enddo
        enddo

        IF(thsmall.eq.1.and.imoist.eq.1)THEN
          ! save qv,ql,qi for buoyancy calculation:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            t12(i,j,k)=dum1(i,j,k)
            t13(i,j,k)=max(0.0,dum2(i,j,k))+max(0.0,dum3(i,j,k))
          enddo
          enddo
          enddo
        ENDIF

        if(timestats.ge.1) time_buoyan=time_buoyan+mytime()


!--------------------------------------------------------------------
! Moisture

        if(imoist.eq.1)then
          DO n=numq,1,-1

      ! t33 = dummy

      bflag=0
      if(stat_qsrc.eq.1 .and. nrk.eq.3) bflag=1

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=qten(i,j,k,n)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_misc=time_misc+mytime()


      if(nrk.eq.3)then
        pdef = 1
      else
        pdef = 0
      endif


    IF( (advweno.eq.1) .or. (advweno.eq.2.and.nrk.eq.3) )THEN
        call wenos(bflag,bsq(n),xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,       &
                   rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,   &
                   divx,rru,rrv,rrw,q3d(ib,jb,kb,n),sten,dt)
    ELSE
      if(hadvorder.eq.5.or.advweno.ge.1)then
        call adv5s(bflag,bsq(n),xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,       &
                   rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,   &
                   rru,rrv,rrw,qa(ib,jb,kb,n),q3d(ib,jb,kb,n),sten,pdef,dttmp)
      elseif(hadvorder.eq.6)then
        call adv6s(bflag,bsq(n),xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,       &
                   rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,   &
                   rru,rrv,rrw,qa(ib,jb,kb,n),q3d(ib,jb,kb,n),sten,pdef,dttmp)
      endif
    ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        q3d(i,j,k,n)=qa(i,j,k,n)+dttmp*sten(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()

      IF(nrk.lt.3)THEN
        qcom = .true.
        if(neweqts.eq.2.and.(n.eq.nqv.or.n.eq.2.or.n.eq.4)) qcom = .false.
        if( qcom )then
        call bcs(q3d(ib,jb,kb,n))
#ifdef MPI
        call comm_3s_start(q3d(ib,jb,kb,n)   &
                     ,qw1(1,1,1,n),qw2(1,1,1,n),qe1(1,1,1,n),qe2(1,1,1,n)     &
                     ,qs1(1,1,1,n),qs2(1,1,1,n),qn1(1,1,1,n),qn2(1,1,1,n)     &
                     ,reqs_q(1,n) )
#endif
        endif
      ENDIF

          ENDDO
        endif

!--------------------------------------------------------------------
!  THETA-equation

      ! t23  = theta used for advection  (full theta if thsmall=0)
      ! t33  = dummy

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        thten(i,j,k)=thten1(i,j,k)
      enddo
      enddo
      enddo

    IF(thsmall.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        t23(i,j,k)=th0(i,j,k)+th3d(i,j,k)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        if(abs(th3d(i,j,k)).lt.smeps)then
          t23(i,j,k)=0.0
        else
          t23(i,j,k)=th3d(i,j,k)
        endif
      enddo
      enddo
      enddo

    ENDIF
      if(timestats.ge.1) time_misc=time_misc+mytime()


    IF( (advweno.eq.1) .or. (advweno.eq.2.and.nrk.eq.3) )THEN
        call wenos(0,bfoo,xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,          &
                   rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,   &
                   divx,rru,rrv,rrw,t23,thten,dt)
    ELSE
      if(hadvorder.eq.5.or.advweno.ge.1)then
        call adv5s(0,bfoo,xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,          &
                   rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,   &
                   rru,rrv,rrw,tha,t23,thten,0,dttmp)
      elseif(hadvorder.eq.6)then
        call adv6s(0,bfoo,xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,          &
                   rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,   &
                   rru,rrv,rrw,tha,t23,thten,0,dttmp)
      endif
    ENDIF


    IF(thsmall.eq.0.and.neweqts.ge.1.and.imoist.eq.1)THEN

      ! this section of code only accessed if thsmall=0
      ! (for which t23 = full theta)

      IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          thten(i,j,k)=thten(i,j,k)-t23(i,j,k)*(              &
                      (u3d(i+1,j,k)-u3d(i,j,k))*rdx*uh(i)      &
                     +(v3d(i,j+1,k)-v3d(i,j,k))*rdy*vh(j)      &
                     +(w3d(i,j,k+1)-w3d(i,j,k))*rdz*mh(i,j,k) )    &
                    *thterm(i,j,k)
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          thten(i,j,k)=thten(i,j,k)-t23(i,j,k)*(              &
                      (xf(i+1)*u3d(i+1,j,k)-xf(i)*u3d(i,j,k))*rdx*uh(i)*rxh(i)      &
                     +(w3d(i,j,k+1)-w3d(i,j,k))*rdz*mh(i,j,k) )    &
                    *thterm(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

        IF(terrain_flag)THEN

          ! dum1 = dudz
          ! dum4 = dvdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=2,nk-1
          do j=1,nj
          do i=1,ni+1
            dum1(i,j,k)=gx(i,j,k)*(u3d(i,j,k+1)-u3d(i,j,k-1))*rdz2
          enddo
          enddo
          enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni+1
            dum1(i,j,1 )=gx(i,j,1 )*(u3d(i,j,2 )-u3d(i,j,1   ))*rdz
            dum1(i,j,nk)=gx(i,j,nk)*(u3d(i,j,nk)-u3d(i,j,nk-1))*rdz
          enddo
          enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=2,nk-1
          do j=1,nj+1
          do i=1,ni
            dum4(i,j,k)=gy(i,j,k)*(v3d(i,j,k+1)-v3d(i,j,k-1))*rdz2
          enddo
          enddo
          enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj+1
          do i=1,ni
            dum4(i,j,1 )=gy(i,j,1 )*(v3d(i,j,2 )-v3d(i,j,1   ))*rdz
            dum4(i,j,nk)=gy(i,j,nk)*(v3d(i,j,nk)-v3d(i,j,nk-1))*rdz
          enddo
          enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            thten(i,j,k)=thten(i,j,k)-t23(i,j,k)*0.5*(      &
                           (dum1(i,j,k)+dum1(i+1,j,k))       &
                          +(dum4(i,j,k)+dum4(i,j+1,k))  )    &
                      *thterm(i,j,k)
          enddo
          enddo
          enddo

        ENDIF   ! for terrain

    ENDIF
    if(timestats.ge.1) time_misc=time_misc+mytime()

    IF(thsmall.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        th3d(i,j,k)=tha(i,j,k)+dttmp*thten(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()

      tcom = .true.
      if(imoist.eq.1.and.neweqts.eq.2) tcom = .false. ! dont start comm for neweqts=2
      if(imoist.eq.1.and.nrk.eq.3) tcom = .false.     ! dont start comm if moist & nrk=3
      if(thsmall.eq.1) tcom = .false.                 ! dont start comm if thsmall=1

      IF( tcom )THEN
        call bcs(th3d)
#ifdef MPI
        call comm_3s_start(th3d,sw31,sw32,se31,se32,   &
                                ss31,ss32,sn31,sn32,reqs_s)
#endif
      ENDIF

    ENDIF

!--------------------------------------------------------------------
!  Pressure equation

      ! t22  = ppterm
      ! t33  = dummy

      IF(psolver.le.3)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          sten(i,j,k)=ppten(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()

#ifdef MPI
        IF(nrk.ge.2)THEN
          call comm_3s_end(pp3d,pw31,pw32,pe31,pe32,   &
                                ps31,ps32,pn31,pn32,reqs_p)
        ENDIF
#endif

        if(hadvorder.eq.5)then
          call adv5s(0,bfoo,xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,         &
                     rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,   &
                     rru,rrv,rrw,ppi,pp3d,sten,0,dttmp)
        elseif(hadvorder.eq.6)then
          call adv6s(0,bfoo,xh,rxh,uh,ruh,xf,vh,rvh,gz,mh,rmh,         &
                     rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,   &
                     rru,rrv,rrw,ppi,pp3d,sten,0,dttmp)
        endif

        IF(imoist.eq.1.and.neweqts.eq.2)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            t33(i,j,k)=ppi(i,j,k)+dttmp*sten(i,j,k)
            t23(i,j,k)=t33(i,j,k)
          enddo
          enddo
          enddo

          IF(thsmall.eq.1)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              th3d(i,j,k)=tha(i,j,k)+dttmp*thten(i,j,k)
            enddo
            enddo
            enddo
          ENDIF
          if(timestats.ge.1) time_satadj=time_satadj+mytime()

          call calcprs(pi0,prs,t33)
          call calcrho(pi0,th0,rho,prs,t33,th3d,q3d)
          IF(nrk.ge.3)THEN
            IF(axisymm.eq.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
              do k=1,nk
              do j=1,nj
              do i=1,ni
                dum3(i,j,k)=rho(i,j,k)
              enddo
              enddo
              enddo
            ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
              do k=1,nk
              do j=1,nj
              do i=1,ni
                dum3(i,j,k) = rho(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
              enddo
              enddo
              enddo
            ENDIF
          ENDIF
          if(ptype.eq.1.or.ptype.eq.3.or.ptype.eq.5.or.ptype.eq.6)then
            call satadj(nrk,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                        rho,dum3,t33,prs,th3d,q3d)
          elseif(ptype.eq.2)then
            call satadj_ice(nrk,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                            rho,dum3,t33,prs,th3d,                      &
                q3d(ib,jb,kb,1),q3d(ib,jb,kb,2),q3d(ib,jb,kb,3),   &
                q3d(ib,jb,kb,4),q3d(ib,jb,kb,5),q3d(ib,jb,kb,6))
          endif

          IF(thsmall.eq.1)THEN
            rdt=1.0/dttmp
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              thten(i,j,k)=(th3d(i,j,k)-tha(i,j,k))*rdt
            enddo
            enddo
            enddo
            if(timestats.ge.1) time_satadj=time_satadj+mytime()
          ENDIF

          if(nrk.lt.3)then
            call bcs(q3d(ib,jb,kb,nqv))
#ifdef MPI
            call comm_3s_start(q3d(ib,jb,kb,nqv),                            &
              qw1(1,1,1,nqv),qw2(1,1,1,nqv),qe1(1,1,1,nqv),qe2(1,1,1,nqv),   &
              qs1(1,1,1,nqv),qs2(1,1,1,nqv),qn1(1,1,1,nqv),qn2(1,1,1,nqv),   &
              reqs_q(1,nqv))
#endif
            call bcs(q3d(ib,jb,kb,2))
#ifdef MPI
            call comm_3s_start(q3d(ib,jb,kb,2),                            &
              qw1(1,1,1,2),qw2(1,1,1,2),qe1(1,1,1,2),qe2(1,1,1,2),   &
              qs1(1,1,1,2),qs2(1,1,1,2),qn1(1,1,1,2),qn2(1,1,1,2),   &
              reqs_q(1,2))
#endif
            if(iice.eq.1)then
              call bcs(q3d(ib,jb,kb,4))
#ifdef MPI
              call comm_3s_start(q3d(ib,jb,kb,4),                            &
                qw1(1,1,1,4),qw2(1,1,1,4),qe1(1,1,1,4),qe2(1,1,1,4),   &
                qs1(1,1,1,4),qs2(1,1,1,4),qn1(1,1,1,4),qn2(1,1,1,4),   &
                reqs_q(1,4))
#endif
            endif
            if(thsmall.eq.0)then
              call bcs(th3d)
#ifdef MPI
              call comm_3s_start(th3d,sw31,sw32,se31,se32,   &
                                      ss31,ss32,sn31,sn32,reqs_s)
#endif
            endif
          endif

          rdt=1.0/dttmp

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            sten(i,j,k)=sten(i,j,k)+(t33(i,j,k)-t23(i,j,k))*rdt
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_satadj=time_satadj+mytime()

        ENDIF

      ENDIF

!--------------------------------------------------------------------
!  U-equation

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          uten(i,j,k)=uten1(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()
 
        if(icor.eq.1)then

        if(pertcor.eq.1)then

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj+1
          do i=0,ni+1
            dum1(i,j,k)=v3d(i,j,k)-v0(i,j,k)
          enddo
          enddo
          enddo

        else

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj+1
          do i=0,ni+1
            dum1(i,j,k)=v3d(i,j,k)
          enddo
          enddo
          enddo

        endif

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni+1
            uten(i,j,k)=uten(i,j,k)+fcor*             &
             0.25*( (dum1(i  ,j,k)+dum1(i  ,j+1,k))   &
                   +(dum1(i-1,j,k)+dum1(i-1,j+1,k)) )
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_cor=time_cor+mytime()

        endif

        if(axisymm.eq.1)then

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=2,ni+1
            uten(i,j,k)=uten(i,j,k)+0.5*(   &
                 ( v3d(i-1,j,k)**2)*rxh(i-1)+(v3d(i,j,k)**2)*rxh(i) )
          enddo
          enddo
          enddo

        endif

        if(hadvorder.eq.5)then
          call adv5u(xf,rxf,uf,vh,gz,mh,rho0,rr0,rf0,rrf0,rru0,dum1,dum2,dum3,dum4,divx,  &
                     rru,u3d,uten,rrv,rrw)
        elseif(hadvorder.eq.6)then
          call adv6u(xf,rxf,uf,vh,gz,mh,rho0,rr0,rf0,rrf0,rru0,dum1,dum2,dum3,dum4,divx,  &
                     rru,u3d,uten,rrv,rrw)
        endif

!--------------------------------------------------------------------
!  V-equation
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vten(i,j,k)=vten1(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()
 
        if(icor.eq.1)then

!--------------
          IF(axisymm.eq.0)THEN

        if(pertcor.eq.1)then

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=1,ni+1
            dum1(i,j,k)=u3d(i,j,k)-u0(i,j,k)
          enddo
          enddo
          enddo

        else

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=1,ni+1
            dum1(i,j,k)=u3d(i,j,k)
          enddo
          enddo
          enddo

        endif

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj+1
          do i=1,ni
            vten(i,j,k)=vten(i,j,k)-fcor*             &
             0.25*( (dum1(i,j  ,k)+dum1(i+1,j  ,k))   &
                   +(dum1(i,j-1,k)+dum1(i+1,j-1,k)) )
          enddo
          enddo
          enddo

          ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            vten(i,j,k)=vten(i,j,k)-fcor*0.5*(xf(i)*u3d(i,j,k)+xf(i+1)*u3d(i+1,j,k))*rxh(i)
          enddo
          enddo
          enddo

          ENDIF
!--------------

          if(timestats.ge.1) time_cor=time_cor+mytime()

        endif

        if(axisymm.eq.1)then

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            vten(i,j,k)=vten(i,j,k)-(v3d(i,j,k)*rxh(i))*0.5*(xf(i)*u3d(i,j,k)+xf(i+1)*u3d(i+1,j,k))*rxh(i)
          enddo
          enddo
          enddo

        endif

        if(hadvorder.eq.5)then
          call adv5v(xh,rxh,uh,xf,vf,gz,mh,rho0,rr0,rf0,rrf0,rrv0,dum1,dum2,dum3,dum4,divx,  &
                     rru,rrv,v3d,vten,rrw)
        elseif(hadvorder.eq.6)then
          call adv6v(xh,rxh,uh,xf,vf,gz,mh,rho0,rr0,rf0,rrf0,rrv0,dum1,dum2,dum3,dum4,divx,  &
                     rru,rrv,v3d,vten,rrw)
        endif

!--------------------------------------------------------------------
!  W-equation

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          wten(i,j,k)=wten1(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()
 
        if( (thsmall.eq.0) .or. (thsmall.eq.1.and.imoist.eq.1) )then
          ! buoyancy:
          IF(thsmall.eq.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum2(i,j,k) = t11(i,j,k)/thv0(i,j,k)-1.0
            enddo
            enddo
            enddo
          ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum2(i,j,k) = repsm1*(t12(i,j,k)-qv0(i,j,k)) - (t13(i,j,k)-qc0(i,j,k))
            enddo
            enddo
            enddo
          ENDIF
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            wten(i,j,k)=wten(i,j,k)+govtwo*(dum2(i,j,k-1)+dum2(i,j,k))
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_buoyan=time_buoyan+mytime()
        endif

        if(hadvorder.eq.5)then
          call adv5w(xh,rxh,uh,xf,vh,gz,mf,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,  &
                     rru,rrv,rrw,w3d,wten)
        elseif(hadvorder.eq.6)then
          call adv6w(xh,rxh,uh,xf,vh,gz,mf,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,  &
                     rru,rrv,rrw,w3d,wten)
        endif


!--------------------------------------------------------------------
!  Update v for axisymmetric model simulations:

        IF(axisymm.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            v3d(i,j,k)=va(i,j,k)+dttmp*vten(i,j,k)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_misc=time_misc+mytime()

          call bcv(v3d)

        ENDIF

!--------------------------------------------------------------------
!  call sound

        IF(psolver.eq.1)THEN

          call soundns(xh,rxh,uh,xf,uf,yh,vh,yf,vf,zh,mh,mf,pi0,thv0,       &
                       radbcw,radbce,radbcs,radbcn,                &
                       divx,u0,ua,u3d,uten,v0,va,v3d,vten,wa,w3d,wten,   &
                       ppi,pp3d,sten,t11,t22,dttmp,nrk,rtime,      &
                       reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,         &
                       uw31,uw32,ue31,ue32,us31,us32,un31,un32,    &
                       vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,    &
                       ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,    &
                       sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,    &
                       pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2)

        ELSEIF(psolver.eq.2)THEN

          call sounde(dt,xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,zh,mh,rmh,mf,rmf,   &
                     pi0,thv0,rho0,rr0,rf0,th0,dzdx,dzdy,            &
                     radbcw,radbce,radbcs,radbcn,                    &
                     dum1,dum2,dum3,dum4,divx,                       &
                     gx,u0,ua,u3d,uten,gy,v0,va,v3d,vten,wa,w3d,wten,      &
                     ppi,pp3d,sten,tha,th3d,thten,thterm,tk,         &
                     t11,t22,nrk,rtime,                              &
                     reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,             &
                     uw31,uw32,ue31,ue32,us31,us32,un31,un32,        &
                     vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,        &
                     ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,        &
                     sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,        &
                     pw31,pw32,pe31,pe32,ps31,ps32,pn31,pn32,        &
                     pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2)

        ELSEIF(psolver.eq.3)THEN

          ! rho,prs,divx are used as dummy arrays by sound

          call sound(dt,xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,zh,mh,rmh,mf,rmf, &
                     pi0,thv0,rho0,rr0,rf0,th0,dzdx,dzdy,            &
                     radbcw,radbce,radbcs,radbcn,                    &
                     dum1,dum2,dum3,dum4,t12,t13,t23,t33,            &
                     gx,u0,ua,u3d,uten,gy,v0,va,v3d,vten,wa,w3d,wten,      &
                     ppi,pp3d,sten,tha,th3d,thten,thterm,tk,         &
                     t11,t22,rho,prs,divx,nrk,rtime,                 &
                     reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,             &
                     uw31,uw32,ue31,ue32,us31,us32,un31,un32,        &
                     vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,        &
                     ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,        &
                     sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,        &
                     pw31,pw32,pe31,pe32,ps31,ps32,pn31,pn32,        &
                     pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2)

        ELSEIF(psolver.eq.4.or.psolver.eq.5)THEN
          ! anelastic/incompressible solver:

          call anelp(xh,uh,xf,uf,yh,vh,yf,vf,                     &
                     zh,mh,rmh,mf,rmf,pi0,thv0,rho0,prs0,rf0,     &
                     radbcw,radbce,radbcs,radbcn,dum1,divx,       &
                     u0,ua,u3d,uten,v0,va,v3d,vten,wa,w3d,wten,   &
                     ppi,pp3d,t11,cfb,cfa,cfc,ad1,ad2,pdt,deft,rhs,trans,dttmp,nrk,rtime)

        ENDIF

!--------------------------------------------------------------------
!  radbc

        if(irbc.eq.4)then

          if(ibw.eq.1 .or. ibe.eq.1)then
            call radbcew4(ruf,radbcw,radbce,ua,u3d,dttmp)
          endif

          if(ibs.eq.1 .or. ibn.eq.1)then
            call radbcns4(rvf,radbcs,radbcn,va,v3d,dttmp)
          endif

        endif

!--------------------------------------------------------------------
! RK loop end

      ENDDO


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   End of RK section   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!--------------------------------------------------------------------
!  Get pressure

    IF(psolver.eq.4.or.psolver.eq.5)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        prs(i,j,k)=prs0(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_misc=time_misc+mytime()

    ELSE

      call calcprs(pi0,prs,pp3d)

    ENDIF

!--------------------------------------------------------------------
!  Get density

    IF(psolver.eq.4.or.psolver.eq.5)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        rho(i,j,k)=rho0(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_misc=time_misc+mytime()

    ELSE

      call calcrho(pi0,th0,rho,prs,pp3d,th3d,q3d)

    ENDIF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   BEGIN microphysics   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(imoist.eq.1)THEN

        getdbz = .false.
        IF(output_dbz.eq.1)THEN
          rtime=sngl(mtime+dt)
          if( (rtime.ge.sngl(taptim)).or.stopit )then
            getdbz = .true.
          endif
          if(getdbz)then
            write(outfile,*) '  Getting dbz ... '
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              sten(i,j,k)=0.0
            enddo
            enddo
            enddo
          endif
        ENDIF

!-----------------------------------------------------------------------
!  store t in dum1

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  prep for efall calculation:  store cvm in dum2

        IF(efall.eq.1)THEN
          if(neweqts.ge.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum2(i,j,k)=q3d(i,j,k,nqv)
            enddo
            enddo
            enddo
            call getqli(0,q3d,dum3,dum4)
          else
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum2(i,j,k)=0.0
              dum3(i,j,k)=0.0
              dum4(i,j,k)=0.0
            enddo
            enddo
            enddo
          endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,ql)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum2(i,j,k)=cv+cvv*dum2(i,j,k)+cpl*dum3(i,j,k)+cpi*dum4(i,j,k)
          enddo
          enddo
          enddo
        ENDIF

!-----------------------------------------------------------------------
!  store appropriate rho for budget calculations in dum3

        IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum3(i,j,k)=rho(i,j,k)
          enddo
          enddo
          enddo

        ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum3(i,j,k) = rho(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
          enddo
          enddo
          enddo

        ENDIF


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  NOTES:
!           sten       is used for     dbz
!
!           dum1   is   T
!           dum2   is   cvm
!           dum3   is   rho for budget calculations
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Kessler scheme   cccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        IF(ptype.eq.1)THEN
          simple_comm = .false.
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq(1.0e-14,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq(1.0e-14,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call k_fallout(rho,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call geterain(dt,cpl,lv1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          if(efall.ge.1)then
            call getefall(dt,cpl,ruh,rvh,mf,pi0,th0,dum1,dum2,dum3,   &
                          pp3d,th3d,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call kessler(dt,qbudget(3),qbudget(4),qbudget(5),ruh,rvh,rmh,pi0,th0,dum1,   &
                       rho,dum3,pp3d,th3d,prs,                            &
                       q3d(ib,jb,kb,nqv),q3d(ib,jb,kb,2),q3d(ib,jb,kb,3))
          call bcs(q3d(ib,jb,kb,3))
#ifdef MPI
          call comm_3s_start(q3d(ib,jb,kb,3),                               &
                     qw1(1,1,1,3),qw2(1,1,1,3),qe1(1,1,1,3),qe2(1,1,1,3),   &
                     qs1(1,1,1,3),qs2(1,1,1,3),qn1(1,1,1,3),qn2(1,1,1,3),   &
                     reqs_q(1,3))
#endif
          call satadj(4,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                      rho,dum3,pp3d,prs,th3d,q3d)
          call bcs(q3d(ib,jb,kb,1))
#ifdef MPI
          call comm_3s_start(q3d(ib,jb,kb,1),                               &
                     qw1(1,1,1,1),qw2(1,1,1,1),qe1(1,1,1,1),qe2(1,1,1,1),   &
                     qs1(1,1,1,1),qs2(1,1,1,1),qn1(1,1,1,1),qn2(1,1,1,1),   &
                     reqs_q(1,1))
#endif
          call bcs(q3d(ib,jb,kb,2))
#ifdef MPI
          call comm_3s_start(q3d(ib,jb,kb,2),                               &
                     qw1(1,1,1,2),qw2(1,1,1,2),qe1(1,1,1,2),qe2(1,1,1,2),   &
                     qs1(1,1,1,2),qs2(1,1,1,2),qn1(1,1,1,2),qn2(1,1,1,2),   &
                     reqs_q(1,2))
#endif
          call bcs(th3d)
#ifdef MPI
          call comm_3s_start(th3d,sw31,sw32,se31,se32,   &
                                  ss31,ss32,sn31,sn32,reqs_s)
#endif
          call bcs(pp3d)
#ifdef MPI
          call comm_3s_start(pp3d,pw31,pw32,pe31,pe32,   &
                                  ps31,ps32,pn31,pn32,reqs_p)
#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Goddard LFO scheme   cccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ELSEIF(ptype.eq.2)THEN
          simple_comm = .false.
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq(1.0e-14,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq(1.0e-14,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq(1.0e-14,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq(1.0e-14,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq(1.0e-14,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
          call goddard(dt,qbudget(3),qbudget(4),qbudget(5),ruh,rvh,rmh,pi0,th0,             &
                       rho,dum3,prs,pp3d,th3d,                            &
     q3d(ib,jb,kb,1), q3d(ib,jb,kb,2),q3d(ib,jb,kb,3),qten(ib,jb,kb,3),   &
     q3d(ib,jb,kb,4),qten(ib,jb,kb,4),q3d(ib,jb,kb,5),qten(ib,jb,kb,5),   &
     q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          call satadj_ice(4,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,     &
                          rho,dum3,pp3d,prs,th3d,                     &
              q3d(ib,jb,kb,1),q3d(ib,jb,kb,2),q3d(ib,jb,kb,3),   &
              q3d(ib,jb,kb,4),q3d(ib,jb,kb,5),q3d(ib,jb,kb,6))
          call bcs(q3d(ib,jb,kb,1))
#ifdef MPI
          call comm_3s_start(q3d(ib,jb,kb,1),                               &
                     qw1(1,1,1,1),qw2(1,1,1,1),qe1(1,1,1,1),qe2(1,1,1,1),   &
                     qs1(1,1,1,1),qs2(1,1,1,1),qn1(1,1,1,1),qn2(1,1,1,1),   &
                     reqs_q(1,1))
#endif
          call bcs(q3d(ib,jb,kb,2))
#ifdef MPI
          call comm_3s_start(q3d(ib,jb,kb,2),                               &
                     qw1(1,1,1,2),qw2(1,1,1,2),qe1(1,1,1,2),qe2(1,1,1,2),   &
                     qs1(1,1,1,2),qs2(1,1,1,2),qn1(1,1,1,2),qn2(1,1,1,2),   &
                     reqs_q(1,2))
#endif
          call geterain(dt,cpl,lv1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call geterain(dt,cpi,ls1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,4),qten(ib,jb,kb,4))
          call geterain(dt,cpi,ls1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,5),qten(ib,jb,kb,5))
          call geterain(dt,cpi,ls1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          if(efall.ge.1)then
            call getefall(dt,cpl,ruh,rvh,mf,pi0,th0,dum1,dum2,dum3,   &
                          pp3d,th3d,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
            call getefall(dt,cpi,ruh,rvh,mf,pi0,th0,dum1,dum2,dum3,   &
                          pp3d,th3d,q3d(ib,jb,kb,4),qten(ib,jb,kb,4))
            call getefall(dt,cpi,ruh,rvh,mf,pi0,th0,dum1,dum2,dum3,   &
                          pp3d,th3d,q3d(ib,jb,kb,5),qten(ib,jb,kb,5))
            call getefall(dt,cpi,ruh,rvh,mf,pi0,th0,dum1,dum2,dum3,   &
                          pp3d,th3d,q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          endif
          call bcs(th3d)
#ifdef MPI
          call comm_3s_start(th3d,sw31,sw32,se31,se32,   &
                                  ss31,ss32,sn31,sn32,reqs_s)
#endif
          call bcs(pp3d)
#ifdef MPI
          call comm_3s_start(pp3d,pw31,pw32,pe31,pe32,   &
                                  ps31,ps32,pn31,pn32,reqs_p)
#endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call bcs(q3d(ib,jb,kb,3))
#ifdef MPI
          call comm_3s_start(q3d(ib,jb,kb,3),                               &
                     qw1(1,1,1,3),qw2(1,1,1,3),qe1(1,1,1,3),qe2(1,1,1,3),   &
                     qs1(1,1,1,3),qs2(1,1,1,3),qn1(1,1,1,3),qn2(1,1,1,3),   &
                     reqs_q(1,3))
#endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,4),qten(ib,jb,kb,4))
          call bcs(q3d(ib,jb,kb,4))
#ifdef MPI
          call comm_3s_start(q3d(ib,jb,kb,4),                               &
                     qw1(1,1,1,4),qw2(1,1,1,4),qe1(1,1,1,4),qe2(1,1,1,4),   &
                     qs1(1,1,1,4),qs2(1,1,1,4),qn1(1,1,1,4),qn2(1,1,1,4),   &
                     reqs_q(1,4))
#endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,5),qten(ib,jb,kb,5))
          call bcs(q3d(ib,jb,kb,5))
#ifdef MPI
          call comm_3s_start(q3d(ib,jb,kb,5),                               &
                     qw1(1,1,1,5),qw2(1,1,1,5),qe1(1,1,1,5),qe2(1,1,1,5),   &
                     qs1(1,1,1,5),qs2(1,1,1,5),qn1(1,1,1,5),qn2(1,1,1,5),   &
                     reqs_q(1,5))
#endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          call bcs(q3d(ib,jb,kb,6))
#ifdef MPI
          call comm_3s_start(q3d(ib,jb,kb,6),                               &
                     qw1(1,1,1,6),qw2(1,1,1,6),qe1(1,1,1,6),qe2(1,1,1,6),   &
                     qs1(1,1,1,6),qs2(1,1,1,6),qn1(1,1,1,6),qn2(1,1,1,6),   &
                     reqs_q(1,6))
#endif
          if(getdbz) call calcdbz(rho,q3d(ib,jb,kb,3),q3d(ib,jb,kb,5),q3d(ib,jb,kb,6),sten)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Thompson scheme   ccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ELSEIF(ptype.eq.3)THEN
          simple_comm = .true.
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq(1.0e-12,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq(1.0e-12,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq(1.0e-12,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq(1.0e-12,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq(1.0e-12,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
!!!          call pdefq(    1.0,asq(7),ruh,rvh,rmh,rho,q3d(ib,jb,kb,7))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! dum1 = pi
            ! dum2 = dz
            ! dum4 = T
            dum1(i,j,k)=pi0(i,j,k)+pp3d(i,j,k)
            dum2(i,j,k)=dz*rmh(i,j,k)
            dum4(i,j,k)=(th0(i,j,k)+th3d(i,j,k))*dum1(i,j,k)
            ! store old T in thten array:
            thten(i,j,k)=dum4(i,j,k)
          enddo
          enddo
          enddo
          call mp_gt_driver(q3d(ib,jb,kb,1),q3d(ib,jb,kb,2),q3d(ib,jb,kb,3), &
                            q3d(ib,jb,kb,4),q3d(ib,jb,kb,5),q3d(ib,jb,kb,6), &
                            q3d(ib,jb,kb,7),q3d(ib,jb,kb,8),                 &
                            th0,dum4,dum1,prs,dum2,dt,rain,                 &
                            qbudget(5),qbudget(6),                           &
                            ruh,rvh,rmh,dum3,sten,getdbz)
        IF(neweqts.ge.1)THEN
          ! for mass conservation:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( abs(dum4(i,j,k)-thten(i,j,k)).ge.tsmall )then
!!!              th3d(i,j,k)=th3d(i,j,k)+(dum4(i,j,k)-thten(i,j,k))
!!!              pp3d(i,j,k)=((rho(i,j,k)*(rd+rv*q3d(i,j,k,nqv))   &
!!!                                      *(th0(i,j,k)+th3d(i,j,k))*rp00)**rddcv)-pi0(i,j,k)
!!!              prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
              prs(i,j,k)=rho(i,j,k)*rd*dum4(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps)
              pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
              th3d(i,j,k)=dum4(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
            endif
          enddo
          enddo
          enddo
        ELSE
          ! traditional thermodynamics:  p,pi remain unchanged
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( abs(dum4(i,j,k)-thten(i,j,k)).ge.tsmall )then
              th3d(i,j,k)=th3d(i,j,k)+(dum4(i,j,k)-thten(i,j,k))/dum1(i,j,k)
              rho(i,j,k)=prs(i,j,k)/(rd*dum4(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps))
            endif
          enddo
          enddo
          enddo
        ENDIF
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
          call satadj(4,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                      rho,dum3,pp3d,prs,th3d,q3d)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   GSR LFO scheme   cccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.4)THEN
          simple_comm = .true.
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq(1.0e-14,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq(1.0e-14,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq(1.0e-14,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq(1.0e-14,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq(1.0e-14,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
          call lfo_ice_drive(dt, mf, pi0, prs0, pp3d, prs, th0, th3d,    &
                             qv0, rho0, q3d, qten, dum1)
          do n=2,numq
            call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                         q3d(ib,jb,kb,n),qten(ib,jb,kb,n))
          enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Morrison scheme   cccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.5)THEN
          simple_comm = .true.
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq(1.0e-12,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq(1.0e-12,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq(1.0e-12,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq(1.0e-12,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq(1.0e-12,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
!!!          call pdefq(    1.0,asq(7),ruh,rvh,rmh,rho,q3d(ib,jb,kb,7))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! dum1 = T  (this should have been calculated already)
            ! dum4 = dz
            dum4(i,j,k)=dz*rmh(i,j,k)
            ! store old T in thten array:
            thten(i,j,k)=dum1(i,j,k)
          enddo
          enddo
          enddo
          call MP_MORR_TWO_MOMENT(nstep,dum1,                                 &
                          q3d(ib,jb,kb, 1),q3d(ib,jb,kb, 2),q3d(ib,jb,kb, 3), &
                          q3d(ib,jb,kb, 4),q3d(ib,jb,kb, 5),q3d(ib,jb,kb, 6), &
                          q3d(ib,jb,kb, 7),q3d(ib,jb,kb, 8),q3d(ib,jb,kb, 9), &
                          q3d(ib,jb,kb,10),                                   &
                               prs,dt,dum4,w3d,rain,                         &
                          qbudget(1),qbudget(2),qbudget(5),qbudget(6),        &
                          ruh,rvh,rmh,dum3,sten,getdbz)
        IF(neweqts.ge.1)THEN
          ! for mass conservation:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall )then
!!!              th3d(i,j,k)=th3d(i,j,k)+(dum1(i,j,k)-thten(i,j,k))
!!!              pp3d(i,j,k)=((rho(i,j,k)*(rd+rv*q3d(i,j,k,nqv))   &
!!!                                      *(th0(i,j,k)+th3d(i,j,k))*rp00)**rddcv)-pi0(i,j,k)
!!!              prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
              prs(i,j,k)=rho(i,j,k)*rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps)
              pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
              th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
            endif
          enddo
          enddo
          enddo
        ELSE
          ! traditional thermodynamics:  p,pi remain unchanged
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall )then
              th3d(i,j,k)=th3d(i,j,k)+(dum1(i,j,k)-thten(i,j,k))/(pi0(i,j,k)+pp3d(i,j,k))
              rho(i,j,k)=prs(i,j,k)/(rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps))
            endif
          enddo
          enddo
          enddo
        ENDIF
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
          call satadj(4,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                      rho,dum3,pp3d,prs,th3d,q3d)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   RE87 scheme   ccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.6)THEN
          simple_comm = .true.
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq(1.0e-14,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if(q3d(i,j,k,2).gt.0.001)then
              qten(i,j,k,2) = v_t
            else
              qten(i,j,k,2) = 0.0
            endif
          enddo
          enddo
          enddo
          call geterain(dt,cpl,lv1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,2),qten(ib,jb,kb,2))
          if(efall.ge.1)then
            call getefall(dt,cpl,ruh,rvh,mf,pi0,th0,dum1,dum2,dum3,   &
                          pp3d,th3d,q3d(ib,jb,kb,2),qten(ib,jb,kb,2))
          endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,2),qten(ib,jb,kb,2))
          call satadj(4,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                      rho,dum3,pp3d,prs,th3d,q3d)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Milbrandt & Yao scheme
!
!        ELSEIF(ptype.eq.7)THEN
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Ziegler/Mansell two-moment scheme
!
!        ELSEIF(ptype.ge.26)THEN
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  insert new microphysics schemes here
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        ELSEIF(ptype.eq.8)THEN
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! otherwise, stop for undefined ptype
        ELSE
          print *,'  Undefined ptype!'
          call stopcm1
        ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Begin:  message passing for simple_comm
        IF(simple_comm)THEN
          call bcs(th3d)
#ifdef MPI
          call comm_3s_start(th3d,sw31,sw32,se31,se32,   &
                                  ss31,ss32,sn31,sn32,reqs_s)
#endif
          call bcs(pp3d)
#ifdef MPI
          call comm_3s_start(pp3d,pw31,pw32,pe31,pe32,   &
                                  ps31,ps32,pn31,pn32,reqs_p)
#endif
          DO n=1,numq
            call bcs(q3d(ib,jb,kb,n))
#ifdef MPI
            call comm_3s_start(q3d(ib,jb,kb,n),                               &
                       qw1(1,1,1,n),qw2(1,1,1,n),qe1(1,1,1,n),qe2(1,1,1,n),   &
                       qs1(1,1,1,n),qs2(1,1,1,n),qn1(1,1,1,n),qn2(1,1,1,n),   &
                       reqs_q(1,n))
#endif
          ENDDO
        ENDIF
!Done:  message passing for simple_comm

      ENDIF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   END microphysics   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!-----------------------------------------------------------------
!  Equate the two arrays

      if(iturb.eq.1)then
#ifdef MPI
          call comm_3t_end(tke3d,tkw1,tkw2,tke1,tke2,   &
                                 tks1,tks2,tkn1,tkn2,reqs_tk)
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=0,nk+2
        do j=0,nj+1
        do i=0,ni+1
          tkea(i,j,k)=tke3d(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_integ=time_integ+mytime()
      endif

      if(iptra.eq.1)then
        do n=1,npt
#ifdef MPI
            call comm_3s_end(pt3d(ib,jb,kb,n),                           &
                  tw1(1,1,1,n),tw2(1,1,1,n),te1(1,1,1,n),te2(1,1,1,n),   &
                  ts1(1,1,1,n),ts2(1,1,1,n),tn1(1,1,1,n),tn2(1,1,1,n),   &
                  reqs_t(1,n))
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=0,nk+1
          do j=0,nj+1
          do i=0,ni+1
            pta(i,j,k,n)=pt3d(i,j,k,n)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_integ=time_integ+mytime()
        enddo
      endif

#ifdef MPI
      call comm_3u_end(u3d,uw31,uw32,ue31,ue32,   &
                           us31,us32,un31,un32,reqs_u)
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+2
        ua(i,j,k)=u3d(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()
 
#ifdef MPI
      call comm_3v_end(v3d,vw31,vw32,ve31,ve32,   &
                           vs31,vs32,vn31,vn32,reqs_v)
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+2
      do i=0,ni+1
        va(i,j,k)=v3d(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()
 
#ifdef MPI
      call comm_3w_end(w3d,ww31,ww32,we31,we32,   &
                           ws31,ws32,wn31,wn32,reqs_w)
#endif
      if(terrain_flag)then
        call bcwsfc(dzdx,dzdy,u3d,v3d,w3d)
        call bc2d(w3d(ib,jb,1))
      endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+2
      do j=0,nj+1
      do i=0,ni+1
        wa(i,j,k)=w3d(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()
 
#ifndef MPI
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        ppi(i,j,k)=pp3d(i,j,k)
      enddo
      enddo
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        tha(i,j,k)=th3d(i,j,k)
      enddo
      enddo
      enddo

      if(imoist.eq.1)then
 
        do n=1,numq

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=0,nk+1
          do j=0,nj+1
          do i=0,ni+1
            qa(i,j,k,n)=q3d(i,j,k,n)
          enddo
          enddo
          enddo

        enddo
 
      endif

      if(timestats.ge.1) time_integ=time_integ+mytime()
#else
    if(imoist.eq.0)then
      call fincomm(tha,th3d,sw31,sw32,se31,se32,   &
                            ss31,ss32,sn31,sn32,reqs_s)
      call fincomm(ppi,pp3d,pw31,pw32,pe31,pe32,   &
                            ps31,ps32,pn31,pn32,reqs_p)
    endif
#endif

!---- finish communication of moisture ----

#ifdef MPI
  IF(imoist.eq.1)THEN
    IF(simple_comm)THEN
      call fincomm(tha,th3d,sw31,sw32,se31,se32,   &
                            ss31,ss32,sn31,sn32,reqs_s)
      call fincomm(ppi,pp3d,pw31,pw32,pe31,pe32,   &
                            ps31,ps32,pn31,pn32,reqs_p)
      do n=1,numq
        call fincomm(qa(ib,jb,kb,n),q3d(ib,jb,kb,n),                        &
                     qw1(1,1,1,n),qw2(1,1,1,n),qe1(1,1,1,n),qe2(1,1,1,n),   &
                     qs1(1,1,1,n),qs2(1,1,1,n),qn1(1,1,1,n),qn2(1,1,1,n),   &
                     reqs_q(1,n))
      enddo
    ELSE
      if(ptype.eq.1)then
        call fincomm(qa(ib,jb,kb,3),q3d(ib,jb,kb,3),                        &
                     qw1(1,1,1,3),qw2(1,1,1,3),qe1(1,1,1,3),qe2(1,1,1,3),   &
                     qs1(1,1,1,3),qs2(1,1,1,3),qn1(1,1,1,3),qn2(1,1,1,3),   &
                     reqs_q(1,3))
        call fincomm(qa(ib,jb,kb,1),q3d(ib,jb,kb,1),                        &
                     qw1(1,1,1,1),qw2(1,1,1,1),qe1(1,1,1,1),qe2(1,1,1,1),   &
                     qs1(1,1,1,1),qs2(1,1,1,1),qn1(1,1,1,1),qn2(1,1,1,1),   &
                     reqs_q(1,1))
        call fincomm(qa(ib,jb,kb,2),q3d(ib,jb,kb,2),                        &
                     qw1(1,1,1,2),qw2(1,1,1,2),qe1(1,1,1,2),qe2(1,1,1,2),   &
                     qs1(1,1,1,2),qs2(1,1,1,2),qn1(1,1,1,2),qn2(1,1,1,2),   &
                     reqs_q(1,2))
        call fincomm(tha,th3d,sw31,sw32,se31,se32,   &
                              ss31,ss32,sn31,sn32,reqs_s)
        call fincomm(ppi,pp3d,pw31,pw32,pe31,pe32,   &
                              ps31,ps32,pn31,pn32,reqs_p)
      elseif(ptype.eq.2)then
        call fincomm(qa(ib,jb,kb,1),q3d(ib,jb,kb,1),                        &
                     qw1(1,1,1,1),qw2(1,1,1,1),qe1(1,1,1,1),qe2(1,1,1,1),   &
                     qs1(1,1,1,1),qs2(1,1,1,1),qn1(1,1,1,1),qn2(1,1,1,1),   &
                     reqs_q(1,1))
        call fincomm(qa(ib,jb,kb,2),q3d(ib,jb,kb,2),                        &
                     qw1(1,1,1,2),qw2(1,1,1,2),qe1(1,1,1,2),qe2(1,1,1,2),   &
                     qs1(1,1,1,2),qs2(1,1,1,2),qn1(1,1,1,2),qn2(1,1,1,2),   &
                     reqs_q(1,2))
        call fincomm(tha,th3d,sw31,sw32,se31,se32,   &
                              ss31,ss32,sn31,sn32,reqs_s)
        call fincomm(ppi,pp3d,pw31,pw32,pe31,pe32,   &
                              ps31,ps32,pn31,pn32,reqs_p)
        call fincomm(qa(ib,jb,kb,3),q3d(ib,jb,kb,3),                        &
                     qw1(1,1,1,3),qw2(1,1,1,3),qe1(1,1,1,3),qe2(1,1,1,3),   &
                     qs1(1,1,1,3),qs2(1,1,1,3),qn1(1,1,1,3),qn2(1,1,1,3),   &
                     reqs_q(1,3))
        call fincomm(qa(ib,jb,kb,4),q3d(ib,jb,kb,4),                        &
                     qw1(1,1,1,4),qw2(1,1,1,4),qe1(1,1,1,4),qe2(1,1,1,4),   &
                     qs1(1,1,1,4),qs2(1,1,1,4),qn1(1,1,1,4),qn2(1,1,1,4),   &
                     reqs_q(1,4))
        call fincomm(qa(ib,jb,kb,5),q3d(ib,jb,kb,5),                        &
                     qw1(1,1,1,5),qw2(1,1,1,5),qe1(1,1,1,5),qe2(1,1,1,5),   &
                     qs1(1,1,1,5),qs2(1,1,1,5),qn1(1,1,1,5),qn2(1,1,1,5),   &
                     reqs_q(1,5))
        call fincomm(qa(ib,jb,kb,6),q3d(ib,jb,kb,6),                        &
                     qw1(1,1,1,6),qw2(1,1,1,6),qe1(1,1,1,6),qe2(1,1,1,6),   &
                     qs1(1,1,1,6),qs2(1,1,1,6),qn1(1,1,1,6),qn2(1,1,1,6),   &
                     reqs_q(1,6))
      endif
    ENDIF

  ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        ppi(i,j,k)=pp3d(i,j,k)
      enddo
      enddo
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        tha(i,j,k)=th3d(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()

  ENDIF
#endif

      if(imove.eq.1.and.imoist.eq.1)then
        call movesfc(0.0,dt,uh,vh,rain(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3))
        if(timestats.ge.1) time_swath=time_swath+mytime()
      endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   All done   cccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!!#ifdef MPI
!!!      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!!!      if(timestats.ge.1) time_mpb=time_mpb+mytime()
!!!#endif

!  Calculate surface "swaths."  Move surface (if necessary). 

    IF( output_sws.eq.1 )THEN

!--------------------------------------------------------------------
! Maximum horizontal wind speed at lowest model level: 
! (include domain movement in calculation)

!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = sqrt( (umove+0.5*(ua(i,j,1)+ua(i+1,j,1)))**2    &
                   +(vmove+0.5*(va(i,j,1)+va(i,j+1,1)))**2 ) 
        do n=1,nrain
          sws(i,j,n)=max(sws(i,j,n),tem)
        enddo
      enddo
      enddo

      if(imove.eq.1)then
        call movesfc(0.0,dt,uh,vh,sws(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3))
      endif

!--------------------------------------------------------------------
!  Maximum vertical vorticity at lowest model level:

!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj+1
      do i=1,ni+1
        tem = (va(i,j,1)-va(i-1,j,1))*rdx*uf(i)   &
             -(ua(i,j,1)-ua(i,j-1,1))*rdy*vf(j)
        do n=1,nrain
          svs(i,j,n)=max(svs(i,j,n),tem)
        enddo
      enddo
      enddo

      if(imove.eq.1)then
        call movesfc(-200000.0,dt,uh,vh,svs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3))
      endif

!--------------------------------------------------------------------
!  Minimum pressure perturbation at lowest model level:

!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = prs(i,j,1)-prs0(i,j,1)
        do n=1,nrain
          sps(i,j,n)=min(sps(i,j,n),tem)
        enddo
      enddo
      enddo

      if(imove.eq.1)then
        call movesfc(-200000.0,dt,uh,vh,sps(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3))
      endif

!--------------------------------------------------------------------
!  Maximum rainwater mixing ratio (qr) at lowest model level:

    IF(imoist.eq.1.and.nqr.ne.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = qa(i,j,1,nqr)
        do n=1,nrain
          srs(i,j,n)=max(srs(i,j,n),tem)
        enddo
      enddo
      enddo

      if(imove.eq.1)then
        call movesfc(0.0,dt,uh,vh,srs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3))
      endif
    ENDIF

!--------------------------------------------------------------------
!  Maximum graupel/hail mixing ratio (qg) at lowest model level:

    IF(imoist.eq.1.and.nqg.ne.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = qa(i,j,1,nqg)
        do n=1,nrain
          sgs(i,j,n)=max(sgs(i,j,n),tem)
        enddo
      enddo
      enddo

      if(imove.eq.1)then
        call movesfc(0.0,dt,uh,vh,sgs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3))
      endif
    ENDIF

!--------------------------------------------------------------------

      ! get height AGL:
      if( terrain_flag )then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dum3(i,j,k) = zh(i,j,k)-zs(i,j)
          wten(i,j,k) = zf(i,j,k)-zs(i,j)
        enddo
        enddo
        enddo
      else
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dum3(i,j,k) = zh(i,j,k)
          wten(i,j,k) = zf(i,j,k)
        enddo
        enddo
        enddo
      endif

!--------------------------------------------------------------------
!  Maximum updraft velocity (w) at 5 km AGL:

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,tem)
      do j=1,nj
      do i=1,ni
        k = 2
        ! wten is height AGL:
        do while( wten(i,j,k).lt.5000.0 )
          k = k + 1
        enddo
        tem = w3d(i,j,k)
        do n=1,nrain
          sus(i,j,n)=max(sus(i,j,n),tem)
        enddo
      enddo
      enddo

      if(imove.eq.1)then
        call movesfc(-200000.0,dt,uh,vh,sus(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3))
      endif

!--------------------------------------------------------------------
!  Maximum integrated updraft helicity:

      ! dum3 is zh (agl), wten is zf (agl)
      call calcuh(uf,vf,dum3,wten,ua,va,wa,dum1(ib,jb,1),dum2)
!$omp parallel do default(shared)  &
!$omp private(i,j,n)
      do j=1,nj
      do i=1,ni
        do n=1,nrain
          shs(i,j,n)=max(shs(i,j,n),dum1(i,j,1))
        enddo
      enddo
      enddo

      if(imove.eq.1)then
        call movesfc(0.0,dt,uh,vh,shs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3))
      endif

      if(timestats.ge.1) time_swath=time_swath+mytime()
    ENDIF

!  Done with "swaths"
!--------------------------------------------------------------------
!  Step time forward, Get statistics

      mtime = mtime + dt

      if( convinit.eq.1 )then
        if( mtime.gt.convtime ) convinit = 0
      endif

      rtime=sngl(mtime)
      if( rtime.ge.sngl(stattim) .or. statfrq.le.0.0 )then
        IF(axisymm.eq.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ppten(i,j,k)=rho(i,j,k)
          enddo
          enddo
          enddo
        ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ppten(i,j,k) = rho(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
          enddo
          enddo
          enddo
        ENDIF
        call statpack(nrec,ndt,dt,rtime,adt,acfl,cloudvar,qname,budname,qbudget,asq,bsq, &
                      xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,vf,zh,mh,rmh,mf,    &
                      rstat,pi0,rho0,thv0,th0,qv0,u0,v0,                &
                      dum1,dum2,dum3,dum4,divx,ppten,prs,               &
                      ua,va,wa,ppi,tha,qa,qten,kmh,kmv,khh,khv,tkea,pta,u10,v10)
        stattim=stattim+statfrq
      else
        if( adapt_dt.eq.1 ) call calccfl(1,rstat,dt,acfl,uf,vf,mf,ua,va,wa,0)
      endif

!--------------------------------------------------------------------
!  Writeout and stuff

      rtime=sngl(mtime)
      if(myid.eq.0)then
        if(timeformat.eq.1)then
          write(6,110) nstep,rtime,' sec '
        elseif(timeformat.eq.2)then
          write(6,110) nstep,rtime/60.0,' min '
        elseif(timeformat.eq.3)then
          write(6,110) nstep,rtime/3600.0,' hour'
        elseif(timeformat.eq.4)then
          write(6,110) nstep,rtime/86400.0,' day '
        else
          write(6,110) nstep,rtime,' sec'
        endif
110     format(2x,i12,4x,f18.6,a5)
      endif
      if(timestats.ge.1) time_misc=time_misc+mytime()

      if( (rtime.ge.sngl(taptim)).or.stopit )then
        nwrite=nwrite+1
      IF(output_format.eq.1)THEN
        nn = 1
        if(terrain_flag .and. output_interp.eq.1) nn = 2
        DO n=1,nn
          if(n.eq.1)then
            fnum = 51
          else
            fnum = 71
          endif
          call writeout(fnum,nwrite,qname,xh,xf,uf,vf,sigma,zh,zf,mf,pi0,prs0,rho0,th0,qv0,u0,v0,  &
                      zs,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,cdu,cdv,ce,dum1,dum2,dum3,dum4, &
                      rho,prs,sten,ua,uten,va,vten,wa,wten,ppi,tha,          &
                   dissten,thpten,qvpten,qcpten,qipten,upten,vpten,          &
                      lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,   &
                      qa,kmh,kmv,khh,khv,tkea,swten,lwten,radsw,rnflx,radswnet,radlwin,pta,   &
                      num_soil_layers,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br)
        ENDDO
#ifdef NETCDF
      ELSEIF(output_format.eq.2)THEN
        call writeout_cdf(nwrite,qname,sigma,sigmaf,xh,xf,uf,yh,yf,vf,mh,zh,mf,zf, &
                      pi0,prs0,rho0,th0,qv0,u0,v0,                     &
                      zs,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,cdu,cdv,ce,dum1,dum2,dum3,dum4,  &
                      rho,prs,sten,ua,uten,va,vten,wa,wten,ppi,tha,    &
                      qa,kmh,kmv,khh,khv,tkea,pta,num_soil_layers,   &
                      lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,   &
                      radsw,rnflx,radswnet,radlwin,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br,   &
                      dissten,thpten,qvpten,qcpten,qipten,upten,vpten,swten,lwten)
#endif
#ifdef HDFOUT
      ELSEIF(output_format.ge.3)THEN
        call writeout_mult_hdf5(rtime,qname,sigma,sigmaf,xh,xf,uf,yh,yf,vf,mh,zh,mf,zf, &
                      pi0,prs0,rho0,th0,qv0,u0,v0,                     &
                      zs,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,cdu,cdv,ce,dum1,dum2,dum3,dum4,  &
                      rho,prs,sten,ua,uten,va,vten,wa,wten,ppi,tha,    &
                      qa,kmh,kmv,khh,khv,tkea,pta,num_soil_layers,   &
                      lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,   &
                      radsw,rnflx,radswnet,radlwin,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br,   &
                      dissten,thpten,qvpten,qcpten,qipten,upten,vpten,swten,lwten)
! See above note about sten being used for dbz calculation
#endif
      ENDIF
        taptim=taptim+tapfrq
        if(timestats.ge.1) time_write=time_write+mytime()
      endif

      if(rtime.ge.rsttim .and. rstfrq.gt.0)then
        nrst=nrst+1
        call write_restart(nstep,nrec,prec,nwrite,nrst,nrad2d,num_soil_layers, &
                               dt,mtime,radtim,qbudget,asq,bsq,                &
                               rain,sws,svs,sps,srs,sgs,sus,shs,tsk,radbcw,radbce,radbcs,radbcn,     &
                               ua,va,wa,ppi,tha,qa,tkea,swten,lwten,   &
                               radsw,rnflx,radswnet,radlwin,rad2d,   &
                               lu_index,kpbl2d,psfc,u10,v10,hfx,qfx,xland,znt,ust, &
                               hpbl,wspd,psim,psih,gz1oz0,br,                      &
                               CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                      &
                               MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                   &
                               CK,CKA,CD,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                               f2d,gsw,glw,chklowq,capg,snowc,tslb,                &
                               tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,              &
                               pta,pdata,rtime)
        rsttim=rsttim+rstfrq
        if(timestats.ge.1) time_write=time_write+mytime()
      endif

!-------------------------------------------------------------------
!  Parcel update (final time step)

      if( (iprcl.eq.1).and.(rtime.ge.timax) )then
        write(outfile,*) '  Calling parcel driver for last step'
        call parcel_driver(prec,dt,xh,uh,ruh,yh,vh,rvh,zh,mh,rmh,mf,        &
                           pi0,thv0,th0,dum1,dum2,dum3,dum4,divx,prs,    &
                           ua,va,wa,ppi,thten,tha,qa,khv,pdata,rtime,    &
                           ploc,packet,reqs_p,                           &
                           pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2)
      endif

!-------------------------------------------------------------------

      if(stopit)then
        if(myid.eq.0)then
          print *
          print *,' Courant number has exceeded 1.5 '
          print *
          print *,' Stopping model .... '
          print *
        endif
        call stopcm1
      endif

!--------------------------------------------------------------------

      return
      end


