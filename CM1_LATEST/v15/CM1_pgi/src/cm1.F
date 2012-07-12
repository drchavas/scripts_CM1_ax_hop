

      program cm1
      implicit none

!-----------------------------------------------------------------------
!  CM1 Numerical Model, Release 15  (cm1r15)
!  13 January 2011
!  http://www.mmm.ucar.edu/people/bryan/cm1/
!-----------------------------------------------------------------------
!
!  Please see documentation at the top of the "solve.F" file.
!
!  See also documentation at the cm1 website, such as:
!
!    "The governing equations for CM1"
!        http://www.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf
!
!-----------------------------------------------------------------------

      include 'input.incl'
      include 'radcst.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      integer :: nstep
      integer :: nrec,prec,nwrite,nrst
      integer :: rbufsz,num_soil_layers,ndt
      real :: dt,dtlast
      real*8 :: mtime,stattim,taptim,rsttim,radtim,adt,acfl
      logical :: dodrag,dosfcflx
      logical, dimension(maxq) :: cloudvar,rhovar
      character*15 :: tdef
      character*3, dimension(maxq) :: qname
      character*6, dimension(maxq) :: budname
      real*8, dimension(:), allocatable :: bud,bud2
      real*8, dimension(:), allocatable :: qbudget
      real*8, dimension(:), allocatable :: asq,bsq
      real, dimension(:), allocatable :: xh,rxh,uh,ruh
      real, dimension(:), allocatable :: xf,rxf,uf,ruf
      real, dimension(:), allocatable :: yh,vh,rvh
      real, dimension(:), allocatable :: yf,vf,rvf
      real, dimension(:), allocatable :: xfref,yfref
      real, dimension(:), allocatable :: sigma,sigmaf
      real, dimension(:,:,:), allocatable :: tauh,taus,zh,mh,rmh
      real, dimension(:,:,:), allocatable :: tauf,zf,mf,rmf
      real, dimension(:), allocatable :: rstat
      real, dimension(:,:), allocatable :: rho0s,pi0s,prs0s,rth0s
      real, dimension(:,:,:), allocatable :: pi0,rho0,prs0,thv0,th0,qv0
      real, dimension(:,:,:), allocatable :: ql0,rr0,rf0,rrf0,rru0,rrv0,u0,v0
      real, dimension(:,:,:), allocatable :: t0,rh0,qc0
      real, dimension(:,:), allocatable :: zs,gz,dzdx,dzdy
      real, dimension(:,:,:), allocatable :: gx,gy
      real, dimension(:,:,:), allocatable :: rain,sws,svs,sps,srs,sgs,sus,shs
      logical, dimension(:,:), allocatable :: doimpl
      real, dimension(:,:), allocatable :: tsk,thflux,qvflux,cdu,cdv,ce,u1,v1,w1
      real, dimension(:,:), allocatable :: radbcw,radbce
      real, dimension(:,:), allocatable :: radbcs,radbcn
      real, dimension(:,:,:), allocatable :: dum1,dum2,dum3,dum4
      real, dimension(:,:,:), allocatable :: divx,rho,prs
      real, dimension(:,:,:), allocatable :: t11,t12,t13,t22,t23,t33
      real, dimension(:,:,:), allocatable :: rru,ua,u3d,uten,uten1
      real, dimension(:,:,:), allocatable :: rrv,va,v3d,vten,vten1
      real, dimension(:,:,:), allocatable :: rrw,wa,w3d,wten,wten1
      real, dimension(:,:,:), allocatable :: ppi,pp3d,ppten,sten
      real, dimension(:,:,:), allocatable :: tha,th3d,thten,thten1
      real, dimension(:,:,:), allocatable :: thterm,tk
      real, dimension(:,:,:,:), allocatable :: qa,q3d,qten
      real, dimension(:,:,:,:), allocatable :: zvdarray
      real, dimension(:,:,:), allocatable :: kmh,kmv,khh,khv
      real, dimension(:,:,:), allocatable :: tkea,tke3d,tketen
      real, dimension(:,:,:), allocatable :: dissten
      real, dimension(:,:,:), allocatable :: thpten,qvpten,qcpten,qipten,upten,vpten
      real, dimension(:,:,:), allocatable :: swten,lwten,o30
      real, dimension(:,:), allocatable :: radsw,rnflx,radswnet,radlwin
      real, dimension(:,:,:), allocatable :: rad2d
      real, dimension(:), allocatable :: x,y,z,za
      real, dimension(:,:,:), allocatable :: zp
      integer, dimension(:,:), allocatable :: lu_index,kpbl2d
      real, dimension(:,:), allocatable :: psfc,u10,v10,hfx,qfx,xland,znt,ust, &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,       &
                                      CK,CKA,CD,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      f2d,gsw,glw,chklowq,capg,snowc,dsxy
      real, dimension(:), allocatable :: slab_zs,slab_dzs
      real, dimension(:,:,:), allocatable :: tslb
      real, dimension(:,:), allocatable :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, dimension(:,:,:,:),  allocatable :: pta,pt3d,ptten
      real, dimension(:,:),      allocatable :: pdata
      real, dimension(:,:,:),    allocatable :: cfb
      real, dimension(:),        allocatable :: cfa,cfc,d1,d2
      complex, dimension(:,:,:), allocatable :: pdt,deft
      complex, dimension(:,:),   allocatable :: rhs,trans

!--- arrays for MPI ---
      integer, dimension(:), allocatable :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_tk
      integer, dimension(:,:),  allocatable :: reqs_q,reqs_t
      real, dimension(:,:), allocatable :: ww1,ww2,we1,we2
      real, dimension(:,:), allocatable :: ws1,ws2,wn1,wn2
      real, dimension(:,:), allocatable :: pw1,pw2,pe1,pe2
      real, dimension(:,:), allocatable :: ps1,ps2,pn1,pn2
      real, dimension(:,:), allocatable :: vw1,vw2,ve1,ve2
      real, dimension(:,:), allocatable :: vs1,vs2,vn1,vn2
      real, dimension(:,:,:), allocatable :: uw31,uw32,ue31,ue32
      real, dimension(:,:,:), allocatable :: us31,us32,un31,un32
      real, dimension(:,:,:), allocatable :: vw31,vw32,ve31,ve32
      real, dimension(:,:,:), allocatable :: vs31,vs32,vn31,vn32
      real, dimension(:,:,:), allocatable :: ww31,ww32,we31,we32
      real, dimension(:,:,:), allocatable :: ws31,ws32,wn31,wn32
      real, dimension(:,:,:), allocatable :: sw31,sw32,se31,se32
      real, dimension(:,:,:), allocatable :: ss31,ss32,sn31,sn32
      real, dimension(:,:,:), allocatable :: pw31,pw32,pe31,pe32
      real, dimension(:,:,:), allocatable :: ps31,ps32,pn31,pn32
      real, dimension(:,:,:), allocatable :: tkw1,tkw2,tke1,tke2
      real, dimension(:,:,:), allocatable :: tks1,tks2,tkn1,tkn2
      real, dimension(:,:,:,:), allocatable :: qw1,qw2,qe1,qe2
      real, dimension(:,:,:,:), allocatable :: qs1,qs2,qn1,qn2
      real, dimension(:,:,:,:), allocatable :: tw1,tw2,te1,te2
      real, dimension(:,:,:,:), allocatable :: ts1,ts2,tn1,tn2
      real, dimension(:,:), allocatable :: ploc,packet

!-----

      integer count,rate,maxr
      real rtime,xtime,time_solve
      real steptime1,steptime2
      integer :: i,j,k,n,nn,fnum
#ifdef MPI
      integer rc
      real mp_total
      real*8 tstart,tend
#endif

      namelist /param0/ nx,ny,nz,nodex,nodey,timeformat,timestats,terrain_flag

!----------------------------------------------------------------------

      nstep = 0
      mtime = 0.0d0
      nrec=1
      prec=1
      nwrite=1
      nrst=0
      outfile=6
      stopit = .false.
      smeps = 1.0e-30
      tsmall = 0.0001
#ifdef DP
      smeps = 1.0e-60
      tsmall = 0.000000000001
#endif

!----------------------------------------------------------------------
!  Initialize MPI

      myid=0
      numprocs=1

#ifdef MPI
      outfile=10
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
#endif

!----------------------------------------------------------------------
!  Get domain dimensions, allocate some arrays, then call PARAM

      open(unit=20,file='namelist.input',form='formatted',status='old',    &
           access='sequential')
      read(20,nml=param0)
      close(unit=20)

      ni = nx / nodex
      nj = ny / nodey
      nk = nz
      nkp1 = nk+1

      ! (The following are needed by ZVD, but are also included for future 
      !  development, e.g., possible distributed-memory decomposition in 
      !  z direction)
      !
      ! number of 'ghost' points in the horizontal directions:
      ngxy  = 3
      ! number of 'ghost' points in the vertical direction:
      ngz   = 1

!---------------------------------------------------------------------
!      For ZVD:
!      ngz   = 3
!      IF( ngz.eq.3 )THEN
!        kb =  1 - ngz
!        ke = nk + ngz
!      ENDIF
!---------------------------------------------------------------------

      ib =  1 - ngxy
      ie = ni + ngxy
      jb =  1 - ngxy
      je = nj + ngxy
      kb =  1 - ngz
      ke = nk + ngz

      allocate(    xh(ib:ie) )
      allocate(   rxh(ib:ie) )
      allocate(    uh(ib:ie) )
      allocate(   ruh(ib:ie) )
      allocate(    xf(ib:ie+1) )
      allocate(   rxf(ib:ie+1) )
      allocate(    uf(ib:ie+1) )
      allocate(   ruf(ib:ie+1) )
      allocate(    yh(jb:je) )
      allocate(    vh(jb:je) )
      allocate(   rvh(jb:je) )
      allocate(    yf(jb:je+1) )
      allocate(    vf(jb:je+1) )
      allocate(   rvf(jb:je+1) )
      allocate( xfref(-2:nx+4) )
      allocate( yfref(-2:ny+4) )
      allocate( sigma(kb:ke) )
      allocate( sigmaf(kb:ke+1) )
      allocate(  tauh(ib:ie,jb:je,kb:ke) )
      allocate(  taus(ib:ie,jb:je,kb:ke) )
      allocate(    zh(ib:ie,jb:je,kb:ke) )
      allocate(    mh(ib:ie,jb:je,kb:ke) )
      allocate(   rmh(ib:ie,jb:je,kb:ke) )
      allocate(  tauf(ib:ie,jb:je,kb:ke+1) )
      allocate(    mf(ib:ie,jb:je,kb:ke+1) )
      allocate(   rmf(ib:ie,jb:je,kb:ke+1) )

      if(terrain_flag)then
        itb=ib
        ite=ie
        jtb=jb
        jte=je
        ktb=kb
        kte=ke
      else
        itb=1
        ite=1
        jtb=1
        jte=1
        ktb=1
        kte=1
      endif

      allocate(   zs(itb:ite,jtb:jte) )
      allocate(   gz(itb:ite,jtb:jte) )
      allocate( dzdx(itb:ite,jtb:jte) )
      allocate( dzdy(itb:ite,jtb:jte) )
      allocate(   gx(itb:ite+1,jtb:jte,ktb:kte) )
      allocate(   gy(itb:ite,jtb:jte+1,ktb:kte) )
      allocate(   zf(ib:ie,jb:je,kb:ke+1) )

      call param(dt,dtlast,stattim,taptim,rsttim,radtim,          &
                 cloudvar,rhovar,qname,budname,                   &
                 xh,rxh,uh,ruh,xf,rxf,uf,ruf,yh,vh,rvh,yf,vf,rvf, &
                 xfref,yfref,                                     &
                 sigma,sigmaf,tauh,taus,zh,mh,rmh,tauf,zf,mf,rmf, &
                 zs,gz,dzdx,dzdy,gx,gy)

!----------------------------------------------------------------------
!  allocate the base state arrays, then call BASE

      allocate( rstat(stat_out) )
      allocate( rho0s(ib:ie,jb:je) )
      allocate(  pi0s(ib:ie,jb:je) )
      allocate( prs0s(ib:ie,jb:je) )
      allocate( rth0s(ib:ie,jb:je) )
      allocate(  pi0(ib:ie,jb:je,kb:ke) )
      allocate( rho0(ib:ie,jb:je,kb:ke) )
      allocate( prs0(ib:ie,jb:je,kb:ke) )
      allocate( thv0(ib:ie,jb:je,kb:ke) )
      allocate(  th0(ib:ie,jb:je,kb:ke) )
      allocate(  qv0(ib:ie,jb:je,kb:ke) )
      allocate(  ql0(ib:ie,jb:je,kb:ke) )
      allocate(  rr0(ib:ie,jb:je,kb:ke) )
      allocate(  rf0(ib:ie,jb:je,kb:ke) )
      allocate( rrf0(ib:ie,jb:je,kb:ke) )
      allocate( rru0(ib:ie,jb:je,kb:ke) )
      allocate( rrv0(ib:ie,jb:je,kb:ke) )
      allocate(   u0(ib:ie+1,jb:je,kb:ke) )
      allocate(   v0(ib:ie,jb:je+1,kb:ke) )

      allocate(   t0(ib:ie,jb:je,kb:ke) )
      allocate(  rh0(ib:ie,jb:je,kb:ke) )
      allocate(  qc0(ib:ie,jb:je,kb:ke) )

      call base(zh,mh,zf,mf,rho0s,pi0s,prs0s,rth0s,pi0,prs0,rho0,thv0,th0,t0,qv0,u0,v0,rh0,    &
                qc0,ql0,rr0,rf0,rrf0,rru0,rrv0)

!----------------------------------------------------------------------
!  Now, allocate the mother lode, then call INIT3D

      allocate(   rain(ib:ie,jb:je,nrain) )
      allocate(    sws(ib:ie,jb:je,nrain) )
      allocate(    svs(ib:ie,jb:je,nrain) )
      allocate(    sps(ib:ie,jb:je,nrain) )
      allocate(    srs(ib:ie,jb:je,nrain) )
      allocate(    sgs(ib:ie,jb:je,nrain) )
      allocate(    sus(ib:ie,jb:je,nrain) )
      allocate(    shs(ib:ie,jb:je,nrain) )

      allocate( doimpl(ib:ie,jb:je) )

      allocate(    tsk(ib:ie,jb:je) )
      allocate( thflux(ib:ie,jb:je) )
      allocate( qvflux(ib:ie,jb:je) )
      allocate(    cdu(ib:ie,jb:je) )
      allocate(    cdv(ib:ie,jb:je) )
      allocate(     ce(ib:ie,jb:je) )
      allocate(     u1(ib:ie,jb:je) )
      allocate(     v1(ib:ie,jb:je) )
      allocate(     w1(ib:ie,jb:je) )

      allocate( radbcw(jb:je,kb:ke) )
      allocate( radbce(jb:je,kb:ke) )
      allocate( radbcs(ib:ie,kb:ke) )
      allocate( radbcn(ib:ie,kb:ke) )

      allocate( dum1(ib:ie,jb:je,kb:ke) )
      allocate( dum2(ib:ie,jb:je,kb:ke) )
      allocate( dum3(ib:ie,jb:je,kb:ke) )
      allocate( dum4(ib:ie,jb:je,kb:ke) )
      allocate( divx(ib:ie,jb:je,kb:ke) )
      allocate(  rho(ib:ie,jb:je,kb:ke) )
      allocate(  prs(ib:ie,jb:je,kb:ke) )
      allocate(  t11(ib:ie,jb:je,kb:ke) )
      allocate(  t12(ib:ie,jb:je,kb:ke) )
      allocate(  t13(ib:ie,jb:je,kb:ke) )
      allocate(  t22(ib:ie,jb:je,kb:ke) )
      allocate(  t23(ib:ie,jb:je,kb:ke) )
      allocate(  t33(ib:ie,jb:je,kb:ke) )

      allocate(   rru(ib:ie+1,jb:je,kb:ke) )
      allocate(    ua(ib:ie+1,jb:je,kb:ke) )
      allocate(   u3d(ib:ie+1,jb:je,kb:ke) )
      allocate(  uten(ib:ie+1,jb:je,kb:ke) )
      allocate( uten1(ib:ie+1,jb:je,kb:ke) )
      allocate(   rrv(ib:ie,jb:je+1,kb:ke) )
      allocate(    va(ib:ie,jb:je+1,kb:ke) )
      allocate(   v3d(ib:ie,jb:je+1,kb:ke) )
      allocate(  vten(ib:ie,jb:je+1,kb:ke) )
      allocate( vten1(ib:ie,jb:je+1,kb:ke) )
      allocate(   rrw(ib:ie,jb:je,kb:ke+1) )
      allocate(    wa(ib:ie,jb:je,kb:ke+1) )
      allocate(   w3d(ib:ie,jb:je,kb:ke+1) )
      allocate(  wten(ib:ie,jb:je,kb:ke+1) )
      allocate( wten1(ib:ie,jb:je,kb:ke+1) )

      allocate(   ppi(ib:ie,jb:je,kb:ke) )
      allocate(  pp3d(ib:ie,jb:je,kb:ke) )
      allocate( ppten(ib:ie,jb:je,kb:ke) )
      allocate(  sten(ib:ie,jb:je,kb:ke) )
      allocate(   tha(ib:ie,jb:je,kb:ke) )
      allocate(  th3d(ib:ie,jb:je,kb:ke) )
      allocate( thten(ib:ie,jb:je,kb:ke) )
      allocate(thten1(ib:ie,jb:je,kb:ke) )
      allocate(thterm(ib:ie,jb:je,kb:ke) )
      allocate(    tk(ib:ie,jb:je,kb:ke) )

      allocate(   bud(nk) )
      allocate(  bud2(nj) )
      allocate( qbudget(nbudget) )
      allocate(    asq(numq) )
      allocate(    bsq(numq) )
      allocate(     qa(ibm:iem,jbm:jem,kbm:kem,numq) )
      allocate(    q3d(ibm:iem,jbm:jem,kbm:kem,numq) )
      allocate(   qten(ibm:iem,jbm:jem,kbm:kem,numq) )
      allocate( zvdarray(ibzvd:iezvd,jbzvd:jezvd,kbzvd:kezvd,nqzvd) )
      allocate(    kmh(ibc:iec,jbc:jec,kbc:kec) )
      allocate(    kmv(ibc:iec,jbc:jec,kbc:kec) )
      allocate(    khh(ibc:iec,jbc:jec,kbc:kec) )
      allocate(    khv(ibc:iec,jbc:jec,kbc:kec) )
      allocate(   tkea(ibt:iet,jbt:jet,kbt:ket) )
      allocate(  tke3d(ibt:iet,jbt:jet,kbt:ket) )
      allocate( tketen(ibt:iet,jbt:jet,kbt:ket) )

      allocate( dissten(ib:ie,jb:je,kb:ke) )

      allocate( thpten(ibb:ieb,jbb:jeb,kbb:keb) )
      allocate( qvpten(ibb:ieb,jbb:jeb,kbb:keb) )
      allocate( qcpten(ibb:ieb,jbb:jeb,kbb:keb) )
      allocate( qipten(ibb:ieb,jbb:jeb,kbb:keb) )
      allocate(  upten(ibb:ieb,jbb:jeb,kbb:keb) )
      allocate(  vpten(ibb:ieb,jbb:jeb,kbb:keb) )

      allocate( swten(ibr:ier,jbr:jer,kbr:ker) )
      allocate( lwten(ibr:ier,jbr:jer,kbr:ker) )
      allocate(   o30(ibr:ier,jbr:jer,kbr:ker) )

      nir = 1
      njr = 1
      nkr = nk+3

      IF( radopt .eq. 1 )THEN
        rbufsz = n2d_radiat*nir*njr + n3d_radiat*nir*njr*nkr
      ELSE
        rbufsz = 1
      ENDIF

      allocate(    rad2d(ni,nj,nrad2d) )
      allocate(    radsw(ni,nj) )
      allocate(    rnflx(ni,nj) )
      allocate( radswnet(ni,nj) )
      allocate(  radlwin(ni,nj) )

      rad2d = 0.0
      radsw = 0.0
      rnflx = 0.0
      radswnet = 0.0
      radlwin = 0.0

      write(outfile,*) '  rbufsz,nrad2d = ',rbufsz,nrad2d

      allocate( x(ni+1) )
      allocate( y(nj+1) )
      allocate( z(nk+3) )
      allocate( za(nk+3) )
      allocate( zp(ni,nj,nk+3) )

      allocate( lu_index(ibl:iel,jbl:jel) )
      allocate(   kpbl2d(ibl:iel,jbl:jel) )
      allocate(     psfc(ibl:iel,jbl:jel) )
      allocate(      u10(ibl:iel,jbl:jel) )
      allocate(      v10(ibl:iel,jbl:jel) )
      allocate(      hfx(ibl:iel,jbl:jel) )
      allocate(      qfx(ibl:iel,jbl:jel) )
      allocate(    xland(ibl:iel,jbl:jel) )
      allocate(      znt(ibl:iel,jbl:jel) )
      allocate(      ust(ibl:iel,jbl:jel) )
      allocate(     hpbl(ibl:iel,jbl:jel) )
      allocate(     wspd(ibl:iel,jbl:jel) )
      allocate(     psim(ibl:iel,jbl:jel) )
      allocate(     psih(ibl:iel,jbl:jel) )
      allocate(   gz1oz0(ibl:iel,jbl:jel) )
      allocate(       br(ibl:iel,jbl:jel) )
      allocate(      chs(ibl:iel,jbl:jel) )
      allocate(     chs2(ibl:iel,jbl:jel) )
      allocate(     cqs2(ibl:iel,jbl:jel) )
      allocate(     cpmm(ibl:iel,jbl:jel) )
      allocate(      zol(ibl:iel,jbl:jel) )
      allocate(   mavail(ibl:iel,jbl:jel) )
      allocate(      mol(ibl:iel,jbl:jel) )
      allocate(     rmol(ibl:iel,jbl:jel) )
      allocate(   regime(ibl:iel,jbl:jel) )
      allocate(       lh(ibl:iel,jbl:jel) )
      allocate(     flhc(ibl:iel,jbl:jel) )
      allocate(     flqc(ibl:iel,jbl:jel) )
      allocate(      qgh(ibl:iel,jbl:jel) )
      allocate(       ck(ibl:iel,jbl:jel) )
      allocate(      cka(ibl:iel,jbl:jel) )
      allocate(       cd(ibl:iel,jbl:jel) )
      allocate(      cda(ibl:iel,jbl:jel) )
      allocate(     ustm(ibl:iel,jbl:jel) )
      allocate(     qsfc(ibl:iel,jbl:jel) )
      allocate(       t2(ibl:iel,jbl:jel) )
      allocate(       q2(ibl:iel,jbl:jel) )
      allocate(      th2(ibl:iel,jbl:jel) )
      allocate(    emiss(ibl:iel,jbl:jel) )
      allocate(      thc(ibl:iel,jbl:jel) )
      allocate(     albd(ibl:iel,jbl:jel) )
      allocate(      f2d(ibl:iel,jbl:jel) )
      allocate(      gsw(ibl:iel,jbl:jel) )
      allocate(      glw(ibl:iel,jbl:jel) )
      allocate(  chklowq(ibl:iel,jbl:jel) )
      allocate(     capg(ibl:iel,jbl:jel) )
      allocate(    snowc(ibl:iel,jbl:jel) )
      allocate(     dsxy(ibl:iel,jbl:jel) )

      ! start with very small, but non-zero, numbers:
      znt = 1.0e-6
      ust = 1.0e-6

      ! start assuming neutral sfclayer:
      mol = 0.0
      zol = 0.0

      num_soil_layers = 5
      allocate(  slab_zs(num_soil_layers) )
      allocate( slab_dzs(num_soil_layers) )
      allocate(  tslb(ibl:iel,jbl:jel,num_soil_layers) )
      allocate(   tmn(ibl:iel,jbl:jel) )

      ! arrays for oml model:
      allocate(   tml(ibl:iel,jbl:jel) )
      allocate(  t0ml(ibl:iel,jbl:jel) )
      allocate(   hml(ibl:iel,jbl:jel) )
      allocate(  h0ml(ibl:iel,jbl:jel) )
      allocate(  huml(ibl:iel,jbl:jel) )
      allocate(  hvml(ibl:iel,jbl:jel) )
      allocate( tmoml(ibl:iel,jbl:jel) )

      allocate(    pta(ibp:iep,jbp:jep,kbp:kep,npt) )
      allocate(   pt3d(ibp:iep,jbp:jep,kbp:kep,npt) )
      allocate(  ptten(ibp:iep,jbp:jep,kbp:kep,npt) )

      allocate(  pdata(npvals,nparcels) )

      allocate(    cfb(ipb:ipe,jpb:jpe,kpb:kpe) )
      allocate(    cfa(kpb:kpe) )
      allocate(    cfc(kpb:kpe) )
      allocate(     d1(kpb:kpe) )
      allocate(     d2(kpb:kpe) )
      allocate(    pdt(ipb:ipe,jpb:jpe,kpb:kpe) )
      allocate(   deft(ipb:ipe,jpb:jpe,kpb:kpe) )
      allocate(    rhs(ipb:ipe,jpb:jpe) )
      allocate(  trans(ipb:ipe,jpb:jpe) )

      call init3d(num_soil_layers,qbudget,asq,bsq,                  &
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
                  pdata,cfb,cfa,cfc,d1,d2,pdt,deft,rhs,trans)

!----------------------------------------------------------------------
!  Now, allocate the MPI arrays (if necessary)

      deallocate( t0 )
      deallocate( rh0 )
!!!      deallocate( qc0 )
      if(ibalance.eq.2 .and.  psolver.ne.4.and.psolver.ne.5 )then
        deallocate( cfb )
        deallocate( cfa )
        deallocate( cfc )
        deallocate( d1 )
        deallocate( d2 )
        deallocate( pdt )
        deallocate( deft )
        deallocate( rhs )
        deallocate( trans )
      endif

#ifdef MPI
      imp = max(1,ni)
      jmp = max(1,nj)
      kmp = max(2,nk)
      kmt = max(2,nk+1)
      rmp = 8
      cmp = 3
#else
      imp = 1
      jmp = 1
      kmp = 2
      rmp = 1
      cmp = 1
#endif
      allocate( reqs_u(rmp) )
      allocate( reqs_v(rmp) )
      allocate( reqs_w(rmp) )
      allocate( reqs_s(rmp) )
      allocate( reqs_p(rmp) )
      allocate( reqs_tk(rmp) )
      allocate( reqs_q(rmp,numq) )
      allocate( reqs_t(rmp,npt) )

      allocate( ww1(jmp,kmp-1) )
      allocate( ww2(jmp,kmp-1) )
      allocate( we1(jmp,kmp-1) )
      allocate( we2(jmp,kmp-1) )
      allocate( ws1(imp,kmp-1) )
      allocate( ws2(imp,kmp-1) )
      allocate( wn1(imp,kmp-1) )
      allocate( wn2(imp,kmp-1) )

      allocate( pw1(jmp,kmp) )
      allocate( pw2(jmp,kmp) )
      allocate( pe1(jmp,kmp) )
      allocate( pe2(jmp,kmp) )
      allocate( ps1(imp,kmp) )
      allocate( ps2(imp,kmp) )
      allocate( pn1(imp,kmp) )
      allocate( pn2(imp,kmp) )

      allocate( vw1(jmp,kmp) )
      allocate( vw2(jmp,kmp) )
      allocate( ve1(jmp,kmp) )
      allocate( ve2(jmp,kmp) )
      allocate( vs1(imp,kmp) )
      allocate( vs2(imp,kmp) )
      allocate( vn1(imp,kmp) )
      allocate( vn2(imp,kmp) )

      allocate( uw31(cmp,jmp,kmp) )
      allocate( uw32(cmp,jmp,kmp) )
      allocate( ue31(cmp,jmp,kmp) )
      allocate( ue32(cmp,jmp,kmp) )
      allocate( us31(imp+1,cmp,kmp) )
      allocate( us32(imp+1,cmp,kmp) )
      allocate( un31(imp+1,cmp,kmp) )
      allocate( un32(imp+1,cmp,kmp) )

      allocate( vw31(cmp,jmp+1,kmp) )
      allocate( vw32(cmp,jmp+1,kmp) )
      allocate( ve31(cmp,jmp+1,kmp) )
      allocate( ve32(cmp,jmp+1,kmp) )
      allocate( vs31(imp,cmp,kmp) )
      allocate( vs32(imp,cmp,kmp) )
      allocate( vn31(imp,cmp,kmp) )
      allocate( vn32(imp,cmp,kmp) )

      allocate( ww31(cmp,jmp,kmp-1) )
      allocate( ww32(cmp,jmp,kmp-1) )
      allocate( we31(cmp,jmp,kmp-1) )
      allocate( we32(cmp,jmp,kmp-1) )
      allocate( ws31(imp,cmp,kmp-1) )
      allocate( ws32(imp,cmp,kmp-1) )
      allocate( wn31(imp,cmp,kmp-1) )
      allocate( wn32(imp,cmp,kmp-1) )

      allocate( sw31(cmp,jmp,kmp) )
      allocate( sw32(cmp,jmp,kmp) )
      allocate( se31(cmp,jmp,kmp) )
      allocate( se32(cmp,jmp,kmp) )
      allocate( ss31(imp,cmp,kmp) )
      allocate( ss32(imp,cmp,kmp) )
      allocate( sn31(imp,cmp,kmp) )
      allocate( sn32(imp,cmp,kmp) )

      allocate( pw31(cmp,jmp,kmp) )
      allocate( pw32(cmp,jmp,kmp) )
      allocate( pe31(cmp,jmp,kmp) )
      allocate( pe32(cmp,jmp,kmp) )
      allocate( ps31(imp,cmp,kmp) )
      allocate( ps32(imp,cmp,kmp) )
      allocate( pn31(imp,cmp,kmp) )
      allocate( pn32(imp,cmp,kmp) )

      allocate( tkw1(cmp,jmp,kmt) )
      allocate( tkw2(cmp,jmp,kmt) )
      allocate( tke1(cmp,jmp,kmt) )
      allocate( tke2(cmp,jmp,kmt) )
      allocate( tks1(imp,cmp,kmt) )
      allocate( tks2(imp,cmp,kmt) )
      allocate( tkn1(imp,cmp,kmt) )
      allocate( tkn2(imp,cmp,kmt) )

      allocate( qw1(cmp,jmp,kmp,numq) )
      allocate( qw2(cmp,jmp,kmp,numq) )
      allocate( qe1(cmp,jmp,kmp,numq) )
      allocate( qe2(cmp,jmp,kmp,numq) )
      allocate( qs1(imp,cmp,kmp,numq) )
      allocate( qs2(imp,cmp,kmp,numq) )
      allocate( qn1(imp,cmp,kmp,numq) )
      allocate( qn2(imp,cmp,kmp,numq) )

      allocate( tw1(cmp,jmp,kmp,npt) )
      allocate( tw2(cmp,jmp,kmp,npt) )
      allocate( te1(cmp,jmp,kmp,npt) )
      allocate( te2(cmp,jmp,kmp,npt) )
      allocate( ts1(imp,cmp,kmp,npt) )
      allocate( ts2(imp,cmp,kmp,npt) )
      allocate( tn1(imp,cmp,kmp,npt) )
      allocate( tn2(imp,cmp,kmp,npt) )

      allocate(          ploc(3,nparcels) )
      allocate( packet(npvals+1,nparcels) )

!----------------------------------------------------------------------

      call setup_output(tdef,qname,budname,xh,xf,yh,yf,xfref,yfref,zh,zf)

      call init_physics(prs0,rf0,cdu,cdv,ce,dum1,dum2,dum3,u0,ua,v0,va,o30,   &
                             lu_index,xland,emiss,thc,albd,znt,mavail,f2d,tsk,u1,v1,w1)

      call init_surface(num_soil_layers,   &
                        dodrag,dosfcflx,xh,ruh,xf,yh,rvh,yf,   &
                        lu_index,xland,tsk,slab_zs,slab_dzs,tslb, &
                        emiss,thc,albd,znt,mavail,dsxy,prs0s,prs0,   &
                        tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml)

      if(irst.eq.1)then
        call read_restart(nstep,nrec,prec,nwrite,nrst,nrad2d,num_soil_layers, &
                              stattim,taptim,rsttim,  &
                              dt,mtime,radtim,qbudget,asq,bsq,              &
                              rain,sws,svs,sps,srs,sgs,sus,shs,tsk,radbcw,radbce,radbcs,radbcn,  &
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
        dtlast = 0.0
      endif

      call getset(dzdx,dzdy,pi0,th0,rho0,prs0,rho,prs,                &
                  ua,u3d,va,v3d,wa,w3d,ppi,pp3d,                      &
                  tha,th3d,qa,q3d,tkea,tke3d,pta,pt3d,                &
                  reqs_u,reqs_v,reqs_w,reqs_s,reqs_tk,                &
                  uw31,uw32,ue31,ue32,us31,us32,un31,un32,            &
                  vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,            &
                  ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,            &
                  sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,            &
                  tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2)

!----------------------------------------------------------------------
!  All done with initialization.  A few more odds and ends ....

      if( adapt_dt.eq.1 )then
        call calccfl(1,rstat,dt,acfl,uf,vf,mf,ua,va,wa,0)
        ndt = 1
        adt = dt
        acfl = cflmax
      endif

      if(irst.ne.1)then
        write(outfile,*)
        write(outfile,*) '  initial conditions:'
        write(outfile,*)
      endif

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
      rtime=sngl(mtime)
      call statpack(nrec,ndt,dt,rtime,adt,acfl,cloudvar,qname,budname,qbudget,asq,bsq, &
                    xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,vf,zh,mh,rmh,mf,     &
                    rstat,pi0,rho0,thv0,th0,qv0,u0,v0,                 &
                    dum1,dum2,dum3,dum4,divx,ppten,prs,                &
                    ua,va,wa,ppi,tha,qa,qten,kmh,kmv,khh,khv,tkea,pta,u10,v10)

    if(irst.ne.1)then
      IF(output_format.eq.1)THEN
        sten = 0.0
        nn = 1
        if(terrain_flag .and. output_interp.eq.1) nn = 2
        DO n=1,nn
          if(n.eq.1)then
            fnum = 51
          else
            fnum = 71
          endif
          call writeout(fnum,1,qname,xh,xf,uf,vf,sigma,zh,zf,mf,pi0,prs0,rho0,th0,qv0,u0,v0,   &
                      zs,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,cdu,cdv,ce,dum1,dum2,dum3,dum4, &
                      rho,prs,sten,ua,uten,va,vten,wa,wten,ppi,tha,      &
                   dissten,thpten,qvpten,qcpten,qipten,upten,vpten,           &
                      lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,     &
                      qa,kmh,kmv,khh,khv,tkea,swten,lwten,radsw,rnflx,radswnet,radlwin,pta,   &
                      num_soil_layers,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br)
        ENDDO
#ifdef NETCDF
      ELSEIF(output_format.eq.2)THEN
        sten = 0.0
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
       sten = 0.0
        call writeout_mult_hdf5(0.0,qname,sigma,sigmaf,xh,xf,uf,yh,yf,vf,mh,zh,mf,zf, &
                      pi0,prs0,rho0,th0,qv0,u0,v0,                     &
                      zs,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,cdu,cdv,ce,dum1,dum2,dum3,dum4,  &
                      rho,prs,sten,ua,uten,va,vten,wa,wten,ppi,tha,    &
                      qa,kmh,kmv,khh,khv,tkea,pta,num_soil_layers,   &
                      lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hnvml,hfx,qfx,gsw,glw,tslb,   &
                      radsw,rnflx,radswnet,radlwin,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br,   &
                      dissten,thpten,qvpten,qcpten,qipten,upten,vpten,swten,lwten)
#endif
      ENDIF
    endif

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

      write(outfile,*)
      write(outfile,*) '-------------Done with Preprocessors-----------'
      write(outfile,*)

      if(iconly.eq.1)then
        write(outfile,*)
        write(outfile,*) '  User has requested initial conditions only'
        write(outfile,*) '     (iconly = 1)'
        write(outfile,*) '  ... stopping ... '
        write(outfile,*)
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        stop 55555
      endif

!----------------------------------------------------------------------

      time_sound=0.
      time_poiss=0.
      time_advs=0.
      time_advu=0.
      time_advv=0.
      time_advw=0.
      time_buoyan=0.
      time_turb=0.
      time_diffu=0.
      time_microphy=0.
      time_stat=0.
      time_bc=0.
      time_misc=0.
      time_integ=0.
      time_rdamp=0.
      time_divx=0.
      time_write=0.
      time_tmix=0.
      time_cor=0.
      time_fall=0.
      time_satadj=0.
      time_sfcphys=0.
      time_parcels=0.0
      time_rad=0.
      time_pbl=0.
      time_swath=0.
      time_pdef=0.
#ifdef MPI
      time_mpu1=0.
      time_mpv1=0.
      time_mpw1=0.
      time_mpp1=0.
      time_mpu2=0.
      time_mpv2=0.
      time_mpw2=0.
      time_mpp2=0.
      time_mps1=0.
      time_mpq1=0.
      time_mptk1=0.
      time_mptk2=0.
      time_mps2=0.
      time_mpq2=0.
      time_mpb=0.

      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

      if(myid.eq.0)then
        tstart=mpi_wtime()
      endif
#endif

      ! This initializes timer
      if(timestats.ge.1)then
        call system_clock(count,rate,maxr)
        clock_rate=1.0/rate
        xtime=mytime()
      endif

!----------------------------------------------------------------------
!  Time loop

      if(timestats.ge.1)then
        steptime1 = 0.0
        steptime2 = 0.0
      endif

      do while( mtime.lt.timax )
        nstep = nstep + 1
        call solve(nstep,nrec,prec,nwrite,nrst,rbufsz,num_soil_layers,ndt,     &
                   dt,dtlast,mtime,stattim,taptim,rsttim,radtim,adt,acfl,  &
                   dodrag,dosfcflx,cloudvar,rhovar,qname,budname,bud,bud2,qbudget,asq,bsq, &
                   xh,rxh,uh,ruh,xf,rxf,uf,ruf,yh,vh,rvh,yf,vf,rvf,   &
                   sigma,sigmaf,tauh,taus,zh,mh,rmh,tauf,zf,mf,rmf,   &
                   rstat,rho0s,pi0s,prs0s,rth0s,pi0,rho0,prs0,thv0,th0,qv0,qc0,  &
                   ql0,rr0,rf0,rrf0,rru0,rrv0,                        &
                   zs,gz,dzdx,dzdy,rain,sws,svs,sps,srs,sgs,sus,shs,  &
                   doimpl,tsk,thflux,qvflux,cdu,cdv,ce,u1,v1,w1,      &
                   radbcw,radbce,radbcs,radbcn,                       &
                   dum1,dum2,dum3,dum4,divx,rho,prs,                  &
                   t11,t12,t13,t22,t23,t33,                           &
                   gx,u0,rru,ua,u3d,uten,uten1,                       &
                   gy,v0,rrv,va,v3d,vten,vten1,                       &
                   rrw,wa,w3d,wten,wten1,ppi,pp3d,ppten,sten,         &
                   tha,th3d,thten,thten1,thterm,tk,qa,q3d,qten,zvdarray, &
                   kmh,kmv,khh,khv,tkea,tke3d,tketen,                 &
                   dissten,thpten,qvpten,qcpten,qipten,upten,vpten,   &
                   swten,lwten,o30,radsw,rnflx,radswnet,radlwin,rad2d, &
                   x,y,z,za,zp,                                       &
                   lu_index,kpbl2d,psfc,u10,v10,hfx,qfx,xland,znt,ust,   &
                   hpbl,wspd,psim,psih,gz1oz0,br,                     &
                   CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                     &
                   MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                  &
                   CK,CKA,CD,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,  &
                   f2d,gsw,glw,chklowq,capg,snowc,dsxy,               &
                   slab_zs,slab_dzs,tslb,tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,        &
                   pta,pt3d,ptten,                                    &
                   pdata,cfb,cfa,cfc,d1,d2,pdt,deft,rhs,trans,        &
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
        if(timestats.eq.2)then
          steptime2=time_sound+time_poiss+time_buoyan+time_turb+            &
                    time_diffu+time_microphy+time_stat+                     &
                    time_bc+time_misc+time_integ+time_rdamp+time_divx+      &
                    time_write+time_tmix+time_cor+time_fall+                &
                    time_satadj+time_sfcphys+time_parcels+                  &
                    time_rad+time_pbl+time_swath+time_pdef+                 &
#ifdef MPI
                    time_mpu1+time_mpv1+time_mpw1+time_mpp1+                &
                    time_mpu2+time_mpv2+time_mpw2+time_mpp2+                &
                    time_mps1+time_mpq1+time_mptk1+                         &
                    time_mps2+time_mpq2+time_mptk2+time_mpb+                &
#endif
                    time_advs+time_advu+time_advv+time_advw
          write(6,157) nstep,steptime2-steptime1
157       format('    timing for time step ',i12,':',f12.4,' s')
          steptime1 = steptime2
        endif
      enddo

!----------------------------------------------------------------------
!  write new stats descriptor file, if necessary:

      IF( output_format.eq.1 .and. myid.eq.0 )THEN
        IF( adapt_dt.eq.1 .and. statfrq.lt.0.0 )THEN
          print *,'  re-writing GrADS stats descriptor file .... '
          call write_statsctl(tdef,qname,budname,nstep+1)
        ENDIF
      ENDIF

!----------------------------------------------------------------------

#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      if(timestats.ge.1) time_mpb=time_mpb+mytime()

      if(myid.eq.0)then
        tend=mpi_wtime()
        print *
        print *,'Total time (s): ',tend-tstart
        print *
      endif

#endif
!----------------------------------------------------------------------

    IF(timestats.ge.1)THEN

      time_solve=time_sound+time_poiss+time_buoyan+time_turb+             &
                  time_diffu+time_microphy+time_stat+                     &
                  time_bc+time_misc+time_integ+time_rdamp+time_divx+      &
                  time_write+time_tmix+time_cor+time_fall+                &
                  time_satadj+time_sfcphys+time_parcels+                  &
                  time_rad+time_pbl+time_swath+time_pdef+                 &
#ifdef MPI
                  time_mpu1+time_mpv1+time_mpw1+time_mpp1+                &
                  time_mpu2+time_mpv2+time_mpw2+time_mpp2+                &
                  time_mps1+time_mpq1+time_mptk1+                         &
                  time_mps2+time_mpq2+time_mptk2+time_mpb+                &
#endif
                  time_advs+time_advu+time_advv+time_advw

#ifdef MPI
      mp_total=time_mpu1+time_mpv1+time_mpw1+time_mpp1+                   &
               time_mpu2+time_mpv2+time_mpw2+time_mpp2+                   &
               time_mps1+time_mpq1+time_mptk1+                            &
               time_mps2+time_mpq2+time_mptk2+time_mpb
#endif

      write(outfile,*)
      write(outfile,*) 'Total time: ',time_solve
      write(outfile,*)
      time_solve=0.01*time_solve
      if(time_solve.lt.0.0001) time_solve=1.

      write(outfile,100) 'sound   ',time_sound,time_sound/time_solve
      write(outfile,100) 'poiss   ',time_poiss,time_poiss/time_solve
      write(outfile,100) 'advs    ',time_advs,time_advs/time_solve
      write(outfile,100) 'advu    ',time_advu,time_advu/time_solve
      write(outfile,100) 'advv    ',time_advv,time_advv/time_solve
      write(outfile,100) 'advw    ',time_advw,time_advw/time_solve
      write(outfile,100) 'divx    ',time_divx,time_divx/time_solve
      write(outfile,100) 'buoyan  ',time_buoyan,time_buoyan/time_solve
      write(outfile,100) 'turb    ',time_turb,time_turb/time_solve
      write(outfile,100) 'sfcphys ',time_sfcphys,time_sfcphys/time_solve
      write(outfile,100) 'tmix    ',time_tmix,time_tmix/time_solve
      write(outfile,100) 'cor     ',time_cor,time_cor/time_solve
      write(outfile,100) 'diffu   ',time_diffu,time_diffu/time_solve
      write(outfile,100) 'rdamp   ',time_rdamp,time_rdamp/time_solve
      write(outfile,100) 'microphy',time_microphy,time_microphy/time_solve
      write(outfile,100) 'satadj  ',time_satadj,time_satadj/time_solve
      write(outfile,100) 'fallout ',time_fall,time_fall/time_solve
      write(outfile,100) 'radiatio',time_rad,time_rad/time_solve
      write(outfile,100) 'pbl     ',time_pbl,time_pbl/time_solve
      write(outfile,100) 'stat    ',time_stat,time_stat/time_solve
      write(outfile,100) 'bc      ',time_bc,time_bc/time_solve
      write(outfile,100) 'integ   ',time_integ,time_integ/time_solve
      write(outfile,100) 'write   ',time_write,time_write/time_solve
      write(outfile,100) 'misc    ',time_misc,time_misc/time_solve
      write(outfile,100) 'swaths  ',time_swath,time_swath/time_solve
      write(outfile,100) 'pdef    ',time_pdef,time_pdef/time_solve
      write(outfile,100) 'parcels ',time_parcels,time_parcels/time_solve
#ifdef MPI
      write(outfile,100) 'mp_total',mp_total,mp_total/time_solve
      write(outfile,*)
      write(outfile,100) 'mpu1    ',time_mpu1,time_mpu1/time_solve
      write(outfile,100) 'mpv1    ',time_mpv1,time_mpv1/time_solve
      write(outfile,100) 'mpw1    ',time_mpw1,time_mpw1/time_solve
      write(outfile,100) 'mpp1    ',time_mpp1,time_mpp1/time_solve
      write(outfile,100) 'mpu2    ',time_mpu2,time_mpu2/time_solve
      write(outfile,100) 'mpv2    ',time_mpv2,time_mpv2/time_solve
      write(outfile,100) 'mpw2    ',time_mpw2,time_mpw2/time_solve
      write(outfile,100) 'mpp2    ',time_mpp2,time_mpp2/time_solve
      write(outfile,100) 'mps1    ',time_mps1,time_mps1/time_solve
      write(outfile,100) 'mpq1    ',time_mpq1,time_mpq1/time_solve
      write(outfile,100) 'mps2    ',time_mps2,time_mps2/time_solve
      write(outfile,100) 'mpq2    ',time_mpq2,time_mpq2/time_solve
      write(outfile,100) 'mptk1   ',time_mptk1,time_mptk1/time_solve
      write(outfile,100) 'mptk2   ',time_mptk2,time_mptk2/time_solve
      write(outfile,100) 'mpb     ',time_mpb,time_mpb/time_solve
#endif
      write(outfile,*)

100   format(3x,a8,' :  ',f10.2,2x,f6.2,'%')

    ENDIF

!  End time loop
!----------------------------------------------------------------------

      close(unit=51)
      close(unit=52)
      close(unit=53)
      close(unit=54)
      close(unit=60)

!----------------------------------------------------------------------

#ifdef MPI
      call MPI_FINALIZE(rc)
      print *,'Program terminated normally:  myid=',myid
#else
      print *,'Program terminated normally'
#endif

      stop
      end


