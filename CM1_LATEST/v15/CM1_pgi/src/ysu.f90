

!WRf:model_layer:physics
!
!
!
!
!
!
!
module module_bl_ysu
contains
!
!
!-------------------------------------------------------------------
!
   subroutine ysu(u3d,v3d,th3d,t3d,qv3d,qc3d,qi3d,p3d,p3di,pi3d,               &
                  rublten,rvblten,rthblten,                                    &
                  rqvblten,rqcblten,rqiblten,flag_qi,                          &
                  cp,g,rovcp,rd,rovg,ep1,ep2,karman,xlv,rv,                    &
                  dz8w,psfc,                                                   &
!!!                  znu,znw,mut,p_top,                                        &
                  znt,ust,hpbl,psim,psih,                                      &
                  xland,hfx,qfx,gz1oz0,wspd,br,                                &
                  dt,kpbl2d,                                                   &
                  exch_h,exch_m,                                               &
                  u10,v10,                                                     &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte,                                   &
                !optional
                  regime                                           )
!-------------------------------------------------------------------
      implicit none
!-------------------------------------------------------------------
!-- u3d         3d u-velocity interpolated to theta points (m/s)
!-- v3d         3d v-velocity interpolated to theta points (m/s)
!-- th3d	3d potential temperature (k)
!-- t3d         temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- qc3d        3d cloud mixing ratio (kg/kg)
!-- qi3d        3d ice mixing ratio (kg/kg)
!               (note: if P_QI<PARAM_FIRST_SCALAR this should be zero filled)
!-- p3d         3d pressure (pa)
!-- p3di        3d pressure (pa) at interface level
!-- pi3d	3d exner function (dimensionless)
!-- rr3d	3d dry air density (kg/m^3)
!-- rublten     u tendency due to
!               pbl parameterization (m/s/s)
!-- rvblten     v tendency due to
!               pbl parameterization (m/s/s)
!-- rthblten    theta tendency due to
!               pbl parameterization (K/s)
!-- rqvblten    qv tendency due to
!               pbl parameterization (kg/kg/s)
!-- rqcblten    qc tendency due to
!               pbl parameterization (kg/kg/s)
!-- rqiblten    qi tendency due to
!               pbl parameterization (kg/kg/s)
!-- cp          heat capacity at constant pressure for dry air (j/kg/k)
!-- g           acceleration due to gravity (m/s^2)
!-- rovcp       r/cp
!-- rd          gas constant for dry air (j/kg/k)
!-- rovg 	r/g
!-- dz8w	dz between full levels (m)
!-- xlv         latent heat of vaporization (j/kg)
!-- rv 		gas constant for water vapor (j/kg/k)
!-- psfc        pressure at the surface (pa)
!-- znt		roughness length (m)
!-- ust		u* in similarity theory (m/s)
!-- hpbl	pbl height (m)
!-- psim        similarity stability function for momentum
!-- psih        similarity stability function for heat
!-- xland	land mask (1 for land, 2 for water)
!-- hfx		upward heat flux at the surface (w/m^2)
!-- qfx		upward moisture flux at the surface (kg/m^2/s)
!-- gz1oz0      log(z/z0) where z0 is roughness length
!-- wspd        wind speed at lowest model level (m/s)
!-- u10         u-wind speed at 10 m (m/s)
!-- v10         v-wind speed at 10 m (m/s)
!-- br          bulk richardson number in surface layer
!-- dt		time step (s)
!-- rvovrd      r_v divided by r_d (dimensionless)
!-- ep1         constant for virtual temperature (r_v/r_d - 1) (dimensionless)
!-- ep2         constant for specific humidity calculation
!-- karman      von karman constant
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-------------------------------------------------------------------
!
   integer,parameter ::  ndiff = 3
   real,parameter    ::  rcl = 1.0
!
   integer,  intent(in   )   ::      ids,ide, jds,jde, kds,kde,                &
                                     ims,ime, jms,jme, kms,kme,                &
                                     its,ite, jts,jte, kts,kte
!
   real,     intent(in   )   ::      dt,cp,g,rovcp,rovg,rd,xlv,rv
!
   real,     intent(in )     ::      ep1,ep2,karman
!
   real,     dimension( ims:ime, jms:jme , kms:kme )                          , &
             intent(in   )   ::                                          qv3d, &
                                         			         qc3d, &
 				                                         qi3d, &
                                   		                          p3d, &
            		                                                 pi3d, &
		                      	                                 th3d, &
		                                		          t3d, &
				                                         dz8w
   real,     dimension( ims:ime, jms:jme , kms:kme )                          , &
             intent(in   )   ::                                          p3di
!
   real,     dimension( ims:ime, jms:jme , kms:kme )                          , &
             intent(inout)   ::                                       rublten, &
                                              			      rvblten, &
                                               			     rthblten, &
                                                 	             rqvblten, &
                                                                     rqcblten
!
   real,     dimension( ims:ime, jms:jme , kms:kme )                         , &
             intent(inout)   ::                                        exch_h, &
                                                                       exch_m
   real,     dimension( ims:ime, jms:jme )                                   , &
             intent(in   )   ::                                           u10, &
                                                                          v10
!
   real,     dimension( ims:ime, jms:jme )                                   , &
             intent(in   )   ::                                         xland, &
     			                                                  hfx, &
                                                                          qfx, &
                                                                           br, &
                                                                         psfc
   real,     dimension( ims:ime, jms:jme )                                   , &
             intent(in   )   ::                                                &
                                                                         psim, &
                                                                         psih, &
                                                                       gz1oz0
   real,     dimension( ims:ime, jms:jme )                                   , &
             intent(inout)   ::                                           znt, &
                                                                          ust, &
                                                                         hpbl, &
                                                                          wspd
!
  real,      dimension( ims:ime, jms:jme , kms:kme )                          , &
             intent(in   )   ::                                           u3d, &
                                                                          v3d
!
  integer,   dimension( ims:ime, jms:jme )                                   , &
             intent(out  )   ::                                        kpbl2d
  logical, intent(in)        ::                                       flag_qi
!optional
   real,     dimension( ims:ime, jms:jme )                                   , &
             optional                                                        , &
             intent(inout)   ::                                        regime
!
   real,     dimension( ims:ime, jms:jme , kms:kme )                          , &
             optional                                                        , &
             intent(inout)   ::                                       rqiblten
!
!!!   real,     dimension( kms:kme )                                            , &
!!!             optional                                                        , &
!!!             intent(in   )   ::                                           znu, &
!!!                                                                          znw
!
!!!   real,     dimension( ims:ime, jms:jme )                                   , &
!!!             optional                                                        , &
!!!             intent(in   )   ::                                           mut
!
!!!   real,     optional, intent(in   )   ::                               p_top
!
!local
   integer ::  i,j,k
   real,     dimension( its:ite, kts:kte*ndiff )  ::                 rqvbl2dt, &
                                                                         qv2d
   real,     dimension( its:ite, kts:kte )  ::                            pdh
   real,     dimension( ims:ime, kms:kme )  ::                                 &
                         u2d,v2d,t2d,pi2d,dz2d,uten,vten,tten,km2d,kh2d
   real,     dimension( its:ite, kts:kte+1 )  ::                         pdhi
   real,     dimension( its:ite )  ::                                          &
                                                                        dusfc, &
                                                                        dvsfc, &
                                                                        dtsfc, &
                                                                        dqsfc
!
   qv2d(:,:) = 0.0

!$omp parallel do default(shared)  &
!$omp private(i,j,k,pdh,pdhi,u2d,v2d,t2d,pi2d,dz2d,qv2d,tten,uten,vten,   &
!$omp rqvbl2dt,dusfc,dvsfc,dtsfc,dqsfc,km2d,kh2d)
   do j = jts,jte
!      if(present(mut))then
!! For ARW we will replace p and p8w with dry hydrostatic pressure
!        do k = kts,kte+1
!          do i = its,ite
!             if(k.le.kte)pdh(i,k) = mut(i,j)*znu(k) + p_top
!             pdhi(i,k) = mut(i,j)*znw(k) + p_top
!          enddo
!        enddo
!      else
        do k = kts,kte+1
          do i = its,ite
             if(k.le.kte)pdh(i,k) = p3d(i,j,k)
             pdhi(i,k) = p3di(i,j,k)
          enddo
        enddo
!      endif
      do k = kts,kte
        do i = its,ite
          u2d(i,k) = u3d(i,j,k)
          v2d(i,k) = v3d(i,j,k)
          t2d(i,k) = t3d(i,j,k)
          pi2d(i,k) = pi3d(i,j,k)
          dz2d(i,k) = dz8w(i,j,k)
          qv2d(i,k) = qv3d(i,j,k)
          qv2d(i,k+kte) = qc3d(i,j,k)
          if(present(rqiblten)) qv2d(i,k+kte+kte) = qi3d(i,j,k)
        enddo
      enddo
      km2d = 0.0
      kh2d = 0.0
!
      call ysu2d(J=j,ux=u2d(ims,kms),vx=v2d(ims,kms)                       &
              ,tx=t2d(ims,kms)                                               &
              ,qx=qv2d(its,kts)                                                &
              ,p2d=pdh(its,kts),p2di=pdhi(its,kts)                             &
              ,pi2d=pi2d(ims,kms)                                            &
              ,utnp=uten(ims,kms),vtnp=vten(ims,kms)                 &
              ,ttnp=tten(ims,kms),qtnp=rqvbl2dt(its,kts),ndiff=ndiff     &
              ,cp=cp,g=g,rovcp=rovcp,rd=rd,rovg=rovg                           &
              ,xlv=xlv,rv=rv                                                   &
              ,ep1=ep1,ep2=ep2,karman=karman                                   &
              ,dz8w2d=dz2d(ims,kms)                                          &
              ,psfcpa=psfc(ims,j),znt=znt(ims,j),ust=ust(ims,j)                &
              ,hpbl=hpbl(ims,j)                                                &
              ,regime=regime(ims,j),psim=psim(ims,j)                           &
              ,psih=psih(ims,j),xland=xland(ims,j)                             &
              ,hfx=hfx(ims,j),qfx=qfx(ims,j)                                   &
              ,wspd=wspd(ims,j),br=br(ims,j)                                   &
              ,dusfc=dusfc,dvsfc=dvsfc,dtsfc=dtsfc,dqsfc=dqsfc                 &
              ,dt=dt,rcl=1.0,kpbl1d=kpbl2d(ims,j)                              &
              ,exch_hx=kh2d(ims,kms),exch_mx=km2d(ims,kms)                     &
              ,u10=u10(ims,j),v10=v10(ims,j)                                   &
              ,gz1oz0=gz1oz0(ims,j)                                            &
              ,ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde               &
              ,ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme               &
              ,its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte   )
!
      do k = kts,kte
        do i = its,ite
          rthblten(i,j,k) = tten(i,k)/pi3d(i,j,k)
          rublten(i,j,k) = uten(i,k)
          rvblten(i,j,k) = vten(i,k)
          rqvblten(i,j,k) = rqvbl2dt(i,k)
          rqcblten(i,j,k) = rqvbl2dt(i,k+kte)
          if(present(rqiblten)) rqiblten(i,j,k) = rqvbl2dt(i,k+kte+kte)
          exch_h(i,j,k) = kh2d(i,k)
          exch_m(i,j,k) = km2d(i,k)
        enddo
      enddo
   enddo
!
   end subroutine ysu
!
!-------------------------------------------------------------------
!
   subroutine ysu2d(j,ux,vx,tx,qx,p2d,p2di,pi2d,                               &
                  utnp,vtnp,ttnp,qtnp,ndiff,                                   &
                  cp,g,rovcp,rd,rovg,ep1,ep2,karman,xlv,rv,                    &
                  dz8w2d,psfcpa,                                               &
                  znt,ust,hpbl,psim,psih,                                      &
                  xland,hfx,qfx,wspd,br,                                       &
                  dusfc,dvsfc,dtsfc,dqsfc,                                     &
                  dt,rcl,kpbl1d,                                               &
                  exch_hx,exch_mx,                                             &
                  u10,v10,                                                     &
                  gz1oz0,                                                      &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte,                                   &
                !optional
                  regime                                           )
!-------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------
!
!     this code is a revised vertical diffusion package ("ysupbl")
!     with a nonlocal turbulent mixing in the pbl after "mrfpbl".
!     the ysupbl (hong et al. 2006) is based on the study of noh 
!     et al.(2003) and accumulated realism of the behavior of the 
!     troen and mahrt (1986) concept implemented by hong and pan(1996). 
!     the major ingredient of the ysupbl is the inclusion of an explicit
!     treatment of the entrainment processes at the entrainment layer.
!     this routine uses an implicit approach for vertical flux
!     divergence and does not require "miter" timesteps.
!     it includes vertical diffusion in the stable atmosphere
!     and moist vertical diffusion in clouds.
!
!     mrfpbl:
!     coded by song-you hong (ncep), implemented by jimy dudhia (ncar)
!              fall 1996
!
!     ysupbl:
!     coded by song-you hong (yonsei university) and implemented by
!              song-you hong (yonsei university) and jimy dudhia (ncar)
!              summer 2002
!
!     further modifications :
!              an enhanced stable layer mixing, april 2008
!               ==> increase pbl height when sfc is stable (hong 2010)
!              pressure-level diffusion, april 2009
!               ==> negligible differences
!              implicit forcing for momentum with clean up, july 2009
!               ==> prevents model blowup when sfc layer is too low
!              incresea of lamda, 30 < 0.1 x del z < 300, feb 2010
!               ==> prevents model blowup when delz is extremely large
!              revised prandtl number at surface, peggy lemone, feb 2010
!               ==> increase kh, decrease mixing due to counter-gradient term
!
!     references:
!
!        hong (2010) quart. j. roy. met. soc
!        hong, noh, and dudhia (2006), mon. wea. rev. 
!        hong and pan (1996), mon. wea. rev.
!        noh, chun, hong, and raasch (2003), boundary layer met.
!        troen and mahrt (1986), boundary layer met.
!
!-------------------------------------------------------------------
!
   real,parameter    ::  xkzmin = 0.01,xkzmax = 1000.,rimin = -100.
   real,parameter    ::  rlam = 30.,prmin = 0.25,prmax = 4.
   real,parameter    ::  brcr_ub = 0.0,brcr_sb = 0.25,cori = 1.e-4
   real,parameter    ::  afac = 6.8,bfac = 6.8,pfac = 2.0,pfac_q = 2.0
   real,parameter    ::  phifac = 8.,sfcfrac = 0.1
   real,parameter    ::  d1 = 0.02, d2 = 0.05, d3 = 0.001
   real,parameter    ::  h1 = 0.33333335, h2 = 0.6666667
   real,parameter    ::  ckz = 0.001,zfmin = 1.e-8,aphi5 = 5.,aphi16 = 16.
   real,parameter    ::  tmin=1.e-2
   real,parameter    ::  gamcrt = 3.,gamcrq = 2.e-3
   real,parameter    ::  xka = 2.4e-5
   integer,parameter ::  imvdif = 1
!
   integer,  intent(in   )   ::     ids,ide, jds,jde, kds,kde,                 &
                                    ims,ime, jms,jme, kms,kme,                 &
                                    its,ite, jts,jte, kts,kte,                 &
                                    j,ndiff
!
   real,     intent(in   )   ::     dt,rcl,cp,g,rovcp,rovg,rd,xlv,rv
!
   real,     intent(in )     ::     ep1,ep2,karman
!
   real,     dimension( ims:ime, kms:kme ),                                    &
             intent(in)      ::                                        dz8w2d, &
                                                                         pi2d
!
   real,     dimension( ims:ime, kms:kme )                                   , &
             intent(in   )   ::                                            tx
   real,     dimension( its:ite, kts:kte*ndiff )                             , &
             intent(in   )   ::                                            qx
!
   real,     dimension( ims:ime, kms:kme )                                   , &
             intent(inout)   ::                                          utnp, &
                                                                         vtnp, &
                                                                         ttnp
   real,     dimension( its:ite, kts:kte*ndiff )                             , &
             intent(inout)   ::                                          qtnp
!
   real,     dimension( its:ite, kts:kte+1 )                                 , &
             intent(in   )   ::                                          p2di
!
   real,     dimension( its:ite, kts:kte )                                   , &
             intent(in   )   ::                                           p2d
!
!
   real,     dimension( ims:ime )                                            , &
             intent(inout)   ::                                           ust, &
                                                                         hpbl, &
                                                                          znt
   real,     dimension( ims:ime )                                            , &
             intent(in   )   ::                                         xland, &
                                                                          hfx, &
                                                                          qfx
!
   real,     dimension( ims:ime ), intent(inout)   ::                    wspd
   real,     dimension( ims:ime ), intent(in  )    ::                      br
!
   real,     dimension( ims:ime ), intent(in   )   ::                    psim, &
                                                                         psih
   real,     dimension( ims:ime ), intent(in   )   ::                  gz1oz0
!
   real,     dimension( ims:ime ), intent(in   )   ::                  psfcpa
   integer,  dimension( ims:ime ), intent(out  )   ::                  kpbl1d
!
   real,     dimension( ims:ime, kms:kme )                                   , &
             intent(in   )   ::                                            ux, &
                                                                           vx
!optional
   real,     dimension( ims:ime )                                            , &
             optional                                                        , &
             intent(inout)   ::                                        regime
!
! local vars
!
   real,     dimension( its:ite )            ::                           hol
   real,     dimension( its:ite, kts:kte+1 ) ::                            zq
!
   real,     dimension( its:ite, kts:kte )   ::                        	       &
                                                                     thx,thvx, &
                                                                          del, &
                                                                          dza, &
                                                                          dzq, &
                                                                           za
!
   real,    dimension( its:ite )             ::                                &
                                                                         rhox, &
                                                                       govrth, &
                                                                  zl1,thermal, &
                                                                 wscale,hgamt, &
                                                                   hgamq,brdn, &
                                                                    brup,phim, &
                                                                         phih, &
                                                                  dusfc,dvsfc, &
                                                                  dtsfc,dqsfc, &
                                                                        prpbl, &
                                                                        wspd1
!
   real,    dimension( its:ite, kts:kte )    ::                     xkzm,xkzh, &
                                                     		        f1,f2, &
                                                     		        r1,r2, &
                                                                        ad,au, &
                                                                           cu, &
                                                                           al, &
                                                                         xkzq, &
                                                                         zfac
!
!jdf added exch_hx
   real,    dimension( ims:ime, kms:kme )                                    , &
            intent(inout)   ::                                        exch_hx, &
                                                                      exch_mx
!
   real,    dimension( ims:ime )                                             , & 
            intent(in  )    ::                                            u10, &
                                                                          v10
   real,    dimension( its:ite )    ::                                         &
                                                                         brcr, &
                                                                        sflux, &
                                                                    brcr_sbro
!
   real,    dimension( its:ite, kts:kte, ndiff)  ::                     r3,f3
   integer, dimension( its:ite )             ::                          kpbl
!
   logical, dimension( its:ite )             ::                        pblflg, &
                                                                       sfcflg, &
                                                                       stable
!
   integer ::  n,i,k,l,ic,is
   integer ::  klpbl, ktrace1, ktrace2, ktrace3
!
!
   real    ::  dt2,rdt,spdk2,fm,fh,hol1,gamfac,vpert,prnum,prnum0
   real    ::  xkzo,ss,ri,qmean,tmean,alph,chi,zk,rl2,dk,sri
   real    ::  brint,dtodsd,dtodsu,rdz,dsdzt,dsdzq,dsdz2,rlamdz
   real    ::  utend,vtend,ttend,qtend
   real    ::  dtstep,govrthv
   real    ::  cont, conq, conw, conwrc 
!
   real, dimension( its:ite, kts:kte )     ::                         wscalek, &
                                                                  xkzml,xkzhl, &
                                                               zfacent,entfac
   real, dimension( its:ite )              ::                            ust3, &
                                                                 wstar3,wstar, &
                                                                  hgamu,hgamv, &
                                                                      wm2, we, &
                                                                       bfxpbl, &
                                                                hfxpbl,qfxpbl, &
                                                                ufxpbl,vfxpbl, &
                                                                  delta,dthvx
   real    ::  prnumfac,bfx0,hfx0,qfx0,delb,dux,dvx,                           &
               dsdzu,dsdzv,wm3,dthx,dqx,wspd10,ross,tem1,dsig,tvcon,conpr,     &
               prfac,phim8z
!
!----------------------------------------------------------------------
!
   klpbl = kte
!
   cont=cp/g
   conq=xlv/g
   conw=1./g
   conwrc = conw*sqrt(rcl)
   conpr = bfac*karman*sfcfrac
!
!  k-start index for tracer diffusion
!
   ktrace1 = 0
   ktrace2 = 0 + kte
   ktrace3 = 0 + kte*2
!
   do k = kts,kte
     do i = its,ite
       thx(i,k) = tx(i,k)/pi2d(i,k)
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       tvcon = (1.+ep1*qx(i,k))
       thvx(i,k) = thx(i,k)*tvcon
     enddo
   enddo
!
   do i = its,ite
     tvcon = (1.+ep1*qx(i,1))
     rhox(i) = psfcpa(i)/(rd*tx(i,1)*tvcon)
     govrth(i) = g/thx(i,1)
   enddo
!
!-----compute the height of full- and half-sigma levels above ground
!     level, and the layer thicknesses.
!
   do i = its,ite
     zq(i,1) = 0.
   enddo
!
   do k = kts,kte
     do i = its,ite
       zq(i,k+1) = dz8w2d(i,k)+zq(i,k)
     enddo
   enddo
! 
   do k = kts,kte
     do i = its,ite
       za(i,k) = 0.5*(zq(i,k)+zq(i,k+1))
       dzq(i,k) = zq(i,k+1)-zq(i,k)
       del(i,k) = p2di(i,k)-p2di(i,k+1)
     enddo
   enddo
!
   do i = its,ite
     dza(i,1) = za(i,1)
   enddo
!
   do k = kts+1,kte
     do i = its,ite
       dza(i,k) = za(i,k)-za(i,k-1)
     enddo
   enddo
!
!
!-----initialize vertical tendencies and
!
   utnp(:,:) = 0.
   vtnp(:,:) = 0.
   ttnp(:,:) = 0.
   qtnp(:,:) = 0.
!
   do i = its,ite
     wspd1(i) = sqrt(ux(i,1)*ux(i,1)+vx(i,1)*vx(i,1))+1.e-9
   enddo
!
!---- compute vertical diffusion
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     compute preliminary variables
!
   dtstep = dt
   dt2 = 2.*dtstep
   rdt = 1./dt2
!
   do i = its,ite
     bfxpbl(i) = 0.0
     hfxpbl(i) = 0.0
     qfxpbl(i) = 0.0
     ufxpbl(i) = 0.0
     vfxpbl(i) = 0.0
     hgamu(i)  = 0.0
     hgamv(i)  = 0.0
     delta(i)  = 0.0
   enddo
!
   do k = kts,klpbl
     do i = its,ite
       wscalek(i,k) = 0.0
     enddo
   enddo
!
   do k = kts,klpbl
     do i = its,ite
       zfac(i,k) = 0.0
     enddo
   enddo
!
   do i = its,ite
     dusfc(i) = 0.
     dvsfc(i) = 0.
     dtsfc(i) = 0.
     dqsfc(i) = 0.
   enddo
!
   do i = its,ite
     hgamt(i)  = 0.
     hgamq(i)  = 0.
     wscale(i) = 0.
     kpbl(i)   = 1
     hpbl(i)   = zq(i,1)
     zl1(i)    = za(i,1)
     thermal(i)= thvx(i,1)
     pblflg(i) = .true.
     sfcflg(i) = .true.
     sflux(i) = hfx(i)/rhox(i)/cp + qfx(i)/rhox(i)*ep1*thx(i,1)
     if(br(i).gt.0.0) sfcflg(i) = .false.
   enddo
!
!     compute the first guess of pbl height
!
   do i = its,ite
     stable(i) = .false.
     brup(i) = br(i)
     brcr(i) = brcr_ub
   enddo
!
   do k = 2,klpbl
     do i = its,ite
       if(.not.stable(i))then
         brdn(i) = brup(i)
         spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
         brup(i) = (thvx(i,k)-thermal(i))*(g*za(i,k)/thvx(i,1))/spdk2
         kpbl(i) = k
         stable(i) = brup(i).gt.brcr(i)
       endif
     enddo
   enddo
!
   do i = its,ite
     k = kpbl(i)
     if(brdn(i).ge.brcr(i))then
       brint = 0.
     elseif(brup(i).le.brcr(i))then
       brint = 1.
     else
       brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
     endif
     hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
     if(hpbl(i).lt.zq(i,2)) kpbl(i) = 1
     if(kpbl(i).le.1) pblflg(i) = .false.
   enddo
!
   do i = its,ite
     fm = gz1oz0(i)-psim(i)
     fh = gz1oz0(i)-psih(i)
     hol(i) = max(br(i)*fm*fm/fh,rimin)
     if(sfcflg(i))then
       hol(i) = min(hol(i),-zfmin)
     else
       hol(i) = max(hol(i),zfmin)
     endif
     hol1 = hol(i)*hpbl(i)/zl1(i)*sfcfrac
     hol(i) = -hol(i)*hpbl(i)/zl1(i)
     if(sfcflg(i))then
       phim(i) = (1.-aphi16*hol1)**(-1./4.)
       phih(i) = (1.-aphi16*hol1)**(-1./2.)
       bfx0 = max(sflux(i),0.)
       hfx0 = max(hfx(i)/rhox(i)/cp,0.)
       qfx0 = max(ep1*thx(i,1)*qfx(i)/rhox(i),0.)
       wstar3(i) = (govrth(i)*bfx0*hpbl(i))
       wstar(i) = (wstar3(i))**h1
     else
       phim(i) = (1.+aphi5*hol1)
       phih(i) = phim(i)
       wstar(i)  = 0.
       wstar3(i) = 0.
     endif
     ust3(i)   = ust(i)**3.
     wscale(i) = (ust3(i)+phifac*karman*wstar3(i)*0.5)**h1
     wscale(i) = min(wscale(i),ust(i)*aphi16)
     wscale(i) = max(wscale(i),ust(i)/aphi5)
   enddo
!
!     compute the surface variables for pbl height estimation
!     under unstable conditions
!
   do i = its,ite
     if(sfcflg(i).and.sflux(i).gt.0.0)then
       gamfac   = bfac/rhox(i)/wscale(i)
       hgamt(i) = min(gamfac*hfx(i)/cp,gamcrt)
       hgamq(i) = min(gamfac*qfx(i),gamcrq)
       vpert = (hgamt(i)+ep1*thx(i,1)*hgamq(i))/bfac*afac
       thermal(i) = thermal(i)+max(vpert,0.)
       hgamt(i) = max(hgamt(i),0.0)
       hgamq(i) = max(hgamq(i),0.0)
       brint    = -15.9*ust(i)*ust(i)/wspd(i)*wstar3(i)/(wscale(i)**4.)
       hgamu(i) = brint*ux(i,1)
       hgamv(i) = brint*vx(i,1)
     else
       pblflg(i) = .false.
     endif
   enddo
!
!     enhance the pbl height by considering the thermal
!
   do i = its,ite
     if(pblflg(i))then
       kpbl(i) = 1
       hpbl(i) = zq(i,1)
     endif
   enddo
!
   do i = its,ite
     if(pblflg(i))then
       stable(i) = .false.
       brup(i) = br(i)
       brcr(i) = brcr_ub
     endif
   enddo
!
   do k = 2,klpbl
     do i = its,ite
       if(.not.stable(i).and.pblflg(i))then
         brdn(i) = brup(i)
         spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
         brup(i) = (thvx(i,k)-thermal(i))*(g*za(i,k)/thvx(i,1))/spdk2
         kpbl(i) = k
         stable(i) = brup(i).gt.brcr(i)
       endif
     enddo
   enddo
!
   do i = its,ite
     if(pblflg(i)) then
       k = kpbl(i)
       if(brdn(i).ge.brcr(i))then
         brint = 0.
       elseif(brup(i).le.brcr(i))then
         brint = 1.
       else
         brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
       endif
       hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
       if(hpbl(i).lt.zq(i,2)) kpbl(i) = 1
       if(kpbl(i).le.1) pblflg(i) = .false.
     endif
   enddo
!
!     stable boundary layer
!
   do i = its,ite
     if((.not.sfcflg(i)).and.hpbl(i).lt.zq(i,2)) then
       brup(i) = br(i)
       stable(i) = .false.
     else
       stable(i) = .true.
     endif
   enddo
!
   do i = its,ite
     if((.not.stable(i)).and.((xland(i)-1.5).ge.0))then
       wspd10 = u10(i)*u10(i) + v10(i)*v10(i)
       wspd10 = sqrt(wspd10)
       ross = wspd10 / (cori*znt(i))
       brcr_sbro(i) = min(0.16*(1.e-7*ross)**(-0.18),.3)
     endif
   enddo
!
   do i = its,ite
     if(.not.stable(i))then
       if((xland(i)-1.5).ge.0)then
         brcr(i) = brcr_sbro(i)
       else
         brcr(i) = brcr_sb
       endif
     endif
   enddo
   do k = 2,klpbl
     do i = its,ite
       if(.not.stable(i))then
         brdn(i) = brup(i)
         spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
         brup(i) = (thvx(i,k)-thermal(i))*(g*za(i,k)/thvx(i,1))/spdk2
         kpbl(i) = k
         stable(i) = brup(i).gt.brcr(i)
       endif
     enddo
   enddo
!
   do i = its,ite
     if((.not.sfcflg(i)).and.hpbl(i).lt.zq(i,2)) then
       k = kpbl(i)
       if(brdn(i).ge.brcr(i))then
         brint = 0.
       elseif(brup(i).le.brcr(i))then
         brint = 1.
       else
         brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
       endif
       hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
       if(hpbl(i).lt.zq(i,2)) kpbl(i) = 1
       if(kpbl(i).le.1) pblflg(i) = .false.
     endif
   enddo
!
!     estimate the entrainment parameters
!
   do i = its,ite
     if(pblflg(i)) then
       k = kpbl(i) - 1
       prpbl(i) = 1.0
       wm3       = wstar3(i) + 5. * ust3(i)
       wm2(i)    = wm3**h2
       bfxpbl(i) = -0.15*thvx(i,1)/g*wm3/hpbl(i)
       dthvx(i)  = max(thvx(i,k+1)-thvx(i,k),tmin)
       dthx  = max(thx(i,k+1)-thx(i,k),tmin)
       dqx   = min(qx(i,k+1)-qx(i,k),0.0)
       we(i) = max(bfxpbl(i)/dthvx(i),-sqrt(wm2(i)))
       hfxpbl(i) = we(i)*dthx
       qfxpbl(i) = we(i)*dqx
!
       dux = ux(i,k+1)-ux(i,k)
       dvx = vx(i,k+1)-vx(i,k)
       if(dux.gt.tmin) then
         ufxpbl(i) = max(prpbl(i)*we(i)*dux,-ust(i)*ust(i))
       elseif(dux.lt.-tmin) then
         ufxpbl(i) = min(prpbl(i)*we(i)*dux,ust(i)*ust(i))
       else
         ufxpbl(i) = 0.0
       endif
       if(dvx.gt.tmin) then
         vfxpbl(i) = max(prpbl(i)*we(i)*dvx,-ust(i)*ust(i))
       elseif(dvx.lt.-tmin) then
         vfxpbl(i) = min(prpbl(i)*we(i)*dvx,ust(i)*ust(i))
       else
         vfxpbl(i) = 0.0
       endif
       delb  = govrth(i)*d3*hpbl(i)
       delta(i) = min(d1*hpbl(i) + d2*wm2(i)/delb,100.)
     endif
   enddo
!
   do k = kts,klpbl
     do i = its,ite
       if(pblflg(i).and.k.ge.kpbl(i))then
         entfac(i,k) = ((zq(i,k+1)-hpbl(i))/delta(i))**2.
       else
         entfac(i,k) = 1.e30
       endif
     enddo
   enddo
!
!     compute diffusion coefficients below pbl
!
   do k = kts,klpbl
     do i = its,ite
       if(k.lt.kpbl(i)) then
         zfac(i,k) = min(max((1.-(zq(i,k+1)-zl1(i))/(hpbl(i)-zl1(i))),zfmin),1.)
         xkzo = ckz*dza(i,k+1)
         zfacent(i,k) = (1.-zfac(i,k))**3.
         wscalek(i,k) = (ust3(i)+phifac*karman*wstar3(i)*(1.-zfac(i,k)))**h1
         if(sfcflg(i)) then 
           prfac = conpr/phim(i)/((1.+4.*karman*wstar3(i)/ust3(i)))**h1
           prnumfac = -3.*(max(zq(i,k+1)-sfcfrac*hpbl(i),0.))**2./hpbl(i)**2.
         else
           prfac = 0.
           prnumfac = 0.
           phim8z = 1.+(phim(i)-1.)*zq(i,k+1)/sfcfrac/hpbl(i)
           wscalek(i,k) = ust(i)/phim8z
           wscalek(i,k) = max(wscalek(i,k),ust(i)/aphi5)
         endif
         prnum0 = (phih(i)/phim(i)+prfac)
         prnum0 = min(prnum0,prmax)
         prnum0 = max(prnum0,prmin)
         xkzm(i,k) = xkzo+wscalek(i,k)*karman*zq(i,k+1)*zfac(i,k)**pfac
         prnum =  1. + (prnum0-1.)*exp(prnumfac)
         xkzq(i,k) = xkzm(i,k)/prnum*zfac(i,k)**(pfac_q-pfac)
         prnum0 = prnum0/(1.+prfac)
         prnum =  1. + (prnum0-1.)*exp(prnumfac)
         xkzh(i,k) = xkzm(i,k)/prnum
         xkzm(i,k) = min(xkzm(i,k),xkzmax)
         xkzm(i,k) = max(xkzm(i,k),xkzmin)
         xkzh(i,k) = min(xkzh(i,k),xkzmax)
         xkzh(i,k) = max(xkzh(i,k),xkzmin)
         xkzq(i,k) = min(xkzq(i,k),xkzmax)
         xkzq(i,k) = max(xkzq(i,k),xkzmin)
       endif
     enddo
   enddo
!
!     compute diffusion coefficients over pbl (free atmosphere)
!
   do k = kts,kte-1
     do i = its,ite
       xkzo = ckz*dza(i,k+1)
       if(k.ge.kpbl(i)) then
         ss = ((ux(i,k+1)-ux(i,k))*(ux(i,k+1)-ux(i,k))                         &
              +(vx(i,k+1)-vx(i,k))*(vx(i,k+1)-vx(i,k)))                        &
              /(dza(i,k+1)*dza(i,k+1))+1.e-9
         govrthv = g/(0.5*(thvx(i,k+1)+thvx(i,k)))
         ri = govrthv*(thvx(i,k+1)-thvx(i,k))/(ss*dza(i,k+1))
         if(imvdif.eq.1.and.ndiff.ge.3)then
           if((qx(i,ktrace2+k)+qx(i,ktrace3+k)).gt.0.01e-3.and.(qx(i           &
             ,ktrace2+k+1)+qx(i,ktrace3+k+1)).gt.0.01e-3)then
!      in cloud
             qmean = 0.5*(qx(i,k)+qx(i,k+1))
             tmean = 0.5*(tx(i,k)+tx(i,k+1))
             alph  = xlv*qmean/rd/tmean
             chi   = xlv*xlv*qmean/cp/rv/tmean/tmean
             ri    = (1.+alph)*(ri-g*g/ss/tmean/cp*((chi-alph)/(1.+chi)))
           endif
         endif
         zk = karman*zq(i,k+1)
         rlamdz = min(max(0.1*dza(i,k+1),rlam),300.)
         rl2 = (zk*rlamdz/(rlamdz+zk))**2
         dk = rl2*sqrt(ss)
         if(ri.lt.0.)then
! unstable regime
           sri = sqrt(-ri)
           xkzm(i,k) = xkzo+dk*(1+8.*(-ri)/(1+1.746*sri))
           xkzh(i,k) = xkzo+dk*(1+8.*(-ri)/(1+1.286*sri))
         else
! stable regime
           xkzh(i,k) = xkzo+dk/(1+5.*ri)**2
           prnum = 1.0+2.1*ri
           prnum = min(prnum,prmax)
           xkzm(i,k) = (xkzh(i,k)-xkzo)*prnum+xkzo
         endif
!
         xkzm(i,k) = min(xkzm(i,k),xkzmax)
         xkzm(i,k) = max(xkzm(i,k),xkzmin)
         xkzh(i,k) = min(xkzh(i,k),xkzmax)
         xkzh(i,k) = max(xkzh(i,k),xkzmin)
         xkzml(i,k) = xkzm(i,k)
         xkzhl(i,k) = xkzh(i,k)
       endif
     enddo
   enddo
!
!     compute tridiagonal matrix elements for heat
!
   do k = kts,kte
     do i = its,ite
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
       f1(i,k) = 0.
     enddo
   enddo
!
   do i = its,ite
     ad(i,1) = 1.
     f1(i,1) = thx(i,1)-300.+hfx(i)/cont/del(i,1)*dt2
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       dtodsd = dt2/del(i,k)
       dtodsu = dt2/del(i,k+1)
       dsig   = p2d(i,k)-p2d(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*xkzh(i,k)*rdz
       if(pblflg(i).and.k.lt.kpbl(i)) then
         dsdzt = tem1*(-hgamt(i)/hpbl(i)-hfxpbl(i)*zfacent(i,k)/xkzh(i,k))
         f1(i,k)   = f1(i,k)+dtodsd*dsdzt
         f1(i,k+1) = thx(i,k+1)-300.-dtodsu*dsdzt
       elseif(pblflg(i).and.k.ge.kpbl(i).and.entfac(i,k).lt.4.6) then
         xkzh(i,k) = -we(i)*dza(i,kpbl(i))*exp(-entfac(i,k))
         xkzh(i,k) = sqrt(xkzh(i,k)*xkzhl(i,k))
         xkzh(i,k) = min(xkzh(i,k),xkzmax)
         xkzh(i,k) = max(xkzh(i,k),xkzmin)
         f1(i,k+1) = thx(i,k+1)-300.
       else
         f1(i,k+1) = thx(i,k+1)-300.
       endif
       tem1   = dsig*xkzh(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2
       al(i,k)   = -dtodsu*dsdz2
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
       exch_hx(i,k+1) = xkzh(i,k)
     enddo
   enddo
!
! copies here to avoid duplicate input args for tridin
!
   do k = kts,kte
     do i = its,ite
       cu(i,k) = au(i,k)
       r1(i,k) = f1(i,k)
     enddo
   enddo
!
   call tridin_ysu(al,ad,cu,r1,au,f1,its,ite,kts,kte,1)
!
!     recover tendencies of heat
!
   do k = kte,kts,-1
     do i = its,ite
       ttend = (f1(i,k)-thx(i,k)+300.)*rdt*pi2d(i,k)
       ttnp(i,k) = ttnp(i,k)+ttend
       dtsfc(i) = dtsfc(i)+ttend*cont*del(i,k)/pi2d(i,k)
     enddo
   enddo
!
!     compute tridiagonal matrix elements for moisture, clouds, and tracers
!
   do k = kts,kte
     do i = its,ite
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
     enddo
   enddo
!
   do ic = 1,ndiff
     do i = its,ite
       do k = kts,kte
         f3(i,k,ic) = 0.
       enddo
     enddo
   enddo
!
   do i = its,ite
     ad(i,1) = 1.
     f3(i,1,1) = qx(i,1)+qfx(i)*g/del(i,1)*dt2
   enddo
!
   if(ndiff.ge.2) then
     do ic = 2,ndiff
       is = (ic-1) * kte
       do i = its,ite
         f3(i,1,ic) = qx(i,1+is)
       enddo
     enddo
   endif
!
   do k = kts,kte
     do i = its,ite
       if(k.ge.kpbl(i)) then
         xkzq(i,k) = xkzh(i,k)
       endif
     enddo
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       dtodsd = dt2/del(i,k)
       dtodsu = dt2/del(i,k+1)
       dsig   = p2d(i,k)-p2d(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*xkzq(i,k)*rdz
       if(pblflg(i).and.k.lt.kpbl(i)) then
         dsdzq = tem1*(-qfxpbl(i)*zfacent(i,k)/xkzq(i,k))
         f3(i,k,1) = f3(i,k,1)+dtodsd*dsdzq
         f3(i,k+1,1) = qx(i,k+1)-dtodsu*dsdzq
       elseif(pblflg(i).and.k.ge.kpbl(i).and.entfac(i,k).lt.4.6) then
         xkzq(i,k) = -we(i)*dza(i,kpbl(i))*exp(-entfac(i,k))
         xkzq(i,k) = sqrt(xkzq(i,k)*xkzhl(i,k))
         xkzq(i,k) = min(xkzq(i,k),xkzmax)
         xkzq(i,k) = max(xkzq(i,k),xkzmin)
         f3(i,k+1,1) = qx(i,k+1)
       else
         f3(i,k+1,1) = qx(i,k+1)
       endif
       tem1   = dsig*xkzq(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2
       al(i,k)   = -dtodsu*dsdz2
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
!      exch_hx(i,k+1) = xkzh(i,k)
     enddo
   enddo
!
   if(ndiff.ge.2) then
     do ic = 2,ndiff
       is = (ic-1) * kte
       do k = kts,kte-1
         do i = its,ite
           f3(i,k+1,ic) = qx(i,k+1+is)
         enddo
       enddo
     enddo
   endif
!
! copies here to avoid duplicate input args for tridin
!
   do k = kts,kte
     do i = its,ite
       cu(i,k) = au(i,k)
     enddo
   enddo
!
   do ic = 1,ndiff
     do k = kts,kte
       do i = its,ite
         r3(i,k,ic) = f3(i,k,ic)
       enddo
     enddo
   enddo
!
!     solve tridiagonal problem for moisture, clouds, and tracers
!
   call tridin_ysu(al,ad,cu,r3,au,f3,its,ite,kts,kte,ndiff)
!
!     recover tendencies of heat and moisture
!
   do k = kte,kts,-1
     do i = its,ite
       qtend = (f3(i,k,1)-qx(i,k))*rdt
       qtnp(i,k) = qtnp(i,k)+qtend
       dqsfc(i) = dqsfc(i)+qtend*conq*del(i,k)
     enddo
   enddo
!
   if(ndiff.ge.2) then
     do ic = 2,ndiff
       is = (ic-1) * kte
       do k = kte,kts,-1
         do i = its,ite
           qtend = (f3(i,k,ic)-qx(i,k+is))*rdt
           qtnp(i,k+is) = qtnp(i,k+is)+qtend
         enddo
       enddo
     enddo
   endif
!
!     compute tridiagonal matrix elements for momentum
!
   do i = its,ite
     do k = kts,kte
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
       f1(i,k) = 0.
       f2(i,k) = 0.
     enddo
   enddo
!
   do i = its,ite
     ad(i,1) = 1.+ust(i)**2/wspd1(i)*rhox(i)*g/del(i,1)*dt2                    &
              *(wspd1(i)/wspd(i))**2
     f1(i,1) = ux(i,1)
     f2(i,1) = vx(i,1)
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       dtodsd = dt2/del(i,k)
       dtodsu = dt2/del(i,k+1)
       dsig   = p2d(i,k)-p2d(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*xkzm(i,k)*rdz
     if(pblflg(i).and.k.lt.kpbl(i))then
       dsdzu     = tem1*(-hgamu(i)/hpbl(i)-ufxpbl(i)*zfacent(i,k)/xkzm(i,k))
       dsdzv     = tem1*(-hgamv(i)/hpbl(i)-vfxpbl(i)*zfacent(i,k)/xkzm(i,k))
       f1(i,k)   = f1(i,k)+dtodsd*dsdzu
       f1(i,k+1) = ux(i,k+1)-dtodsu*dsdzu
       f2(i,k)   = f2(i,k)+dtodsd*dsdzv
       f2(i,k+1) = vx(i,k+1)-dtodsu*dsdzv
     elseif(pblflg(i).and.k.ge.kpbl(i).and.entfac(i,k).lt.4.6) then
       xkzm(i,k) = prpbl(i)*xkzh(i,k)
       xkzm(i,k) = sqrt(xkzm(i,k)*xkzml(i,k))
       xkzm(i,k) = min(xkzm(i,k),xkzmax)
       xkzm(i,k) = max(xkzm(i,k),xkzmin)
       f1(i,k+1) = ux(i,k+1)
       f2(i,k+1) = vx(i,k+1)
     else
       f1(i,k+1) = ux(i,k+1)
       f2(i,k+1) = vx(i,k+1)
     endif
       tem1   = dsig*xkzm(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2
       al(i,k)   = -dtodsu*dsdz2
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
       exch_mx(i,k+1) = xkzm(i,k)
     enddo
   enddo
!
! copies here to avoid duplicate input args for tridin
!
   do k = kts,kte
     do i = its,ite
       cu(i,k) = au(i,k)
       r1(i,k) = f1(i,k)
       r2(i,k) = f2(i,k)
     enddo
   enddo
!
!     solve tridiagonal problem for momentum
!
   call tridi1n(al,ad,cu,r1,r2,au,f1,f2,its,ite,kts,kte,1)
!
!     recover tendencies of momentum
!
   do k = kte,kts,-1
     do i = its,ite
       utend = (f1(i,k)-ux(i,k))*rdt
       vtend = (f2(i,k)-vx(i,k))*rdt
       utnp(i,k) = utnp(i,k)+utend
       vtnp(i,k) = vtnp(i,k)+vtend
       dusfc(i) = dusfc(i) + utend*conwrc*del(i,k)
       dvsfc(i) = dvsfc(i) + vtend*conwrc*del(i,k)
     enddo
   enddo
!
!---- end of vertical diffusion
!
   do i = its,ite
     kpbl1d(i) = kpbl(i)
   enddo
!
   end subroutine ysu2d
!
   subroutine tridi1n(cl,cm,cu,r1,r2,au,f1,f2,its,ite,kts,kte,nt)
!----------------------------------------------------------------
   implicit none
!----------------------------------------------------------------
!
   integer, intent(in )      ::     its,ite, kts,kte, nt
!
   real, dimension( its:ite, kts+1:kte+1 )                                   , &
         intent(in   )  ::                                                 cl
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(in   )  ::                                                 cm, &
						                           r1
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(in   )  ::                                                 r2
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(inout)  ::                                                 au, &
						                           cu, &
							                   f1
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(inout)  ::                                                 f2
!
   real    :: fk
   integer :: i,k,l,n,it
!
!----------------------------------------------------------------
!
   l = ite
   n = kte
!
   do i = its,l
     fk = 1./cm(i,1)
     au(i,1) = fk*cu(i,1)
     f1(i,1) = fk*r1(i,1)
   enddo
   do it = 1,nt
     do i = its,l
       fk = 1./cm(i,1)
       f2(i,1,it) = fk*r2(i,1,it)
     enddo
   enddo
   do k = kts+1,n-1
     do i = its,l
       fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
       au(i,k) = fk*cu(i,k)
       f1(i,k) = fk*(r1(i,k)-cl(i,k)*f1(i,k-1))
     enddo
   enddo
   do it = 1,nt
   do k = kts+1,n-1
     do i = its,l
       fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
       f2(i,k,it) = fk*(r2(i,k,it)-cl(i,k)*f2(i,k-1,it))
     enddo
   enddo
   enddo
   do i = its,l
     fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
     f1(i,n) = fk*(r1(i,n)-cl(i,n)*f1(i,n-1))
   enddo
   do it = 1,nt
   do i = its,l
     fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
     f2(i,n,it) = fk*(r2(i,n,it)-cl(i,n)*f2(i,n-1,it))
   enddo
   enddo
   do k = n-1,kts,-1
     do i = its,l
       f1(i,k) = f1(i,k)-au(i,k)*f1(i,k+1)
     enddo
   enddo
   do it = 1,nt
   do k = n-1,kts,-1
     do i = its,l
       f2(i,k,it) = f2(i,k,it)-au(i,k)*f2(i,k+1,it)
     enddo
   enddo
   enddo
!
   end subroutine tridi1n
!
   subroutine tridin_ysu(cl,cm,cu,r2,au,f2,its,ite,kts,kte,nt)
!----------------------------------------------------------------
   implicit none
!----------------------------------------------------------------
!
   integer, intent(in )      ::     its,ite, kts,kte, nt
!
   real, dimension( its:ite, kts+1:kte+1 )                                   , &
         intent(in   )  ::                                                 cl
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(in   )  ::                                                 cm
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(in   )  ::                                                 r2
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(inout)  ::                                                 au, &
						                           cu
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(inout)  ::                                                 f2
!
   real    :: fk
   integer :: i,k,l,n,it
!
!----------------------------------------------------------------
!
   l = ite
   n = kte
!
   do it = 1,nt
     do i = its,l
       fk = 1./cm(i,1)
       au(i,1) = fk*cu(i,1)
       f2(i,1,it) = fk*r2(i,1,it)
     enddo
   enddo
   do it = 1,nt
   do k = kts+1,n-1
     do i = its,l
       fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
       au(i,k) = fk*cu(i,k)
       f2(i,k,it) = fk*(r2(i,k,it)-cl(i,k)*f2(i,k-1,it))
     enddo
   enddo
   enddo
   do it = 1,nt
   do i = its,l
     fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
     f2(i,n,it) = fk*(r2(i,n,it)-cl(i,n)*f2(i,n-1,it))
   enddo
   enddo
   do it = 1,nt
   do k = n-1,kts,-1
     do i = its,l
       f2(i,k,it) = f2(i,k,it)-au(i,k)*f2(i,k+1,it)
     enddo
   enddo
   enddo
!
   end subroutine tridin_ysu
!
   subroutine ysuinit(rublten,rvblten,rthblten,rqvblten,                       &
                      rqcblten,rqiblten,p_qi,p_first_scalar,                   &
                      restart, allowed_to_read,                                &
                      ids, ide, jds, jde, kds, kde,                            &
                      ims, ime, jms, jme, kms, kme,                            &
                      its, ite, jts, jte, kts, kte                 )
!-------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------
!
   logical , intent(in)          :: restart, allowed_to_read
   integer , intent(in)          ::  ids, ide, jds, jde, kds, kde,             &
                                     ims, ime, jms, jme, kms, kme,             &
                                     its, ite, jts, jte, kts, kte
   integer , intent(in)          ::  p_qi,p_first_scalar
   real , dimension( ims:ime , kms:kme , jms:jme ), intent(out) ::             &
                                                                      rublten, &
                                                                      rvblten, &
                                                                     rthblten, &
                                                                     rqvblten, &
                                                                     rqcblten, &
                                                                     rqiblten
   integer :: i, j, k, itf, jtf, ktf
!
   jtf = min0(jte,jde-1)
   ktf = min0(kte,kde-1)
   itf = min0(ite,ide-1)
!
   if(.not.restart)then
     do j = jts,jtf
     do k = kts,ktf
     do i = its,itf
        rublten(i,k,j) = 0.
        rvblten(i,k,j) = 0.
        rthblten(i,k,j) = 0.
        rqvblten(i,k,j) = 0.
        rqcblten(i,k,j) = 0.
     enddo
     enddo
     enddo
   endif
!
   if (p_qi .ge. p_first_scalar .and. .not.restart) then
      do j = jts,jtf
      do k = kts,ktf
      do i = its,itf
         rqiblten(i,k,j) = 0.
      enddo
      enddo
      enddo
   endif
!
   end subroutine ysuinit
!-------------------------------------------------------------------
end module module_bl_ysu
