!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE IRRAD                      ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE irrad(idim,jdim,m,n,np,                                      &
           taucl,ccld,pl,ta,wa,oa,co2,ts,                               &
           high,radlwin,flx,flc,dfdts,st4,                              &
           fclr,dbs,trant,th2o,tcon,tco2,                               &
           pa,dt,sh2o,swpre,swtem,sco3,scopre,scotem,                   &
           dh2o,dcont,dco2,do3,flxu,flxd,clr,blayer,                    &
           h2oexp,conexp,co2exp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate IR fluxes due to water vapor, co2, and o3. Clouds in
!  different layers are assumed randomly overlapped.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Adopted the original code and formatted it in accordance with the
!  ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!******************** CLIRAD IR1  Date: Oct. 17, 1994 ****************
!*********************************************************************
!
! This routine computes ir fluxes due to water vapor, co2, and o3.
!   Clouds in different layers are assumed randomly overlapped.
!
! This is a vectorized code.  It computes fluxes simultaneously for
!   (m x n) soundings, which is a subset of (m x ndim) soundings.
!   In a global climate model, m and ndim correspond to the numbers of
!   grid boxes in the zonal and meridional directions, respectively.
!
! Detailed description of the radiation routine is given in
!   Chou and Suarez (1994).
!
! There are two options for computing cooling rate profiles.
!
!   if high = .true., transmission functions in the co2, o3, and the
!   three water vapor bands with strong absorption are computed using
!   table look-up.  cooling rates are computed accurately from the
!   surface up to 0.01 mb.
!   if high = .false., transmission functions are computed using the
!   k-distribution method with linear pressure scaling.  cooling rates
!   are not calculated accurately for pressures less than 20 mb.
!   the computation is faster with high=.false. than with high=.true.
!
! The IR spectrum is divided into eight bands:
!
!   bnad     wavenumber (/cm)   absorber         method
!
!    1           0 - 340           h2o            K/T
!    2         340 - 540           h2o            K/T
!    3         540 - 800       h2o,cont,co2       K,S,K/T
!    4         800 - 980       h2o,cont           K,S
!    5         980 - 1100      h2o,cont,o3        K,S,T
!    6        1100 - 1380      h2o,cont           K,S
!    7        1380 - 1900          h2o            K/T
!    8        1900 - 3000          h2o            K
!
! Note : "h2o" for h2o line absorption
!     "cont" for h2o continuum absorption
!     "K" for k-distribution method
!     "S" for one-parameter temperature scaling
!     "T" for table look-up
!
! The 15 micrometer region (540-800/cm) is further divided into
!   3 sub-bands :
!
!   subbnad   wavenumber (/cm)
!
!    1          540 - 620
!    2          620 - 720
!    3          720 - 800
!
!---- Input parameters                               units    size
!
!   number of soundings in zonal direction (m)        n/d      1
!   number of soundings in meridional direction (n)   n/d      1
!   maximum number of soundings in
!              meridional direction (ndim)         n/d      1
!   number of atmospheric layers (np)                 n/d      1
!   cloud optical thickness (taucl)                   n/d     m*ndim*np
!   cloud cover (ccld)                              fraction  m*ndim*np
!   level pressure (pl)                               mb      m*ndim*(np+1)
!   layer temperature (ta)                            k       m*ndim*np
!   layer specific humidity (wa)                      g/g     m*ndim*np
!   layer ozone mixing ratio by mass (oa)             g/g     m*ndim*np
!   surface temperature (ts)                          k       m*ndim
!   co2 mixing ratio by volumn (co2)                  pppv     1
!   high                                                       1
!
! pre-computed tables used in table look-up for transmittance calculations:
!
!   c1 , c2, c3: for co2 (band 3)
!   o1 , o2, o3: for  o3 (band 5)
!   h11,h12,h13: for h2o (band 1)
!   h21,h22,h23: for h2o (band 2)
!   h71,h72,h73: for h2o (band 7)
!
!---- output parameters
!
!   net downward flux, all-sky   (flx)             w/m**2     m*ndim*(np+1)
!   net downward flux, clear-sky (flc)             w/m**2     m*ndim*(np+1)
!   sensitivity of net downward flux
!    to surface temperature (dfdts)             w/m**2/k   m*ndim*(np+1)
!   emission by the surface (st4)                  w/m**2     m*ndim
!
! Notes:
!
!   (1)  Water vapor continuum absorption is included in 540-1380 /cm.
!   (2)  Scattering by clouds is not included.
!   (3)  Clouds are assumed "gray" bodies.
!   (4)  The diffuse cloud transmission is computed to be exp(-1.66*taucl).
!   (5)  If there are no clouds, flx=flc.
!   (6)  plevel(1) is the pressure at the top of the model atmosphere, and
!     plevel(np+1) is the surface pressure.
!
!    ARPS note: pl was replaced by pa at scalar points (layers)
!
!   (7)  Downward flux is positive, and upward flux is negative.
!   (8)  dfdts is always negative because upward flux is defined as negative.
!   (9)  For questions and coding errors, please contact with Ming-Dah Chou,
!     Code 913, NASA/Goddard Space Flight Center, Greenbelt, MD 20771.
!     Phone: 301-286-4012, Fax: 301-286-1759,
!     e-mail: chou@climate.gsfc.nasa.gov
!
!-----parameters defining the size of the pre-computed tables for transmittance
!  calculations using table look-up.
!
!  "nx" is the number of intervals in pressure
!  "no" is the number of intervals in o3 amount
!  "nc" is the number of intervals in co2 amount
!  "nh" is the number of intervals in h2o amount
!  "nt" is the number of copies to be made from the pre-computed
!       transmittance tables to reduce "memory-bank conflict"
!       in parallel machines and, hence, enhancing the speed of
!       computations using table look-up.
!       If such advantage does not exist, "nt" can be set to 1.
!***************************************************************************
!
!fpp$ expand (expmn)
!!dir$ inline always expmn
!*$*  inline routine (expmn)
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: idim,jdim

!!!  INTEGER :: nx,no,nc,nh,nt,m,n,np
  INTEGER :: m,n,np
  INTEGER :: i,j,k,ip,iw,it,ib,ik,iq,isb,k1,k2
!!!  PARAMETER (nx=26,no=21,nc=24,nh=31,nt=7)

!---- input parameters ------

  REAL :: taucl(idim,jdim,np)
  REAL :: ccld(idim,jdim,np)
  REAL :: pl(idim,jdim,np+1)
  REAL :: ta(idim,jdim,np)
  REAL :: wa(idim,jdim,np)
  REAL :: oa(idim,jdim,np)
  REAL :: ts(idim,jdim)
  REAL :: co2

  LOGICAL :: high

!---- output parameters ------

  REAL :: radlwin(idim,jdim) 
  REAL :: flx(idim,jdim,np+1)
  REAL :: flc(idim,jdim,np+1)
  REAL :: dfdts(idim,jdim,np+1)
  REAL :: st4(idim,jdim)
!
!-----------------------------------------------------------------------
!
!  Temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL :: fclr(m,n)
  REAL :: dbs(m,n)
  REAL :: trant(m,n)
  REAL :: th2o(m,n,6)
  REAL :: tcon(m,n,3)
  REAL :: tco2(m,n,6,2)

  REAL :: pa(m,n,np)
  REAL :: dt(m,n,np)
  REAL :: sh2o(m,n,np+1)
  REAL :: swpre(m,n,np+1)
  REAL :: swtem(m,n,np+1)
  REAL :: sco3(m,n,np+1)
  REAL :: scopre(m,n,np+1)
  REAL :: scotem(m,n,np+1)
  REAL :: dh2o(m,n,np)
  REAL :: dcont(m,n,np)
  REAL :: dco2(m,n,np)
  REAL :: do3(m,n,np)
  REAL :: flxu(m,n,np+1)
  REAL :: flxd(m,n,np+1)
  REAL :: clr(m,n,0:np+1)
  REAL :: blayer(m,n,0:np+1)

  REAL :: h2oexp(m,n,np,6)
  REAL :: conexp(m,n,np,3)

  REAL :: co2exp(m,n,np,6,2)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  LOGICAL :: oznbnd
  LOGICAL :: co2bnd
  LOGICAL :: h2otbl
  LOGICAL :: conbnd

!!!  REAL :: c1 (nx,nc,nt),c2 (nx,nc,nt),c3 (nx,nc,nt)
!!!  REAL :: o1 (nx,no,nt),o2 (nx,no,nt),o3 (nx,no,nt)
!!!  REAL :: h11(nx,nh,nt),h12(nx,nh,nt),h13(nx,nh,nt)
!!!  REAL :: h21(nx,nh,nt),h22(nx,nh,nt),h23(nx,nh,nt)
!!!  REAL :: h71(nx,nh,nt),h72(nx,nh,nt),h73(nx,nh,nt)
  include 'irrad.incl'

  REAL :: xx, w1,p1,dwe,dpe

!---- static data -----

  REAL :: cb(5,8)

!-----the following coefficients (table 2 of chou and suarez, 1995)
!  are for computing spectrally integtrated planck fluxes of
!  the 8 bands using eq. (22)

  DATA cb/                                                              &
      -2.6844E-1,-8.8994E-2, 1.5676E-3,-2.9349E-6, 2.2233E-9,           &
      3.7315E+1,-7.4758E-1, 4.6151E-3,-6.3260E-6, 3.5647E-9,            &
      3.7187E+1,-3.9085E-1,-6.1072E-4, 1.4534E-5,-1.6863E-8,            &
      -4.1928E+1, 1.0027E+0,-8.5789E-3, 2.9199E-5,-2.5654E-8,           &
      -4.9163E+1, 9.8457E-1,-7.0968E-3, 2.0478E-5,-1.5514E-8,           &
      -1.0345E+2, 1.8636E+0,-1.1753E-2, 2.7864E-5,-1.1998E-8,           &
      -6.9233E+0,-1.5878E-1, 3.9160E-3,-2.4496E-5, 4.9301E-8,           &
      1.1483E+2,-2.2376E+0, 1.6394E-2,-5.3672E-5, 6.6456E-8/

!-----copy tables to enhance the speed of co2 (band 3), o3 (band5),
!  and h2o (bands 1, 2, and 7 only) transmission calculations
!  using table look-up.

!!!  LOGICAL :: first
!!!  SAVE first
!!!  DATA first /.true./
!
!-----------------------------------------------------------------------
!
!  Functions:
!
!-----------------------------------------------------------------------
!
  REAL :: expmn
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
!  include "h2o.tran3"
!  include "co2.tran3"
!  include "o3.tran3"

!!!  COMMON /radtab001/ h11,h12,h13,h21,h22,h23,h71,h72,h73
!!!  COMMON /radtab002/ c1,c2,c3
!!!  COMMON /radtab003/ o1,o2,o3
!
!-----------------------------------------------------------------------
!
!  Save variables:
!
!-----------------------------------------------------------------------
!
!  save c1,c2,c3,o1,o2,o3
!  save h11,h12,h13,h21,h22,h23,h71,h72,h73
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!----- moved this section of code to new initrad subroutine -----!

!-----compute layer pressure (pa) and layer temperature minus 250K (dt)

  DO k=1,np
    DO j=1,n
      DO i=1,m
        dt(i,j,k)=ta(i,j,k)-250.0
!!!     GHB, 101129:  pa is now given directly through radtrans3d
!!!        pa(i,j,k)=0.5*(pl(i,j,k)+pl(i,j,k+1))
!!!        print *,'  old,new,error:',k,0.5*(pl(i,j,k)+pl(i,j,k+1)),pa(i,j,k),   &
!!!                                     0.5*(pl(i,j,k)+pl(i,j,k+1))-pa(i,j,k)
      END DO
    END DO
  END DO

!-----compute layer absorber amount

!  dh2o : water vapor amount (g/cm**2)
!  dcont: scaled water vapor amount for continuum absorption (g/cm**2)
!  dco2 : co2 amount (cm-atm)stp
!  do3  : o3 amount (cm-atm)stp
!  the factor 1.02 is equal to 1000/980
!  factors 789 and 476 are for unit conversion
!  the factor 0.001618 is equal to 1.02/(.622*1013.25)
!  the factor 6.081 is equal to 1800/296

  DO k=1,np
    DO j=1,n
      DO i=1,m
        dh2o(i,j,k) = 1.02*wa(i,j,k)*(pl(i,j,k+1)-pl(i,j,k))+1.e-10
        dco2(i,j,k) = 789.*co2*(pl(i,j,k+1)-pl(i,j,k))+1.e-10
        do3 (i,j,k) = 476.0*oa(i,j,k)*(pl(i,j,k+1)-pl(i,j,k))+1.e-10

!-----compute scaled water vapor amount for h2o continuum absorption
!  following eq. (43).

        xx=pa(i,j,k)*0.001618*wa(i,j,k)*wa(i,j,k)                       &
            *(pl(i,j,k+1)-pl(i,j,k))
        dcont(i,j,k) = xx*expmn(1800./ta(i,j,k)-6.081)+1.e-10

!-----compute effective cloud-free fraction, clr, for each layer.
!  the cloud diffuse transmittance is approximated by using a
!  diffusivity factor of 1.66.

        clr(i,j,k)=1.0-(ccld(i,j,k)*(1.-expmn(-1.66*taucl(i,j,k))))

      END DO
    END DO
  END DO

!-----compute column-integrated h2o amoumt, h2o-weighted pressure
!  and temperature.  it follows eqs. (37) and (38).

  IF (high) THEN

    CALL column(m,n,np,pa,dt,dh2o,sh2o,swpre,swtem)

  END IF

!-----the surface (with an index np+1) is treated as a layer filled with
!  black clouds.

  DO j=1,n
    DO i=1,m
      clr(i,j,0)    = 1.0
      clr(i,j,np+1) = 0.0
      st4(i,j)      = 0.0
    END DO
  END DO

!-----initialize fluxes

  DO k=1,np+1
    DO j=1,n
      DO i=1,m
        flx(i,j,k)  = 0.0
        flc(i,j,k)  = 0.0
        dfdts(i,j,k)= 0.0
        flxu(i,j,k) = 0.0
        flxd(i,j,k) = 0.0
      END DO
    END DO
  END DO

!-----integration over spectral bands

  DO ib=1,8

!-----if h2otbl, compute h2o (line) transmittance using table look-up.
!  if conbnd, compute h2o (continuum) transmittance in bands 3, 4, 5 and 6.
!  if co2bnd, compute co2 transmittance in band 3.
!  if oznbnd, compute  o3 transmittance in band 5.

    h2otbl=high.AND.(ib == 1.OR.ib == 2.OR.ib == 7)
    conbnd=ib >= 3.AND.ib <= 6
    co2bnd=ib == 3
    oznbnd=ib == 5

!-----blayer is the spectrally integrated planck flux of the mean layer
!  temperature derived from eq. (22)
!  the fitting for the planck flux is valid in the range 160-345 K.

    DO k=1,np
      DO j=1,n
        DO i=1,m
          blayer(i,j,k)=ta(i,j,k)*(ta(i,j,k)*(ta(i,j,k)                 &
                       *(ta(i,j,k)*cb(5,ib)+cb(4,ib))+cb(3,ib))         &
                       +cb(2,ib))+cb(1,ib)
        END DO
      END DO
    END DO

!-----the earth's surface, with an index "np+1", is treated as a layer

    DO j=1,n
      DO i=1,m
        blayer(i,j,0)   = 0.0
        blayer(i,j,np+1)=ts(i,j)*(ts(i,j)*(ts(i,j)                      &
                        *(ts(i,j)*cb(5,ib)+cb(4,ib))+cb(3,ib))          &
                        +cb(2,ib))+cb(1,ib)

!-----dbs is the derivative of the surface planck flux with respect to
!  surface temperature (eq. 59).

        dbs(i,j)=ts(i,j)*(ts(i,j)*(ts(i,j)                              &
                *4.*cb(5,ib)+3.*cb(4,ib))+2.*cb(3,ib))+cb(2,ib)

      END DO
    END DO

!-----compute column-integrated absorber amoumt, absorber-weighted
!  pressure and temperature for co2 (band 3) and o3 (band 5).
!  it follows eqs. (37) and (38).

!-----this is in the band loop to save storage

    IF( high .AND. co2bnd) THEN

      CALL column(m,n,np,pa,dt,dco2,sco3,scopre,scotem)

    END IF

    IF(oznbnd) THEN

      CALL column(m,n,np,pa,dt,do3,sco3,scopre,scotem)

    END IF

!-----compute the exponential terms (eq. 32) at each layer for
!  water vapor line absorption when k-distribution is used

    IF( .NOT. h2otbl) THEN

      CALL h2oexps(ib,m,n,np,dh2o,pa,dt,h2oexp)

    END IF

!-----compute the exponential terms (eq. 46) at each layer for
!  water vapor continuum absorption

    IF( conbnd) THEN

      CALL conexps(ib,m,n,np,dcont,conexp)

    END IF


!-----compute the  exponential terms (eq. 32) at each layer for
!  co2 absorption

    IF( .NOT.high .AND. co2bnd) THEN

      CALL co2exps(m,n,np,dco2,pa,dt,co2exp)

    END IF

!-----compute transmittances for regions between levels k1 and k2
!  and update the fluxes at the two levels.

    DO k1=1,np

!-----initialize fclr, th2o, tcon, and tco2

      DO j=1,n
        DO i=1,m
          fclr(i,j)=1.0
        END DO
      END DO

!-----for h2o line absorption

      IF(.NOT. h2otbl) THEN
        DO ik=1,6
          DO j=1,n
            DO i=1,m
              th2o(i,j,ik)=1.0
            END DO
          END DO
        END DO
      END IF

!-----for h2o continuum absorption

      IF (conbnd) THEN
        DO iq=1,3
          DO j=1,n
            DO i=1,m
              tcon(i,j,iq)=1.0
            END DO
          END DO
        END DO
      END IF

!-----for co2 absorption when using k-distribution method.
!  band 3 is divided into 3 sub-bands, but sub-bands 3a and 3c
!  are combined in computing the co2 transmittance.

      IF (.NOT. high .AND. co2bnd) THEN
        DO isb=1,2
          DO ik=1,6
            DO j=1,n
              DO i=1,m
                tco2(i,j,ik,isb)=1.0
              END DO
            END DO
          END DO
        END DO
      END IF

!-----loop over the bottom level of the region (k2)

      DO k2=k1+1,np+1

        DO j=1,n
          DO i=1,m
            trant(i,j)=1.0
          END DO
        END DO

        IF(h2otbl) THEN

          w1=-8.0
          p1=-2.0
          dwe=0.3
          dpe=0.2

!-----compute water vapor transmittance using table look-up

          IF (ib == 1 ) THEN

            CALL tablup(k1,k2,m,n,np,nx,nh,nt,sh2o,swpre,swtem,         &
                        w1,p1,dwe,dpe,h11,h12,h13,trant)

          END IF
          IF (ib == 2 ) THEN

            CALL tablup(k1,k2,m,n,np,nx,nh,nt,sh2o,swpre,swtem,         &
                        w1,p1,dwe,dpe,h21,h22,h23,trant)

          END IF
          IF (ib == 7 ) THEN

            CALL tablup(k1,k2,m,n,np,nx,nh,nt,sh2o,swpre,swtem,         &
                        w1,p1,dwe,dpe,h71,h72,h73,trant)

          END IF

        ELSE

!-----compute water vapor transmittance using k-distribution.

          CALL wvkdis(ib,m,n,np,k2-1,h2oexp,conexp,th2o,tcon,trant)

        END IF

        IF(co2bnd) THEN

          IF( high ) THEN

!-----compute co2 transmittance using table look-up method

            w1=-4.0
            p1=-2.0
            dwe=0.3
            dpe=0.2

            CALL tablup(k1,k2,m,n,np,nx,nc,nt,sco3,scopre,scotem,       &
                        w1,p1,dwe,dpe,c1,c2,c3,trant)

          ELSE

!-----compute co2 transmittance using k-distribution method

            CALL co2kdis(m,n,np,k2-1,co2exp,tco2,trant)

          END IF

        END IF

!-----compute o3 transmittance using table look-up

        IF (oznbnd) THEN

          w1=-6.0
          p1=-2.0
          dwe=0.3
          dpe=0.2

          CALL tablup(k1,k2,m,n,np,nx,no,nt,sco3,scopre,scotem,         &
                      w1,p1,dwe,dpe,o1,o2,o3,trant)

        END IF

!-----fclr is the clear line-of-sight between levels k1 and k2.
!  in computing fclr, clouds are assumed randomly overlapped
!  using eq. (10).

        DO j=1,n
          DO i=1,m
            fclr(i,j) = fclr(i,j)*clr(i,j,k2-1)
          END DO
        END DO

!-----compute upward and downward fluxes


!-----add "boundary" terms to the net downward flux.
!  these are the first terms on the right-hand-side of
!  eqs. (56a) and (56b).
!  downward fluxes are positive.

        IF (k2 == k1+1) THEN
          DO j=1,n
            DO i=1,m
              flc(i,j,k1)=flc(i,j,k1)-blayer(i,j,k1)
              flc(i,j,k2)=flc(i,j,k2)+blayer(i,j,k1)
            END DO
          END DO
        END IF

!-----add flux components involving the four layers above and below
!  the levels k1 and k2.  it follows eqs. (56a) and (56b).

        DO j=1,n
          DO i=1,m
            xx=trant(i,j)*(blayer(i,j,k2-1)-blayer(i,j,k2))
            flc(i,j,k1) =flc(i,j,k1)+xx
            xx=trant(i,j)*(blayer(i,j,k1-1)-blayer(i,j,k1))
            flc(i,j,k2) =flc(i,j,k2)+xx
          END DO
        END DO

!-----compute upward and downward fluxes for all-sky situation

        IF (k2 == k1+1) THEN
          DO j=1,n
            DO i=1,m
              flxu(i,j,k1)=flxu(i,j,k1)-blayer(i,j,k1)
              flxd(i,j,k2)=flxd(i,j,k2)+blayer(i,j,k1)
            END DO
          END DO
        END IF

        DO j=1,n
          DO i=1,m
            xx=trant(i,j)*(blayer(i,j,k2-1)-blayer(i,j,k2))
            flxu(i,j,k1) =flxu(i,j,k1)+xx*fclr(i,j)
            xx=trant(i,j)*(blayer(i,j,k1-1)-blayer(i,j,k1))
            flxd(i,j,k2) =flxd(i,j,k2)+xx*fclr(i,j)
          END DO
        END DO


      END DO

!-----compute the partial derivative of fluxes with respect to
!  surface temperature (eq. 59).

      DO j=1,n
        DO i=1,m
          dfdts(i,j,k1) =dfdts(i,j,k1)-dbs(i,j)*trant(i,j)*fclr(i,j)
        END DO
      END DO

    END DO

!-----add contribution from the surface to the flux terms at the surface.

    DO j=1,n
      DO i=1,m
        dfdts(i,j,np+1) =dfdts(i,j,np+1)-dbs(i,j)
      END DO
    END DO

    DO j=1,n
      DO i=1,m
        flc(i,j,np+1)=flc(i,j,np+1)-blayer(i,j,np+1)
        flxu(i,j,np+1)=flxu(i,j,np+1)-blayer(i,j,np+1)
        st4(i,j)=st4(i,j)-blayer(i,j,np+1)
      END DO
    END DO


!  write(7,3211) ib, flxd(1,1,52),flxu(1,1,52)
!  write(7,3211) ib, flxd(1,1,np+1),flxu(1,1,np+1)
!    3211 FORMAT ('ib, fluxd, fluxu=', i3,2F12.3)

  END DO

  DO k=1,np+1
    DO j=1,n
      DO i=1,m
        flx(i,j,k)=flxd(i,j,k)+flxu(i,j,k)
        radlwin(i,j) = flxd(i,j,np+1) 
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE irrad
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE H2OEXPS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE h2oexps(ib,m,n,np,dh2o,pa,dt,h2oexp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate exponentials for water vapor line absorption in
!  individual layers.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION:
!
!  03/11/1996 (Yuhe Liu)
!  Formatted code to ARPS standard format
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!   compute exponentials for water vapor line absorption
!   in individual layers.
!
!---- input parameters
!  spectral band (ib)
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of layers (np)
!  layer water vapor amount for line absorption (dh2o)
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!
!---- output parameters
!  6 exponentials for each layer  (h2oexp)
!
!**********************************************************************
!
!fpp$ expand (expmn)
!!dir$ inline always expmn
!*$*  inline routine (expmn)
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: ib,m,n,np

!---- input parameters ------

  REAL :: dh2o(m,n,np)
  REAL :: pa(m,n,np)
  REAL :: dt(m,n,np)

!---- output parameters -----

  REAL :: h2oexp(m,n,np,6)

!---- static data -----

  INTEGER :: mw(8)

  REAL :: xkw(8)
  REAL :: aw(8)
  REAL :: bw(8)

!---- temporary arrays -----

  REAL :: xh

!---- local misc. variables

  INTEGER :: i,j,k,ik

!-----xkw  are the absorption coefficients for the first
!  k-distribution function due to water vapor line absorption
!  (tables 4 and 7).  units are cm**2/g

  DATA xkw / 29.55  , 4.167E-1, 1.328E-2, 5.250E-4,                     &
              5.25E-4, 2.340E-3, 1.320E-0, 5.250E-4/

!-----mw are the ratios between neighboring absorption coefficients
!  for water vapor line absorption (tables 4 and 7).

  DATA mw /6,6,8,6,6,8,6,16/

!-----aw and bw (table 3) are the coefficients for temperature scaling
!  in eq. (25).

  DATA aw/ 0.0021, 0.0140, 0.0167, 0.0302,                              &
           0.0307, 0.0154, 0.0008, 0.0096/
  DATA bw/ -1.01E-5, 5.57E-5, 8.54E-5, 2.96E-4,                         &
            2.86E-4, 7.53E-5,-3.52E-6, 1.64E-5/

!-----expmn is an external function

  REAL :: expmn
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!**********************************************************************
!    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
!    and bw.  therefore, h2oexp for these sub-bands are identical.
!**********************************************************************

  DO k=1,np
    DO j=1,n
      DO i=1,m

!-----xh is   the scaled water vapor amount for line absorption
!  computed from (27).

        xh = dh2o(i,j,k)*(pa(i,j,k)*0.002)                              &
            * ( 1.+(aw(ib)+bw(ib)* dt(i,j,k))*dt(i,j,k) )

!-----h2oexp is the water vapor transmittance of the layer (k2-1)
!  due to line absorption

        h2oexp(i,j,k,1) = expmn(-xh*xkw(ib))

      END DO
    END DO
  END DO

  DO ik=2,6

    IF(mw(ib) == 6) THEN

      DO k=1,np
        DO j=1,n
          DO i=1,m
            xh = h2oexp(i,j,k,ik-1)*h2oexp(i,j,k,ik-1)
            h2oexp(i,j,k,ik) = xh*xh*xh
          END DO
        END DO
      END DO

    ELSE IF(mw(ib) == 8) THEN

      DO k=1,np
        DO j=1,n
          DO i=1,m
            xh = h2oexp(i,j,k,ik-1)*h2oexp(i,j,k,ik-1)
            xh = xh*xh
            h2oexp(i,j,k,ik) = xh*xh
          END DO
        END DO
      END DO

    ELSE

      DO k=1,np
        DO j=1,n
          DO i=1,m
            xh = h2oexp(i,j,k,ik-1)*h2oexp(i,j,k,ik-1)
            xh = xh*xh
            xh = xh*xh
            h2oexp(i,j,k,ik) = xh*xh
          END DO
        END DO
      END DO

    END IF
  END DO

  RETURN
END SUBROUTINE h2oexps
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CONEXPS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE conexps(ib,m,n,np,dcont,conexp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate exponentials for continuum absorption in individual
!  layers.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Adopted the original code and formatted it in accordance with the
!  ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!**********************************************************************
!   compute exponentials for continuum absorption in individual layers.
!
!---- input parameters
!  spectral band (ib)
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of layers (np)
!  layer scaled water vapor amount for continuum absorption (dcont)
!
!---- output parameters
!  1 or 3 exponentials for each layer (conexp)
!
!**********************************************************************
!
!fpp$ expand (expmn)
!!dir$ inline always expmn
!*$*  inline routine (expmn)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: ib,m,n,np,i,j,k,iq

!---- input parameters ------

  REAL :: dcont(m,n,np)

!---- updated parameters -----

  REAL :: conexp(m,n,np,3)

!---- static data -----

  INTEGER :: NE(8)
  REAL :: xke(8)


!-----xke are the absorption coefficients for the first
!  k-distribution function due to water vapor continuum absorption
!  (table 6).  units are cm**2/g

  DATA xke /  0.00,   0.00,   27.40,   15.8,                            &
              9.40,   7.75,     0.0,    0.0/


!-----ne is the number of terms in computing water vapor
!  continuum transmittance (Table 6).
!  band 3 is divided into 3 sub-bands.

  DATA NE /0,0,3,1,1,1,0,0/
!
!-----------------------------------------------------------------------
!
!  Functions:
!
!-----------------------------------------------------------------------
!

!-----expmn is an external function
  REAL :: expmn
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,np
    DO j=1,n
      DO i=1,m
        conexp(i,j,k,1) = expmn(-dcont(i,j,k)*xke(ib))
      END DO
    END DO
  END DO

  IF (ib == 3) THEN

!-----the absorption coefficients for sub-bands 3b (iq=2) and 3a (iq=3)
!  are, respectively, double and quadruple that for sub-band 3c (iq=1)
!  (table 6).

    DO iq=2,3
      DO k=1,np
        DO j=1,n
          DO i=1,m
            conexp(i,j,k,iq) = conexp(i,j,k,iq-1) *conexp(i,j,k,iq-1)
          END DO
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE conexps
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CO2EXPS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE co2exps(m,n,np,dco2,pa,dt,co2exp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate co2 exponentials for individual layers.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Adopted the original code and formatted it in accordance with the
!  ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!
!**********************************************************************
!   compute co2 exponentials for individual layers.
!
!---- input parameters
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of layers (np)
!  layer co2 amount (dco2)
!  layer pressure (pa)
!  layer temperature minus 250K (dt)
!
!---- output parameters
!  6 exponentials for each layer (co2exp)
!**********************************************************************
!
!fpp$ expand (expmn)
!!dir$ inline always expmn
!*$*  inline routine (expmn)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: m,n,np,i,j,k

!---- input parameters -----

  REAL :: dco2(m,n,np)
  REAL :: pa(m,n,np)
  REAL :: dt(m,n,np)

!---- output parameters -----

  REAL :: co2exp(m,n,np,6,2)

!---- static data -----

  REAL :: xkc(2),ac(2),bc(2),pm(2),prc(2)

!---- temporary arrays -----

  REAL :: xc

!-----xkc is the absorption coefficients for the
!  first k-distribution function due to co2 (table 7).
!  units are 1/(cm-atm)stp.

  DATA xkc/2.656E-5,2.656E-3/

!-----parameters (table 3) for computing the scaled co2 amount
!  using (27).

  DATA prc/  300.0,   30.0/
  DATA pm /    0.5,   0.85/
  DATA ac / 0.0182, 0.0042/
  DATA bc /1.07E-4,2.00E-5/
!
!-----------------------------------------------------------------------
!
!  Functions:
!
!-----------------------------------------------------------------------
!

!-----function expmn

  REAL :: expmn
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,np
    DO j=1,n
      DO i=1,m

!-----compute the scaled co2 amount from eq. (27) for band-wings
!  (sub-bands 3a and 3c).

        xc = dco2(i,j,k)*(pa(i,j,k)/prc(1))**pm(1)                      &
                *(1.+(ac(1)+bc(1)*dt(i,j,k))*dt(i,j,k))

!-----six exponential by powers of 8 (table 7).

        co2exp(i,j,k,1,1)=expmn(-xc*xkc(1))

        xc=co2exp(i,j,k,1,1)*co2exp(i,j,k,1,1)
        xc=xc*xc
        co2exp(i,j,k,2,1)=xc*xc

        xc=co2exp(i,j,k,2,1)*co2exp(i,j,k,2,1)
        xc=xc*xc
        co2exp(i,j,k,3,1)=xc*xc

        xc=co2exp(i,j,k,3,1)*co2exp(i,j,k,3,1)
        xc=xc*xc
        co2exp(i,j,k,4,1)=xc*xc

        xc=co2exp(i,j,k,4,1)*co2exp(i,j,k,4,1)
        xc=xc*xc
        co2exp(i,j,k,5,1)=xc*xc

        xc=co2exp(i,j,k,5,1)*co2exp(i,j,k,5,1)
        xc=xc*xc
        co2exp(i,j,k,6,1)=xc*xc

!-----compute the scaled co2 amount from eq. (27) for band-center
!  region (sub-band 3b).

        xc = dco2(i,j,k)*(pa(i,j,k)/prc(2))**pm(2)                      &
                *(1.+(ac(2)+bc(2)*dt(i,j,k))*dt(i,j,k))

        co2exp(i,j,k,1,2)=expmn(-xc*xkc(2))

        xc=co2exp(i,j,k,1,2)*co2exp(i,j,k,1,2)
        xc=xc*xc
        co2exp(i,j,k,2,2)=xc*xc

        xc=co2exp(i,j,k,2,2)*co2exp(i,j,k,2,2)
        xc=xc*xc
        co2exp(i,j,k,3,2)=xc*xc

        xc=co2exp(i,j,k,3,2)*co2exp(i,j,k,3,2)
        xc=xc*xc
        co2exp(i,j,k,4,2)=xc*xc

        xc=co2exp(i,j,k,4,2)*co2exp(i,j,k,4,2)
        xc=xc*xc
        co2exp(i,j,k,5,2)=xc*xc

        xc=co2exp(i,j,k,5,2)*co2exp(i,j,k,5,2)
        xc=xc*xc
        co2exp(i,j,k,6,2)=xc*xc

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE co2exps
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE COLUMN                     ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE column(m,n,np,pa,dt,sabs0,sabs,spre,stem)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate column-integrated (from top of the model atmosphere)
!  absorber amount, absorber-weighted pressure and temperature.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Adopted the original code and formatted it in accordance with the
!  ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!**************************************************************************
!-----compute column-integrated (from top of the model atmosphere)
!  absorber amount (sabs), absorber-weighted pressure (spre) and
!  temperature (stem).
!  computations of spre and stem follows eqs. (37) and (38).
!
!--- input parameters
!   number of soundings in zonal direction (m)
!   number of soundings in meridional direction (n)
!   number of atmospheric layers (np)
!   layer pressure (pa)
!   layer temperature minus 250K (dt)
!   layer absorber amount (sabs0)
!
!--- output parameters
!   column-integrated absorber amount (sabs)
!   column absorber-weighted pressure (spre)
!   column absorber-weighted temperature (stem)
!
!--- units of pa and dt are mb and k, respectively.
!    units of sabs are g/cm**2 for water vapor and (cm-atm)stp for co2 and o3
!**************************************************************************
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: m,n,np,i,j,k

!---- input parameters -----

  REAL :: pa(m,n,np)
  REAL :: dt(m,n,np)
  REAL :: sabs0(m,n,np)

!---- output parameters -----

  REAL :: sabs(m,n,np+1)
  REAL :: spre(m,n,np+1)
  REAL :: stem(m,n,np+1)

!*********************************************************************
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,n
    DO i=1,m
      sabs(i,j,1)=0.0
      spre(i,j,1)=0.0
      stem(i,j,1)=0.0
    END DO
  END DO

  DO k=1,np
    DO j=1,n
      DO i=1,m
        sabs(i,j,k+1)=sabs(i,j,k)+sabs0(i,j,k)
        spre(i,j,k+1)=spre(i,j,k)+pa(i,j,k)*sabs0(i,j,k)
        stem(i,j,k+1)=stem(i,j,k)+dt(i,j,k)*sabs0(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE column
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TABLUP                     ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE tablup(k1,k2,m,n,np,nx,nh,nt,sabs,spre,stem,w1,p1,           &
           dwe,dpe,coef1,coef2,coef3,tran)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate water vapor, co2, and o3 transmittances using table
!  look-up. rlwopt = 0, or high=.false.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Adopted the original code and formatted it in accordance with the
!  ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!**********************************************************************
!   compute water vapor, co2, and o3 transmittances between levels k1 and k2
!   using table look-up for m x n soundings.
!
!   Calculations follow Eq. (40) of Chou and Suarez (1995)
!
!---- input ---------------------
!  indices for pressure levels (k1 and k2)
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of atmospheric layers (np)
!  number of pressure intervals in the table (nx)
!  number of absorber amount intervals in the table (nh)
!  number of tables copied (nt)
!  column-integrated absorber amount (sabs)
!  column absorber amount-weighted pressure (spre)
!  column absorber amount-weighted temperature (stem)
!  first value of absorber amount (log10) in the table (w1)
!  first value of pressure (log10) in the table (p1)
!  size of the interval of absorber amount (log10) in the table (dwe)
!  size of the interval of pressure (log10) in the table (dpe)
!  pre-computed coefficients (coef1, coef2, and coef3)
!
!---- updated ---------------------
!  transmittance (tran)
!
!  Note:
!   (1) units of sabs are g/cm**2 for water vapor and (cm-atm)stp for co2 and o3.
!   (2) units of spre and stem are, respectively, mb and K.
!   (3) there are nt identical copies of the tables (coef1, coef2, and
!    coef3).  the prupose of using the multiple copies of tables is
!    to increase the speed in parallel (vectorized) computations.
!    if such advantage does not exist, nt can be set to 1.
!
!**********************************************************************
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: k1,k2,m,n,np,nx,nh,nt,i,j

!---- input parameters -----

  REAL :: w1,p1,dwe,dpe
  REAL :: sabs(m,n,np+1)
  REAL :: spre(m,n,np+1)
  REAL :: stem(m,n,np+1)

  REAL :: coef1(nx,nh,nt)
  REAL :: coef2(nx,nh,nt)
  REAL :: coef3(nx,nh,nt)

!---- update parameter -----

  REAL :: tran(m,n)

!---- temporary variables -----

  REAL :: x1,x2,x3,we,pe,fw,fp,pa,pb,pc,ax,ba,bb,t1,ca,cb,t2
  INTEGER :: iw,ip,nn

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,n
    DO i=1,m

      nn=MOD(i,nt)+1

      x1=sabs(i,j,k2)-sabs(i,j,k1)
      x2=(spre(i,j,k2)-spre(i,j,k1))/x1
      x3=(stem(i,j,k2)-stem(i,j,k1))/x1

      we=(LOG10(x1)-w1)/dwe
      pe=(LOG10(x2)-p1)/dpe

      we=MAX(we,w1-2.*dwe)
      pe=MAX(pe,p1)

      iw=INT(we+1.5)
      ip=INT(pe+1.5)

      iw=MIN(iw,nh-1)
      iw=MAX(iw, 2)

      ip=MIN(ip,nx-1)
      ip=MAX(ip, 1)

      fw=we-FLOAT(iw-1)
      fp=pe-FLOAT(ip-1)

!-----linear interpolation in pressure

      pa = coef1(ip,iw-1,nn)*(1.-fp)+coef1(ip+1,iw-1,nn)*fp
      pb = coef1(ip,iw,  nn)*(1.-fp)+coef1(ip+1,iw,  nn)*fp
      pc = coef1(ip,iw+1,nn)*(1.-fp)+coef1(ip+1,iw+1,nn)*fp

!-----quadratic interpolation in absorber amount for coef1

      ax = (-pa*(1.-fw)+pc*(1.+fw)) *fw*0.5 + pb*(1.-fw*fw)

!-----linear interpolation in absorber amount for coef2 and coef3

      ba = coef2(ip,iw,  nn)*(1.-fp)+coef2(ip+1,iw,  nn)*fp
      bb = coef2(ip,iw+1,nn)*(1.-fp)+coef2(ip+1,iw+1,nn)*fp
      t1 = ba*(1.-fw) + bb*fw

      ca = coef3(ip,iw,  nn)*(1.-fp)+coef3(ip+1,iw,  nn)*fp
      cb = coef3(ip,iw+1,nn)*(1.-fp)+coef3(ip+1,iw+1,nn)*fp
      t2 = ca*(1.-fw) + cb*fw

!-----update the total transmittance between levels k1 and k2

      tran(i,j)= (ax + (t1+t2*x3) * x3)*tran(i,j)

    END DO
  END DO

  RETURN
END SUBROUTINE tablup
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WVKDIS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wvkdis(ib,m,n,np,k,h2oexp,conexp,th2o,tcon,tran)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate water vapor transmittance using the k-distribution
!  method.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Adopted the original code and formatted it in accordance with the
!  ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!**********************************************************************
!   compute water vapor transmittance between levels k1 and k2 for
!   m x n soundings using the k-distribution method.
!
!   computations follow eqs. (34), (46), (50) and (52).
!
!---- input parameters
!  spectral band (ib)
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of levels (np)
!  current level (k)
!  exponentials for line absorption (h2oexp)
!  exponentials for continuum absorption (conexp)
!
!---- updated parameters
!  transmittance between levels k1 and k2 due to
!    water vapor line absorption (th2o)
!  transmittance between levels k1 and k2 due to
!    water vapor continuum absorption (tcon)
!  total transmittance (tran)
!
!**********************************************************************
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: ib,m,n,np,k,i,j

!---- input parameters ------

  REAL :: conexp(m,n,np,3)
  REAL :: h2oexp(m,n,np,6)

!---- updated parameters -----

  REAL :: th2o(m,n,6)
  REAL :: tcon(m,n,3)
  REAL :: tran(m,n)

!---- static data -----

  INTEGER :: NE(8)
  REAL :: fkw(6,8),gkw(6,3)

!---- temporary variable -----

  REAL :: trnth2o

!-----fkw is the planck-weighted k-distribution function due to h2o
!  line absorption given in table 4 of Chou and Suarez (1995).
!  the k-distribution function for the third band, fkw(*,3), is not used

  DATA fkw / 0.2747,0.2717,0.2752,0.1177,0.0352,0.0255,                 &
             0.1521,0.3974,0.1778,0.1826,0.0374,0.0527,                 &
             6*1.00,                                                    &
             0.4654,0.2991,0.1343,0.0646,0.0226,0.0140,                 &
             0.5543,0.2723,0.1131,0.0443,0.0160,0.0000,                 &
             0.1846,0.2732,0.2353,0.1613,0.1146,0.0310,                 &
             0.0740,0.1636,0.4174,0.1783,0.1101,0.0566,                 &
             0.1437,0.2197,0.3185,0.2351,0.0647,0.0183/

!-----gkw is the planck-weighted k-distribution function due to h2o
!  line absorption in the 3 subbands (800-720,620-720,540-620 /cm)
!  of band 3 given in table 7.  Note that the order of the sub-bands
!  is reversed.

  DATA gkw/  0.1782,0.0593,0.0215,0.0068,0.0022,0.0000,                 &
             0.0923,0.1675,0.0923,0.0187,0.0178,0.0000,                 &
             0.0000,0.1083,0.1581,0.0455,0.0274,0.0041/

!-----ne is the number of terms used in each band to compute water vapor
!  continuum transmittance (table 6).

  DATA NE /0,0,3,1,1,1,0,0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!-----tco2 are the six exp factors between levels k1 and k2
!  tran is the updated total transmittance between levels k1 and k2


!-----th2o is the 6 exp factors between levels k1 and k2 due to
!  h2o line absorption.

!-----tcon is the 3 exp factors between levels k1 and k2 due to
!  h2o continuum absorption.

!-----trnth2o is the total transmittance between levels k1 and k2 due
!  to both line and continuum absorption computed from eq. (52).

  DO j=1,n
    DO i=1,m
      th2o(i,j,1) = th2o(i,j,1)*h2oexp(i,j,k,1)
      th2o(i,j,2) = th2o(i,j,2)*h2oexp(i,j,k,2)
      th2o(i,j,3) = th2o(i,j,3)*h2oexp(i,j,k,3)
      th2o(i,j,4) = th2o(i,j,4)*h2oexp(i,j,k,4)
      th2o(i,j,5) = th2o(i,j,5)*h2oexp(i,j,k,5)
      th2o(i,j,6) = th2o(i,j,6)*h2oexp(i,j,k,6)
    END DO
  END DO


  IF (NE(ib) == 0) THEN


    DO j=1,n
      DO i=1,m

        trnth2o      =(fkw(1,ib)*th2o(i,j,1)                            &
                     + fkw(2,ib)*th2o(i,j,2)                            &
                     + fkw(3,ib)*th2o(i,j,3)                            &
                     + fkw(4,ib)*th2o(i,j,4)                            &
                     + fkw(5,ib)*th2o(i,j,5)                            &
                     + fkw(6,ib)*th2o(i,j,6))

        tran(i,j)=tran(i,j)*trnth2o

      END DO
    END DO

  ELSE IF (NE(ib) == 1) THEN


    DO j=1,n
      DO i=1,m

        tcon(i,j,1)= tcon(i,j,1)*conexp(i,j,k,1)

        trnth2o      =(fkw(1,ib)*th2o(i,j,1)                            &
                     + fkw(2,ib)*th2o(i,j,2)                            &
                     + fkw(3,ib)*th2o(i,j,3)                            &
                     + fkw(4,ib)*th2o(i,j,4)                            &
                     + fkw(5,ib)*th2o(i,j,5)                            &
                     + fkw(6,ib)*th2o(i,j,6))*tcon(i,j,1)

        tran(i,j)=tran(i,j)*trnth2o

      END DO
    END DO

  ELSE

    DO j=1,n
      DO i=1,m

        tcon(i,j,1)= tcon(i,j,1)*conexp(i,j,k,1)
        tcon(i,j,2)= tcon(i,j,2)*conexp(i,j,k,2)
        tcon(i,j,3)= tcon(i,j,3)*conexp(i,j,k,3)

        trnth2o      = (  gkw(1,1)*th2o(i,j,1)                          &
                        + gkw(2,1)*th2o(i,j,2)                          &
                        + gkw(3,1)*th2o(i,j,3)                          &
                        + gkw(4,1)*th2o(i,j,4)                          &
                        + gkw(5,1)*th2o(i,j,5)                          &
                        + gkw(6,1)*th2o(i,j,6) ) * tcon(i,j,1)          &
                     + (  gkw(1,2)*th2o(i,j,1)                          &
                        + gkw(2,2)*th2o(i,j,2)                          &
                        + gkw(3,2)*th2o(i,j,3)                          &
                        + gkw(4,2)*th2o(i,j,4)                          &
                        + gkw(5,2)*th2o(i,j,5)                          &
                        + gkw(6,2)*th2o(i,j,6) ) * tcon(i,j,2)          &
                     + (  gkw(1,3)*th2o(i,j,1)                          &
                        + gkw(2,3)*th2o(i,j,2)                          &
                        + gkw(3,3)*th2o(i,j,3)                          &
                        + gkw(4,3)*th2o(i,j,4)                          &
                        + gkw(5,3)*th2o(i,j,5)                          &
                        + gkw(6,3)*th2o(i,j,6) ) * tcon(i,j,3)

        tran(i,j)=tran(i,j)*trnth2o

      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE wvkdis
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CO2KDIS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE co2kdis(m,n,np,k,co2exp,tco2,tran)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate co2 transmittances using the k-distribution method
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Adopted the original code and formatted it in accordance with the
!  ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!**********************************************************************
!   compute co2 transmittances between levels k1 and k2 for m x n soundings
!   using the k-distribution method with linear pressure scaling.
!
!   computations follow eq. (34).
!
!---- input parameters
!   number of grid intervals in zonal direction (m)
!   number of grid intervals in meridional direction (n)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to co2 absorption
!  for the various values of the absorption coefficient (tco2)
!   total transmittance (tran)
!
!**********************************************************************
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: m,n,np,k,i,j

!---- input parameters -----

  REAL :: co2exp(m,n,np,6,2)

!---- updated parameters -----

  REAL :: tco2(m,n,6,2)
  REAL :: tran(m,n)

!---- static data -----

  REAL :: gkc(6,2)

!---- temporary variable -----

  REAL :: xc

!-----gkc is the planck-weighted co2 k-distribution function
!  in the band-wing and band-center regions given in table 7.
!  for computing efficiency, sub-bands 3a and 3c are combined.

  DATA gkc/  0.1395,0.1407,0.1549,0.1357,0.0182,0.0220,                 &
             0.0766,0.1372,0.1189,0.0335,0.0169,0.0059/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!-----tco2 is the 6 exp factors between levels k1 and k2.
!  xc is the total co2 transmittance given by eq. (53).

  DO j=1,n
    DO i=1,m

!-----band-wings

      tco2(i,j,1,1)=tco2(i,j,1,1)*co2exp(i,j,k,1,1)
      xc=             gkc(1,1)*tco2(i,j,1,1)

      tco2(i,j,2,1)=tco2(i,j,2,1)*co2exp(i,j,k,2,1)
      xc=xc+gkc(2,1)*tco2(i,j,2,1)

      tco2(i,j,3,1)=tco2(i,j,3,1)*co2exp(i,j,k,3,1)
      xc=xc+gkc(3,1)*tco2(i,j,3,1)

      tco2(i,j,4,1)=tco2(i,j,4,1)*co2exp(i,j,k,4,1)
      xc=xc+gkc(4,1)*tco2(i,j,4,1)

      tco2(i,j,5,1)=tco2(i,j,5,1)*co2exp(i,j,k,5,1)
      xc=xc+gkc(5,1)*tco2(i,j,5,1)

      tco2(i,j,6,1)=tco2(i,j,6,1)*co2exp(i,j,k,6,1)
      xc=xc+gkc(6,1)*tco2(i,j,6,1)

!-----band-center region

      tco2(i,j,1,2)=tco2(i,j,1,2)*co2exp(i,j,k,1,2)
      xc=xc+gkc(1,2)*tco2(i,j,1,2)

      tco2(i,j,2,2)=tco2(i,j,2,2)*co2exp(i,j,k,2,2)
      xc=xc+gkc(2,2)*tco2(i,j,2,2)

      tco2(i,j,3,2)=tco2(i,j,3,2)*co2exp(i,j,k,3,2)
      xc=xc+gkc(3,2)*tco2(i,j,3,2)

      tco2(i,j,4,2)=tco2(i,j,4,2)*co2exp(i,j,k,4,2)
      xc=xc+gkc(4,2)*tco2(i,j,4,2)

      tco2(i,j,5,2)=tco2(i,j,5,2)*co2exp(i,j,k,5,2)
      xc=xc+gkc(5,2)*tco2(i,j,5,2)

      tco2(i,j,6,2)=tco2(i,j,6,2)*co2exp(i,j,k,6,2)
      xc=xc+gkc(6,2)*tco2(i,j,6,2)

      tran(i,j)=tran(i,j)*xc

    END DO
  END DO

  RETURN
END SUBROUTINE co2kdis
