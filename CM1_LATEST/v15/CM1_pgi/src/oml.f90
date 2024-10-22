

!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_sf_oml

CONTAINS

!----------------------------------------------------------------
   SUBROUTINE OCEANML(tml,t0ml,hml,h0ml,huml,hvml,ust,u_phy,v_phy, &
                      tmoml,f,g,oml_gamma,                         &
                      XLAND,HFX,LH,TSK,GSW,GLW,EMISS,              &
                      DELTSM,STBOLT,                               &
                      ids,ide, jds,jde, kds,kde,                   &
                      ims,ime, jms,jme, kms,kme,                   &
                      its,ite, jts,jte, kts,kte                    )

!----------------------------------------------------------------
   IMPLICIT NONE
!----------------------------------------------------------------
!
!  SUBROUTINE OCEANML CALCULATES THE SEA SURFACE TEMPERATURE (TSK)
!  FROM A SIMPLE OCEAN MIXED LAYER MODEL BASED ON
!  (Pollard, Rhines and Thompson (1973).
!
!-- TML         ocean mixed layer temperature (K)
!-- T0ML        ocean mixed layer temperature (K) at initial time
!-- TMOML       top 200 m ocean mean temperature (K) at initial time
!-- HML         ocean mixed layer depth (m)
!-- H0ML        ocean mixed layer depth (m) at initial time
!-- HUML        ocean mixed layer u component of wind
!-- HVML        ocean mixed layer v component of wind
!-- OML_GAMMA   deep water lapse rate (K m-1)
!-- UAIR,VAIR   lowest model level wind component
!-- UST         frictional velocity
!-- HFX         upward heat flux at the surface (W/m^2)
!-- LH          latent heat flux at the surface (W/m^2)
!-- TSK         surface temperature (K)
!-- GSW         downward short wave flux at ground surface (W/m^2)
!-- GLW         downward long wave flux at ground surface (W/m^2)
!-- EMISS       emissivity of the surface
!-- XLAND       land mask (1 for land, 2 for water)
!-- STBOLT      Stefan-Boltzmann constant (W/m^2/K^4)
!-- F           Coriolis parameter
!-- DT          time step (second)
!-- G           acceleration due to gravity

   INTEGER,  INTENT(IN   )   ::     ids,ide, jds,jde, kds,kde,  &
                                    ims,ime, jms,jme, kms,kme,  &
                                    its,ite, jts,jte, kts,kte

   REAL,     INTENT(IN   )   ::     DELTSM, STBOLT

   REAL,    DIMENSION( ims:ime, jms:jme )                     , &
            INTENT(IN   )    ::                          EMISS, &
                                                         XLAND, &
                                                           GSW, &
                                                           GLW, &
                                                           HFX, &
                                                            LH

   REAL,    DIMENSION( ims:ime, jms:jme )                     , &
            INTENT(INOUT)    ::                            TSK

   REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) ::     &
                                    TML,T0ML,HML,H0ML,HUML,HVML

   REAL,    DIMENSION( ims:ime, jms:jme, kms:kme ), INTENT(IN   ) ::     &
                                             U_PHY,V_PHY

   REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) ::     &
                                             UST, F, TMOML

   REAL,    INTENT(IN   )   ::     G
   REAL,    INTENT(IN   )   ::     OML_GAMMA

! LOCAL VARS

   INTEGER ::  I,J

   DO J=jts,jte

!!!      CALL wrf_debug( 100, 'Call OML' )
         DO i=its,ite
            IF (XLAND(I,J).GT.1.5) THEN
               CALL OML1D(I,J,TML(i,j),T0ML(i,j),HML(i,j),H0ML(i,j),           &
                          HUML(i,j),HVML(i,j),TSK(i,j),HFX(i,j),               &
                          LH(i,j),GSW(i,j),GLW(i,j),TMOML(i,j),                &
                          U_PHY(i,j,kts),V_PHY(i,j,kts),UST(i,j),F(i,j),       &
                          EMISS(i,j),STBOLT,G,DELTSM,OML_GAMMA,                &
                          ids,ide, jds,jde, kds,kde,                           &
                          ims,ime, jms,jme, kms,kme,                           &
                          its,ite, jts,jte, kts,kte                            )
            ENDIF
         ENDDO

   ENDDO

   END SUBROUTINE OCEANML

!----------------------------------------------------------------
   SUBROUTINE OML1D(I,J,TML,T0ML,H,H0,HUML,                              &
                    HVML,TSK,HFX,                                        &
                    LH,GSW,GLW,TMOML,                                    &
                    UAIR,VAIR,UST,F,EMISS,STBOLT,G,DT,OML_GAMMA,         &
                    ids,ide, jds,jde, kds,kde,                           &
                    ims,ime, jms,jme, kms,kme,                           &
                    its,ite, jts,jte, kts,kte                            )

!----------------------------------------------------------------
   IMPLICIT NONE
!----------------------------------------------------------------
!
!  SUBROUTINE OCEANML CALCULATES THE SEA SURFACE TEMPERATURE (TSK) 
!  FROM A SIMPLE OCEAN MIXED LAYER MODEL BASED ON 
!  (Pollard, Rhines and Thompson (1973).
!
!-- TML         ocean mixed layer temperature (K)
!-- T0ML        ocean mixed layer temperature (K) at initial time
!-- TMOML       top 200 m ocean mean temperature (K) at initial time
!-- H           ocean mixed layer depth (m)
!-- H0          ocean mixed layer depth (m) at initial time
!-- HUML        ocean mixed layer u component of wind
!-- HVML        ocean mixed layer v component of wind
!-- OML_GAMMA   deep water lapse rate (K m-1)
!-- OMLCALL     whether to call oml model
!-- UAIR,VAIR   lowest model level wind component
!-- UST         frictional velocity
!-- HFX         upward heat flux at the surface (W/m^2)
!-- LH          latent heat flux at the surface (W/m^2)
!-- TSK         surface temperature (K)
!-- GSW         downward short wave flux at ground surface (W/m^2)
!-- GLW         downward long wave flux at ground surface (W/m^2)
!-- EMISS       emissivity of the surface
!-- STBOLT      Stefan-Boltzmann constant (W/m^2/K^4)
!-- F           Coriolis parameter
!-- DT          time step (second)
!-- G           acceleration due to gravity
!
!----------------------------------------------------------------
   INTEGER, INTENT(IN   )    ::      I, J
   INTEGER, INTENT(IN   )    ::      ids,ide, jds,jde, kds,kde, &
                                     ims,ime, jms,jme, kms,kme, &
                                     its,ite, jts,jte, kts,kte

   REAL,    INTENT(INOUT)    :: TML, H, H0, HUML, HVML, TSK

   REAL,    INTENT(IN   )    :: T0ML, HFX, LH, GSW, GLW,        &
                                UAIR, VAIR, UST, F, EMISS, TMOML

   REAL,    INTENT(IN) :: STBOLT, G, DT, OML_GAMMA

! Local
   REAL :: rhoair, rhowater, Gam, alp, BV2, A1, A2, B2, u, v, wspd, &
           hu1, hv1, hu2, hv2, taux, tauy, tauxair, tauyair, q, hold, &
           hsqrd, thp, cwater, ust2
   CHARACTER(LEN=120) :: time_series

      hu1=huml
      hv1=hvml
      rhoair=1.
      rhowater=1000.
      cwater=4200.
! Deep ocean lapse rate (K/m) - from Rich
      Gam=oml_gamma
!     if(i.eq.1 .and. j.eq.1 .or. i.eq.105.and.j.eq.105) print *, 'gamma = ', gam
!     Gam=0.14
!     Gam=5.6/40.
!     Gam=5./100.
! Thermal expansion coeff (/K)
!     alp=.0002
!     temp dependence (/K)
      alp=max((tml-273.15)*1.e-5, 1.e-6)
      BV2=alp*g*Gam
      thp=t0ml-Gam*(h-h0)
      A1=(tml-thp)*h - 0.5*Gam*h*h
      if(h.ne.0.)then
        u=hu1/h
        v=hv1/h
      else
        u=0.
        v=0.
      endif

!  time step

        q=(-hfx-lh+gsw+glw-stbolt*emiss*tml*tml*tml*tml)/(rhowater*cwater)
!       wspd=max(sqrt(uair*uair+vair*vair),0.1)
        wspd=sqrt(uair*uair+vair*vair)
        if (wspd .lt. 1.e-10 ) then
!          print *, 'i,j,wspd are ', i,j,wspd
           wspd = 1.e-10
        endif
! limit ust to 1.6 to give a value of ust for water of 0.05
!       ust2=min(ust, 1.6)
! new limit for ust: reduce atmospheric ust by half for ocean
        ust2=0.5*ust
        tauxair=ust2*ust2*uair/wspd
        taux=rhoair/rhowater*tauxair
        tauyair=ust2*ust2*vair/wspd
        tauy=rhoair/rhowater*tauyair
! note: forward-backward coriolis force for effective time-centering
        hu2=hu1+dt*( f*hv1 + taux)
        hv2=hv1+dt*(-f*hu2 + tauy)
! consider the flux effect
        A2=A1+q*dt

        huml=hu2
        hvml=hv2

        hold=h
        B2=hu2*hu2+hv2*hv2
        hsqrd=-A2/Gam + sqrt(A2*A2/(Gam*Gam) + 2.*B2/BV2)
        h=sqrt(max(hsqrd,0.0))
! limit to positive h change
        if(h.lt.hold)h=hold

        if(h.ne.0.)then
          tml=max(t0ml - Gam*(h-h0) + 0.5*Gam*h + A2/h, tmoml)
          u=hu2/h
          v=hv2/h
        else
          tml=t0ml
          u=0.
          v=0.
        endif
        tsk=tml
!        if(h.gt.100.)print *,i,j,h,tml,' h,tml'

! ww: output point data
!     if( (i.eq.190 .and. j.eq.115) .or. (i.eq.170 .and. j.eq.125) ) then
!        write(jtime,fmt='("TS ",f10.0)') float(itimestep)
!        CALL wrf_message ( TRIM(jtime) )
!        write(time_series,fmt='("OML",2I4,2F9.5,2F8.2,2E15.5,F8.3)') &
!              i,j,u,v,tml,h,taux,tauy,a2
!        CALL wrf_message ( TRIM(time_series) )
!     end if

   END SUBROUTINE OML1D

!================================================================
   SUBROUTINE omlinit(outfile,oml_hml0, tsk,                   &
                      tml,t0ml,hml,h0ml,huml,hvml,tmoml,       &
                      allowed_to_read, start_of_simulation,    &
                      ids,ide, jds,jde, kds,kde,               &
                      ims,ime, jms,jme, kms,kme,               &
                      its,ite, jts,jte, kts,kte                )
!----------------------------------------------------------------
   IMPLICIT NONE
!----------------------------------------------------------------
   INTEGER, INTENT(IN   )    ::      outfile
   LOGICAL , INTENT(IN)      ::      allowed_to_read
   LOGICAL , INTENT(IN)      ::      start_of_simulation
   INTEGER, INTENT(IN   )    ::      ids,ide, jds,jde, kds,kde, &
                                     ims,ime, jms,jme, kms,kme, &
                                     its,ite, jts,jte, kts,kte

   REAL,    DIMENSION( ims:ime, jms:jme )                     , &
            INTENT(IN)    ::                               TSK

   REAL,    DIMENSION( ims:ime, jms:jme )                     , &
            INTENT(INOUT)    ::     TML, T0ML, HML, H0ML, HUML, HVML, TMOML
   REAL   , INTENT(IN   )    ::     oml_hml0

!  LOCAR VAR

   INTEGER                   ::      L,J,I,itf,jtf
   CHARACTER*1024 message

!----------------------------------------------------------------
 
   itf=min0(ite,ide-1)
   jtf=min0(jte,jde-1)

   IF(start_of_simulation) THEN
     DO J=jts,jtf
     DO I=its,itf
       TML(I,J)=TSK(I,J)
       T0ML(I,J)=TSK(I,J)
     ENDDO
     ENDDO
     IF (oml_hml0 .gt. 0.) THEN
!!!        WRITE(message,*)'Initializing OML with HML0 = ', oml_hml0
!!!        CALL wrf_debug (0, TRIM(message))
        write(outfile,*) 'Initializing OML with HML0 = ', oml_hml0
        DO J=jts,jtf
        DO I=its,itf
          HML(I,J)=oml_hml0
          H0ML(I,J)=HML(I,J)
          HUML(I,J)=0.
          HVML(I,J)=0.
          TMOML(I,J)=TSK(I,J)-5.
        ENDDO
        ENDDO
     ELSE
!!!        WRITE(message,*)'Initializing OML with real HML0, h(1,1) = ', h0ml(1,1)
!!!        CALL wrf_debug (0, TRIM(message))
        write(outfile,*) 'Initializing OML with real HML0, h(1,1) = ', h0ml(1,1)
        stop 34543
        DO J=jts,jtf
        DO I=its,itf
          HML(I,J)=H0ML(I,J)
! fill in near coast area with SST: 200 K was set as missing value in ocean pre-processing code
          IF(TMOML(I,J).GT.200. .and. TMOML(I,J).LE.201.) TMOML(I,J)=TSK(I,J)
        ENDDO
        ENDDO
     ENDIF
   ENDIF

   END SUBROUTINE omlinit
!-------------------------------------------------------------------          
END MODULE module_sf_oml
