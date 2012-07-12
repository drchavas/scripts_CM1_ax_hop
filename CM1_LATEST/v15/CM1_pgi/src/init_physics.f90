

      subroutine init_physics(prs0,rf0,cdu,cdv,ce,dum1,dum2,dum3,u0,ua,v0,va,o30,   &
                             lu_index,xland,emiss,thc,albd,znt,mavail,f2d,tsk,u1,v1,w1)
      use module_sf_slab
      use module_sf_sfclay
      implicit none

      include 'input.incl'
      include 'radcst.incl'

      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: prs0,rf0
      real, intent(inout), dimension(ib:ie,jb:je) :: cdu,cdv,ce
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0,ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0,va
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: o30
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: lu_index
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: xland,emiss,thc,albd,znt,mavail,f2d
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,u1,v1,w1

      real :: foo1,foo2,foo3,foo4,foo5,foo6

!-----------------------------------------------------------------------
!-----  USERS SHOULD NOT NEED TO CHANGE ANYTHING IN THIS FILE ----------
!-----  (unless you really, really know what you are doing -------------
!-----------------------------------------------------------------------

      if(radopt.eq.1)then
        ! initialize radiation code:
        call setradwrk(nir,njr,nkr)
        call julday( year, month, day, jday )
        write(outfile,*) '  jday = ',jday
        call initrad(myid,year,month,day,hour,minute,second,jday,nir,njr,nkr)
        o30 = 1.0e-6
        call fito3(nir,njr,1,1,nkr,dum1(1,1,1),dum2(1,1,1),prs0,o30,ib,ie,jb,je,kb,ke,nk)
        ! Settings from Goddard scheme:
        call getgoddardvars(foo1,foo2,foo3,foo4,foo5,foo6)
        roqr = foo1
        tnw  = foo2
        roqs = foo3
        tns  = foo4
        roqg = foo5
        tng  = foo6
      endif

      IF( idrag.eq.1 .or. isfcflx.eq.1 )THEN
        if( sfcmodel.eq.1 )then
          call getcecd(cdu,cdv,ce,u0,v0,rf0,u1,v1,w1,ua,va)
        elseif( sfcmodel.eq.2 )then
          call sfclayinit
        endif
      ENDIF

      f2d = fcor

!-----------------------------------------------------------------------

      end subroutine init_physics


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getgoddardvars(foo1,foo2,foo3,foo4,foo5,foo6)
      implicit none
      include 'goddard.incl'

      real, intent(inout) :: foo1,foo2,foo3,foo4,foo5,foo6

      foo1 = roqr
      foo2 = tnw
      foo3 = roqs
      foo4 = tnss
      foo5 = roqg
      foo6 = tng

      end subroutine getgoddardvars


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine initrad(myid,year,month,day,hour,minute,second,jday,   &
                         nir,njr,nkr)
      implicit none

      include 'irrad.incl'
      include 'radzen.incl'
      include 'radmore.incl'

      integer, intent(in) :: myid,year,month,day,hour,minute,second,jday,   &
                             nir,njr,nkr

      integer :: ip,iw,it
      logical :: high

!----------------------------------------------------------------------

  IF ( rlwopt == 0 ) THEN
    high = .false.
  ELSE
    high = .true.
  END IF

  if(myid.eq.0) print *,'  high = ',high

!----------------------------------------------------------------------
!  from zenangl:

    pi2 = 2.0 * pi
    deg2rad = pi/180.0
    rad2deg = 1./deg2rad

    hour0 = FLOAT(hour)                                                 &
          + FLOAT(minute)/60.0                                          &
          + FLOAT(second)/3600.0

    IF ( MOD(year, 4) == 0 ) THEN
      yrday = 366.
    ELSE
      yrday = 365.
    END IF

!!! not using arps MPI code:  GHB, 100720
!!! hard-wire these in, just in case:
    nxmid = 1
    nymid = 1
    source = 0

!----------------------------------------------------------------------
!  from irrad:

!-----tables co2 and h2o are only used with 'high' option

    IF (high) THEN

      DO iw=1,nh
        DO ip=1,nx
          h11(ip,iw,1)=1.0-h11(ip,iw,1)
          h21(ip,iw,1)=1.0-h21(ip,iw,1)
          h71(ip,iw,1)=1.0-h71(ip,iw,1)
        END DO
      END DO

      DO iw=1,nc
        DO ip=1,nx
          c1(ip,iw,1)=1.0-c1(ip,iw,1)
        END DO
      END DO

!-----tables are replicated to avoid memory bank conflicts

      DO it=2,nt
        DO iw=1,nc
          DO ip=1,nx
            c1 (ip,iw,it)= c1(ip,iw,1)
            c2 (ip,iw,it)= c2(ip,iw,1)
            c3 (ip,iw,it)= c3(ip,iw,1)
          END DO
        END DO
        DO iw=1,nh
          DO ip=1,nx
            h11(ip,iw,it)=h11(ip,iw,1)
            h12(ip,iw,it)=h12(ip,iw,1)
            h13(ip,iw,it)=h13(ip,iw,1)
            h21(ip,iw,it)=h21(ip,iw,1)
            h22(ip,iw,it)=h22(ip,iw,1)
            h23(ip,iw,it)=h23(ip,iw,1)
            h71(ip,iw,it)=h71(ip,iw,1)
            h72(ip,iw,it)=h72(ip,iw,1)
            h73(ip,iw,it)=h73(ip,iw,1)
          END DO
        END DO
      END DO

    END IF

!-----always use table look-up for ozone transmittance

    DO iw=1,no
      DO ip=1,nx
        o1(ip,iw,1)=1.0-o1(ip,iw,1)
      END DO
    END DO

    DO it=2,nt
      DO iw=1,no
        DO ip=1,nx
          o1 (ip,iw,it)= o1(ip,iw,1)
          o2 (ip,iw,it)= o2(ip,iw,1)
          o3 (ip,iw,it)= o3(ip,iw,1)
        END DO
      END DO
    END DO

      return
      end subroutine initrad


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getlanduse(season,myid,ib,ie,jb,je,ibl,iel,jbl,jel,   &
                            lu_index,xland,emiss,thc,albedo,znt,mavail)
      implicit none





      integer, intent(in) :: season,myid,ib,ie,jb,je,ibl,iel,jbl,jel
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: lu_index
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: xland,emiss,thc,albedo,znt,mavail

      integer, parameter :: max_cats = 33    ! max categories
      integer, parameter :: max_seas =  2    ! max seasons

      integer, parameter :: iswater = 16  ! water must be category 16

      real, dimension(:,:), allocatable :: albd,slmo,sfem,sfz0,therin,scfx,sfhc
      integer :: cat,sea

      integer :: i,j,is,isn,ierr

      allocate(   albd(max_cats,max_seas) )
      allocate(   slmo(max_cats,max_seas) )
      allocate(   sfem(max_cats,max_seas) )
      allocate(   sfz0(max_cats,max_seas) )
      allocate( therin(max_cats,max_seas) )
      allocate(   scfx(max_cats,max_seas) )
      allocate(   sfhc(max_cats,max_seas) )

!-----------------------------------------------------------------------

      IF( myid.eq.0 )THEN
        open(unit=11,file='LANDUSE.TBL',status='old',err=888)
        read(11,*)
        read(11,*)
        do isn=1,max_seas
          read(11,*)
          do is=1,max_cats
            read(11,*) i,albd(is,isn),slmo(is,isn),sfem(is,isn),sfz0(is,isn), &
                       therin(is,isn),scfx(is,isn),sfhc(is,isn)
          enddo
        enddo
!        print *
!        print *,'  summer:'
!        print *,'  albd = ',albd(:,1)
!        print *,'  slmo = ',slmo(:,1)
!        print *,'  sfem = ',sfem(:,1)
!        print *,'  sfz0 = ',sfz0(:,1)
!        print *,'  ther = ',therin(:,1)
!        print *,'  scfx = ',scfx(:,1)
!        print *,'  sfhc = ',sfhc(:,1)
!        print *
!        print *,'  winter:'
!        print *,'  albd = ',albd(:,2)
!        print *,'  slmo = ',slmo(:,2)
!        print *,'  sfem = ',sfem(:,2)
!        print *,'  sfz0 = ',sfz0(:,2)
!        print *,'  ther = ',therin(:,2)
!        print *,'  scfx = ',scfx(:,2)
!        print *,'  sfhc = ',sfhc(:,2)
      ENDIF


!-----------------------------------------------------------------------

      ! ISN = season:    summer = 1    winter = 2
      ISN = season

      do j=jb,je
      do i=ib,ie
        IS = lu_index(i,j)
        ! SET NO-DATA POINTS (IS=0) TO WATER
        IF(IS.EQ.0)THEN
          IS=ISWATER
        ENDIF
        if( albd(is,isn).le.0.0 )then
          print *,'  category not found '
          call stopcm1
        endif
        ALBEDO(I,J) = ALBD(IS,ISN)/100.
        THC(i,j) = THERIN(IS,ISN)/100.
        EMISS(I,J) = SFEM(IS,ISN)
        ZNT(I,J) = SFZ0(IS,ISN)/100.
        MAVAIL(I,J) = SLMO(IS,ISN)
        IF(IS.NE.ISWATER)THEN
          XLAND(I,J)=1.0
        ELSE
          XLAND(I,J)=2.0
        ENDIF
      enddo
      enddo

!-----------------------------------------------------------------------

      deallocate(   albd )
      deallocate(   slmo )
      deallocate(   sfem )
      deallocate(   sfz0 )
      deallocate( therin )
      deallocate(   scfx )
      deallocate(   sfhc )


      return

!-----------------------------------------------------------------------

888   print *
      print *,'  There was an error opening the LANDUSE.TBL file '
      print *
      print *,'  Please make sure that LANDUSE.TBL is in the same directory '
      print *,'  as cm1.exe ... it is distributed with CM1 in the "run" directory '
      print *
      call stopcm1
      end subroutine getlanduse
