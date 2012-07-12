

      subroutine setup_output(tdef,qname,budname,xh,xf,yh,yf,xfref,yfref,zh,zf)
      implicit none

      include 'input.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      character*15 :: tdef
      character*3, dimension(maxq) :: qname
      character*6, dimension(maxq) :: budname
      real, dimension(ib:ie) :: xh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je) :: yh
      real, dimension(jb:je+1) :: yf
      real, dimension(-2:nx+4) :: xfref
      real, dimension(-2:ny+4) :: yfref
      real, dimension(ib:ie,jb:je,kb:ke) :: zh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: zf

!-----------------------------------------------------------------------

      integer :: i,j,k,n,flag
      character*8 text1
      character*30 text2
      character*50 fname
      character*8,  dimension(:), allocatable :: varname
      character*30, dimension(:), allocatable :: vardesc

!-----------------------------------------------------------------------
! get length of output_path string

    flag=0
    n=0
    do while( flag.eq.0 .and. n.le.70 )
      n=n+1
      if( output_path(n:n).eq.' ' .or. output_path(n:n).eq.'.' ) flag=1
    enddo

    strlen=n-1

!--------------------------------------
! get length of output_basename string

    flag=0
    n=0
    do while( flag.eq.0 .and. n.le.70 )
      n=n+1
      if( output_basename(n:n).eq.' ' .or. output_basename(n:n).eq.'.' ) flag=1
    enddo

    baselen=n-1

!------

    totlen = strlen + baselen

      string = '                                                                      '
    statfile = '                                                                      '
     sstring = '                                                                      '

  if(strlen.gt.0)then
      string(1:strlen) = output_path(1:strlen)
    statfile(1:strlen) = output_path(1:strlen)
  endif

      string(strlen+1:strlen+baselen) = output_basename(1:baselen)
    statfile(strlen+1:strlen+baselen) = output_basename(1:baselen)
     sstring(1:baselen) = output_basename(1:baselen)

    statfile(totlen+1:totlen+1+12) = '_stats.dat  '

    write(outfile,*)
    write(outfile,*) '  writing ctl files ... '
    write(outfile,*)
    write(outfile,*) '  strlen          = ',strlen
    write(outfile,*) '  baselen         = ',baselen
    write(outfile,*) '  totlen          = ',totlen
  if(strlen.gt.0)then
    write(outfile,*) '  output_path     = ',output_path(1:strlen)
  endif
    write(outfile,*) '  output_basename = ',output_basename(1:baselen)
    write(outfile,*) '  statfile        = ',statfile
    write(outfile,*)

      IF( myid.eq.0 )THEN
        if(output_filetype.eq.2)then
          tdef = '00:00Z03JUL0001'
        else
          tdef = '00:00Z03JUL2000'
        endif
        IF( radopt.ge.1 )THEN
          write(tdef( 1: 2),237) hour
          write(tdef( 4: 5),237) minute
          write(tdef( 7: 8),237) day
        if(output_filetype.eq.2)then
          write(tdef(12:15),238) 1
        else
          write(tdef(12:15),238) year
        endif
237       format(i2.2)
238       format(i4.4)
          IF( month.eq.1 )THEN
            write(tdef(9:11),239) 'JAN'
          ELSEIF( month.eq.2 )THEN
            write(tdef(9:11),239) 'FEB'
          ELSEIF( month.eq.3 )THEN
            write(tdef(9:11),239) 'MAR'
          ELSEIF( month.eq.4 )THEN
            write(tdef(9:11),239) 'APR'
          ELSEIF( month.eq.5 )THEN
            write(tdef(9:11),239) 'MAY'
          ELSEIF( month.eq.6 )THEN
            write(tdef(9:11),239) 'JUN'
          ELSEIF( month.eq.7 )THEN
            write(tdef(9:11),239) 'JUL'
          ELSEIF( month.eq.8 )THEN
            write(tdef(9:11),239) 'AUG'
          ELSEIF( month.eq.9 )THEN
            write(tdef(9:11),239) 'SEP'
          ELSEIF( month.eq.10 )THEN
            write(tdef(9:11),239) 'OCT'
          ELSEIF( month.eq.11 )THEN
            write(tdef(9:11),239) 'NOV'
          ELSEIF( month.eq.12 )THEN
            write(tdef(9:11),239) 'DEC'
          ELSE
            print *
            print *,'  Invalid value for MONTH '
            print *
            print *,'  Stopping CM1 .... '
            print *
            call stopcm1
          ENDIF
239       format(a3)
        ENDIF
      ENDIF

!-----------------------------------------------------------------------
!  GrADS descriptor files
!-----------------------------------------------------------------------

  grads_descriptors: IF( output_format.eq.1 )THEN

      IF(myid.eq.0)THEN

        allocate( varname(1000) )
        allocate( vardesc(1000) )

!----------------------------
! s file:
! accounts for both 2d and 3d variables:

    sout2d = 0
    s_out = 0

    ! all 2d variables MUST be listed first:

    if(output_rain   .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'rn      '
      vardesc(s_out) = 'accumulated rainfall (cm)     '
    endif
    if(output_sws    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'sws     '
      vardesc(s_out) = 'max wind speed lwst lvl (m/s) '
      s_out = s_out + 1
      varname(s_out) = 'svs     '
      vardesc(s_out) = 'max vert vort lwst lvl (s-1)  '
      s_out = s_out + 1
      varname(s_out) = 'sps     '
      vardesc(s_out) = 'min pressure lowest level (Pa)'
      s_out = s_out + 1
      varname(s_out) = 'srs     '
      vardesc(s_out) = 'max sfc rainwater (kg/kg)     '
      s_out = s_out + 1
      varname(s_out) = 'sgs     '
      vardesc(s_out) = 'max sfc graupel/hail (kg/kg)  '
      s_out = s_out + 1
      varname(s_out) = 'sus     '
      vardesc(s_out) = 'max w at 5 km AGL (m/s)       '
      s_out = s_out + 1
      varname(s_out) = 'shs     '
      vardesc(s_out) = 'max integrated uh (m2/s2)     '
    endif
    if(nrain.eq.2)then
      if(output_rain   .eq.1)then
        s_out = s_out + 1
        varname(s_out) = 'rn2     '
        vardesc(s_out) = 'translated rainfall (cm)      '
      endif
      if(output_sws    .eq.1)then
        s_out = s_out + 1
        varname(s_out) = 'sws2    '
        vardesc(s_out) = 'translated max wind (m/s)     '
        s_out = s_out + 1
        varname(s_out) = 'svs2    '
        vardesc(s_out) = 'translated max vorticity (s-1)'
        s_out = s_out + 1
        varname(s_out) = 'sps2    '
        vardesc(s_out) = 'translated min pressure (Pa)  '
        s_out = s_out + 1
        varname(s_out) = 'srs2    '
        vardesc(s_out) = 'translated max rainwater      '
        s_out = s_out + 1
        varname(s_out) = 'sgs2    '
        vardesc(s_out) = 'translated max graupel/hail   '
        s_out = s_out + 1
        varname(s_out) = 'sus2    '
        vardesc(s_out) = 'translated max w at 5 km (m/s)'
        s_out = s_out + 1
        varname(s_out) = 'shs2    '
        vardesc(s_out) = 'translated max integrated uh  '
      endif
    endif
    if(output_uh.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'uh      '
      vardesc(s_out) = 'integ. updraft helicity (m2/s2'
    endif
    if(output_coldpool.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'cpc     '
      vardesc(s_out) = 'cold pool intensity C (m/s)   '
      s_out = s_out + 1
      varname(s_out) = 'cph     '
      vardesc(s_out) = 'cold pool depth h (m AGL)     '
    endif
    if(output_sfcflx .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'thflux  '
      vardesc(s_out) = 'sfc theta flux (K m/s)        '
      s_out = s_out + 1
      varname(s_out) = 'qvflux  '
      vardesc(s_out) = 'sfc water vapor flux (g/g m/s)'
      s_out = s_out + 1
      varname(s_out) = 'cd      '
      vardesc(s_out) = 'exchange coeff for momentum   '
      s_out = s_out + 1
      varname(s_out) = 'ce      '
      vardesc(s_out) = 'exchange coeff for enthalpy   '
      s_out = s_out + 1
      varname(s_out) = 'tsk     '
      vardesc(s_out) = 'soil/ocean temperature (K)    '
    endif
    if(output_zs     .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'zs      '
      vardesc(s_out) = 'terrain height (m)            '
    endif
    if(output_dbz    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'cref    '
      vardesc(s_out) = 'composite reflectivity (dBZ)  '
    endif
    if(output_sfcparams.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'xland   '
      vardesc(s_out) = 'land/water flag (1=land,2=wtr)'
      s_out = s_out + 1
      varname(s_out) = 'lu      '
      vardesc(s_out) = 'land use index                '
      s_out = s_out + 1
      varname(s_out) = 'mavail  '
      vardesc(s_out) = 'surface moisture availability '
    endif
    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.oceanmodel.eq.2))then
      s_out = s_out + 1
      varname(s_out) = 'tmn     '
      vardesc(s_out) = 'deep-layer soil temperature (K'
      s_out = s_out + 1
      varname(s_out) = 'hfx     '
      vardesc(s_out) = 'heat flux at surface (W/m^2)  '
      s_out = s_out + 1
      varname(s_out) = 'qfx     '
      vardesc(s_out) = 'moisture flux at sfc (kg/m^2/s'
      s_out = s_out + 1
      varname(s_out) = 'gsw     '
      vardesc(s_out) = 'downward SW flux at sfc (W/m2)'
      s_out = s_out + 1
      varname(s_out) = 'glw     '
      vardesc(s_out) = 'downward LW flux at sfc (W/m2)'
    endif
    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2))then
      s_out = s_out + 1
      varname(s_out) = 'tslb1   '
      vardesc(s_out) = 'soil temp, layer 1 (K)        '
      s_out = s_out + 1
      varname(s_out) = 'tslb2   '
      vardesc(s_out) = 'soil temp, layer 2 (K)        '
      s_out = s_out + 1
      varname(s_out) = 'tslb3   '
      vardesc(s_out) = 'soil temp, layer 3 (K)        '
      s_out = s_out + 1
      varname(s_out) = 'tslb4   '
      vardesc(s_out) = 'soil temp, layer 4 (K)        '
      s_out = s_out + 1
      varname(s_out) = 'tslb5   '
      vardesc(s_out) = 'soil temp, layer 5 (K)        '
    endif
    if(output_sfcparams.eq.1.and.oceanmodel.eq.2)then
      s_out = s_out + 1
      varname(s_out) = 'tml     '
      vardesc(s_out) = 'ocean mixed layer temp (K)    '
      s_out = s_out + 1
      varname(s_out) = 'hml     '
      vardesc(s_out) = 'ocean mixed layer depth (m)   '
      s_out = s_out + 1
      varname(s_out) = 'huml    '
      vardesc(s_out) = 'ocean mixed layer u vel. (m/s)'
      s_out = s_out + 1
      varname(s_out) = 'hvml    '
      vardesc(s_out) = 'ocean mixed layer v vel. (m/s)'
    endif
    if(output_radten.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'radsw   '
      vardesc(s_out) = 'solar radiation at surface    '
      s_out = s_out + 1
      varname(s_out) = 'rnflx   '
      vardesc(s_out) = 'net radiation absorbed by sfc '
      s_out = s_out + 1
      varname(s_out) = 'radswnet'
      vardesc(s_out) = 'net solar radiation           '
      s_out = s_out + 1
      varname(s_out) = 'radlwin '
      vardesc(s_out) = 'incoming longwave radiation   '
    endif
    IF(output_sfcdiags.eq.1)THEN
      s_out = s_out + 1
      varname(s_out) = 'u10     '
      vardesc(s_out) = 'diagnostic 10m u wind (m/s)   '
      s_out = s_out + 1
      varname(s_out) = 'v10     '
      vardesc(s_out) = 'diagnostic 10m v wind (m/s)   '
      s_out = s_out + 1
      varname(s_out) = 't2      '
      vardesc(s_out) = 'diagnostic 2m temperature (K) '
      s_out = s_out + 1
      varname(s_out) = 'q2      '
      vardesc(s_out) = 'diagnostic 2m mixing ratio g/g'
      s_out = s_out + 1
      varname(s_out) = 'znt     '
      vardesc(s_out) = 'roughness length (m)          '
      s_out = s_out + 1
      varname(s_out) = 'ust     '
      vardesc(s_out) = 'u* in similarity theory (m/s) '
      s_out = s_out + 1
      varname(s_out) = 'hpbl    '
    if(ipbl.eq.1)then
      vardesc(s_out) = 'PBL height (m) (from PBL schem'
    else
      vardesc(s_out) = 'rough estimate of PBL hght (m)'
    endif
      s_out = s_out + 1
      varname(s_out) = 'zol     '
      vardesc(s_out) = 'z/L (z over Monin-Obukhov len)'
      s_out = s_out + 1
      varname(s_out) = 'mol     '
      vardesc(s_out) = 'T* (similarity theory) (K)    '
      s_out = s_out + 1
      varname(s_out) = 'br      '
      vardesc(s_out) = 'bulk Richardson No in sfc lay.'
    ENDIF

    ! done with 2d variables

    sout2d = s_out

    ! Now, all 3d variables:

    if(output_zh     .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'zh      '
      vardesc(s_out) = 'height on model levels (m)    '
    endif
    if(output_th     .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'th      '
      vardesc(s_out) = 'potential temp. (K)           '
    endif
    if(output_thpert .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'thpert  '
      vardesc(s_out) = 'potential temp. pert. (K)     '
    endif
    if(output_prs    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'prs     '
      vardesc(s_out) = 'pressure (Pa)                 '
    endif
    if(output_prspert.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'prspert '
      vardesc(s_out) = 'pressure pert. (Pa)           '
    endif
    if(output_pi     .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'pi      '
      vardesc(s_out) = 'nondimensional pressure       '
    endif
    if(output_pipert .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'pipert  '
      vardesc(s_out) = 'nondimensional pressure pert. '
    endif
    if(output_rho    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'rho     '
      vardesc(s_out) = 'density (kg/m^3)              '
    endif
    if(output_rhopert.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'rhopert '
      vardesc(s_out) = 'density pert. (kg/m^3)        '
    endif
    if(iptra         .eq.1)then
      do n=1,npt
        text1='pt      '
        if(n.lt.10)then
          write(text1(3:3),155) n
155       format(i1.1)
        else
          write(text1(3:4),154) n
154       format(i2.2)
        endif
        s_out = s_out + 1
        varname(s_out) = text1
        vardesc(s_out) = 'passive tracer                '
      enddo
    endif
    if(output_qv     .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'qv      '
      vardesc(s_out) = 'water vapor mixing ratio      '
    endif
    if(output_qvpert .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'qvpert  '
      vardesc(s_out) = 'qv pert                       '
    endif
    if(output_q      .eq.1)then
      do n=1,numq
        if(n.ne.nqv)then
          text1='        '
          text2='                              '
          write(text1(1:3),156) qname(n)
          write(text2(1:3),156) qname(n)
156       format(a3)
          s_out = s_out + 1
          varname(s_out) = text1
          vardesc(s_out) = text2
        endif
      enddo
    endif
    if(output_dbz    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'dbz     '
      vardesc(s_out) = 'reflectivity (dBZ)            '
    endif
    if(output_uinterp.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'uinterp '
      vardesc(s_out) = 'u interp. to scalar points    '
    endif
    if(output_vinterp.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'vinterp '
      vardesc(s_out) = 'v interp. to scalar points    '
    endif
    if(output_winterp.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'winterp '
      vardesc(s_out) = 'w interp. to scalar points    '
    endif
    if(output_vort.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'xvort   '
      vardesc(s_out) = 'horiz vorticity (x) (s^-1)    '
      s_out = s_out + 1
      varname(s_out) = 'yvort   '
      vardesc(s_out) = 'horiz vorticity (y) (s^-1)    '
      s_out = s_out + 1
      varname(s_out) = 'zvort   '
      vardesc(s_out) = 'vertical vorticity (s^-1)     '
    endif
    if(output_basestate.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'pi0     '
      vardesc(s_out) = 'base-state nondim. pressure   '
      s_out = s_out + 1
      varname(s_out) = 'th0     '
      vardesc(s_out) = 'base-state potential temp (K) '
      s_out = s_out + 1
      varname(s_out) = 'prs0    '
      vardesc(s_out) = 'base-state pressure (Pa)      '
      s_out = s_out + 1
      varname(s_out) = 'qv0     '
      vardesc(s_out) = 'base-state qv (kg/kg)         '
    endif
    if(output_dissten.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'dissten '
      vardesc(s_out) = 'dissipative heating tendency  '
    endif
    if(output_pblten.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'thpten  '
      vardesc(s_out) = 'pbl tendency:  theta          '
      s_out = s_out + 1
      varname(s_out) = 'qvpten  '
      vardesc(s_out) = 'pbl tendency:  qv             '
      s_out = s_out + 1
      varname(s_out) = 'qcpten  '
      vardesc(s_out) = 'pbl tendency:  qc             '
      s_out = s_out + 1
      varname(s_out) = 'qipten  '
      vardesc(s_out) = 'pbl tendency:  qi             '
      s_out = s_out + 1
      varname(s_out) = 'upten   '
      vardesc(s_out) = 'pbl tendency:  u              '
      s_out = s_out + 1
      varname(s_out) = 'vpten   '
      vardesc(s_out) = 'pbl tendency:  v              '
    endif
    if(output_radten.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'swten   '
      vardesc(s_out) = 'pot temp tendency, sw rad (K/s'
      s_out = s_out + 1
      varname(s_out) = 'lwten   '
      vardesc(s_out) = 'pot temp tendency, lw rad (K/s'
    endif

    sout3d = s_out - sout2d

!----------------------------
!  ready to write GrADS descriptor file:

  IF(s_out.ge.1)THEN
    string(totlen+1:totlen+1+12) = '_s.ctl'
    write(outfile,*) string
    open(unit=50,file=string,status='unknown')

  if(output_filetype.eq.1)then
    sstring(baselen+1:baselen+1+12) = '_s.dat'
  elseif(output_filetype.eq.2)then
    sstring(baselen+1:baselen+1+12) = '_%y4_s.dat'
  endif
    write(50,201) sstring
!!!    write(50,222)
    if(output_filetype.eq.2) write(50,221)
    write(50,202)
    write(50,203)
    if(stretch_x.ge.1)then
      write(50,214) nx
      do i=1,nx
        write(50,217) 0.001*0.5*(xfref(i)+xfref(i+1))
      enddo
    else
      write(50,204) nx,xh(1)/1000.0,dx/1000.0
    endif
    if(stretch_y.ge.1)then
      write(50,215) ny
      do j=1,ny
        write(50,217) 0.001*0.5*(yfref(j)+yfref(j+1))
      enddo
    else
      write(50,205) ny,yh(1)/1000.0,dy/1000.0
    endif
    if(stretch_z.eq.0)then
      write(50,206) nz,0.5*dz/1000.0,dz/1000.0
    else
      write(50,216) nz
      do k=1,nz
        write(50,217) 0.001*zh(1,1,k)
      enddo
    endif
  if(output_filetype.eq.1)then
    write(50,207) int(1+timax/tapfrq),tdef,max(1,int(tapfrq/60.0))
  elseif(output_filetype.eq.2)then
    write(50,227) int(1+timax/tapfrq),tdef
  endif
    write(50,208) s_out
    ! account for both 2d and 3d output files:
    do n=1,sout2d
      write(50,209) varname(n), 0,vardesc(n)
    enddo
    do n=sout2d+1,s_out
      write(50,209) varname(n),nk,vardesc(n)
    enddo
    write(50,210)
    close(unit=50)
  ENDIF

!-----------------------------------
! i file:  (for interpolated output when using terrain)
!   follows s file very closely:
!   no need to re-define varname,vardesc...

  IF(s_out.ge.1 .and. terrain_flag .and. output_interp.eq.1)THEN
    string(totlen+1:totlen+1+12) = '_i.ctl'
    write(outfile,*) string
    open(unit=50,file=string,status='unknown')

  if(output_filetype.eq.1)then
    sstring(baselen+1:baselen+1+12) = '_i.dat'
  elseif(output_filetype.eq.2)then
    sstring(baselen+1:baselen+1+12) = '_%y4_i.dat'
  endif

    write(50,201) sstring
!!!    write(50,222)
    if(output_filetype.eq.2) write(50,221)
    write(50,202)
    write(50,203)
    if(stretch_x.ge.1)then
      write(50,214) nx
      do i=1,nx
        write(50,217) 0.001*0.5*(xfref(i)+xfref(i+1))
      enddo
    else
      write(50,204) nx,xh(1)/1000.0,dx/1000.0
    endif
    if(stretch_y.ge.1)then
      write(50,215) ny
      do j=1,ny
        write(50,217) 0.001*0.5*(yfref(j)+yfref(j+1))
      enddo
    else
      write(50,205) ny,yh(1)/1000.0,dy/1000.0
    endif
    if(stretch_z.eq.0)then
      write(50,206) nz,0.5*dz/1000.0,dz/1000.0
    else
      write(50,216) nz
      do k=1,nz
        write(50,217) 0.001*zh(1,1,k)
      enddo
    endif
  if(output_filetype.eq.1)then
    write(50,207) int(1+timax/tapfrq),tdef,max(1,int(tapfrq/60.0))
  elseif(output_filetype.eq.2)then
    write(50,227) int(1+timax/tapfrq),tdef
  endif
    write(50,208) s_out
    ! account for both 2d and 3d output files:
    do n=1,sout2d
      write(50,209) varname(n), 0,vardesc(n)
    enddo
    do n=sout2d+1,s_out
      write(50,209) varname(n),nk,vardesc(n)
    enddo
    write(50,210)
    close(unit=50)
  ENDIF

!-----------------------------------
! u file:
! I have assumed that all variables are 3d for this file.

    u_out = 0

    if(output_u    .eq.1)then
      u_out = u_out + 1
      varname(u_out) = 'u       '
      vardesc(u_out) = 'E-W velocity (m/s)            '
    endif
    if(output_upert.eq.1)then
      u_out = u_out + 1
      varname(u_out) = 'upert   '
      vardesc(u_out) = 'u pert. (m/s)                 '
    endif
    if(output_basestate.eq.1)then
      u_out = u_out + 1
      varname(u_out) = 'u0      '
      vardesc(u_out) = 'base-state u (m/s)            '
    endif

  IF(u_out.ge.1)THEN
    string(totlen+1:totlen+1+12) = '_u.ctl'
    write(outfile,*) string
    open(unit=50,file=string,status='unknown')

  if(output_filetype.eq.1)then
    sstring(baselen+1:baselen+1+12) = '_u.dat'
  elseif(output_filetype.eq.2)then
    sstring(baselen+1:baselen+1+12) = '_%y4_u.dat'
  endif

    write(50,201) sstring
    if(output_filetype.eq.2) write(50,221)
    write(50,202)
    write(50,203)
    if(stretch_x.ge.1)then
      write(50,214) nx+1
      do i=1,nx+1
        write(50,217) 0.001*xfref(i)
      enddo
    else
      write(50,204) nx+1,xf(1)/1000.0,dx/1000.0
    endif
    if(stretch_y.ge.1)then
      write(50,215) ny
      do j=1,ny
        write(50,217) 0.001*0.5*(yfref(j)+yfref(j+1))
      enddo
    else
      write(50,205) ny,yh(1)/1000.0,dy/1000.0
    endif
    if(stretch_z.eq.0)then
      write(50,206) nz,0.5*dz/1000.0,dz/1000.0
    else
      write(50,216) nz
      do k=1,nz
        write(50,217) 0.001*zh(1,1,k)
      enddo
    endif
  if(output_filetype.eq.1)then
    write(50,207) int(1+timax/tapfrq),tdef,max(1,int(tapfrq/60.0))
  elseif(output_filetype.eq.2)then
    write(50,227) int(1+timax/tapfrq),tdef
  endif
    write(50,208) u_out
    ! assumes all variables are 3d:
    do n=1,u_out
      write(50,209) varname(n),nk,vardesc(n)
    enddo
    write(50,210)
    close(unit=50)
  ENDIF

!-----------------------------------
! v file:
! I have assumed that all variables are 3d for this file.

    v_out = 0

    if(output_v    .eq.1)then
      v_out = v_out + 1
      varname(v_out) = 'v       '
      vardesc(v_out) = 'N-S velocity (m/s)            '
    endif
    if(output_vpert.eq.1)then
      v_out = v_out + 1
      varname(v_out) = 'vpert   '
      vardesc(v_out) = 'v pert (m/s)                  '
    endif
    if(output_basestate.eq.1)then
      v_out = v_out + 1
      varname(v_out) = 'v0      '
      vardesc(v_out) = 'base-state v (m/s)            '
    endif

  IF(v_out.ge.1)THEN
    string(totlen+1:totlen+1+12) = '_v.ctl'
    write(outfile,*) string
    open(unit=50,file=string,status='unknown')

  if(output_filetype.eq.1)then
    sstring(baselen+1:baselen+1+12) = '_v.dat'
  elseif(output_filetype.eq.2)then
    sstring(baselen+1:baselen+1+12) = '_%y4_v.dat'
  endif

    write(50,201) sstring
    if(output_filetype.eq.2) write(50,221)
    write(50,202)
    write(50,203)
    if(stretch_x.ge.1)then
      write(50,214) nx
      do i=1,nx
        write(50,217) 0.001*0.5*(xfref(i)+xfref(i+1))
      enddo
    else
      write(50,204) nx,xh(1)/1000.0,dx/1000.0
    endif
    if(stretch_y.ge.1)then
      write(50,215) ny+1
      do j=1,ny+1
        write(50,217) 0.001*yfref(j)
      enddo
    else
      write(50,205) ny+1,yf(1)/1000.0,dy/1000.0
    endif
    if(stretch_z.eq.0)then
      write(50,206) nz,0.5*dz/1000.0,dz/1000.0
    else
      write(50,216) nz
      do k=1,nz
        write(50,217) 0.001*zh(1,1,k)
      enddo
    endif
  if(output_filetype.eq.1)then
    write(50,207) int(1+timax/tapfrq),tdef,max(1,int(tapfrq/60.0))
  elseif(output_filetype.eq.2)then
    write(50,227) int(1+timax/tapfrq),tdef
  endif
    write(50,208) v_out
    ! assumes all variables are 3d:
    do n=1,v_out
      write(50,209) varname(n),nk,vardesc(n)
    enddo
    write(50,210)
    close(unit=50)
  ENDIF

!-----------------------------------
! w file:
! I have assumed that all variables are 3d for this file.

    w_out = 0

    if(output_w  .eq.1)then
      w_out = w_out + 1
      varname(w_out) = 'w       '
      vardesc(w_out) = 'vertical velocity (m/s)       '
    endif
    if(output_tke.eq.1)then
      w_out = w_out + 1
      varname(w_out) = 'tke     '
      vardesc(w_out) = 'turb. kinetic energy (m^2/s^2)'
    endif
    if(output_km .eq.1)then
      w_out = w_out + 1
      varname(w_out) = 'kmh     '
      IF( ipbl.eq.1 )THEN
        vardesc(w_out) = 'turb. coef. for mo. (2D Smag.)'
      ELSE
        vardesc(w_out) = 'turb. coef. for mo. (m^2/s)   '
      ENDIF
      w_out = w_out + 1
      varname(w_out) = 'kmv     '
      IF( ipbl.eq.1 )THEN
        vardesc(w_out) = 'turb. coef. for mo. (from YSU)'
      ELSE
        vardesc(w_out) = 'turb. coef. for mo. (m^2/s)   '
      ENDIF
    endif
    if(output_kh .eq.1)then
      w_out = w_out + 1
      varname(w_out) = 'khh     '
      IF( ipbl.eq.1 )THEN
        vardesc(w_out) = 'turb. coef. for scalar (2D Sm)'
      ELSE
        vardesc(w_out) = 'turb. coef. for scalar (m^2/s)'
      ENDIF
      w_out = w_out + 1
      varname(w_out) = 'khv     '
      IF( ipbl.eq.1 )THEN
        vardesc(w_out) = 'turb. coef. for scalar (YSU)  '
      ELSE
        vardesc(w_out) = 'turb. coef. for scalar (m^2/s)'
      ENDIF
    endif

  IF(w_out.ge.1)THEN
    string(totlen+1:totlen+1+12) = '_w.ctl'
    write(outfile,*) string
    open(unit=50,file=string,status='unknown')

  if(output_filetype.eq.1)then
    sstring(baselen+1:baselen+1+12) = '_w.dat'
  elseif(output_filetype.eq.2)then
    sstring(baselen+1:baselen+1+12) = '_%y4_w.dat'
  endif

    write(50,201) sstring
    if(output_filetype.eq.2) write(50,221)
    write(50,202)
    write(50,203)
    if(stretch_x.ge.1)then
      write(50,214) nx
      do i=1,nx
        write(50,217) 0.001*0.5*(xfref(i)+xfref(i+1))
      enddo
    else
      write(50,204) nx,xh(1)/1000.0,dx/1000.0
    endif
    if(stretch_y.ge.1)then
      write(50,215) ny
      do j=1,ny
        write(50,217) 0.001*0.5*(yfref(j)+yfref(j+1))
      enddo
    else
      write(50,205) ny,yh(1)/1000.0,dy/1000.0
    endif
    if(stretch_z.eq.0)then
      write(50,206) nz+1,0.0,dz/1000.0
    else
      write(50,216) nz+1
      do k=1,nz+1
        write(50,217) 0.001*zf(1,1,k)
      enddo
    endif
  if(output_filetype.eq.1)then
    write(50,207) int(1+timax/tapfrq),tdef,max(1,int(tapfrq/60.0))
  elseif(output_filetype.eq.2)then
    write(50,227) int(1+timax/tapfrq),tdef
  endif
    write(50,208) w_out
    ! assumes all variables are 3d:
    do n=1,w_out
      write(50,209) varname(n),nk+1,vardesc(n)
    enddo
    write(50,210)
    close(unit=50)
  ENDIF

!-----------------------------------

    write(outfile,*)

201   format('dset ^',a70)
202   format('title CM1 output')
221   format('options template')
222   format('byteswapped')
203   format('undef -99999999.')
204   format('xdef ',i6,' linear ',f13.6,1x,f13.6)
214   format('xdef ',i6,' levels ')
205   format('ydef ',i6,' linear ',f13.6,1x,f13.6)
215   format('ydef ',i6,' levels ')
206   format('zdef ',i6,' linear ',f13.6,1x,f13.6)
216   format('zdef ',i6,' levels ')
217   format(2x,f13.6)
207   format('tdef ',i10,' linear ',a15,' ',i5,'MN')
227   format('tdef ',i10,' linear ',a15,' 1YR')
208   format('vars ',i4)
209   format(a8,2x,i6,'  99  ',a30)
210   format('endvars')

211   format(2x,f7.3)

!-----------------------------------------------------------------------

      call write_statsctl(tdef,qname,budname,1+nint(timax/max(statfrq,dtl)))

!-----------------------------------------------------------------------
!  Parcel data file:

      if(iprcl.eq.1.and.myid.eq.0)then

        string(totlen+1:totlen+1+12) = '_pdata.ctl  '
        write(outfile,*) string
        open(unit=50,file=string,status='unknown')

        sstring(baselen+1:baselen+1+12) = '_pdata.dat  '

        write(50,401) sstring
        write(50,402)
        write(50,403)
        write(50,404) nparcels
        write(50,405)
        write(50,406)
        write(50,407) 1+int(timax/prclfrq),tdef,max(1,int(prclfrq/60.0))
        write(50,408) npvals - 3
        write(50,409) 'x       ','x (m)                         '
        write(50,409) 'y       ','y (m)                         '
        write(50,409) 'z       ','z (m)                         '
        write(50,409) 'qv      ','water vapor mixing ratio      '
        write(50,409) 'qc      ','cloud water mixing ratio      '
        write(50,409) 'qr      ','rain water mixing ratio       '
        write(50,409) 'nm      ','squared Brunt-Vaisala frqncy  '
        write(50,409) 'u       ','u (m/s)                       '
        write(50,409) 'v       ','v (m/s)                       '
        write(50,409) 'w       ','w (m/s)                       '
        write(50,409) 'kh      ','turb. coef. for scalar (m^2/s)'
        write(50,409) 'the     ','theta-e (K)                   '
        write(50,409) 'b       ','buoyancy (m/s^2)              '
        write(50,409) 'dpdz    ','dpdz tendency (m/s^2)         '
        write(50,410)

401     format('dset ^',a70)
402     format('undef -9999.')
403     format('title ctl file for pdata.dat')
404     format('xdef ',i10,' linear 1 1')
405     format('ydef          1 linear 1 1')
406     format('zdef          1 linear 1 1')
407     format('tdef ',i10,' linear ',a15,' ',i5,'MN')
408     format('vars ',i6)
409     format(a8,' 1 99 ',a30)
410     format('endvars')

        close(unit=50)

      endif

!-----------------------------------------------------------------------

        deallocate( varname )
        deallocate( vardesc )

      ENDIF     ! endif for myid=0

#ifdef MPI
      call MPI_BCAST(sout2d,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sout3d,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(s_out ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(u_out ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(v_out ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(w_out ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(z_out ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

      write(outfile,*)
      write(outfile,*) '  sout2d = ',sout2d
      write(outfile,*) '  sout3d = ',sout3d
      write(outfile,*) '  s_out  = ',s_out
      write(outfile,*) '  u_out  = ',u_out
      write(outfile,*) '  v_out  = ',v_out
      write(outfile,*) '  w_out  = ',w_out
      write(outfile,*) '  z_out  = ',z_out

  ENDIF grads_descriptors

      write(outfile,*)

!-----------------------------------------------------------------------

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write_statsctl(tdef,qname,budname,numt)
      implicit none
      include 'input.incl'

      character*15, intent(in) :: tdef
      character*3, intent(in), dimension(maxq) :: qname
      character*6, intent(in), dimension(maxq) :: budname
      integer, intent(in) :: numt

      integer :: n
      character*8 text1
      character*30 text2
      character*50 fname

!  Subroutine to write GrADS stats descriptor file:
!-----------------------------------------------------------------------
!  write descriptors for stats file:

    string(totlen+1:totlen+1+12) = '_stats.ctl  '
    write(outfile,*) string
    open(unit=50,file=string,status='unknown')

    sstring(baselen+1:baselen+1+12) = '_stats.dat  '

    write(50,301) sstring
    write(50,302)
    write(50,303)
    write(50,304)
    write(50,305)
    write(50,306)
    write(50,307) numt,tdef,max(1,int(max(statfrq,60.0)/60.0))
    write(50,308) stat_out
    IF( adapt_dt.eq.1 )   write(50,309) 'dt      ','average timestep dt (s)       '
    if(stat_w      .eq.1) write(50,309) 'wmax    ','max vertical velocity (m/s)   '
    if(stat_w      .eq.1) write(50,309) 'wmin    ','min vertical velocity (m/s)   '
    if(stat_u      .eq.1) write(50,309) 'umax    ','max E-W velocity (m/s)        '
    if(stat_u      .eq.1) write(50,309) 'umin    ','min E-W velocity (m/s)        '
    if(stat_u      .eq.1) write(50,309) 'sumax   ','max E-W velocity lwst lvl(m/s)'
    if(stat_u      .eq.1) write(50,309) 'sumin   ','min E-W velocity lwst lvl(m/s)'
    if(stat_v      .eq.1) write(50,309) 'vmax    ','max N-S velocity (m/s)        '
    if(stat_v      .eq.1) write(50,309) 'vmin    ','min N-S velocity (m/s)        '
    if(stat_v      .eq.1) write(50,309) 'svmax   ','max N-S velocity lwst lvl(m/s)'
    if(stat_v      .eq.1) write(50,309) 'svmin   ','min N-S velocity lwst lvl(m/s)'
    if(stat_rmw    .eq.1) write(50,309) 'rmw     ','radius of maximum V (m)       '
    if(stat_pipert .eq.1) write(50,309) 'ppimax  ','max pi pert.                  '
    if(stat_pipert .eq.1) write(50,309) 'ppimin  ','min pi pert.                  '
    if(stat_prspert.eq.1) write(50,309) 'ppmax   ','max prs pert.(Pa)             '
    if(stat_prspert.eq.1) write(50,309) 'ppmin   ','min prs pert.(Pa)             '
    if(stat_thpert .eq.1) write(50,309) 'thpmax  ','max potential temp. pert. (K) '
    if(stat_thpert .eq.1) write(50,309) 'thpmin  ','min potential temp. pert. (K) '
    if(stat_thpert .eq.1) write(50,309) 'sthpmax ','max pot temp pert lwst lvl (K)'
    if(stat_thpert .eq.1) write(50,309) 'sthpmin ','min pot temp pert lwst lvl (K)'
    if(stat_q      .eq.1)then
      do n=1,numq
        text1='max     '
        text2='max                           '
        write(text1(4:6),156) qname(n)
        write(text2(5:7),156) qname(n)
        write(50,309) text1,text2
        text1='min     '
        text2='min                           '
        write(text1(4:6),156) qname(n)
        write(text2(5:7),156) qname(n)
        write(50,309) text1,text2
      enddo
    endif
    if(stat_tke    .eq.1) write(50,309) 'tkemax  ','max tke (m^2/s^2)             '
    if(stat_tke    .eq.1) write(50,309) 'tkemin  ','min tke (m^2/s^2)             '
    if(stat_km     .eq.1) write(50,309) 'kmhmax  ','max kmh (m^2/s)               '
    if(stat_km     .eq.1) write(50,309) 'kmhmin  ','min kmh (m^2/s)               '
    if(stat_km     .eq.1) write(50,309) 'kmvmax  ','max kmv (m^2/s)               '
    if(stat_km     .eq.1) write(50,309) 'kmvmin  ','min kmv (m^2/s)               '
    if(stat_kh     .eq.1) write(50,309) 'khhmax  ','max khh (m^2/s)               '
    if(stat_kh     .eq.1) write(50,309) 'khhmin  ','min khh (m^2/s)               '
    if(stat_kh     .eq.1) write(50,309) 'khvmax  ','max khv (m^2/s)               '
    if(stat_kh     .eq.1) write(50,309) 'khvmin  ','min khv (m^2/s)               '
    if(stat_div    .eq.1) write(50,309) 'divmax  ','max 3d divergence             '
    if(stat_div    .eq.1) write(50,309) 'divmin  ','min 3d divergence             '
    if(stat_rh     .eq.1) write(50,309) 'rhmax   ','max relative humidity         '
    if(stat_rh     .eq.1) write(50,309) 'rhmin   ','min relative humidity         '
    if(stat_rhi    .eq.1) write(50,309) 'rhimax  ','max relative humidity wrt ice '
    if(stat_rhi    .eq.1) write(50,309) 'rhimin  ','min relative humidity wrt ice '
    if(iptra       .eq.1)then
      do n=1,npt
        text1='maxpt   '
        text2='max pt                        '
        write(text1(6:6),157) n
        write(text2(7:7),157) n
157     format(i1)
        write(50,309) text1,text2
        text1='minpt   '
        text2='min pt                        '
        write(text1(6:6),157) n
        write(text2(7:7),157) n
        write(50,309) text1,text2
      enddo
    endif
    if(stat_the    .eq.1) write(50,309) 'themax  ','max theta-e below 10 km       '
    if(stat_the    .eq.1) write(50,309) 'themin  ','min theta-e below 10 km       '
    if(stat_the    .eq.1) write(50,309) 'sthemax ','max theta-e at lowest level   '
    if(stat_the    .eq.1) write(50,309) 'sthemin ','min theta-e at lowest level   '
    if(stat_cloud  .eq.1) write(50,309) 'qctop   ','max cloud top height (m)      '
    if(stat_cloud  .eq.1) write(50,309) 'qcbot   ','min cloud base height (m)     '
    if(stat_sfcprs .eq.1) write(50,309) 'sprsmax ','max pressure at lowest lvl (Pa'
    if(stat_sfcprs .eq.1) write(50,309) 'sprsmin ','min pressure at lowest lvl (Pa'
    if(stat_wsp    .eq.1) write(50,309) 'wspmax  ','max wind speed (m/s)          '
    if(stat_wsp    .eq.1) write(50,309) 'wspmin  ','min wind speed (m/s)          '
    if(stat_wsp    .eq.1) write(50,309) 'swspmax ','max wind speed lowst lvl (m/s)'
    if(stat_wsp    .eq.1) write(50,309) 'swspmin ','min wind speed lowst lvl (m/s)'
  IF(idrag.eq.1)THEN
    if(stat_wsp    .eq.1) write(50,309) 'wsp10max','max 10 m wind speed (m/s)     '
    if(stat_wsp    .eq.1) write(50,309) 'wsp10min','min 10 m wind speed (m/s)     '
  ENDIF
  IF( adapt_dt.eq.1 )THEN
    if(stat_cfl    .eq.1) write(50,309) 'cflmax  ','max Courant number (average)  '
  ELSE
    if(stat_cfl    .eq.1) write(50,309) 'cflmax  ','max Courant number            '
  ENDIF
    if(stat_cfl    .eq.1) write(50,309) 'kshmax  ','max horiz K stability factor  '
    if(stat_cfl    .eq.1) write(50,309) 'ksvmax  ','max vert K stability factor   '
    if(stat_vort   .eq.1) write(50,309) 'vortsfc ','max vert. vort. lwst lvl (1/s)'
    if(stat_vort   .eq.1) write(50,309) 'vort1km ','max vert. vort. at 1 km (1/s) '
    if(stat_vort   .eq.1) write(50,309) 'vort2km ','max vert. vort. at 2 km (1/s) '
    if(stat_vort   .eq.1) write(50,309) 'vort3km ','max vert. vort. at 3 km (1/s) '
    if(stat_vort   .eq.1) write(50,309) 'vort4km ','max vert. vort. at 4 km (1/s) '
    if(stat_vort   .eq.1) write(50,309) 'vort5km ','max vert. vort. at 5 km (1/s) '
    if(stat_tmass  .eq.1) write(50,309) 'tmass   ','total mass of (dry) air       '
    if(stat_tmois  .eq.1) write(50,309) 'tmois   ','total moisture                '
    if(stat_qmass  .eq.1)then
      do n=1,numq
        IF( (n.eq.nqv) .or.                                 &
            (n.ge.nql1.and.n.le.nql2) .or.                  &
            (n.ge.nqs1.and.n.le.nqs2.and.iice.eq.1) )THEN
          text1='mass    '
          text2='total mass of                 '
          write(text1( 5: 7),156) qname(n)
          write(text2(15:17),156) qname(n)
          write(50,309) text1,text2
        ENDIF
      enddo
    endif
    if(stat_tenerg .eq.1) write(50,309) 'ek      ','total kinetic energy          '
    if(stat_tenerg .eq.1) write(50,309) 'ei      ','total internal energy         '
    if(stat_tenerg .eq.1) write(50,309) 'ep      ','total potential energy        '
    if(stat_tenerg .eq.1) write(50,309) 'le      ','total latent energy (sort of) '
    if(stat_tenerg .eq.1) write(50,309) 'et      ','total energy                  '
    if(stat_mo     .eq.1) write(50,309) 'tmu     ','total E-W momentum            '
    if(stat_mo     .eq.1) write(50,309) 'tmv     ','total N-S momentum            '
    if(stat_mo     .eq.1) write(50,309) 'tmw     ','total vertical momentum       '
    if(stat_tmf    .eq.1) write(50,309) 'tmfu    ','total upward mass flux        '
    if(stat_tmf    .eq.1) write(50,309) 'tmfd    ','total downward mass flux      '
    if(stat_pcn    .eq.1)then
      do n=1,nbudget
        text1='        '
        text2='                              '
        write(text1(1:6),158) budname(n)
        write(text2(1:6),158) budname(n)
158     format(a6)
        write(50,309) text1,text2
      enddo
    endif
    if(stat_qsrc   .eq.1)then
      do n=1,numq
        text1='as      '
        text2='artificial source of          '
        write(text1( 3: 5),156) qname(n)
        write(text2(22:24),156) qname(n)
        write(50,309) text1,text2
      enddo
      do n=1,numq
        text1='bs      '
        text2='bndry source/sink of          '
        write(text1( 3: 5),156) qname(n)
        write(text2(22:24),156) qname(n)
        write(50,309) text1,text2
      enddo
    endif
    write(50,310)

156   format(a3)
301   format('dset ^',a70)
302   format('undef -9999.')
303   format('title ctl file for stats.dat')
304   format('xdef 1 linear 1 1')
305   format('ydef 1 linear 1 1')
306   format('zdef 1 linear 1 1')
307   format('tdef ',i10,' linear ',a15,' ',i5,'MN')
308   format('vars ',i6)
309   format(a8,' 1 99 ',a30)
310   format('endvars')

      close(unit=50)

      return
      end subroutine write_statsctl


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine writeout(fnum,nwrite,qname,xh,xf,uf,vf,sigma,zh,zf,mf,pi0,prs0,rho0,th0,qv0,u0,v0,  &
                          zs,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,cdu,cdv,ce,dum1,dum2,dum3,dum4, &
                          rho,prs,dbz,ua,dumu,va,dumv,wa,dumw,ppi,tha,        &
                   dissten,thpten,qvpten,qcpten,qipten,upten,vpten,           &
                          lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,   &
                          qa,kmh,kmv,khh,khv,tkea,swten,lwten,radsw,rnflx,radswnet,radlwin,pta,   &
                          num_soil_layers,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br)
      implicit none

      include 'input.incl'
      include 'constants.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      integer, intent(in) :: fnum,nwrite
      character*3, dimension(maxq), intent(in) :: qname
      real, intent(in), dimension(ib:ie) :: xh
      real, dimension(ib:ie+1), intent(in) :: xf,uf
      real, dimension(jb:je+1), intent(in) :: vf
      real, dimension(kb:ke), intent(in) :: sigma
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: zh
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: zf,mf
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: pi0,prs0,rho0,th0,qv0
      real, dimension(itb:ite,jtb:jte), intent(in) :: zs
      real, dimension(ib:ie,jb:je,nrain), intent(in) :: rain,sws,svs,sps,srs,sgs,sus,shs
      real, dimension(ib:ie,jb:je), intent(in) :: thflux,qvflux,cdu,cdv,ce
      real, dimension(ib:ie,jb:je,kb:ke), intent(inout) :: dum1,dum2,dum3,dum4
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: rho,prs,dbz
      real, dimension(ib:ie+1,jb:je,kb:ke), intent(in) :: u0,ua
      real, dimension(ib:ie+1,jb:je,kb:ke), intent(inout) :: dumu
      real, dimension(ib:ie,jb:je+1,kb:ke), intent(in) :: v0,va
      real, dimension(ib:ie,jb:je+1,kb:ke), intent(inout) :: dumv
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: wa
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(inout) :: dumw
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: ppi,tha
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: dissten
      real, dimension(ibb:ieb,jbb:jeb,kbb:keb), intent(in) :: thpten,qvpten,qcpten,qipten,upten,vpten
      integer, dimension(ibl:iel,jbl:jel), intent(in) :: lu_index
      real, dimension(ib:ie,jb:je), intent(in) :: tsk
      real, dimension(ibl:iel,jbl:jel), intent(in) :: xland,mavail,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw
      real, dimension(ibl:iel,jbl:jel,num_soil_layers), intent(in) :: tslb
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq), intent(in) :: qa
      real, dimension(ibc:iec,jbc:jec,kbc:kec), intent(in) :: kmh,kmv,khh,khv
      real, dimension(ibt:iet,jbt:jet,kbt:ket), intent(in) :: tkea
      real, dimension(ibr:ier,jbr:jer,kbr:ker), intent(in) :: swten,lwten
      real, dimension(ni,nj), intent(in) :: radsw,rnflx,radswnet,radlwin
      real, dimension(ibp:iep,jbp:jep,kbp:kep,npt), intent(in) :: pta
      integer, intent(in) :: num_soil_layers
      real, dimension(ibl:iel,jbl:jel), intent(in) :: u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br

      integer i,j,k,n,irec
#ifdef MPI
      character*80 sname,uname,vname,wname,zname
#endif

!--------------------------------------------------------------
!  writeout data on scalar-points

#ifndef MPI
    write(outfile,*)
  IF(output_filetype.eq.1)THEN
    if(s_out.ge.1)then
      if(fnum.eq.51)then
        string(totlen+1:totlen+1+12) = '_s.dat'
      elseif(fnum.eq.71)then
        string(totlen+1:totlen+1+12) = '_i.dat'
      endif
      write(outfile,*) string
      open(unit=fnum,file=string,form='unformatted',access='direct',   &
           recl=(ni*nj*4),status='unknown')
!!!      irec=1+(nwrite-1)*( nrain*output_rain + nrain*output_sws + output_zs + 4*output_sfcflx +   &
!!!                          nk*(s_out-nrain*output_rain-nrain*output_sws-output_zs-4*output_sfcflx) )
      irec=1+(nwrite-1)*( sout2d + nk*sout3d )
    endif
    if(u_out.ge.1)then
      string(totlen+1:totlen+1+12) = '_u.dat'
      write(outfile,*) string
      open(unit=52,file=string,form='unformatted',access='direct',   &
           recl=((ni+1)*nj*4),status='unknown')
    endif
    if(v_out.ge.1)then
      string(totlen+1:totlen+1+12) = '_v.dat'
      write(outfile,*) string
      open(unit=53,file=string,form='unformatted',access='direct',   &
           recl=(ni*(nj+1)*4),status='unknown')
    endif
    if(w_out.ge.1)then
      string(totlen+1:totlen+1+12) = '_w.dat'
      write(outfile,*) string
      open(unit=54,file=string,form='unformatted',access='direct',   &
           recl=(ni*nj*4),status='unknown')
    endif
  ELSEIF(output_filetype.eq.2)THEN
    if(s_out.ge.1)then
      if(fnum.eq.51)then
        string(totlen+1:totlen+1+12) = '_XXXX_s.dat'
      elseif(fnum.eq.71)then
        string(totlen+1:totlen+1+12) = '_XXXX_i.dat'
      endif
      write(string(totlen+2:totlen+5),102) nwrite
102   format(i4.4)
      write(outfile,*) string
      open(unit=fnum,file=string,form='unformatted',access='direct',   &
           recl=(ni*nj*4),status='unknown')
      irec=1
    endif
    if(u_out.ge.1)then
      string(totlen+1:totlen+1+12) = '_XXXX_u.dat'
      write(string(totlen+2:totlen+5),102) nwrite
      write(outfile,*) string
      open(unit=52,file=string,form='unformatted',access='direct',   &
           recl=((ni+1)*nj*4),status='unknown')
    endif
    if(v_out.ge.1)then
      string(totlen+1:totlen+1+12) = '_XXXX_v.dat'
      write(string(totlen+2:totlen+5),102) nwrite
      write(outfile,*) string
      open(unit=53,file=string,form='unformatted',access='direct',   &
           recl=(ni*(nj+1)*4),status='unknown')
    endif
    if(w_out.ge.1)then
      string(totlen+1:totlen+1+12) = '_XXXX_w.dat'
      write(string(totlen+2:totlen+5),102) nwrite
      write(outfile,*) string
      open(unit=54,file=string,form='unformatted',access='direct',   &
           recl=(ni*nj*4),status='unknown')
    endif
  ELSE
    write(outfile,*) '  Invalid option for output_filetype'
    call stopcm1
  ENDIF
#else
      irec=1

      sname = '                                                                                '
      uname = '                                                                                '
      vname = '                                                                                '
      wname = '                                                                                '
      zname = '                                                                                '

    if(strlen.gt.0)then
      sname(1:strlen) = output_path(1:strlen)
      uname(1:strlen) = output_path(1:strlen)
      vname(1:strlen) = output_path(1:strlen)
      wname(1:strlen) = output_path(1:strlen)
      zname(1:strlen) = output_path(1:strlen)
    endif

      sname(strlen+1:strlen+baselen) = output_basename(1:baselen)
      uname(strlen+1:strlen+baselen) = output_basename(1:baselen)
      vname(strlen+1:strlen+baselen) = output_basename(1:baselen)
      wname(strlen+1:strlen+baselen) = output_basename(1:baselen)
      zname(strlen+1:strlen+baselen) = output_basename(1:baselen)

      if(fnum.eq.51)then
        sname(totlen+1:totlen+1+16) = '_XXXX_YYYY_s.dat'
      elseif(fnum.eq.71)then
        sname(totlen+1:totlen+1+16) = '_XXXX_YYYY_i.dat'
      endif
      uname(totlen+1:totlen+1+16) = '_XXXX_YYYY_u.dat'
      vname(totlen+1:totlen+1+16) = '_XXXX_YYYY_v.dat'
      wname(totlen+1:totlen+1+16) = '_XXXX_YYYY_w.dat'
      zname(totlen+1:totlen+1+16) = '_XXXX_YYYY_z.dat'

      write(sname(totlen+2:totlen+ 5),100) myid
      write(sname(totlen+7:totlen+10),100) nwrite

      write(uname(totlen+2:totlen+ 5),100) myid
      write(uname(totlen+7:totlen+10),100) nwrite

      write(vname(totlen+2:totlen+ 5),100) myid
      write(vname(totlen+7:totlen+10),100) nwrite

      write(wname(totlen+2:totlen+ 5),100) myid
      write(wname(totlen+7:totlen+10),100) nwrite

      write(zname(totlen+2:totlen+ 5),100) myid
      write(zname(totlen+7:totlen+10),100) nwrite

100   format(i4.4)

      write(outfile,*)
      write(outfile,*) '  myid,sname=',myid,'   ',sname
      open(unit=fnum,file=sname,                   &
           form='unformatted',access='direct',   &
           recl=(ni*nj*4),status='unknown')

      if(u_out.ge.1)then
        write(outfile,*) '  myid,uname=',myid,'   ',uname
        open(unit=52,file=uname,                   &
             form='unformatted',access='direct',   &
             recl=((ni+1)*nj*4),status='unknown')
      endif

      if(v_out.ge.1)then
        write(outfile,*) '  myid,vname=',myid,'   ',vname
        open(unit=53,file=vname,                   &
             form='unformatted',access='direct',   &
             recl=(ni*(nj+1)*4),status='unknown')
      endif

      if(w_out.ge.1)then
        write(outfile,*) '  myid,wname=',myid,'   ',wname
        open(unit=54,file=wname,                   &
             form='unformatted',access='direct',   &
             recl=(ni*nj*4),status='unknown')
      endif
#endif

      if(output_rain.eq.1) call write2d(fnum,ni,nj,ngxy,irec,rain(ib,jb,1))
      if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,sws(ib,jb,1))
      if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,svs(ib,jb,1))
      if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,sps(ib,jb,1))
      if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,srs(ib,jb,1))
      if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,sgs(ib,jb,1))
      if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,sus(ib,jb,1))
      if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,shs(ib,jb,1))
      if(nrain.eq.2)then
        if(output_rain.eq.1) call write2d(fnum,ni,nj,ngxy,irec,rain(ib,jb,2))
        if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,sws(ib,jb,2))
        if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,svs(ib,jb,2))
        if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,sps(ib,jb,2))
        if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,srs(ib,jb,2))
        if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,sgs(ib,jb,2))
        if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,sus(ib,jb,2))
        if(output_sws .eq.1) call write2d(fnum,ni,nj,ngxy,irec,shs(ib,jb,2))
      endif
      if(output_uh.eq.1)then
        ! get height AGL:
        if( terrain_flag )then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            dum3(i,j,k) = zh(i,j,k)-zs(i,j)
            dumw(i,j,k) = zf(i,j,k)-zs(i,j)
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
            dumw(i,j,k) = zf(i,j,k)
          enddo
          enddo
          enddo
        endif
        call calcuh(uf,vf,dum3,dumw,ua,va,wa,dum1(ib,jb,1),dum2)
        call write2d(fnum,ni,nj,ngxy,irec,dum1(ib,jb,1))
      endif
      if(output_coldpool.eq.1)then
        call calccpch(zf,th0,qv0,dum1(ib,jb,1),dum1(ib,jb,2),tha,qa)
        call write2d(fnum,ni,nj,ngxy,irec,dum1(ib,jb,1))
        call write2d(fnum,ni,nj,ngxy,irec,dum1(ib,jb,2))
      endif
      if(output_sfcflx.eq.1) call write2d(fnum,ni,nj,ngxy,irec,thflux)
      if(output_sfcflx.eq.1) call write2d(fnum,ni,nj,ngxy,irec,qvflux)
      if(output_sfcflx.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1)=0.25*( (cdu(i,j)+cdu(i+1,j))   &
                            +(cdv(i,j)+cdv(i,j+1)) )
        enddo
        enddo
        call write2d(fnum,ni,nj,ngxy,irec,dum1(ib,jb,1))
      endif
      if(output_sfcflx.eq.1) call write2d(fnum,ni,nj,ngxy,irec,ce)
      if(output_sfcflx.eq.1) call write2d(fnum,ni,nj,ngxy,irec,tsk(ib,jb))
      if(output_zs  .eq.1) call write2d(fnum,ni,nj,ngxy,irec,zs)
      if(output_dbz   .eq.1)then
        call calccref(dum1(ib,jb,1),dbz)
        call write2d(fnum,ni,nj,ngxy,irec,dum1(ib,jb,1))
      endif

    if(output_sfcparams.eq.1)then
      call write2d(fnum,ni,nj,ngxy,irec,xland(ib,jb))
!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        dum1(i,j,1)=float(lu_index(i,j))
      enddo
      enddo
      call write2d(fnum,ni,nj,ngxy,irec,dum1(ib,jb,1))
      call write2d(fnum,ni,nj,ngxy,irec,mavail(ib,jb))
    endif
    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.oceanmodel.eq.2))then
      call write2d(fnum,ni,nj,ngxy,irec,tmn(ib,jb))
      call write2d(fnum,ni,nj,ngxy,irec,hfx(ib,jb))
      call write2d(fnum,ni,nj,ngxy,irec,qfx(ib,jb))
      call write2d(fnum,ni,nj,ngxy,irec,gsw(ib,jb))
      call write2d(fnum,ni,nj,ngxy,irec,glw(ib,jb))
    endif

    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2))then
      call write2d(fnum,ni,nj,ngxy,irec,tslb(ib,jb,1))
      call write2d(fnum,ni,nj,ngxy,irec,tslb(ib,jb,2))
      call write2d(fnum,ni,nj,ngxy,irec,tslb(ib,jb,3))
      call write2d(fnum,ni,nj,ngxy,irec,tslb(ib,jb,4))
      call write2d(fnum,ni,nj,ngxy,irec,tslb(ib,jb,5))
    endif

      if(output_sfcparams.eq.1.and.oceanmodel.eq.2)then
        call write2d(fnum,ni,nj,ngxy,irec,tml(ib,jb))
        call write2d(fnum,ni,nj,ngxy,irec,hml(ib,jb))
        call write2d(fnum,ni,nj,ngxy,irec,huml(ib,jb))
        call write2d(fnum,ni,nj,ngxy,irec,hvml(ib,jb))
      endif

      if( output_radten.eq.1 )then
        write(fnum,rec=irec) ((radsw(i,j),i=1,ni),j=1,nj)
        irec=irec+1
        write(fnum,rec=irec) ((rnflx(i,j),i=1,ni),j=1,nj)
        irec=irec+1
        write(fnum,rec=irec) ((radswnet(i,j),i=1,ni),j=1,nj)
        irec=irec+1
        write(fnum,rec=irec) ((radlwin(i,j),i=1,ni),j=1,nj)
        irec=irec+1
      endif

      IF(output_sfcdiags.eq.1)THEN
        call write2d(fnum,ni,nj,ngxy,irec,u10)
        call write2d(fnum,ni,nj,ngxy,irec,v10)
        call write2d(fnum,ni,nj,ngxy,irec, t2)
        call write2d(fnum,ni,nj,ngxy,irec, q2)
        call write2d(fnum,ni,nj,ngxy,irec,znt)
        call write2d(fnum,ni,nj,ngxy,irec,ust)
        call write2d(fnum,ni,nj,ngxy,irec,hpbl)
        call write2d(fnum,ni,nj,ngxy,irec,zol)
        call write2d(fnum,ni,nj,ngxy,irec,mol)
        call write2d(fnum,ni,nj,ngxy,irec,br)
      ENDIF

!--- 3D vars below here:

      dum1=zh
      if(fnum.eq.71)then
        do k=1,nk
        do j=1,nj
        do i=1,ni
!!!          dum1(i,j,k)=(k*dz-0.5*dz)-zs(i,j)
          dum1(i,j,k)=sigma(k)-zs(i,j)
        enddo
        enddo
        enddo
      endif
      if(output_zh  .eq.1) call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      if(output_th  .eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=th0(i,j,k)+tha(i,j,k)
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      endif
      dum1=tha
      if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
      if(output_thpert .eq.1) call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      dum1=prs
      if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
      if(output_prs    .eq.1) call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      if(output_prspert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=prs(i,j,k)-p00*(pi0(i,j,k)**cpdrd)
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      endif
      if(output_pi.eq.1)then  
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=pi0(i,j,k)+ppi(i,j,k)
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      endif
      dum1=ppi
      if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
      if(output_pipert .eq.1) call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      dum1=rho
      if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
      if(output_rho    .eq.1) call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      if(output_rhopert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=rho(i,j,k)-rho0(i,j,k)
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      endif
      if(iptra.eq.1)then
        do n=1,npt
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k)=pta(i,j,k,n)
          enddo
          enddo
          enddo
          if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
          call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        enddo
      endif
      if(imoist.eq.1)then
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=qa(i,j,k,nqv)
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        if(output_qv    .eq.1) call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        if(output_qvpert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k)=qa(i,j,k,nqv)-qv0(i,j,k)
          enddo
          enddo
          enddo
          if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
          call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        endif
        if(output_q.eq.1)then
          do n=1,numq
            if(n.ne.nqv)then
              do k=1,nk
              do j=1,nj
              do i=1,ni
                dum1(i,j,k)=qa(i,j,k,n)
              enddo
              enddo
              enddo
              if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
              call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
            endif
          enddo
        endif
        dum1=dbz
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        if(output_dbz   .eq.1) call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      endif
      if(output_uinterp.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=0.5*(ua(i,j,k)+ua(i+1,j,k))
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      endif
      if(output_vinterp.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=0.5*(va(i,j,k)+va(i,j+1,k))
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      endif
      if(output_winterp.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=0.5*(wa(i,j,k)+wa(i,j,k+1))
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      endif

      if(output_vort.eq.1)then
        call calcvort(xh,xf,uf,vf,zh,mf,zf,ua,va,wa,dum2,dum3,dum4,dum1)
        dum1=dum2
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        dum1=dum3
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        dum1=dum4
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      endif

      if(output_basestate.eq.1)then
        dum1=pi0
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        dum1=th0
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        dum1=prs0
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        dum1=qv0
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      endif

      if(output_dissten.eq.1)then
        dum1=dissten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      endif

      if(output_pblten.eq.1)then
        dum1=thpten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        dum1=qvpten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        dum1=qcpten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        dum1=qipten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        dum1=upten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        dum1=vpten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      endif

      if( output_radten.eq.1 )then
        dum1=swten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
        dum1=lwten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call write3d(fnum,ni,nj,nk,ngxy,ngz,irec,dum1)
      endif


!--------------------------------------------------------------
!  writeout data on u-points

      irec=1+(nwrite-1)*nk*u_out
      if(output_filetype.eq.2) irec=1
#ifdef MPI
      irec=1
#endif

      if(output_u    .eq.1) call write3d(52,ni+1,nj,nk,ngxy,ngz,irec,ua)

      if(output_upert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          dumu(i,j,k)=ua(i,j,k)-u0(i,j,k)
        enddo
        enddo
        enddo
        call write3d(52,ni+1,nj,nk,ngxy,ngz,irec,dumu)
      endif

      if(output_basestate.eq.1) call write3d(52,ni+1,nj,nk,ngxy,ngz,irec,u0)

!--------------------------------------------------------------
!  writeout data on v-points

      irec=1+(nwrite-1)*nk*v_out
      if(output_filetype.eq.2) irec=1
#ifdef MPI
      irec=1
#endif

      if(output_v    .eq.1) call write3d(53,ni,nj+1,nk,ngxy,ngz,irec,va)

      if(output_vpert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          dumv(i,j,k)=va(i,j,k)-v0(i,j,k)
        enddo
        enddo
        enddo
        call write3d(53,ni,nj+1,nk,ngxy,ngz,irec,dumv)
      endif

      if(output_basestate.eq.1) call write3d(53,ni,nj+1,nk,ngxy,ngz,irec,v0)

!--------------------------------------------------------------
!  writeout data on w-points

      irec=1+(nwrite-1)*(nk+1)*w_out
      if(output_filetype.eq.2) irec=1
#ifdef MPI
      irec=1
#endif

      if(output_w  .eq.1)                call write3d(54,ni,nj,nk+1,ngxy,ngz,irec,wa)
      if(output_tke.eq.1.and.iturb.eq.1) call write3d(54,ni,nj,nk+1,ngxy,ngz,irec,tkea)
      if(output_km .eq.1)                call write3d(54,ni,nj,nk+1,ngxy,ngz,irec,kmh)
      if(output_km .eq.1)                call write3d(54,ni,nj,nk+1,ngxy,ngz,irec,kmv)
      if(output_kh .eq.1)                call write3d(54,ni,nj,nk+1,ngxy,ngz,irec,khh)
      if(output_kh .eq.1)                call write3d(54,ni,nj,nk+1,ngxy,ngz,irec,khv)

!--------------------------------------------------------------

      write(outfile,*)
      write(outfile,*) 'Done Writing Data to File: nwrite=',nwrite
      write(outfile,*)

      close(unit=fnum)
      if(u_out.ge.1)then
        close(unit=52)
      endif
      if(v_out.ge.1)then
        close(unit=53)
      endif
      if(w_out.ge.1)then
        close(unit=54)
      endif
#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write2d(fileunit,numi,numj,ngxy,irec,var)
      implicit none

      integer, intent(in) :: fileunit,numi,numj,ngxy
      integer, intent(inout) :: irec
      real, intent(in), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy) :: var

      integer i,j

!----------------------------------

#ifdef DP
      write(fileunit,rec=irec) ((sngl(var(i,j)),i=1,numi),j=1,numj)
#else
      write(fileunit,rec=irec) ((var(i,j),i=1,numi),j=1,numj)
#endif
      irec=irec+1

!----------------------------------

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write3d(fileunit,numi,numj,numk,ngxy,ngz,irec,var)
      implicit none

      integer, intent(in) :: fileunit,numi,numj,numk,ngxy,ngz
      integer, intent(inout) :: irec
      real, intent(in), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,1-ngz:numk+ngz) :: var

      integer i,j,k

!----------------------------------

      do k=1,numk
#ifdef DP
        write(fileunit,rec=irec) ((sngl(var(i,j,k)),i=1,numi),j=1,numj)
#else
        write(fileunit,rec=irec) ((var(i,j,k),i=1,numi),j=1,numj)
#endif
        irec=irec+1
      enddo

!----------------------------------

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write_restart(nstep,nrec,prec,nwrite,nrst,nrad2d,num_soil_layers,      &
                               dt,mtime,radtim,qbudget,asq,bsq,                &
                               rain,sws,svs,sps,srs,sgs,sus,shs,tsk,radbcw,radbce,radbcs,radbcn,     &
                               ua,va,wa,ppi,tha,qa,tkea,swten,lwten,               &
                               radsw,rnflx,radswnet,radlwin,rad2d,                 &
                               lu_index,kpbl2d,psfc,u10,v10,hfx,qfx,xland,znt,ust, &
                               hpbl,wspd,psim,psih,gz1oz0,br,                      &
                               CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                      &
                               MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                   &
                               CK,CKA,CD,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                               f2d,gsw,glw,chklowq,capg,snowc,tslb,                &
                               tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,              &
                               pta,pdata,rtime)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      integer :: nstep,nrec,prec,nwrite,nrst
      integer :: nrad2d,num_soil_layers
      real :: dt
      real*8 :: mtime,radtim
      real*8, dimension(nbudget) :: qbudget
      real*8, dimension(numq) :: asq,bsq
      real, dimension(ib:ie,jb:je,nrain) :: rain,sws,svs,sps,srs,sgs,sus,shs
      real, dimension(ib:ie,jb:je) :: tsk
      real, dimension(jb:je,kb:ke) :: radbcw,radbce
      real, dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten
      real, dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin
      real, dimension(ni,nj,nrad2d) :: rad2d
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: lu_index
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: kpbl2d
      real, intent(in), dimension(ibl:iel,jbl:jel) :: psfc,u10,v10,hfx,qfx,xland,znt,ust, &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,   &
                                      CK,CKA,CD,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      f2d,gsw,glw,chklowq,capg,snowc
      real, intent(in), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(in), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      real, dimension(npvals,nparcels) :: pdata
      real rtime

      character*80 fname

!-----------------------------------------------------------------------

      fname = '                                                                                '
    if(strlen.gt.0)then
      fname(1:strlen) = output_path(1:strlen)
    endif
      fname(strlen+1:strlen+baselen) = output_basename(1:baselen)
      fname(totlen+1:totlen+1+18) = '_rst_XXXX_YYYY.dat'

      write(fname(totlen+ 6:totlen+ 9),101) myid
      write(fname(totlen+11:totlen+14),101) nrst
101   format(i4.4)

      write(outfile,*)
      write(outfile,*) '  Writing to restart file!'
      write(outfile,*) '  fname=',fname
      write(outfile,*)

      open(unit=50,file=fname,form='unformatted',status='unknown')

      write(50) nstep
      write(50) nrec
      write(50) prec
      write(50) nwrite
      write(50) nrst
      write(50) rtime
      write(50) dt
      write(50) mtime
      write(50) radtim
      write(50) qbudget
      write(50) asq
      write(50) bsq
      write(50) rain
      write(50) sws
      write(50) svs
      write(50) sps
      write(50) srs
      write(50) sgs
      write(50) sus
      write(50) shs
      write(50) tsk
      write(50) ua
      write(50) va
      write(50) wa
      write(50) ppi
      write(50) tha
      if(imoist.eq.1) write(50) qa
      if(iturb.eq.1) write(50) tkea
      if(radopt.eq.1)then
        write(50) swten
        write(50) lwten
        write(50) radsw
        write(50) rnflx
        write(50) radswnet
        write(50) radlwin
        write(50) rad2d
      endif
      if((oceanmodel.eq.2).or.(ipbl.eq.1).or.(sfcmodel.ge.1))then
        ! I don't know how many of these are really needed in restart
        ! files, but let's include them all for now ... just to be safe
        write(50) lu_index
        write(50) kpbl2d
        write(50) psfc
        write(50) u10
        write(50) v10
        write(50) hfx
        write(50) qfx
        write(50) xland
        write(50) znt
        write(50) ust
        write(50) hpbl
        write(50) wspd
        write(50) psim
        write(50) psih
        write(50) gz1oz0
        write(50) br
        write(50) CHS
        write(50) CHS2
        write(50) CQS2
        write(50) CPMM
        write(50) ZOL
        write(50) MAVAIL
        write(50) MOL
        write(50) RMOL
        write(50) REGIME
        write(50) LH
        write(50) TMN
        write(50) FLHC
        write(50) FLQC
        write(50) QGH
        write(50) CK
        write(50) CKA
        write(50) CD
        write(50) CDA
        write(50) USTM
        write(50) QSFC
        write(50) T2
        write(50) Q2
        write(50) TH2
        write(50) EMISS
        write(50) THC
        write(50) ALBD
        write(50) gsw
        write(50) glw
        write(50) chklowq
        write(50) capg
        write(50) snowc
        write(50) tslb
      endif
      if(oceanmodel.eq.2)then
        write(50) tml
        write(50) t0ml
        write(50) hml
        write(50) h0ml
        write(50) huml
        write(50) hvml
        write(50) tmoml
      endif
      if(iptra.eq.1) write(50) pta
      if(iprcl.eq.1) write(50) pdata
      if(irbc.eq.4.and.ibw.eq.1) write(50) radbcw
      if(irbc.eq.4.and.ibe.eq.1) write(50) radbce
      if(irbc.eq.4.and.ibs.eq.1) write(50) radbcs
      if(irbc.eq.4.and.ibn.eq.1) write(50) radbcn

      close(unit=50)

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine read_restart(nstep,nrec,prec,nwrite,nrst,nrad2d,num_soil_layers,   &
                              stattim,taptim,rsttim,  &
                              dt,mtime,radtim,qbudget,asq,bsq,              &
                              rain,sws,svs,sps,srs,sgs,sus,shs,tsk,radbcw,radbce,radbcs,radbcn,         &
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
      implicit none

      include 'input.incl'
      include 'constants.incl'

      integer :: nstep,nrec,prec,nwrite,nrst
      integer :: nrad2d,num_soil_layers
      real*8 :: stattim,taptim,rsttim
      real :: dt
      real*8 :: mtime,radtim
      real*8, dimension(nbudget) :: qbudget
      real*8, dimension(numq) :: asq,bsq
      real, dimension(ib:ie,jb:je,nrain) :: rain,sws,svs,sps,srs,sgs,sus,shs
      real, dimension(ib:ie,jb:je) :: tsk
      real, dimension(jb:je,kb:ke) :: radbcw,radbce
      real, dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten
      real, dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin
      real, dimension(ni,nj,nrad2d) :: rad2d
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: lu_index
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: kpbl2d
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: psfc,u10,v10,hfx,qfx,xland,znt,ust, &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,   &
                                      CK,CKA,CD,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      f2d,gsw,glw,chklowq,capg,snowc
      real, intent(inout), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      real, dimension(npvals,nparcels) :: pdata
      real rtime

      character*80 fname

!-----------------------------------------------------------------------

      fname = '                                                                                '
    if(strlen.gt.0)then
      fname(1:strlen) = output_path(1:strlen)
    endif
      fname(strlen+1:strlen+baselen) = output_basename(1:baselen)
      fname(totlen+1:totlen+1+18) = '_rst_XXXX_YYYY.dat'

      write(fname(totlen+ 6:totlen+ 9),101) myid
      write(fname(totlen+11:totlen+14),101) rstnum
101   format(i4.4)

      write(outfile,*)
      write(outfile,*) '  Reading from restart file!'
      write(outfile,*) '  fname=',fname
      write(outfile,*)

      open(unit=50,file=fname,form='unformatted',status='old')

      read(50) nstep
      read(50) nrec
      read(50) prec
      read(50) nwrite
      read(50) nrst
      read(50) rtime
      read(50) dt
      read(50) mtime
      read(50) radtim
      read(50) qbudget
      read(50) asq
      read(50) bsq
      read(50) rain
      read(50) sws
      read(50) svs
      read(50) sps
      read(50) srs
      read(50) sgs
      read(50) sus
      read(50) shs
      read(50) tsk
      read(50) ua
      read(50) va
      read(50) wa
      read(50) ppi
      read(50) tha
      if(imoist.eq.1) read(50) qa
      if(iturb.eq.1) read(50) tkea
      if(radopt.eq.1)then
        read(50) swten
        read(50) lwten
        read(50) radsw
        read(50) rnflx
        read(50) radswnet
        read(50) radlwin
        read(50) rad2d
      endif
      if((oceanmodel.eq.2).or.(ipbl.eq.1).or.(sfcmodel.ge.1))then
        ! I don't know how many of these are really needed in restart
        ! files, but let's include them all for now ... just to be safe
        read(50) lu_index
        read(50) kpbl2d
        read(50) psfc
        read(50) u10
        read(50) v10
        read(50) hfx
        read(50) qfx
        read(50) xland
        read(50) znt
        read(50) ust
        read(50) hpbl
        read(50) wspd
        read(50) psim
        read(50) psih
        read(50) gz1oz0
        read(50) br
        read(50) CHS
        read(50) CHS2
        read(50) CQS2
        read(50) CPMM
        read(50) ZOL
        read(50) MAVAIL
        read(50) MOL
        read(50) RMOL
        read(50) REGIME
        read(50) LH
        read(50) TMN
        read(50) FLHC
        read(50) FLQC
        read(50) QGH
        read(50) CK
        read(50) CKA
        read(50) CD
        read(50) CDA
        read(50) USTM
        read(50) QSFC
        read(50) T2
        read(50) Q2
        read(50) TH2
        read(50) EMISS
        read(50) THC
        read(50) ALBD
        read(50) gsw
        read(50) glw
        read(50) chklowq
        read(50) capg
        read(50) snowc
        read(50) tslb
      endif
      if(oceanmodel.eq.2)then
        read(50) tml
        read(50) t0ml
        read(50) hml
        read(50) h0ml
        read(50) huml
        read(50) hvml
        read(50) tmoml
      endif
      if(iptra.eq.1) read(50) pta
      if(iprcl.eq.1) read(50) pdata
      if(irbc.eq.4.and.ibw.eq.1) read(50) radbcw
      if(irbc.eq.4.and.ibe.eq.1) read(50) radbce
      if(irbc.eq.4.and.ibs.eq.1) read(50) radbcs
      if(irbc.eq.4.and.ibn.eq.1) read(50) radbcn

      close(unit=50)

!---------

      stattim=rtime+statfrq
      taptim=rtime+tapfrq
      rsttim=rtime+rstfrq
      prcltim=rtime+prclfrq
      if( output_format.eq.2 )then
        nrec=nrec-1
      else
        nrec=nrec-stat_out
      endif

      write(outfile,*) nrec,stattim,taptim,rsttim
      write(outfile,*)

!---------

      return
      end


