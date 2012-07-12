


!-------------------------------------------------------------
!
!  This subroutine writes data in NetCDF files.
!
!  Code originally written by Daniel Kirshbaum
!  Code last modified by George Bryan, 110111
!
!-------------------------------------------------------------


      subroutine writeout_cdf(nwrite,qname,sigma,sigmaf,xh,xf,uf,yh,yf,vf,mh,zh,mf,zf,    &
                          pi0,prs0,rho0,th0,qv0,u0,v0,                        &
                          zs,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,cdu,cdv,ce,dum1,dum2,dum3,dum4, &
                          rho,prs,dbz,ua,dumu,va,dumv,wa,dumw,ppi,tha,        &
                          qa,kmh,kmv,khh,khv,tkea,pta,num_soil_layers,   &
                          lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,   &
                          radsw,rnflx,radswnet,radlwin,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br,  &
                          dissten,thpten,qvpten,qcpten,qipten,upten,vpten,swten,lwten)
      implicit none

      include 'input.incl'
      include 'constants.incl'





      integer :: nwrite
      character*3, dimension(maxq) :: qname
      real, dimension(kb:ke) :: sigma
      real, dimension(kb:ke+1) :: sigmaf
      real, dimension(ib:ie) :: xh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je) :: yh
      real, dimension(jb:je+1) :: yf,vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh,zh,pi0,prs0,rho0,th0,qv0
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf,zf
      real, dimension(itb:ite,jtb:jte) :: zs
      real, dimension(ib:ie,jb:je,nrain) :: rain,sws,svs,sps,srs,sgs,sus,shs
      real, dimension(ib:ie,jb:je) :: thflux,qvflux,cdu,cdv,ce
      real, dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,rho,prs,dbz
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u0,ua,dumu
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v0,va,dumv
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa,dumw
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      integer, intent(in) :: num_soil_layers
      integer, dimension(ibl:iel,jbl:jel), intent(in) :: lu_index
      real, dimension(ib:ie,jb:je), intent(in) :: tsk
      real, dimension(ibl:iel,jbl:jel), intent(in) :: xland,mavail,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw
      real, dimension(ibl:iel,jbl:jel,num_soil_layers), intent(in) :: tslb
      real, dimension(ni,nj), intent(in) :: radsw,rnflx,radswnet,radlwin
      real, dimension(ibl:iel,jbl:jel), intent(in) :: u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: dissten
      real, dimension(ibb:ieb,jbb:jeb,kbb:keb), intent(in) :: thpten,qvpten,qcpten,qipten,upten,vpten
      real, dimension(ibr:ier,jbr:jer,kbr:ker), intent(in) :: swten,lwten



      include 'netcdf.inc'
      integer i,j,k,n,irec




      ! Users of GrADS might want to set coards to .true.
      logical, parameter :: coards = .false.

      integer :: cdfid,ncid    ! ID for the netCDF file to be created
      integer :: time_index,status,dimid,varid
      integer :: niid,njid,nkid,nip1id,njp1id,nkp1id,timeid,oneid
      character*80  cdf_out    ! Name of the netCDF output file
      character*8 chid

      real, dimension(ni,nj) :: d2d
      real, dimension(ni,nj,nk) :: ds
      real, dimension(ni+1,nj,nk) :: du
      real, dimension(ni,nj+1,nk) :: dv
      real, dimension(ni,nj,nk+1) :: dw

!-------------------------------------------------------------
! Declare and set integer values for the netCDF dimensions 
!-------------------------------------------------------------

      integer :: nip1,njp1
      real :: actual_time
      real :: x_min, x_max, y_min, y_max, z_min, z_max
      real :: x_mns, x_mxs, y_mns, y_mxs, z_mns, z_mxs

      logical :: allinfo

!--------------------------------------------------------------
! Initializing some things
!--------------------------------------------------------------

      time_index = nwrite
      if(output_filetype.eq.2) time_index = 1

    if(coards)then
      if(tapfrq.lt.60.0)then
        print *
        print *,'  Output frequency cannot be less than 60 s for coards format'
        print *
        call stopcm1
      endif
      actual_time = (nwrite-1)*tapfrq/60.0
    else
      actual_time = (nwrite-1)*tapfrq
    endif

      nip1=ni+1
      njp1=nj+1

!--------------------------------------------------------------
!  Write data to cdf file
!--------------------------------------------------------------


    IF(output_filetype.eq.1)THEN
      cdf_out = 'cm1out.nc'
    ELSEIF(output_filetype.eq.2)THEN
      cdf_out = 'cm1out_XXXX.nc'
      write(cdf_out(8:11),100) nwrite
    ELSE
      write(outfile,*) '  Invalid option for output_filetype'
      call stopcm1
    ENDIF








100   format(i4.4)

!--------------------------------------------------------------
!  if this is the start of a file, then do this stuff:

    allinfo = .false.
    IF(nwrite.eq.1) allinfo=.true.
    IF(output_filetype.eq.2) allinfo=.true.




    ifallinfo: IF(allinfo)THEN

!--------------------------------------------------------------
! Determine start and end locations for data
!--------------------------------------------------------------

      x_min = xf(1)
      x_max = xf(nip1)
      x_mns = xh(1)
      x_mxs = xh(ni)

      y_min = yf(1)
      y_max = yf(njp1)
      y_mns = yh(1)
      y_mxs = yh(nj)

    if(terrain_flag)then
      z_min = sigmaf(1)
      z_max = sigmaf(nkp1)
      z_mns = sigma(1)
      z_mxs = sigma(nk)
    else
      z_min = zf(1,1,1)
      z_max = zf(1,1,nkp1)
      z_mns = zh(1,1,1)
      z_mxs = zh(1,1,nk)
    endif

!-----------------------------------------------------------------------
!  BEGIN NEW:
      call disp_err( nf_create(cdf_out,nf_write,ncid), .true. )

      status = nf_def_dim(ncid,'ni',ni,niid)
      status = nf_def_dim(ncid,'nj',nj,njid)
      status = nf_def_dim(ncid,'nk',nk,nkid)
      status = nf_def_dim(ncid,'nip1',nip1,nip1id)
      status = nf_def_dim(ncid,'njp1',njp1,njp1id)
      status = nf_def_dim(ncid,'nkp1',nkp1,nkp1id)
      status = nf_def_dim(ncid,'time',nf_unlimited,timeid)
      status = nf_def_dim(ncid,'one',1,oneid)

    IF(icor.eq.1)THEN
      status = nf_def_var(ncid,"f_cor",nf_real,1,oneid,varid)
      status = nf_put_att_text(ncid,varid,"def",18,"Coriolis parameter")
      status = nf_put_att_text(ncid,varid,"units",3,"1/s")
    ENDIF

    if(.not.coards)then
      status = nf_def_var(ncid,"ztop",nf_real,1,oneid,varid)
      status = nf_put_att_text(ncid,varid,"units",2,"km")
    endif

    IF(coards)THEN

      status = nf_def_var(ncid,"time",nf_real,1,timeid,varid)
      status = nf_put_att_text(ncid,varid,"def",34,"time since beginning of simulation")
      status = nf_put_att_text(ncid,varid,"units",33,"minutes since 2000-07-03 00:00:00")

      status = nf_def_var(ncid,"ni",nf_real,1,niid,varid)
      status = nf_put_att_text(ncid,varid,"def",40,"west-east location of scalar grid points")
      status = nf_put_att_text(ncid,varid,"units",11,"degree_east")

      status = nf_def_var(ncid,"nip1",nf_real,1,nip1id,varid)
      status = nf_put_att_text(ncid,varid,"def",45,"west-east location of staggered u grid points")
      status = nf_put_att_text(ncid,varid,"units",11,"degree_east")

      status = nf_def_var(ncid,"nj",nf_real,1,njid,varid)
      status = nf_put_att_text(ncid,varid,"def",42,"south-north location of scalar grid points")
      status = nf_put_att_text(ncid,varid,"units",12,"degree_north")

      status = nf_def_var(ncid,"njp1",nf_real,1,njp1id,varid)
      status = nf_put_att_text(ncid,varid,"def",47,"south-north location of staggered v grid points")
      status = nf_put_att_text(ncid,varid,"units",12,"degree_north")

      status = nf_def_var(ncid,"nk",nf_real,1,nkid,varid)
      status = nf_put_att_text(ncid,varid,"def",36,"nominal height of scalar grid points")
      status = nf_put_att_text(ncid,varid,"units",2,"km")

      status = nf_def_var(ncid,"nkp1",nf_real,1,nkp1id,varid)
      status = nf_put_att_text(ncid,varid,"def",41,"nominal height of staggered w grid points")
      status = nf_put_att_text(ncid,varid,"units",2,"km")

    ELSE

      status = nf_def_var(ncid,"time",nf_real,1,timeid,varid)
      status = nf_put_att_text(ncid,varid,"def",34,"time since beginning of simulation")
      status = nf_put_att_text(ncid,varid,"units",33,"seconds since 2000-07-03 00:00:00")

      status = nf_def_var(ncid,"xh",nf_real,1,niid,varid)
      status = nf_put_att_text(ncid,varid,"def",40,"west-east location of scalar grid points")
      status = nf_put_att_text(ncid,varid,"units",2,"km")

      status = nf_def_var(ncid,"xf",nf_real,1,nip1id,varid)
      status = nf_put_att_text(ncid,varid,"def",45,"west-east location of staggered u grid points")
      status = nf_put_att_text(ncid,varid,"units",2,"km")

      status = nf_def_var(ncid,"yh",nf_real,1,njid,varid)
      status = nf_put_att_text(ncid,varid,"def",42,"south-north location of scalar grid points")
      status = nf_put_att_text(ncid,varid,"units",2,"km")

      status = nf_def_var(ncid,"yf",nf_real,1,njp1id,varid)
      status = nf_put_att_text(ncid,varid,"def",47,"south-north location of staggered v grid points")
      status = nf_put_att_text(ncid,varid,"units",2,"km")

      status = nf_def_var(ncid,"z",nf_real,1,nkid,varid)
      status = nf_put_att_text(ncid,varid,"def",36,"nominal height of scalar grid points")
      status = nf_put_att_text(ncid,varid,"units",2,"km")

      status = nf_def_var(ncid,"zf",nf_real,1,nkp1id,varid)
      status = nf_put_att_text(ncid,varid,"def",41,"nominal height of staggered w grid points")
      status = nf_put_att_text(ncid,varid,"units",2,"km")

    ENDIF

!--------------------------------------------------------
!  Just to be sure:

        status = nf_inq_dimid(ncid,'time',timeid)
        status = nf_inq_dimid(ncid,'ni',niid)
        status = nf_inq_dimid(ncid,'nj',njid)
        status = nf_inq_dimid(ncid,'nk',nkid)
        status = nf_inq_dimid(ncid,'nip1',nip1id)
        status = nf_inq_dimid(ncid,'njp1',njp1id)
        status = nf_inq_dimid(ncid,'nkp1',nkp1id)

!--------------------------------------------------------

      if (output_zs.eq.1.and.terrain_flag) then
        status = nf_def_var(ncid,"zs",nf_real,2,(/niid,njid/),varid)
        status = nf_put_att_text(ncid,varid,"def",14,"terrain height")
        status = nf_put_att_text(ncid,varid,"units",1,"m")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
      endif

      if (output_zh.eq.1.and.terrain_flag) then
        status = nf_def_var(ncid,"zh",nf_real,3,(/niid,njid,nkid/),varid)
        status = nf_put_att_text(ncid,varid,"def",43,"height (above ground) of scalar grid points")
        status = nf_put_att_text(ncid,varid,"units",1,"m")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif

      if (output_basestate.eq.1) then

        status = nf_def_var(ncid,"pi0",nf_real,3,(/niid,njid,nkid/),varid)
        status = nf_put_att_text(ncid,varid,"def",34,"base-state nondimensional pressure")
        status = nf_put_att_text(ncid,varid,"units",13,"dimensionless")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

        status = nf_def_var(ncid,"th0",nf_real,3,(/niid,njid,nkid/),varid)
        status = nf_put_att_text(ncid,varid,"def",32,"base-state potential temperature")
        status = nf_put_att_text(ncid,varid,"units",1,"K")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

        status = nf_def_var(ncid,"prs0",nf_real,3,(/niid,njid,nkid/),varid)
        status = nf_put_att_text(ncid,varid,"def",19,"base-state pressure")
        status = nf_put_att_text(ncid,varid,"units",2,"Pa")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

        status = nf_def_var(ncid,"qv0",nf_real,3,(/niid,njid,nkid/),varid)
        status = nf_put_att_text(ncid,varid,"def",35,"base-state water vapor mixing ratio")
        status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

        status = nf_def_var(ncid,"u0",nf_real,3,(/nip1id,njid,nkid/),varid)
        status = nf_put_att_text(ncid,varid,"def",34,"base-state x-component of velocity")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_min)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

        status = nf_def_var(ncid,"v0",nf_real,3,(/niid,njp1id,nkid/),varid)
        status = nf_put_att_text(ncid,varid,"def",34,"base-state y-component of velocity")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_min)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

      endif

! 2D variables

      if(output_rain.eq.1)then
        status = nf_def_var(ncid,"rain",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",28,"accumulated surface rainfall")
        status = nf_put_att_text(ncid,varid,"units",2,"cm")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
      endif
      if(output_sws.eq.1) then
        status = nf_def_var(ncid,"sws",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",29,"max windspeed at lowest level")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_def_var(ncid,"svs",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",34,"max vert vorticity at lowest level")
        status = nf_put_att_text(ncid,varid,"units",3,"s-1")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_def_var(ncid,"sps",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",28,"min pressure at lowest level")
        status = nf_put_att_text(ncid,varid,"units",2,"Pa")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_def_var(ncid,"srs",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",21,"max surface rainwater")
        status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_def_var(ncid,"sgs",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",24,"max surface graupel/hail")
        status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_def_var(ncid,"sus",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",17,"max w at 5 km AGL")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_def_var(ncid,"shs",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",31,"max integrated updraft helicity")
        status = nf_put_att_text(ncid,varid,"units",5,"m2/s2")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
      endif
      IF(nrain.eq.2)THEN
        if(output_rain.eq.1)then
          status = nf_def_var(ncid,"rain2",nf_real,3,(/niid,njid,timeid/),varid)
          status = nf_put_att_text(ncid,varid,"def",59,"accumulated surface rainfall, translated with moving domain")
          status = nf_put_att_text(ncid,varid,"units",2,"cm")
          status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
          status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        endif
        if(output_sws.eq.1) then
          status = nf_def_var(ncid,"sws2",nf_real,3,(/niid,njid,timeid/),varid)
          status = nf_put_att_text(ncid,varid,"def",60,"max windspeed at lowest level, translated with moving domain")
          status = nf_put_att_text(ncid,varid,"units",3,"m/s")
          status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
          status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
          status = nf_def_var(ncid,"svs2",nf_real,3,(/niid,njid,timeid/),varid)
          status = nf_put_att_text(ncid,varid,"def",60,"max vorticity at lowest level, translated with moving domain")
          status = nf_put_att_text(ncid,varid,"units",3,"s-1")
          status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
          status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
          status = nf_def_var(ncid,"sps2",nf_real,3,(/niid,njid,timeid/),varid)
          status = nf_put_att_text(ncid,varid,"def",59,"min pressure at lowest level, translated with moving domain")
          status = nf_put_att_text(ncid,varid,"units",2,"Pa")
          status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
          status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
          status = nf_def_var(ncid,"srs2",nf_real,3,(/niid,njid,timeid/),varid)
          status = nf_put_att_text(ncid,varid,"def",52,"max surface rainwater, translated with moving domain")
          status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")
          status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
          status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
          status = nf_def_var(ncid,"sgs2",nf_real,3,(/niid,njid,timeid/),varid)
          status = nf_put_att_text(ncid,varid,"def",55,"max surface graupel/hail, translated with moving domain")
          status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")
          status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
          status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
          status = nf_def_var(ncid,"sus2",nf_real,3,(/niid,njid,timeid/),varid)
          status = nf_put_att_text(ncid,varid,"def",48,"max w at 5 km AGL, translated with moving domain")
          status = nf_put_att_text(ncid,varid,"units",3,"m/s")
          status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
          status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
          status = nf_def_var(ncid,"shs2",nf_real,3,(/niid,njid,timeid/),varid)
          status = nf_put_att_text(ncid,varid,"def",42,"translated max integrated updraft helicity")
          status = nf_put_att_text(ncid,varid,"units",5,"m2/s2")
          status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
          status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        endif
      ENDIF

      IF(output_uh.eq.1)THEN
        status = nf_def_var(ncid,"uh",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",27,"integrated updraft helicity")
        status = nf_put_att_text(ncid,varid,"units",5,"m2/s2")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
      ENDIF

      IF (output_dbz.eq.1) THEN
        status = nf_def_var(ncid,"cref",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",28,"composite reflectivity (dBZ)")
        status = nf_put_att_text(ncid,varid,"units",3,"dBZ")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
      ENDIF

      IF (output_coldpool.eq.1) THEN
        status = nf_def_var(ncid,"cpc",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",21,"cold pool intensity C")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_def_var(ncid,"cph",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",17,"cold pool depth h")
        status = nf_put_att_text(ncid,varid,"units",5,"m AGL")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
      ENDIF

      if (output_sfcflx.eq.1) then
        status = nf_def_var(ncid,"thflux",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",40,"surface potential temperature flux")
        status = nf_put_att_text(ncid,varid,"units",10,"K m s^{-1}")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"qvflux",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",40,"surface water vapor flux")
        status = nf_put_att_text(ncid,varid,"units",19,"kg kg^{-1} m s^{-1}")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"cd",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",2,"cd")
        status = nf_put_att_text(ncid,varid,"units",14,"nondimensional")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"ce",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",2,"ce")
        status = nf_put_att_text(ncid,varid,"units",14,"nondimensional")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"tsk",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",22,"soil/ocean temperature")
        status = nf_put_att_text(ncid,varid,"units",1,"K")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
      endif

      if(output_sfcparams.eq.1)then
        status = nf_def_var(ncid,"xland",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",32,"land/water flag (1=land,2=water)")
        status = nf_put_att_text(ncid,varid,"units",12,"integer flag")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"lu",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",14,"land use index")
        status = nf_put_att_text(ncid,varid,"units",12,"integer flag")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"mavail",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",29,"surface moisture availability")
        status = nf_put_att_text(ncid,varid,"units",12,"integer flag")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
      endif

      if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.oceanmodel.eq.2))then
        status = nf_def_var(ncid,"tmn",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",27,"deep-layer soil temperature")
        status = nf_put_att_text(ncid,varid,"units",1,"K")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"hfx",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",20,"heat flux at surface")
        status = nf_put_att_text(ncid,varid,"units",5,"W/m^2")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"qfx",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",21,"surface moisture flux")
        status = nf_put_att_text(ncid,varid,"units",5,"W/m^2")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"gsw",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",27,"downward SW flux at surface")
        status = nf_put_att_text(ncid,varid,"units",5,"W/m^2")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"glw",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",27,"downward LW flux at surface")
        status = nf_put_att_text(ncid,varid,"units",5,"W/m^2")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
      endif

      if((output_sfcparams.eq.1).and.(sfcmodel.eq.2))then
        status = nf_def_var(ncid,"tslb1",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",18,"soil temp, layer 1")
        status = nf_put_att_text(ncid,varid,"units",1,"K")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"tslb2",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",18,"soil temp, layer 2")
        status = nf_put_att_text(ncid,varid,"units",1,"K")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"tslb3",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",18,"soil temp, layer 3")
        status = nf_put_att_text(ncid,varid,"units",1,"K")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"tslb4",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",18,"soil temp, layer 4")
        status = nf_put_att_text(ncid,varid,"units",1,"K")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"tslb5",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",18,"soil temp, layer 5")
        status = nf_put_att_text(ncid,varid,"units",1,"K")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
      endif

      if(output_sfcparams.eq.1.and.oceanmodel.eq.2)then
        status = nf_def_var(ncid,"tml",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",29,"ocean mixed layer temperature")
        status = nf_put_att_text(ncid,varid,"units",1,"K")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"hml",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",23,"ocean mixed layer depth")
        status = nf_put_att_text(ncid,varid,"units",1,"m")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"huml",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",24,"ocean mixed layer u vel.")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"hvml",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",24,"ocean mixed layer v vel.")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
      endif

      if(output_radten.eq.1)then
        status = nf_def_var(ncid,"radsw",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",26,"solar radiation at surface")
        status = nf_put_att_text(ncid,varid,"units",5,"w/m^2")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"rnflx",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",33,"net radiation absorbed by surface")
        status = nf_put_att_text(ncid,varid,"units",5,"W/m^2")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"radswnet",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",19,"net solar radiation")
        status = nf_put_att_text(ncid,varid,"units",5,"W/m^2")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"radlwin",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",27,"incoming longwave radiation")
        status = nf_put_att_text(ncid,varid,"units",5,"W/m^2")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
      endif

      IF(output_sfcdiags.eq.1)THEN
        status = nf_def_var(ncid,"u10",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",22,"diagnostic 10 m u wind")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"v10",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",22,"diagnostic 10 m v wind")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"t2",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",26,"diagnostic 2 m temperature")
        status = nf_put_att_text(ncid,varid,"units",1,"K")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"q2",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",27,"diagnostic 2 m mixing ratio")
        status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"znt",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",16,"roughness length")
        status = nf_put_att_text(ncid,varid,"units",1,"m")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"ust",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",23,"u* in similarity theory")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"hpbl",nf_real,3,(/niid,njid,timeid/),varid)
      if(ipbl.eq.1)then
        status = nf_put_att_text(ncid,varid,"def",28,"PBL height (from PBL scheme)")
      else
        status = nf_put_att_text(ncid,varid,"def",28,"rough estimate of PBL height")
      endif
        status = nf_put_att_text(ncid,varid,"units",1,"m")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"zol",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",33,"z/L (z over Monin-Obukhov length)")
        status = nf_put_att_text(ncid,varid,"units",3,"   ")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"mol",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",22,"T* (similarity theory)")
        status = nf_put_att_text(ncid,varid,"units",1,"K")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)

        status = nf_def_var(ncid,"br",nf_real,3,(/niid,njid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",39,"bulk Richardson number in surface layer")
        status = nf_put_att_text(ncid,varid,"units",1," ")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
      ENDIF

! 3D variables

      if(output_u.eq.1)then
        status = nf_def_var(ncid,"u",nf_real,4,(/nip1id,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",23,"velocity in x-direction")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_min)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
      if(output_upert.eq.1)then
        status = nf_def_var(ncid,"upert",nf_real,4,(/nip1id,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",36,"perturbation velocity in x-direction")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_min)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
      if(output_v.eq.1)then
        status = nf_def_var(ncid,"v",nf_real,4,(/niid,njp1id,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",23,"velocity in y-direction")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_min)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
      if(output_vpert.eq.1)then
        status = nf_def_var(ncid,"vpert",nf_real,4,(/niid,njp1id,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",36,"perturbation velocity in y-direction")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_min)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
      if(output_w.eq.1)then
        status = nf_def_var(ncid,"w",nf_real,4,(/niid,njid,nkp1id,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",23,"velocity in z-direction")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_min)
      endif

      if(output_rho.eq.1)then
        status = nf_def_var(ncid,"rho",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",18,"density of dry air")
        status = nf_put_att_text(ncid,varid,"units",6,"kg/m^3")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
      if(output_rhopert.eq.1)then
        status = nf_def_var(ncid,"rhopert",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",31,"perturbation density of dry air")
        status = nf_put_att_text(ncid,varid,"units",6,"kg/m^3")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
      if(output_th.eq.1)then
        status = nf_def_var(ncid,"th",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",21,"potential temperature")
        status = nf_put_att_text(ncid,varid,"units",1,"K")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
      if(output_thpert.eq.1)then
        status = nf_def_var(ncid,"thpert",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",34,"perturbation potential temperature")
        status = nf_put_att_text(ncid,varid,"units",1,"K")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif

      if(output_pi.eq.1)then
        status = nf_def_var(ncid,"pi",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",23,"nondimensional pressure")
        status = nf_put_att_text(ncid,varid,"units",13,"dimensionless")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
      if(output_pipert.eq.1)then
        status = nf_def_var(ncid,"pipert",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",36,"perturbation nondimensional pressure")
        status = nf_put_att_text(ncid,varid,"units",13,"dimensionless")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
      if(output_prs.eq.1)then
        status = nf_def_var(ncid,"prs",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",8,"pressure")
        status = nf_put_att_text(ncid,varid,"units",2,"Pa")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
      if(output_prspert.eq.1)then
        status = nf_def_var(ncid,"prspert",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",21,"perturbation pressure")
        status = nf_put_att_text(ncid,varid,"units",2,"Pa")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
    IF(imoist.eq.1)THEN
      if(output_qv.eq.1)then
        status = nf_def_var(ncid,"qv",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",24,"water vapor mixing ratio")
        status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
      if(output_qvpert.eq.1)then
        status = nf_def_var(ncid,"qvpert",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",37,"perturbation water vapor mixing ratio")
        status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
        if(output_q.eq.1)then
        do n=1,numq
          if(n.ne.nqv)then
            status = nf_def_var(ncid,qname(n),nf_real,4,(/niid,njid,nkid,timeid/),varid)
            if(idm.eq.1.and.n.ge.nnc1.and.n.le.nnc2)then
              status = nf_put_att_text(ncid,varid,"def",20,"number concentration")
              status = nf_put_att_text(ncid,varid,"units",7,"kg^{-1}")
            else
              status = nf_put_att_text(ncid,varid,"def",12,"mixing ratio")
              status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")
            endif
            status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
            status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
            status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
          endif
        enddo
      endif
      if(output_dbz.eq.1)then
        status = nf_def_var(ncid,"dbz",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",12,"reflectivity")
        status = nf_put_att_text(ncid,varid,"units",3,"dBZ")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
    ENDIF
    IF(iptra.eq.1)THEN
      do n=1,npt
        chid = 'pt      '
        write(chid(3:4),111) n
111     format(i2.2)
        status = nf_def_var(ncid,chid,nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",30,"mixing ratio of passive tracer")
        status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      enddo
    ENDIF
    IF(iturb.ge.1)THEN
      IF(output_km.eq.1)THEN
        !----
        status = nf_def_var(ncid,"kmh",nf_real,4,(/niid,njid,nkp1id,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",64,"eddy mixing coefficient for momentum in the horizontal direction")
        status = nf_put_att_text(ncid,varid,"units",5,"m^2/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_min)
        !----
        status = nf_def_var(ncid,"kmv",nf_real,4,(/niid,njid,nkp1id,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",62,"eddy mixing coefficient for momentum in the vertical direction")
        status = nf_put_att_text(ncid,varid,"units",5,"m^2/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_min)
        !----
      ENDIF
      IF(output_kh.eq.1)THEN
        !----
        status = nf_def_var(ncid,"khh",nf_real,4,(/niid,njid,nkp1id,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",63,"eddy mixing coefficient for scalars in the horizontal direction")
        status = nf_put_att_text(ncid,varid,"units",5,"m^2/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_min)
        !----
        status = nf_def_var(ncid,"khv",nf_real,4,(/niid,njid,nkp1id,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",61,"eddy mixing coefficient for scalars in the vertical direction")
        status = nf_put_att_text(ncid,varid,"units",5,"m^2/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_min)
        !----
      ENDIF
    ENDIF
      IF((iturb.eq.1).and.(output_tke.eq.1))THEN
        status = nf_def_var(ncid,"tke",nf_real,4,(/niid,njid,nkp1id,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",33,"subgrid turbulence kinetic energy")
        status = nf_put_att_text(ncid,varid,"units",7,"m^2/s^2")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_min)
      ENDIF

      if(output_uinterp.eq.1)then
        status = nf_def_var(ncid,"uinterp",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",54,"velocity in x-direction, interpolated to scalar points")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
      if(output_vinterp.eq.1)then
        status = nf_def_var(ncid,"vinterp",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",54,"velocity in y-direction, interpolated to scalar points")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif
      if(output_winterp.eq.1)then
        status = nf_def_var(ncid,"winterp",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",54,"velocity in z-direction, interpolated to scalar points")
        status = nf_put_att_text(ncid,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif

      if(output_vort.eq.1)then
        status = nf_def_var(ncid,"xvort",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",24,"horizontal vorticity (x)")
        status = nf_put_att_text(ncid,varid,"units",4,"s^-1")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

        status = nf_def_var(ncid,"yvort",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",24,"horizontal vorticity (y)")
        status = nf_put_att_text(ncid,varid,"units",4,"s^-1")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

        status = nf_def_var(ncid,"zvort",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",18,"vertical vorticity")
        status = nf_put_att_text(ncid,varid,"units",4,"s^-1")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif

      if(output_dissten.eq.1)then
        status = nf_def_var(ncid,"dissten",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",28,"dissipative heating tendency")
        status = nf_put_att_text(ncid,varid,"units",3,"K/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif

      if(output_pblten.eq.1)then
        status = nf_def_var(ncid,"thpten",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",19,"pbl tendency: theta")
        status = nf_put_att_text(ncid,varid,"units",4,"    ")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

        status = nf_def_var(ncid,"qvpten",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",16,"pbl tendency: qv")
        status = nf_put_att_text(ncid,varid,"units",4,"    ")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

        status = nf_def_var(ncid,"qcpten",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",16,"pbl tendency: qc")
        status = nf_put_att_text(ncid,varid,"units",4,"    ")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

        status = nf_def_var(ncid,"qipten",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",16,"pbl tendency: qi")
        status = nf_put_att_text(ncid,varid,"units",4,"    ")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

        status = nf_def_var(ncid,"upten",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",15,"pbl tendency: u")
        status = nf_put_att_text(ncid,varid,"units",4,"    ")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

        status = nf_def_var(ncid,"vpten",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",15,"pbl tendency: v")
        status = nf_put_att_text(ncid,varid,"units",4,"    ")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif

      if(output_radten.eq.1)then
        status = nf_def_var(ncid,"swten",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",25,"pot temp tendency, sw rad")
        status = nf_put_att_text(ncid,varid,"units",4,"    ")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)

        status = nf_def_var(ncid,"lwten",nf_real,4,(/niid,njid,nkid,timeid/),varid)
        status = nf_put_att_text(ncid,varid,"def",25,"pot temp tendency, lw rad")
        status = nf_put_att_text(ncid,varid,"units",3,"K/s")
        status = nf_put_att_real(ncid,varid,"x_min",nf_real,1,0.001*x_mns)
        status = nf_put_att_real(ncid,varid,"y_min",nf_real,1,0.001*y_mns)
        status = nf_put_att_real(ncid,varid,"z_min",nf_real,1,0.001*z_mns)
      endif

      status = nf_put_att_text(ncid,NF_GLOBAL,'Conventions',6,'COARDS')
      status = nf_put_att_real(ncid,NF_GLOBAL,'x_min'  ,nf_real,1,0.001*x_min)
      status = nf_put_att_real(ncid,NF_GLOBAL,'x_max'  ,nf_real,1,0.001*x_max)
      status = nf_put_att_real(ncid,NF_GLOBAL,'x_delta',nf_real,1,0.001*dx)
      status = nf_put_att_text(ncid,NF_GLOBAL,'x_units',2,'km')
      status = nf_put_att_text(ncid,NF_GLOBAL,'x_label',1,'x')
      status = nf_put_att_real(ncid,NF_GLOBAL,'y_min'  ,nf_real,1,0.001*y_min)
      status = nf_put_att_real(ncid,NF_GLOBAL,'y_max'  ,nf_real,1,0.001*y_max)
      status = nf_put_att_real(ncid,NF_GLOBAL,'y_delta',nf_real,1,0.001*dy)
      status = nf_put_att_text(ncid,NF_GLOBAL,'y_units',2,'km')
      status = nf_put_att_text(ncid,NF_GLOBAL,'y_label',1,'y')
      status = nf_put_att_real(ncid,NF_GLOBAL,'z_min'  ,nf_real,1,0.001*z_min)
      status = nf_put_att_real(ncid,NF_GLOBAL,'z_max'  ,nf_real,1,0.001*z_max)
      status = nf_put_att_real(ncid,NF_GLOBAL,'z_delta',nf_real,1,0.001*dz)
      status = nf_put_att_text(ncid,NF_GLOBAL,'z_units',2,'km')
      status = nf_put_att_text(ncid,NF_GLOBAL,'z_label',1,'z')

      status = nf_enddef(ncid)

! ... end of defs
!--------------------------------------------------
! begin data ... initial time ...

    IF(icor.eq.1)THEN
      status = nf_inq_varid(ncid,'f_cor',varid)
      status = nf_put_var_real(ncid,varid,fcor)
    ENDIF

    if(.not.coards)then
      status = nf_inq_varid(ncid,'ztop',varid)
      status = nf_put_var_real(ncid,varid,0.001*ztop)
    endif

      if(coards)then
        status = nf_inq_varid(ncid,'ni',varid)
      else
        status = nf_inq_varid(ncid,'xh',varid)
      endif
      do i=1,ni
        status = nf_put_var1_real(ncid,varid,i,0.001*xh(i))
      enddo

      if(coards)then
        status = nf_inq_varid(ncid,'nip1',varid)
      else
        status = nf_inq_varid(ncid,'xf',varid)
      endif
      do i=1,ni+1
        status = nf_put_var1_real(ncid,varid,i,0.001*xf(i))
      enddo

      if(coards)then
        status = nf_inq_varid(ncid,'nj',varid)
      else
        status = nf_inq_varid(ncid,'yh',varid)
      endif
      do j=1,nj
        status = nf_put_var1_real(ncid,varid,j,0.001*yh(j))
      enddo

      if(coards)then
        status = nf_inq_varid(ncid,'njp1',varid)
      else
        status = nf_inq_varid(ncid,'yf',varid)
      endif
      do j=1,nj+1
        status = nf_put_var1_real(ncid,varid,j,0.001*yf(j))
      enddo

      if(coards)then
        status = nf_inq_varid(ncid,'nk',varid)
      else
        status = nf_inq_varid(ncid,'z',varid)
      endif
      if(terrain_flag)then
        do k=1,nk
          status = nf_put_var1_real(ncid,varid,k,0.001*sigma(k))
        enddo
      else
        do k=1,nk
          status = nf_put_var1_real(ncid,varid,k,0.001*zh(1,1,k))
        enddo
      endif

      if(coards)then
        status = nf_inq_varid(ncid,'nkp1',varid)
      else
        status = nf_inq_varid(ncid,'zf',varid)
      endif
      if(terrain_flag)then
        do k=1,nk+1
          status = nf_put_var1_real(ncid,varid,k,0.001*sigmaf(k))
        enddo
      else
        do k=1,nk+1
          status = nf_put_var1_real(ncid,varid,k,0.001*zf(1,1,k))
        enddo
      endif

!----------------------------------------------------------

      if (output_zs.eq.1.and.terrain_flag) then
        do j=1,nj
        do i=1,ni
          d2d(i,j)=zs(i,j)
        enddo
        enddo
        status = nf_inq_varid(ncid,'zs',varid)
        call checkstatus(status)
        status = nf_put_var_real(ncid,varid,d2d)
        call checkstatus(status)
      endif
      if (output_zh.eq.1.and.terrain_flag) then
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=zh(i,j,k)
        enddo
        enddo
        enddo
        status = nf_inq_varid(ncid,'zh',varid)
        call checkstatus(status)
        status = nf_put_var_real(ncid,varid,ds)
        call checkstatus(status)
      endif

      if (output_basestate.eq.1) then
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=pi0(i,j,k)
        enddo
        enddo
        enddo
        status = nf_inq_varid(ncid,'pi0',varid)
        call checkstatus(status)
        status = nf_put_var_real(ncid,varid,ds)
        call checkstatus(status)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=th0(i,j,k)
        enddo
        enddo
        enddo
        status = nf_inq_varid(ncid,'th0',varid)
        call checkstatus(status)
        status = nf_put_var_real(ncid,varid,ds)
        call checkstatus(status)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=prs0(i,j,k)
        enddo
        enddo
        enddo
        status = nf_inq_varid(ncid,'prs0',varid)
        call checkstatus(status)
        status = nf_put_var_real(ncid,varid,ds)
        call checkstatus(status)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=qv0(i,j,k)
        enddo
        enddo
        enddo
        status = nf_inq_varid(ncid,'qv0',varid)
        call checkstatus(status)
        status = nf_put_var_real(ncid,varid,ds)
        call checkstatus(status)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          du(i,j,k)=u0(i,j,k)
        enddo
        enddo
        enddo
        status = nf_inq_varid(ncid,'u0',varid)
        call checkstatus(status)
        status = nf_put_var_real(ncid,varid,du)
        call checkstatus(status)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          dv(i,j,k)=v0(i,j,k)
        enddo
        enddo
        enddo
        status = nf_inq_varid(ncid,'v0',varid)
        call checkstatus(status)
        status = nf_put_var_real(ncid,varid,dv)
        call checkstatus(status)
      endif

      ! ... end if info at initial time only

    ELSE

      ! Not initial time ... open file ...

      call disp_err( nf_open(cdf_out,nf_write,ncid), .true. )

    ENDIF ifallinfo

!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
! Main data dump:

      status = nf_inq_varid(ncid,'time',varid)
      call checkstatus(status)
      status = nf_put_var1_real(ncid,varid,time_index,actual_time)
      call checkstatus(status)

        if(output_rain.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            d2d(i,j)=rain(i,j,1)
          enddo
          enddo
          call write2d_nc('rain    ',ncid,time_index,ni,nj,d2d)
        endif
        if(output_sws .eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            d2d(i,j)=sws(i,j,1)
          enddo
          enddo
          call write2d_nc('sws     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            d2d(i,j)=svs(i,j,1)
          enddo
          enddo
          call write2d_nc('svs     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            d2d(i,j)=sps(i,j,1)
          enddo
          enddo
          call write2d_nc('sps     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            d2d(i,j)=srs(i,j,1)
          enddo
          enddo
          call write2d_nc('srs     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            d2d(i,j)=sgs(i,j,1)
          enddo
          enddo
          call write2d_nc('sgs     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            d2d(i,j)=sus(i,j,1)
          enddo
          enddo
          call write2d_nc('sus     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            d2d(i,j)=shs(i,j,1)
          enddo
          enddo
          call write2d_nc('shs     ',ncid,time_index,ni,nj,d2d)
        endif
        IF(nrain.eq.2)THEN
          if(output_rain.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              d2d(i,j)=rain(i,j,2)
            enddo
            enddo
            call write2d_nc('rain2   ',ncid,time_index,ni,nj,d2d)
          endif
          if(output_sws .eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              d2d(i,j)=sws(i,j,2)
            enddo
            enddo
            call write2d_nc('sws2    ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              d2d(i,j)=svs(i,j,2)
            enddo
            enddo
            call write2d_nc('svs2    ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              d2d(i,j)=sps(i,j,2)
            enddo
            enddo
            call write2d_nc('sps2    ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              d2d(i,j)=srs(i,j,2)
            enddo
            enddo
            call write2d_nc('srs2    ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              d2d(i,j)=sgs(i,j,2)
            enddo
            enddo
            call write2d_nc('sgs2    ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              d2d(i,j)=sus(i,j,2)
            enddo
            enddo
            call write2d_nc('sus2    ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              d2d(i,j)=shs(i,j,2)
            enddo
            enddo
            call write2d_nc('shs2    ',ncid,time_index,ni,nj,d2d)
          endif
        ENDIF
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
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=dum1(i,j,1)
        enddo
        enddo
        call write2d_nc('uh      ',ncid,time_index,ni,nj,d2d)
      endif
      if(output_dbz.eq.1)then
        call calccref(dum1(ib,jb,1),dbz)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=dum1(i,j,1)
        enddo
        enddo
        call write2d_nc('cref    ',ncid,time_index,ni,nj,d2d)
      endif
      if(output_coldpool.eq.1)then
        call calccpch(zf,th0,qv0,dum1(ib,jb,1),dum1(ib,jb,2),tha,qa)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=dum1(i,j,1)
        enddo
        enddo
        call write2d_nc('cpc     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=dum1(i,j,2)
        enddo
        enddo
        call write2d_nc('cph     ',ncid,time_index,ni,nj,d2d)
      endif
      if(output_sfcflx.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=thflux(i,j)
        enddo
        enddo
        call write2d_nc('thflux  ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=qvflux(i,j)
        enddo
        enddo
        call write2d_nc('qvflux  ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=0.25*( (cdu(i,j)+cdu(i+1,j))   &
                         +(cdv(i,j)+cdv(i,j+1)) )
        enddo
        enddo
        call write2d_nc('cd      ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=ce(i,j)
        enddo
        enddo
        call write2d_nc('ce      ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=tsk(i,j)
        enddo
        enddo
        call write2d_nc('tsk     ',ncid,time_index,ni,nj,d2d)
      endif

      if(output_sfcparams.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=xland(i,j)
        enddo
        enddo
        call write2d_nc('xland   ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=lu_index(i,j)
        enddo
        enddo
        call write2d_nc('lu      ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=mavail(i,j)
        enddo
        enddo
        call write2d_nc('mavail  ',ncid,time_index,ni,nj,d2d)
      endif

      if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.oceanmodel.eq.2))then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=tmn(i,j)
        enddo
        enddo
        call write2d_nc('tmn     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=hfx(i,j)
        enddo
        enddo
        call write2d_nc('hfx     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=qfx(i,j)
        enddo
        enddo
        call write2d_nc('qfx     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=gsw(i,j)
        enddo
        enddo
        call write2d_nc('gsw     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=glw(i,j)
        enddo
        enddo
        call write2d_nc('glw     ',ncid,time_index,ni,nj,d2d)
      endif

      if((output_sfcparams.eq.1).and.(sfcmodel.eq.2))then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=tslb(i,j,1)
        enddo
        enddo
        call write2d_nc('tslb1   ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=tslb(i,j,2)
        enddo
        enddo
        call write2d_nc('tslb2   ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=tslb(i,j,3)
        enddo
        enddo
        call write2d_nc('tslb3   ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=tslb(i,j,4)
        enddo
        enddo
        call write2d_nc('tslb4   ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=tslb(i,j,5)
        enddo
        enddo
        call write2d_nc('tslb5   ',ncid,time_index,ni,nj,d2d)
      endif

      if(output_sfcparams.eq.1.and.oceanmodel.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=tml(i,j)
        enddo
        enddo
        call write2d_nc('tml     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=hml(i,j)
        enddo
        enddo
        call write2d_nc('hml     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=huml(i,j)
        enddo
        enddo
        call write2d_nc('huml    ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=hvml(i,j)
        enddo
        enddo
        call write2d_nc('hvml    ',ncid,time_index,ni,nj,d2d)
      endif

      if( output_radten.eq.1 )then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=radsw(i,j)
        enddo
        enddo
        call write2d_nc('radsw   ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=rnflx(i,j)
        enddo
        enddo
        call write2d_nc('rnflx   ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=radswnet(i,j)
        enddo
        enddo
        call write2d_nc('radswnet',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=radlwin(i,j)
        enddo
        enddo
        call write2d_nc('radlwin ',ncid,time_index,ni,nj,d2d)
      endif

      IF(output_sfcdiags.eq.1)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=u10(i,j)
        enddo
        enddo
        call write2d_nc('u10     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=v10(i,j)
        enddo
        enddo
        call write2d_nc('v10     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=t2(i,j)
        enddo
        enddo
        call write2d_nc('t2      ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=q2(i,j)
        enddo
        enddo
        call write2d_nc('q2      ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=znt(i,j)
        enddo
        enddo
        call write2d_nc('znt     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=ust(i,j)
        enddo
        enddo
        call write2d_nc('ust     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=hpbl(i,j)
        enddo
        enddo
        call write2d_nc('hpbl    ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=zol(i,j)
        enddo
        enddo
        call write2d_nc('zol     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=mol(i,j)
        enddo
        enddo
        call write2d_nc('mol     ',ncid,time_index,ni,nj,d2d)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          d2d(i,j)=br(i,j)
        enddo
        enddo
        call write2d_nc('br      ',ncid,time_index,ni,nj,d2d)
      ENDIF

!---------------------------------------------

      if(output_u.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          du(i,j,k)=ua(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('u       ',ncid,time_index,ni+1,nj,nk,du)
      endif
      if(output_upert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          du(i,j,k)=ua(i,j,k)-u0(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('upert   ',ncid,time_index,ni+1,nj,nk,du)
      endif
      if(output_v.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          dv(i,j,k)=va(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('v       ',ncid,time_index,ni,nj+1,nk,dv)
      endif
      if(output_vpert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          dv(i,j,k)=va(i,j,k)-v0(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('vpert   ',ncid,time_index,ni,nj+1,nk,dv)
      endif
      if(output_w.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dw(i,j,k)=wa(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('w       ',ncid,time_index,ni,nj,nk+1,dw)
      endif

      if(output_rho.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=rho(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('rho     ',ncid,time_index,ni,nj,nk,ds)
      endif
      if(output_rhopert.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=rho(i,j,k)-rho0(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('rhopert ',ncid,time_index,ni,nj,nk,ds)
      endif
      if(output_th.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=tha(i,j,k)+th0(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('th      ',ncid,time_index,ni,nj,nk,ds)
      endif
      if(output_thpert.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=tha(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('thpert  ',ncid,time_index,ni,nj,nk,ds)
      endif
      if(output_pi.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=pi0(i,j,k)+ppi(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('pi      ',ncid,time_index,ni,nj,nk,ds)
      endif
      if(output_pipert.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=ppi(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('pipert  ',ncid,time_index,ni,nj,nk,ds)
      endif
      if(output_prs.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=prs(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('prs     ',ncid,time_index,ni,nj,nk,ds)
      endif
      if(output_prspert.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=prs(i,j,k)-p00*(pi0(i,j,k)**cpdrd)
        enddo
        enddo
        enddo
        call write3d_nc('prspert ',ncid,time_index,ni,nj,nk,ds)
      endif
    IF(imoist.eq.1)THEN
      if(output_qv.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=qa(i,j,k,nqv)
        enddo
        enddo
        enddo
        call write3d_nc('qv      ',ncid,time_index,ni,nj,nk,ds)
      endif
      if(output_qvpert.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=qa(i,j,k,nqv)-qv0(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('qvpert  ',ncid,time_index,ni,nj,nk,ds)
      endif
      if (output_q.eq.1) then
        chid = '        '
        do n=1,numq
          if(n.ne.nqv)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              ds(i,j,k)=qa(i,j,k,n)
            enddo
            enddo
            enddo
            write(chid(1:3),110) qname(n)
110         format(a3)
            call write3d_nc(chid,ncid,time_index,ni,nj,nk,ds)
          endif
        enddo
      endif
      if(output_dbz.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=dbz(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('dbz     ',ncid,time_index,ni,nj,nk,ds)
      endif
    ENDIF  ! endif for imoist=1
      IF(iptra.eq.1)THEN
        do n=1,npt
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=pta(i,j,k,n)
          enddo
          enddo
          enddo
          chid = 'pt      '
          write(chid(3:4),112) n
112       format(i2.2)
          call write3d_nc(chid,ncid,time_index,ni,nj,nk,ds)
        enddo
      ENDIF

    IF(iturb.ge.1)THEN
      if(output_km.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dw(i,j,k)=kmh(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('kmh     ',ncid,time_index,ni,nj,nk+1,dw)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dw(i,j,k)=kmv(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('kmv     ',ncid,time_index,ni,nj,nk+1,dw)
      endif
      if(output_kh.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dw(i,j,k)=khh(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('khh     ',ncid,time_index,ni,nj,nk+1,dw)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dw(i,j,k)=khv(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('khv     ',ncid,time_index,ni,nj,nk+1,dw)
      endif
    ENDIF

      IF((iturb.eq.1).and.(output_tke.eq.1))then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dw(i,j,k)=tkea(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('tke     ',ncid,time_index,ni,nj,nk+1,dw)
      ENDIF

      if(output_uinterp.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=0.5*(ua(i,j,k)+ua(i+1,j,k))
        enddo
        enddo
        enddo
        call write3d_nc('uinterp ',ncid,time_index,ni,nj,nk,ds)
      endif
      if(output_vinterp.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=0.5*(va(i,j,k)+va(i,j+1,k))
        enddo
        enddo
        enddo
        call write3d_nc('vinterp ',ncid,time_index,ni,nj,nk,ds)
      endif
      if(output_winterp.eq.1) then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=0.5*(wa(i,j,k)+wa(i,j,k+1))
        enddo
        enddo
        enddo
        call write3d_nc('winterp ',ncid,time_index,ni,nj,nk,ds)
      endif

      if(output_vort.eq.1)then
        call calcvort(xh,xf,uf,vf,zh,mf,zf,ua,va,wa,dum1,dum2,dum3,dum4)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=dum1(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('xvort   ',ncid,time_index,ni,nj,nk,ds)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=dum2(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('yvort   ',ncid,time_index,ni,nj,nk,ds)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=dum3(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('zvort   ',ncid,time_index,ni,nj,nk,ds)
      endif

      if(output_dissten.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=dissten(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('dissten ',ncid,time_index,ni,nj,nk,ds)
      endif

      if(output_pblten.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=thpten(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('thpten  ',ncid,time_index,ni,nj,nk,ds)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=qvpten(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('qvpten  ',ncid,time_index,ni,nj,nk,ds)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=qcpten(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('qcpten  ',ncid,time_index,ni,nj,nk,ds)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=qipten(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('qipten  ',ncid,time_index,ni,nj,nk,ds)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=upten(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('upten   ',ncid,time_index,ni,nj,nk,ds)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=vpten(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('vpten   ',ncid,time_index,ni,nj,nk,ds)
      endif

      if( output_radten.eq.1 )then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=swten(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('swten   ',ncid,time_index,ni,nj,nk,ds)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ds(i,j,k)=lwten(i,j,k)
        enddo
        enddo
        enddo
        call write3d_nc('lwten   ',ncid,time_index,ni,nj,nk,ds)
      endif

!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
!  All done ... close file ...

      call disp_err( nf_close(ncid) , .true. )



      write(outfile,*)
      write(outfile,*) 'Done Writing Data to File: nwrite=',nwrite
      write(outfile,*)

      return
      end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine checkstatus(status)
      implicit none

      integer :: status

      include 'netcdf.inc'

      if(status.ne.nf_noerr)then
        print *,'  Error ... '
        print *,nf_strerror(status)
        call stopcm1
      endif

      return
      end subroutine checkstatus

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine write2d_nc(chid,ncid,time_index,ni,nj,d2d)
      implicit none

      include 'netcdf.inc'

      character*8, intent(in) :: chid
      integer, intent(in) :: ncid,time_index,ni,nj
      real, dimension(ni,nj), intent(inout) :: d2d

      integer :: varid,status

!----------------------------------

      status = nf_inq_varid(ncid,chid,varid)
      if(status.ne.nf_noerr)then
        print *,'  Error in write2d_nc, chid = ',chid
        print *,nf_strerror(status)
        call stopcm1
      endif

      status = nf_put_vara_real(ncid,varid,(/1,1,time_index/),(/ni,nj,1/),d2d)
      if(status.ne.nf_noerr)then
        print *,'  Error in write2d_nc, chid = ',chid
        print *,nf_strerror(status)
        call stopcm1
      endif

!----------------------------------

      return
      end subroutine write2d_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write3d_nc(chid,ncid,time_index,ni,nj,nk,ds)
      implicit none

      include 'netcdf.inc'

      character*8, intent(in) :: chid
      integer, intent(in) :: ncid,time_index,ni,nj,nk
      real, dimension(ni,nj), intent(inout) :: ds

      integer :: varid,status

!----------------------------------

      status = nf_inq_varid(ncid,chid,varid)
      if(status.ne.nf_noerr)then
        print *,'  Error in write3d_nc, chid = ',chid
        print *,nf_strerror(status)
        call stopcm1
      endif

      status = nf_put_vara_real(ncid,varid,(/1,1,1,time_index/),(/ni,nj,nk,1/),ds)
      if(status.ne.nf_noerr)then
        print *,'  Error in write3d_nc, chid = ',chid
        print *,nf_strerror(status)
        call stopcm1
      endif

!----------------------------------

      return
      end subroutine write3d_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine writestat_nc(nrec,rtime,nstat,rstat,qname,budname)
      implicit none

      include 'input.incl'

      integer, intent(inout) :: nrec
      real,    intent(in)    :: rtime
      integer, intent(in)    :: nstat
      real, dimension(stat_out), intent(in) :: rstat
      character*3, dimension(maxq), intent(in) :: qname
      character*6, dimension(maxq), intent(in) :: budname

      include 'netcdf.inc'
      integer :: n,ncid,status,dimid,varid,time_index
      character*8  :: text1
      character*30 :: text2

  IF(nrec.eq.1)THEN
    ! Definitions/descriptions:

    call disp_err( nf_create('cm1out_stats.nc',nf_write,ncid), .true. )

    status = nf_def_dim(ncid,"xh",1,dimid)
    status = nf_def_dim(ncid,"yh",1,dimid)
    status = nf_def_dim(ncid,"zh",1,dimid)
    status = nf_def_dim(ncid,"time",nf_unlimited,dimid)

    status = nf_def_var(ncid,"xh",nf_real,1,(/1/),varid)
    status = nf_put_att_text(ncid,varid,"def",18,"west-east location")
    status = nf_put_att_text(ncid,varid,"units",11,"degree_east")

    status = nf_def_var(ncid,"yh",nf_real,1,(/2/),varid)
    status = nf_put_att_text(ncid,varid,"def",20,"south-north location")
    status = nf_put_att_text(ncid,varid,"units",12,"degree_north")

    status = nf_def_var(ncid,"zh",nf_real,1,(/3/),varid)
    status = nf_put_att_text(ncid,varid,"def",6,"height")
    status = nf_put_att_text(ncid,varid,"units",1,"m")

    status = nf_def_var(ncid,"time",nf_real,1,(/4/),varid)
    status = nf_put_att_text(ncid,varid,"def",4,"time")
    status = nf_put_att_text(ncid,varid,"units",7,"seconds")

    IF(adapt_dt.eq.1)THEN
      status = nf_def_var(ncid,"dt",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"average timestep dt           ")
      status = nf_put_att_text(ncid,varid,"units",7,"seconds")
    ENDIF
    IF(stat_w.eq.1)THEN
      status = nf_def_var(ncid,"wmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"maximum vertical velocity     ")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")
      status = nf_def_var(ncid,"wmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"minimum vertical velocity     ")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")
    ENDIF
    IF(stat_u.eq.1)THEN
      status = nf_def_var(ncid,"umax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max E-W velocity              ")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")
      status = nf_def_var(ncid,"umin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min E-W velocity              ")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")
      status = nf_def_var(ncid,"sumax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max E-W velocity at lowest lvl")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")
      status = nf_def_var(ncid,"sumin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min E-W velocity at lowest lvl")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")
    ENDIF
    IF(stat_v.eq.1)THEN
      status = nf_def_var(ncid,"vmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max N-S velocity              ")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")
      status = nf_def_var(ncid,"vmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min N-S velocity              ")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")
      status = nf_def_var(ncid,"svmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max N-S velocity at lowest lvl")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")
      status = nf_def_var(ncid,"svmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min N-S velocity at lowest lvl")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")
    ENDIF
    IF(stat_rmw.eq.1)THEN
      status = nf_def_var(ncid,"rmw",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"radius of maximum V           ")
      status = nf_put_att_text(ncid,varid,"units",1,"m")
    ENDIF
    IF(stat_pipert.eq.1)THEN
      status = nf_def_var(ncid,"ppimax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max pi pert.                  ")
      status = nf_put_att_text(ncid,varid,"units",14,"nondimensional")
      status = nf_def_var(ncid,"ppimin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min pi pert.                  ")
      status = nf_put_att_text(ncid,varid,"units",14,"nondimensional")
    ENDIF
    IF(stat_prspert.eq.1)THEN
      status = nf_def_var(ncid,"ppmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max pressure pert.            ")
      status = nf_put_att_text(ncid,varid,"units",2,"Pa")
      status = nf_def_var(ncid,"ppmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min pressure pert.            ")
      status = nf_put_att_text(ncid,varid,"units",2,"Pa")
    ENDIF
    IF(stat_thpert.eq.1)THEN
      status = nf_def_var(ncid,"thpmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max potential temp. pert.     ")
      status = nf_put_att_text(ncid,varid,"units",1,"K")
      status = nf_def_var(ncid,"thpmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min potential temp. pert.     ")
      status = nf_put_att_text(ncid,varid,"units",1,"K")
      status = nf_def_var(ncid,"sthpmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max pot temp pert lowest level")
      status = nf_put_att_text(ncid,varid,"units",1,"K")
      status = nf_def_var(ncid,"sthpmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min pot temp pert lowest level")
      status = nf_put_att_text(ncid,varid,"units",1,"K")
    ENDIF
    IF(stat_q.eq.1)THEN
      do n=1,numq
        text1='max     '
        text2='max                           '
        write(text1(4:6),156) qname(n)
        write(text2(5:7),156) qname(n)
156     format(a3)
        status = nf_def_var(ncid,text1,nf_real,1,4,varid)
        status = nf_put_att_text(ncid,varid,"def",30,text2)
        status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")
        text1='min     '
        text2='min                           '
        write(text1(4:6),156) qname(n)
        write(text2(5:7),156) qname(n)
        status = nf_def_var(ncid,text1,nf_real,1,4,varid)
        status = nf_put_att_text(ncid,varid,"def",30,text2)
        status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")
      enddo
    ENDIF
    IF(stat_tke.eq.1)THEN
      status = nf_def_var(ncid,"tkemax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max tke                       ")
      status = nf_put_att_text(ncid,varid,"units",7,"m^2/s^2")
      status = nf_def_var(ncid,"tkemin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min tke                       ")
      status = nf_put_att_text(ncid,varid,"units",7,"m^2/s^2")
    ENDIF
    IF(stat_km.eq.1)THEN
      status = nf_def_var(ncid,"kmhmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max kmh                       ")
      status = nf_put_att_text(ncid,varid,"units",5,"m^2/s")
      status = nf_def_var(ncid,"kmhmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min kmh                       ")
      status = nf_put_att_text(ncid,varid,"units",5,"m^2/s")
      status = nf_def_var(ncid,"kmvmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max kmv                       ")
      status = nf_put_att_text(ncid,varid,"units",5,"m^2/s")
      status = nf_def_var(ncid,"kmvmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min kmv                       ")
      status = nf_put_att_text(ncid,varid,"units",5,"m^2/s")
    ENDIF
    IF(stat_kh.eq.1)THEN
      status = nf_def_var(ncid,"khhmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max khh                       ")
      status = nf_put_att_text(ncid,varid,"units",5,"m^2/s")
      status = nf_def_var(ncid,"khhmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min khh                       ")
      status = nf_put_att_text(ncid,varid,"units",5,"m^2/s")
      status = nf_def_var(ncid,"khvmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max khv                       ")
      status = nf_put_att_text(ncid,varid,"units",5,"m^2/s")
      status = nf_def_var(ncid,"khvmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min khv                       ")
      status = nf_put_att_text(ncid,varid,"units",5,"m^2/s")
    ENDIF
    IF(stat_div.eq.1)THEN
      status = nf_def_var(ncid,"divmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max 3d divergence             ")
      status = nf_def_var(ncid,"divmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min 3d divergence             ")
    ENDIF
    IF(stat_rh.eq.1)THEN
      status = nf_def_var(ncid,"rhmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max relative humidity         ")
      status = nf_def_var(ncid,"rhmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min relative humidity         ")
    ENDIF
    IF(stat_rhi.eq.1)THEN
      status = nf_def_var(ncid,"rhimax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max relative humidity wrt ice ")
      status = nf_def_var(ncid,"rhimin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min relative humidity wrt ice ")
    ENDIF
    IF(iptra.eq.1)then
      do n=1,npt
        text1='maxpt   '
        text2='max pt                        '
        write(text1(6:6),157) n
        write(text2(7:7),157) n
157     format(i1)
        status = nf_def_var(ncid,text1,nf_real,1,4,varid)
        status = nf_put_att_text(ncid,varid,"def",30,text2)
        text1='minpt   '
        text2='min pt                        '
        write(text1(6:6),157) n
        write(text2(7:7),157) n
        status = nf_def_var(ncid,text1,nf_real,1,4,varid)
        status = nf_put_att_text(ncid,varid,"def",30,text2)
      enddo
    endif
    IF(stat_the.eq.1)THEN
      status = nf_def_var(ncid,"themax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max theta-e below 10 km       ")
      status = nf_put_att_text(ncid,varid,"units",1,"K")
      status = nf_def_var(ncid,"themin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min theta-e below 10 km       ")
      status = nf_put_att_text(ncid,varid,"units",1,"K")
      status = nf_def_var(ncid,"sthemax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max theta-e at lowest level   ")
      status = nf_put_att_text(ncid,varid,"units",1,"K")
      status = nf_def_var(ncid,"sthemin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min theta-e at lowest level   ")
      status = nf_put_att_text(ncid,varid,"units",1,"K")
    ENDIF
    IF(stat_cloud.eq.1)THEN
      status = nf_def_var(ncid,"qctop",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max cloud top height          ")
      status = nf_put_att_text(ncid,varid,"units",1,"m")
      status = nf_def_var(ncid,"qcbot",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min cloud base height         ")
      status = nf_put_att_text(ncid,varid,"units",1,"m")
    ENDIF
    IF(stat_sfcprs.eq.1)THEN
      status = nf_def_var(ncid,"sprsmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max pressure at lowest level  ")
      status = nf_put_att_text(ncid,varid,"units",2,"Pa")
      status = nf_def_var(ncid,"sprsmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min pressure at lowest level  ")
      status = nf_put_att_text(ncid,varid,"units",2,"Pa")
    ENDIF
    IF(stat_wsp.eq.1)THEN
      status = nf_def_var(ncid,"wspmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max wind speed                ")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")

      status = nf_def_var(ncid,"wspmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min wind speed                ")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")

      status = nf_def_var(ncid,"swspmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max wind speed at lowest level")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")

      status = nf_def_var(ncid,"swspmin",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min wind speed at lowest level")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")

    IF(idrag.eq.1)THEN
      status = nf_def_var(ncid,"wsp10max",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max 10 m wind speed           ")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")

      status = nf_def_var(ncid,"wsp10min",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"min 10 m wind speed           ")
      status = nf_put_att_text(ncid,varid,"units",3,"m/s")
    ENDIF
    ENDIF
    IF(stat_cfl.eq.1)THEN
    IF(adapt_dt.eq.1)THEN
      status = nf_def_var(ncid,"cflmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max Courant number (average)  ")
      status = nf_put_att_text(ncid,varid,"units",14,"nondimensional")
    ELSE
      status = nf_def_var(ncid,"cflmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max Courant number            ")
      status = nf_put_att_text(ncid,varid,"units",14,"nondimensional")
    ENDIF
      status = nf_def_var(ncid,"kshmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max horiz K stability factor  ")
      status = nf_put_att_text(ncid,varid,"units",14,"nondimensional")
      status = nf_def_var(ncid,"ksvmax",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max vert K stability factor   ")
      status = nf_put_att_text(ncid,varid,"units",14,"nondimensional")
    ENDIF
    IF(stat_vort.eq.1)THEN
      status = nf_def_var(ncid,"vortsfc",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max vert. vort. at lowest lvl ")
      status = nf_put_att_text(ncid,varid,"units",3,"1/s")
      status = nf_def_var(ncid,"vort1km",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max vert. vort. at z = 1 km   ")
      status = nf_put_att_text(ncid,varid,"units",3,"1/s")
      status = nf_def_var(ncid,"vort2km",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max vert. vort. at z = 2 km   ")
      status = nf_put_att_text(ncid,varid,"units",3,"1/s")
      status = nf_def_var(ncid,"vort3km",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max vert. vort. at z = 3 km   ")
      status = nf_put_att_text(ncid,varid,"units",3,"1/s")
      status = nf_def_var(ncid,"vort4km",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max vert. vort. at z = 4 km   ")
      status = nf_put_att_text(ncid,varid,"units",3,"1/s")
      status = nf_def_var(ncid,"vort5km",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"max vert. vort. at z = 5 km   ")
      status = nf_put_att_text(ncid,varid,"units",3,"1/s")
    ENDIF
    IF(stat_tmass.eq.1)THEN
      status = nf_def_var(ncid,"tmass",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"total mass of (dry) air       ")
      status = nf_put_att_text(ncid,varid,"units",2,"kg")
    ENDIF
    IF(stat_tmois.eq.1)THEN
      status = nf_def_var(ncid,"tmois",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"total moisture                ")
      status = nf_put_att_text(ncid,varid,"units",2,"kg")
    ENDIF
    IF(stat_qmass.eq.1)THEN
      do n=1,numq
        IF( (n.eq.nqv) .or.                                 &
            (n.ge.nql1.and.n.le.nql2) .or.                  &
            (n.ge.nqs1.and.n.le.nqs2.and.iice.eq.1) )THEN
          text1='mass    '
          text2='total mass of                 '
          write(text1( 5: 7),156) qname(n)
          write(text2(15:17),156) qname(n)
          status = nf_def_var(ncid,text1,nf_real,1,4,varid)
          status = nf_put_att_text(ncid,varid,"def",30,text2)
          status = nf_put_att_text(ncid,varid,"units",2,"kg")
        ENDIF
      enddo
    ENDIF
    IF(stat_tenerg.eq.1)THEN
      status = nf_def_var(ncid,"ek",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"total kinetic energy          ")
      status = nf_def_var(ncid,"ei",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"total internal energy         ")
      status = nf_def_var(ncid,"ep",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"total potential energy        ")
      status = nf_def_var(ncid,"le",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"total latent energy           ")
      status = nf_def_var(ncid,"et",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"total energy                  ")
    ENDIF
    IF(stat_mo.eq.1)THEN
      status = nf_def_var(ncid,"tmu",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"total E-W momentum            ")
      status = nf_def_var(ncid,"tmv",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"total N-S momentum            ")
      status = nf_def_var(ncid,"tmw",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"total vertical momentum       ")
    ENDIF
    IF(stat_tmf.eq.1)THEN
      status = nf_def_var(ncid,"tmfu",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"total upward mass flux        ")
      status = nf_def_var(ncid,"tmfd",nf_real,1,4,varid)
      status = nf_put_att_text(ncid,varid,"def",30,"total downward mass flux      ")
    ENDIF
    IF(stat_pcn.eq.1)THEN
      do n=1,nbudget
        text1='        '
        text2='                              '
        write(text1(1:6),158) budname(n)
        write(text2(1:6),158) budname(n)
158     format(a6)
        status = nf_def_var(ncid,text1,nf_real,1,4,varid)
        status = nf_put_att_text(ncid,varid,"def",30,text2)
      enddo
    ENDIF
    IF(stat_qsrc.eq.1)THEN
      do n=1,numq
        text1='as      '
        text2='artificial source of          '
        write(text1( 3: 5),156) qname(n)
        write(text2(22:24),156) qname(n)
        status = nf_def_var(ncid,text1,nf_real,1,4,varid)
        status = nf_put_att_text(ncid,varid,"def",30,text2)
      enddo
      do n=1,numq
        text1='bs      '
        text2='bndry source/sink of          '
        write(text1( 3: 5),156) qname(n)
        write(text2(22:24),156) qname(n)
        status = nf_def_var(ncid,text1,nf_real,1,4,varid)
        status = nf_put_att_text(ncid,varid,"def",30,text2)
      enddo
    ENDIF

    status = nf_put_att_text(ncid, NF_GLOBAL, 'Conventions', 6, 'COARDS')

    status = nf_enddef(ncid)

    status = nf_put_var_real(ncid,1,0.0)
    status = nf_put_var_real(ncid,2,0.0)
    status = nf_put_var_real(ncid,3,0.0)

  ELSE

    ! open file:

    call disp_err( nf_open('cm1out_stats.nc',nf_write,ncid), .true. )

  ENDIF

    ! Write data:

    time_index = nrec

    status = nf_put_var1_real(ncid,4,time_index,rtime)

    DO n=1,nstat
      varid = 4 + n
      status = nf_put_var1_real(ncid,varid,time_index,rstat(n))
    ENDDO

    ! close file

    call disp_err( nf_close(ncid) , .true. )

    nrec = nrec + 1

    ! all done

      return
      end subroutine writestat_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine disp_err( status , stop_on_error )
      implicit none

      integer, intent(in) :: status
      logical, intent(in) :: stop_on_error

      include 'netcdf.inc'

      IF( status.ne.nf_noerr )THEN
        IF( stop_on_error )THEN
          print *,'  netcdf status returned an error: ', status,' ... stopping program'
          call stopcm1
        ENDIF
      ENDIF

      return
      end subroutine disp_err

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine writepdata_nc(prec,rtime,pdata)
      implicit none

      include 'input.incl'

      integer :: prec
      real :: rtime
      real, dimension(npvals,nparcels) :: pdata

      include 'netcdf.inc'

      integer :: ncid,status,dimid,varid,time_index,n,np

!-----------------------------------------------------------------------

  IF(prec.eq.1)THEN
    ! Definitions/descriptions:

    call disp_err( nf_create('cm1out_pdata.nc',nf_write,ncid), .true. )

    status = nf_def_dim(ncid,"xh",nparcels,dimid)
    status = nf_def_dim(ncid,"yh",1,dimid)
    status = nf_def_dim(ncid,"zh",1,dimid)
    status = nf_def_dim(ncid,"time",nf_unlimited,dimid)

    status = nf_def_var(ncid,"xh",nf_real,1,(/1/),varid)
    status = nf_put_att_text(ncid,varid,"def",51,"west-east location ... actually, really parcel info")
    status = nf_put_att_text(ncid,varid,"units",11,"degree_east")

    status = nf_def_var(ncid,"yh",nf_real,1,(/2/),varid)
    status = nf_put_att_text(ncid,varid,"def",20,"south-north location")
    status = nf_put_att_text(ncid,varid,"units",12,"degree_north")

    status = nf_def_var(ncid,"zh",nf_real,1,(/3/),varid)
    status = nf_put_att_text(ncid,varid,"def",6,"height")
    status = nf_put_att_text(ncid,varid,"units",1,"m")

    status = nf_def_var(ncid,"time",nf_real,1,(/4/),varid)
    status = nf_put_att_text(ncid,varid,"def",4,"time")
    status = nf_put_att_text(ncid,varid,"units",7,"seconds")

!------------------------

    status = nf_def_var(ncid,"x",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"x                             ")
    status = nf_put_att_text(ncid,varid,"units",1,"m")

    status = nf_def_var(ncid,"y",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"y                             ")
    status = nf_put_att_text(ncid,varid,"units",1,"m")

    status = nf_def_var(ncid,"z",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"z                             ")
    status = nf_put_att_text(ncid,varid,"units",1,"m")

    status = nf_def_var(ncid,"qv",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"water vapor mixing ratio      ")
    status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")

    status = nf_def_var(ncid,"qc",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"cloud water mixing ratio      ")
    status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")

    status = nf_def_var(ncid,"qr",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"rain water mixing ratio       ")
    status = nf_put_att_text(ncid,varid,"units",5,"kg/kg")

    status = nf_def_var(ncid,"nm",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"squared Brunt-Vaisala frqncy  ")
    status = nf_put_att_text(ncid,varid,"units",3,"s-2")

    status = nf_def_var(ncid,"u",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"u                             ")
    status = nf_put_att_text(ncid,varid,"units",3,"m/s")

    status = nf_def_var(ncid,"v",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"v                             ")
    status = nf_put_att_text(ncid,varid,"units",3,"m/s")

    status = nf_def_var(ncid,"w",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"w                             ")
    status = nf_put_att_text(ncid,varid,"units",3,"m/s")

    status = nf_def_var(ncid,"kh",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"turb. coef. for scalar        ")
    status = nf_put_att_text(ncid,varid,"units",5,"m^2/s")

    status = nf_def_var(ncid,"the",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"theta-e                       ")
    status = nf_put_att_text(ncid,varid,"units",1,"K")

    status = nf_def_var(ncid,"b",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"buoyancy                      ")
    status = nf_put_att_text(ncid,varid,"units",5,"m/s^2")

    status = nf_def_var(ncid,"dpdz",nf_real,2,(/1,4/),varid)
    status = nf_put_att_text(ncid,varid,"def",30,"dpdz tendency                 ")
    status = nf_put_att_text(ncid,varid,"units",5,"m/s^2")

!------------------------

    status = nf_put_att_text(ncid, NF_GLOBAL, 'Conventions', 6, 'COARDS')

    status = nf_enddef(ncid)

  do np=1,nparcels
    status = nf_put_var1_real(ncid,1,np,float(np))
  enddo
    status = nf_put_var_real(ncid,2,0.0)
    status = nf_put_var_real(ncid,3,0.0)

!------------------------

  ELSE

    ! open file:

    call disp_err( nf_open('cm1out_pdata.nc',nf_write,ncid), .true. )

  ENDIF

      ! Write data:

      time_index = prec

      status = nf_put_var1_real(ncid,4,time_index,rtime)

      DO np=1,nparcels
      DO n=1,npvals - 3
        varid = 4 + n
        status = nf_put_var1_real(ncid,varid,(/np,time_index/),pdata(n,np))
      ENDDO
      ENDDO

      ! close file

      call disp_err( nf_close(ncid) , .true. )

      prec = prec + 1

      ! all done

      return
      end subroutine writepdata_nc
