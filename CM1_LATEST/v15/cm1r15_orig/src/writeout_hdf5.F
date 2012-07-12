! HDF5 output written by Leigh Orf <leigh.orf@cmich.edu>
!
! In namelist.input, you can choose which of the 2D and 3D variables are written.
!
! Three output options are available, set in namelist.input:
! output_format = 3 uses scale-offset lossy compression followed by gzip compression. (smallest files)
! output_format = 4 uses gzip compression (level 1) on floating point data and is lossless.
! output_format = 5 uses no compression whatsoever, just floating point (largest files)
!
! These compression options apply to 3D data. See comments in subroutine
! writefloat_hdf_unlimited for how to to add gzip compression to
! unlimited dimension arrays in parcel and stats hdf5 files. For more
! information on compression, please refer to the HDF5 user's guide.
!
!
! The modular nature of HDF5 makes it so that you do not need to worry
! about how the data was compressed when you read it - filters are
! activated automatically and your data appears as 32 bit floating point
! data regardless of the compression choice.
!
! Note that with the scaleoffset filter, you need to a scaling parameter
! which is related to absolute magnitude of your data. The values I have
! chosen work fine for me, but you may wish to tweak them. See hdf5
! docs.
!
! I have written this code to work with MPI only. If you wish to run
! with hdf5 output on a single core, you still need to compile with
! MPI and run with 1 MPI process.
!
! Each MPI rank is matched to one hdf5 file for history data. A
! directory is created which houses all of the hdf5 files at a given
! model dump time. Software has been developed to simplify dealing with
! this format where model domain is split amongst many files. Conversion
! code to create netcdf or vis5d format has also been created and is
! available. If you are interested in exploring additional compression
! options and/or conversion code contact Leigh Orf <leigh.orf@cmich.edu>
! or post to the cm1 CM1 google group.
!
! All hdf5 datasets now contain description and units attributes. New
! code has been added which writes both parcel (trajectory) data to hdf5
! format, as well as statistics data. Both are stored in a single file
! with an unlimited dimension in time, growing as the model progresses.
! This is similar to the netcdf code, which I used as a guide.
!
! Leigh Orf 1/11/11

!--------------------------------------------------------------------------------------------MARK


subroutine writeout_mult_hdf5(model_time,qname,sigma,sigmaf,xh,xf,uf,yh,yf,vf,mh,zh,mf,zf,    &
                pi0,prs0,rho0,th0,qv0,u0,v0,                        &
                zs,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,cdu,cdv,ce,dum1,dum2,dum3,dum4, &
                rho,prs,dbz,ua,dumu,va,dumv,wa,dumw,ppi,tha,        &
                qa,kmh,kmv,khh,khv,tkea,pta,num_soil_layers,   &
                lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,   &
                radsw,rnflx,radswnet,radlwin,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br,  &
                dissten,thpten,qvpten,qcpten,qipten,upten,vpten,swten,lwten)
#ifdef HDFOUT
      use HDF5
#endif
      implicit none

      include 'input.incl'
      include 'constants.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      real :: model_time
      character*3, dimension(maxq) :: qname
      real, dimension(kb:ke) :: sigma
      real, dimension(kb:ke+1) :: sigmaf
      real, dimension(ib:ie)   :: xh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je)   :: yh
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


#ifdef HDFOUT
      real, dimension(ni,nj,nk) :: ds
      real, dimension(ni,nj) :: ds2
      real, dimension(1:nk) :: zhcol,zfcol,u0col,v0col,th0col,pi0col,rho0col,pres0col,qv0col
      real, dimension(1:ni) :: xfcol,xhcol
      real, dimension(1:nj) :: yfcol,yhcol

      integer i,j,k,n
      character(len=120)  description
      character(len=20)   units
      character(len=20)   varname
      character(len=160)  filename
      character(LEN=20) :: qvarname
      character(len=8) chid
      character(len=4) cid
      character(len=5) ctime

      character(len=100) filedir
      character(len=109) dircmd

      integer ierror
      integer(HID_T) :: file_id,dset_id,dspace_id
      integer(HID_T) :: grid_group_id, mesh_group_id, basestate_group_id, twod_group_id, threed_group_id,threed_basestate_group_id
      character(LEN=4),  parameter :: grid_group_name = 'grid'
      character(LEN=4),  parameter :: mesh_group_name = 'mesh'
      character(LEN=9),  parameter :: basestate_group_name = 'basestate'
      character(LEN=12),  parameter :: threed_basestate_group_name = '3d_basestate'
      character(LEN=2),  parameter :: twod_group_name = '2d'
      character(LEN=2),  parameter :: threed_group_name = '3d'
      integer sd_id, sds_id, status, dim_id, data_type
      integer ix,iy
      integer rank
      integer(HSIZE_T), DIMENSION(3) :: dims
      integer itime,lt,flag
      integer x0,x1,y0,y1; ! bounding box for each MPI rank within total domain

      if (len_trim(output_path) .eq. 0) then
            write(outfile,*) 'empty output_path'
            write(outfile,*) 'defaulting to current working directory'
            output_path='.'
      endif
      write(ctime,102)int(model_time)
      write(cid,103)myid
102   format(i5.5)
103   format(i4.4)
      filedir = trim(output_path)//'/'//trim(output_basename)//'.'//ctime//'.cdir/'
      dircmd = 'mkdir -p '// trim(filedir)
      if (myid.eq.0) call system(dircmd)
#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr) ! TODO: create directory structure ahead of time
#endif
      filename = trim(filedir)//trim(output_basename)//'.'//ctime//'_'//cid//'.cm1hdf5'
      write(outfile,*) 'myid,filename=',myid,'   ',trim(filename)

      call h5open_f (ierror)
      call h5fcreate_f(trim(filename),H5F_ACC_TRUNC_F,file_id,ierror)

! Let's create some groups for organizational purposes

      call h5gcreate_f(file_id,grid_group_name,grid_group_id,ierror)      
      call h5gcreate_f(file_id,mesh_group_name,mesh_group_id,ierror)      
      call h5gcreate_f(file_id,basestate_group_name,basestate_group_id,ierror)      
      call h5gcreate_f(file_id,threed_basestate_group_name,threed_basestate_group_id,ierror)      
      call h5gcreate_f(file_id,twod_group_name,twod_group_id,ierror)      
      call h5gcreate_f(file_id,threed_group_name,threed_group_id,ierror)      

      x0 = (myi-1)*ni;
      x1 = myi*ni-1;
      y0 = (myj-1)*nj;
      y1 = myj*nj-1;

      varname='time'
      description='model time'
      units='seconds'
      call h5_write_float(file_id,varname,description,units,model_time)

      varname='dx'
      description='grid spacing in the x direction. If E/W grid stretching, this should be somehwat representative'
      units='meters'
      call h5_write_float(mesh_group_id,varname,description,units,dx)

      varname='dy'
      description='grid spacing in the y direction. If N/S grid stretching, this should be somehwat representative'
      units='meters'
      call h5_write_float(mesh_group_id,varname,description,units,dy)

      varname='x0'
      description='index of westmost point of 3d field in this hdf5 file with respect to full model domain'
      units='index'
      call h5_write_int  (grid_group_id,varname,description,units,x0)

      varname='x1'
      description='index of eastmost point of 3d field in this hdf5 file with respect to full model domain'
      units='index'
      call h5_write_int  (grid_group_id,varname,description,units,x1)

      varname='y0'
      description='index of southmost point of 3d field in this hdf5 file with respect to full model domain'
      units='index'
      call h5_write_int  (grid_group_id,varname,description,units,y0)

      varname='y1'
      description='index of northmost point of 3d field in this hdf5 file with respect to full model domain'
      units='index'
      call h5_write_int  (grid_group_id,varname,description,units,y1)

      varname='myi'
      description='east/west index of hdf5 file in 2d decomposition. (1 < myi <= nodex)'
      units='index'
      call h5_write_int  (grid_group_id,varname,description,units,myi)

      varname='myj'
      description='north/south index of hdf5 file in 2d decomposition. (1 < myj <= nodej)'
      units='index'
      call h5_write_int  (grid_group_id,varname,description,units,myj)

      varname='ni'
      description='number of gridpoints in east/west direction in this hdf5 file for 2d and 3d data'
      units='index'
      call h5_write_int  (grid_group_id,varname,description,units,ni)

      varname='nj'
      description='number of gridpoints in north/south direction in this hdf5 file for 2d and 3d data'
      units='index'
      call h5_write_int  (grid_group_id,varname,description,units,nj)

      varname='nodex'
      description='number of hdf5 files spanning the full east/west extent of the model domain'
      units='index'
      call h5_write_int  (grid_group_id,varname,description,units,nodex)

      varname='nodey'
      description='number of hdf5 files spanning the full north/south extent of the model domain'
      units='index'
      call h5_write_int  (grid_group_id,varname,description,units,nodey)

      varname='nx'
      description='number of gridpoints spanning the east/west direction of full the model domain'
      units='index'
      call h5_write_int  (grid_group_id,varname,description,units,nx)

      varname='ny'
      description='number of gridpoints spanning the north/south direction of full the model domain'
      units='index'
      call h5_write_int  (grid_group_id,varname,description,units,ny)

      varname='nz'
      description='number of gridpoints spanning the vertical extent of the full the model domain'
      units='index'
      call h5_write_int  (grid_group_id,varname,description,units,nz)


! ORF 1/11/11

! Always write out 1D vertical base state column varables. Full 3D
! base state variables are written below if output_basestate.eq.1 for
! you folks who run with terrain. If you have terrain then ignore the
! following 1D basestate variables (i.e., x0col(k) )

! In other words, if you run with terrain and want base state
! variables, set output_basesteate.eq.1 and use vars in group
! /3d_basestate. If you do not run with terrain (flat bottom), you can
! set output_basestate.eq.0 and use vars in /basestate since all base
! state variables will be horizontlaly homogeneous. 1D vars save space
! over 3D vars. Of course nothing is stopping you from writing full 3D
! basestate even if you are not running with terrain!

      do k=1,nk
          zhcol(k)=zh(1,1,k)
          zfcol(k)=zf(1,1,k)
          u0col(k) = u0(1,1,k)
          v0col(k) = v0(1,1,k)
          th0col(k) = th0(1,1,k)
          pi0col(k) = pi0(1,1,k)
          rho0col(k) = rho0(1,1,k)
          pres0col(k) = p00*(pi0(1,1,k))**cpdrd
          qv0col(k) = qv0(1,1,k)
      enddo

      do j = 1,nj
          yfcol(j) = yf(j)
          yhcol(j) = yh(j)
      enddo

      do i = 1,ni
          xfcol(i) = xf(i)
          xhcol(i) = xh(i)
      enddo

      varname='zh'
      description='Height above ground of scalar points (model uses Arakawa C grid)'
      units='meters'
      call h5_write_1d_float(mesh_group_id,varname,description,units,zhcol,nk)

      varname='zf'
      description='Height above ground of W points (model uses Arakawa C grid)'
      units='meters'
      call h5_write_1d_float(mesh_group_id,varname,description,units,zfcol,nk)

      varname='yh'
      description='N/S location of scalar points (model uses Arakawa C grid)'
      units='meters'
      call h5_write_1d_float(mesh_group_id,varname,description,units,yhcol,nj)

      varname='yf'
      description='N/S location of V points (model uses Arakawa C grid)'
      units='meters'
      call h5_write_1d_float(mesh_group_id,varname,description,units,yfcol,nj)

      varname='xh'
      description='E/W location of scalar points (model uses Arakawa C grid)'
      units='meters'
      call h5_write_1d_float(mesh_group_id,varname,description,units,xhcol,ni)

      varname='xf'
      description='E/W location of U points (model uses Arakawa C grid)'
      units='meters'
      call h5_write_1d_float(mesh_group_id,varname,description,units,xfcol,ni)

      varname='u0'
      description='U (E/W wind) base state, (i.e., from sounding)'
      units='m/s'
      call h5_write_1d_float(basestate_group_id,varname,description,units,u0col,nk)

      varname='v0'
      description='V (N/S wind) base state, (i.e., from sounding)'
      units='m/s'
      call h5_write_1d_float(basestate_group_id,varname,description,units,v0col,nk)

      varname='th0'
      description='Potential temperature base state, (i.e., from sounding)'
      units='K'
      call h5_write_1d_float(basestate_group_id,varname,description,units,th0col,nk)

      varname='pi0'
      description='Dimensionless pressure (Exner function) base state, (i.e., from sounding)'
      units='dimensionless'
      call h5_write_1d_float(basestate_group_id,varname,description,units,pi0col,nk)

      varname='pres0'
      description='Pressure base state, (i.e., from sounding)'
      units='Pa'
      call h5_write_1d_float(basestate_group_id,varname,description,units,pres0col,nk)

      varname='rh0'
      description='Density base state, (i.e., from sounding)'
      units='kg/m^3'
      call h5_write_1d_float(basestate_group_id,varname,description,units,rho0col,nk)

      call h5gclose_f(grid_group_id,ierror)
      call h5gclose_f(mesh_group_id,ierror)
      call h5gclose_f(basestate_group_id,ierror)

! terrain

     if (output_zs.eq.1.and.terrain_flag) then
       do j=1,nj
       do i=1,ni
         ds2(i,j)=zs(i,j)
       enddo
       enddo

       varname='zs'
       description='terrain height'
       units='meters'
       call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

      endif

     if (output_zh.eq.1.and.terrain_flag) then
       do k=1,nk
       do j=1,nj
       do i=1,ni
         ds(i,j,k)=zh(i,j,k)
       enddo
       enddo
       enddo

       varname='zh'
       description='height (above ground) of scalar grid points'
       units='meters'
       call h5_write_3d_float(threed_basestate_group_id,varname,description,units,3,ds,ni,nj,nk)

     endif
            
! ORF 1/11/11 Following netcdf example, allow option for output of 3D
! base state. 1D basestate variables in /basestate group are always
! written but are not useful for folks running with terrain who require
! the 3D basestate.

     if (output_basestate.eq.1) then

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=pi0(i,j,k)
          enddo
          enddo
          enddo

          varname='pi0'
          description = '3D base state nondimensional pressure (Exner function)'
          units = 'dimensionless'
          call h5_write_3d_float(threed_basestate_group_id,varname,description,units,3,ds,ni,nj,nk)

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=th0(i,j,k)
          enddo
          enddo
          enddo

          varname='th0'
          description = '3D base state potential temperature'
          units = 'K'
          call h5_write_3d_float(threed_basestate_group_id,varname,description,units,3,ds,ni,nj,nk)

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=prs0(i,j,k)
          enddo
          enddo
          enddo

          varname='p0'
          description = '3D base state pressure'
          units = 'Pa'
          call h5_write_3d_float(threed_basestate_group_id,varname,description,units,3,ds,ni,nj,nk)

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=qv0(i,j,k)
          enddo
          enddo
          enddo

          varname='qv0'
          description = '3D base state water vapor mixing ratio'
          units = 'kg/kg'
          call h5_write_3d_float(threed_basestate_group_id,varname,description,units,3,ds,ni,nj,nk)

          do k=1,nk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k)=u0(i,j,k)
          enddo
          enddo
          enddo

          varname='u0'
          description = '3D base state x component of velocity on native Arakawa C grid'
          units = 'm/s'
          call h5_write_3d_float(threed_basestate_group_id,varname,description,units,3,dumu,ni+1,nj,nk)

          do k=1,nk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k)=v0(i,j,k)
          enddo
          enddo
          enddo

          varname='v0'
          description = '3D base state y component of velocity on native Arakawa C grid'
          units = 'm/s'
          call h5_write_3d_float(threed_basestate_group_id,varname,description,units,3,dumv,ni,nj+1,nk)

      endif !output_basestate.eq.1

      if(output_rain.eq.1) then
          do j=1,nj
          do i=1,ni
            ds2(i,j)=rain(i,j,1)
          enddo
          enddo
          varname='sfcrain'
          description='accumulated surface rainfall'
          units='cm'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)
          if (nrain.eq.2) then
              do j=1,nj
              do i=1,ni
                ds2(i,j)=rain(i,j,2)
              enddo
              enddo
              varname='sfcrain_trans'
              description='accumulated surface rainfall, translated with moving domain'
              units='cm'
              call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)
          endif
      endif

      if(output_sws.eq.1) then
           do j=1,nj
           do i=1,ni
             ds2(i,j)=sws(i,j,1)
           enddo
           enddo
           varname='maxsws'
           description='maximum windspeed at lowest level'
           units='m/s'
           call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

           do j=1,nj
           do i=1,ni
             ds2(i,j)=svs(i,j,1)
           enddo
           enddo
           varname='maxsvs'
           description='maximum vertical vorticity at lowest level'
           units='s^{-1}'
           call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

           do j=1,nj
           do i=1,ni
             ds2(i,j)=sps(i,j,1)
           enddo
           enddo
           varname='maxsps'
           description='minimum pressure at lowest level'
           units='Pa'
           call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

           do j=1,nj
           do i=1,ni
             ds2(i,j)=srs(i,j,1)
           enddo
           enddo
           varname='maxsrs'
           description='maximum surface rainwater mixing ratio'
           units='kg/kg'
           call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

           do j=1,nj
           do i=1,ni
             ds2(i,j)=sgs(i,j,1)
           enddo
           enddo
           varname='maxsgs'
           description='maximum surface graupel/hail mixing ratio'
           units='kg/kg'
           call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

           do j=1,nj
           do i=1,ni
             ds2(i,j)=sus(i,j,1)
           enddo
           enddo
           varname='maxsus'
           description='maximum updraft velocity'
           units='m/s'
           call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

           do j=1,nj
           do i=1,ni
             ds2(i,j)=sus(i,j,1)
           enddo
           enddo
           varname='maxshs'
           description='maximum integrated updraft helicity'
           units='m^2/s^2'
           call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

      endif !output_sws.eq.1

      if (nrain.eq.2) then
         if(output_rain.eq.1) then
            do j=1,nj
            do i=1,ni
              ds2(i,j)=rain(i,j,2)
            enddo
            enddo
            varname='maxrain_trans'
            description='accumulated surface rainfall, translated with moving domain'
            units='cm'
            call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)
          endif

          if(output_sws.eq.1) then
              do j=1,nj
              do i=1,ni
                ds2(i,j)=sws(i,j,2)
              enddo
              enddo
              varname='maxsws_trans'
              description='maximum windspeed at lowest level, translated with moving domain'
              units='m/s'
              call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

              do j=1,nj
              do i=1,ni
                ds2(i,j)=svs(i,j,2)
              enddo
              enddo
              varname='maxsvs_trans'
              description='maximum vertical vorticity at lowest level, translated with moving domain'
              units='s^{-1}'
              call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

              do j=1,nj
              do i=1,ni
                ds2(i,j)=sps(i,j,2)
              enddo
              enddo
              varname='maxsps_trans'
              description='minimum pressure at lowest level, translated with moving domain'
              units='Pa'
              call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

              do j=1,nj
              do i=1,ni
                ds2(i,j)=srs(i,j,2)
              enddo
              enddo
              varname='maxsrs_trans'
              description='maximum surface rainwater mixing ratio, translated with moving domain'
              units='kg/kg'
              call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

              do j=1,nj
              do i=1,ni
                ds2(i,j)=sgs(i,j,2)
              enddo
              enddo
              varname='maxsgs_trans'
              description='maximum surface graupel/hail mixing ratio, translated with moving domain'
              units='kg/kg'
              call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

              do j=1,nj
              do i=1,ni
                ds2(i,j)=sus(i,j,2)
              enddo
              enddo
              varname='maxsus_trans'
              description='maximum updraft velocity, translated with moving domain'
              units='m/s'
              call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

              do j=1,nj
              do i=1,ni
                ds2(i,j)=sus(i,j,2)
              enddo
              enddo
              varname='maxshs_trans'
              description='maximum integrated updraft helicity, translated with moving domain'
              units='m^2/s^2'
              call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)
           endif !output_sws.eq.1
      endif !nrain.eq.2

! Composite reflectivity

      if(output_dbz.eq.1)then
          call calccref(dum1(ib,jb,1),dbz)
          do j=1,nj
          do i=1,ni
            ds2(i,j)=dum1(i,j,1)
          enddo
          enddo
          varname='cref'
          description='composite reflectivity'
          units='dbZ'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)
      endif

! Cold pool C and H (RKW)

      if(output_coldpool.eq.1) then
          call calccpch(zf,th0,qv0,dum1(ib,jb,1),dum1(ib,jb,2),tha,qa)
          do j=1,nj
          do i=1,ni
            ds2(i,j)=dum1(i,j,1)
          enddo
          enddo
          varname='cpc'
          description='cold pool intensity C'
          units='m/s'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=dum1(i,j,2)
          enddo
          enddo
          varname='cph'
          description='cold pool depth h'
          units='meters'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)
      endif !output_coldpool.eq.1

! Surface fluxes

      if(output_sfcflx.eq.1) then
          do j=1,nj
          do i=1,ni
            ds2(i,j)=thflux(i,j)
          enddo
          enddo
          varname='thflux'
          description='surface potential temperature flux'
          units='K m s^{-1}'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)
        
          do j=1,nj
          do i=1,ni
            ds2(i,j)=qvflux(i,j)
          enddo
          enddo

          varname='qvflux'
          description='surface water vapor flux'
          units='kg/kg m s^{-1}'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
            do i=1,ni
             ds2(i,j)=0.25*( (cdu(i,j)+cdu(i+1,j))   &
                           +(cdv(i,j)+cdv(i,j+1)) )
            enddo
          enddo

          varname='cd'
          description='cd' !ORF WTF IS THIS?
          units='nondimensional'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=ce(i,j)
          enddo
          enddo

          varname='ce'
          description='ce' !ORF WTF IS THIS?
          units='nondimensional'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=tsk(i,j)
          enddo
          enddo

          varname='tsk'
          description='soil/ocean temperature'
          units='K'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

      endif !output_sfcflx.eq.1

      if(output_sfcparams.eq.1) then

          do j=1,nj
          do i=1,ni
            ds2(i,j)=xland(i,j)
          enddo
          enddo

          varname='xland'
          description='land/water flag (1=land, 2=water)'
          units='integer flag'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=lu_index(i,j)
          enddo
          enddo

          varname='lu'
          description='land usage index'
          units='integer flag'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=mavail(i,j)
          enddo
          enddo

          varname='mavail'
          description='surface moisture availability'
          units='integer flag'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

      endif !output_sfcparams.eq.1

      if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.oceanmodel.eq.2))then

          do j=1,nj
          do i=1,ni
            ds2(i,j)=tmn(i,j)
          enddo
          enddo
          varname='tmn'
          description='deep-layer soil temperature'
          units='K'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=hfx(i,j)
          enddo
          enddo
          varname='hfx'
          description='heat flux at surface'
          units='W/m^2'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=qfx(i,j)
          enddo
          enddo
          varname='qfx'
          description='surface moisture flux'
          units='W/m^2'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=gsw(i,j)
          enddo
          enddo
          varname='gsw'
          description='downward SW flux at surface'
          units='W/m^2'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=glw(i,j)
          enddo
          enddo
          varname='glw'
          description='downward LW flux at surface'
          units='W/m^2'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)
      endif

      if((output_sfcparams.eq.1).and.(sfcmodel.eq.2))then

          do j=1,nj
          do i=1,ni
            ds2(i,j)=tslb(i,j,1)
          enddo
          enddo
          varname='tslb1'
          description='soil temp, layer 1'
          units='K'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=tslb(i,j,2)
          enddo
          enddo
          varname='tslb2'
          description='soil temp, layer 2'
          units='K'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=tslb(i,j,3)
          enddo
          enddo
          varname='tslb3'
          description='soil temp, layer 3'
          units='K'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=tslb(i,j,4)
          enddo
          enddo
          varname='tslb4'
          description='soil temp, layer 4'
          units='K'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=tslb(i,j,5)
          enddo
          enddo
          varname='tslb5'
          description='soil temp, layer 5'
          units='K'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)
      endif

      if(output_sfcparams.eq.1.and.oceanmodel.eq.2)then

          do j=1,nj
          do i=1,ni
            ds2(i,j)=tml(i,j)
          enddo
          enddo
          varname='tml'
          description='ocean mixed layer temperature'
          units='K'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=hml(i,j)
          enddo
          enddo
          varname='hml'
          description='ocean mixed layer depth'
          units='meters'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=huml(i,j)
          enddo
          enddo
          varname='huml'
          description='ocean mixed layer u velocity'
          units='m/s'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=hvml(i,j)
          enddo
          enddo
          varname='hvml'
          description='ocean mixed layer v velocity'
          units='m/s'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)
      endif

      if(output_radten.eq.1)then

          do j=1,nj
          do i=1,ni
            ds2(i,j)=radsw(i,j)
          enddo
          enddo
          varname='radsw'
          description='solar radiation at surface'
          units='W/m^2'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=rnflx(i,j)
          enddo
          enddo
          varname='rnflx'
          description='net radiation absorbed by surface'
          units='W/m^2'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=radswnet(i,j)
          enddo
          enddo
          varname='radswnet'
          description='net solar radiation'
          units='W/m^2'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=radlwin(i,j)
          enddo
          enddo
          varname='radlwin'
          description='incoming longwave radiation'
          units='W/m^2'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)
      endif

      if(output_sfcdiags.eq.1)then

          do j=1,nj
          do i=1,ni
            ds2(i,j)=u10(i,j)
          enddo
          enddo
          varname='ul0'
          description='diagnostic 10m u wind'
          units='m/s'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=v10(i,j)
          enddo
          enddo
          varname='vl0'
          description='diagnostic 10m v wind'
          units='m/s'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=t2(i,j)
          enddo
          enddo
          varname='t2'
          description='diagnostic 2m temperature'
          units='K'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=q2(i,j)
          enddo
          enddo
          varname='q2'
          description='diagnostic 2m water vapor mixing ratio'
          units='kg/kg'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=znt(i,j)
          enddo
          enddo
          varname='znt'
          description='roughness length'
          units='m'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=ust(i,j)
          enddo
          enddo
          varname='ust'
          description='u* in similarity theory'
          units='m/s'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          if(ipbl.eq.1)then
                description='PBL height (from PBL scheme)'
          else
                description='rough estimate of PBL height'
          endif
          do j=1,nj
          do i=1,ni
            ds2(i,j)=hpbl(i,j)
          enddo
          enddo
          varname='hpbl'
          units='meters'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=zol(i,j)
          enddo
          enddo
          varname='zol'
          description='z/L (z over Monin-Obukhov length)'
          units='dimensionless'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=mol(i,j)
          enddo
          enddo
          varname='mol'
          description='T* (similarity theory)'
          units='K'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

          do j=1,nj
          do i=1,ni
            ds2(i,j)=br(i,j)
          enddo
          enddo
          varname='br'
          description='bulk Richardson number in surface layer'
          units='dimensionless'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)
      endif

! begin 3d group

      if(output_th  .eq.1)then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=th0(i,j,k)+tha(i,j,k)
          enddo
          enddo
          enddo
          varname='th'
          description = 'Potential temperature'
          units = 'K'
          call h5_write_3d_float(threed_group_id,varname,description,units,3,ds,ni,nj,nk)
      endif

      if(output_thpert .eq.1) then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=tha(i,j,k)
          enddo
          enddo
          enddo
          varname='thpert'
          description = 'Perturbation potential temperature'
          units = 'K'
          call h5_write_3d_float(threed_group_id,varname,description,units,3,ds,ni,nj,nk)
      endif

      if(output_prs    .eq.1) then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=prs(i,j,k)
          enddo
          enddo
          enddo
          varname='p'
          description = 'Pressure'
          units = 'Pa'
          call h5_write_3d_float(threed_group_id,varname,description,units,1,ds,ni,nj,nk)
      endif

      if(output_prspert.eq.1)then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=prs(i,j,k)-p00*(pi0(i,j,k)**cpdrd)
          enddo
          enddo
          enddo
          varname="ppert"
          description = 'Perturbation pressure'
          units = 'Pa'
          call h5_write_3d_float(threed_group_id,varname,description,units,3,ds,ni,nj,nk)
      endif
      if(output_pi.eq.1)then  
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=pi0(i,j,k)+ppi(i,j,k)
          enddo
          enddo
          enddo
          varname = "pi"
          description = 'Nondimensional pressure (Exner function)'
          units = 'dimensionless'
          call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)
      endif
      if(output_pipert .eq.1) then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=ppi(i,j,k)
          enddo
          enddo
          enddo
          varname='pipert'
          description = 'Nondimensional pressure (Exner function) perturbation'
          units = 'dimensionless'
          call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)
      endif

      if(output_rho    .eq.1) then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=rho(i,j,k)
          enddo
          enddo
          enddo
          varname = "rho"
          description = 'Density'
          units = 'kg/m^3'
          call h5_write_3d_float(threed_group_id,varname,description,units,3,ds,ni,nj,nk)
      endif
      if(output_rhopert.eq.1)then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=rho(i,j,k)-rho0(i,j,k)
          enddo
          enddo
          enddo
          varname="rhopert"
          description = 'Density perturbation'
          units = 'kg/m^3'
          call h5_write_3d_float(threed_group_id,varname,description,units,4,ds,ni,nj,nk)
      endif
      if(iturb.eq.1)then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=tkea(i,j,k)
          enddo
          enddo
          enddo
          varname="tke"
          description = 'Subgrid turbulence kinetic energy'
          units = 'J/kg'
          call h5_write_3d_float(threed_group_id,varname,description,units,4,ds,ni,nj,nk)
      endif

      if(iturb.ge.1)then

          if(output_km  .eq.1) then
              do k=1,nk
              do j=1,nj
              do i=1,ni
                ds(i,j,k)=kmh(i,j,k)
              enddo
              enddo
              enddo

              varname='kmh'
              description = 'Eddy mixing coefficient for momentum in the horizontal direction'
              units = 'm^2/s'
              call h5_write_3d_float(threed_group_id,varname,description,units,4,ds,ni,nj,nk)

              do k=1,nk
              do j=1,nj
              do i=1,ni
                ds(i,j,k)=kmv(i,j,k)
              enddo
              enddo
              enddo
              varname='kmv'
              description = 'Eddy mixing coefficient for momentum in the vertical direction'
              units = 'm^2/s'
              call h5_write_3d_float(threed_group_id,varname,description,units,4,ds,ni,nj,nk)
          endif

         if(output_kh  .eq.1) then
           do k=1,nk
           do j=1,nj
           do i=1,ni
             ds(i,j,k)=khh(i,j,k)
           enddo
           enddo
           enddo
           varname='khh'
           description = 'Eddy mixing coefficient for scalars in the horizontal direction'
           units = 'm^2/s'
           call h5_write_3d_float(threed_group_id,varname,description,units,4,ds,ni,nj,nk)

           do k=1,nk
           do j=1,nj
           do i=1,ni
             ds(i,j,k)=khv(i,j,k)
           enddo
           enddo
           enddo
           varname='khv'
           description = 'Eddy mixing coefficient for scalars in the vertical direction'
           units = 'm^2/s'
           call h5_write_3d_float(threed_group_id,varname,description,units,4,ds,ni,nj,nk)
         endif
      endif !iturb .ge. 1

      if(imoist.eq.1)then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=qa(i,j,k,nqv)
          enddo
          enddo
          enddo
          varname='qv'
          description = 'Water vapor mixing ratio'
          units = 'kg/kg                         '
          if(output_qv .eq. 1)    call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)
          if(output_qvpert.eq.1)then
              do k=1,nk
              do j=1,nj
              do i=1,ni
                ds(i,j,k)=qa(i,j,k,nqv)-qv0(i,j,k)
              enddo
              enddo
              enddo
              varname="qvpert"
              description = 'Perturbation water vapor mixing ratio'
              units = 'kg/kg                         '
              call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)
          endif

          if(output_q.eq.1)then
              do n=1,numq
                  if(n.ne.nqv)then
                      do k=1,nk
                      do j=1,nj
                      do i=1,ni
                        ds(i,j,k)=qa(i,j,k,n)
                      enddo
                      enddo
                      enddo
                      qvarname='                    '
                      write(qvarname(1:3),156) qname(n)
156                   Format(a3)
                      if(idm.eq.1.and.n.ge.nnc1.and.n.le.nnc2)then
                           description = 'Number concentration'
                           units = 'kg^{-1}'
                      else
                           description = 'Hydrometeor mixing ratio'
                           units = 'kg/kg'
                      endif
                      call h5_write_3d_float(threed_group_id,qvarname,description,units,6,ds,ni,nj,nk)
                  endif
              enddo
          endif !output_q .eq.1
          if(output_dbz.eq.1) then
              do k=1,nk
              do j=1,nj
              do i=1,ni
                ds(i,j,k)=dbz(i,j,k)
              enddo
              enddo
              enddo
              varname="dbz"
              description = 'Reflectivity'
              units = 'dBZ'
              call h5_write_3d_float(threed_group_id,varname,description,units,2,ds,ni,nj,nk)
          endif

      endif !imoist.eq.1

      if(output_uinterp.eq.1)then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=0.5*(ua(i,j,k)+ua(i+1,j,k))
          enddo
          enddo
          enddo
          varname="uinterp"
          description = 'velocity in x-direction, interpolated to scalar points'
          units = 'm/s'
          call h5_write_3d_float(threed_group_id,varname,description,units,3,ds,ni,nj,nk)
      endif
      if(output_vinterp.eq.1)then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=0.5*(va(i,j,k)+va(i,j+1,k))
          enddo
          enddo
          enddo
          varname="vinterp"
          description = 'velocity in y-direction, interpolated to scalar points'
          units = 'm/s'
          call h5_write_3d_float(threed_group_id,varname,description,units,3,ds,ni,nj,nk)
      endif
      if(output_winterp.eq.1)then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=0.5*(wa(i,j,k)+wa(i,j,k+1))
          enddo
          enddo
          enddo
          varname='winterp'
          description = 'velocity in z direction, interpolated to scalar points'
          units = 'm/s'
          call h5_write_3d_float(threed_group_id,varname,description,units,3,ds,ni,nj,nk)
      endif


      if(iptra.eq.1)then
          do n=1,npt
            do k=1,nk
            do j=1,nj
            do i=1,ni
              ds(i,j,k)=pta(i,j,k,n)
            enddo
            enddo
            enddo

            chid = 'pt      '
            write(chid(3:4),111) n
111         format(i2.2)

            varname=chid
            description='mixing ratio of passive tracer'
            units='kg/kg'
            call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)
          enddo
      endif


      if(output_pblten.eq.1)then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=thpten(i,j,k)
          enddo
          enddo
          enddo

          varname='thpten'
          description='pbl tendency: theta'
          units='Not sure, email leigh.orf@cmich.edu if you do know'
          call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=qvpten(i,j,k)
          enddo
          enddo
          enddo

          varname='qvpten'
          description='pbl tendency: qv'
          units='Not sure, email leigh.orf@cmich.edu if you do know'
          call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=qcpten(i,j,k)
          enddo
          enddo
          enddo

          varname='qcpten'
          description='pbl tendency: qc'
          units='Not sure, email leigh.orf@cmich.edu if you do know'
          call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=qipten(i,j,k)
          enddo
          enddo
          enddo

          varname='qipten'
          description='pbl tendency: qi'
          units='Not sure, email leigh.orf@cmich.edu if you do know'
          call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=upten(i,j,k)
          enddo
          enddo
          enddo

          varname='upten'
          description='pbl tendency: u'
          units='Not sure, email leigh.orf@cmich.edu if you do know'
          call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=vpten(i,j,k)
          enddo
          enddo
          enddo

          varname='vpten'
          description='pbl tendency: v'
          units='Not sure, email leigh.orf@cmich.edu if you do know'
          call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)

      endif !output_pblten.eq.1

      if(output_uh.eq.1) then

        ! get height AGL:
          if( terrain_flag )then
            do k=1,nk+1
            do j=1,nj
            do i=1,ni
              dum3(i,j,k) = zh(i,j,k)-zs(i,j)
              dumw(i,j,k) = zf(i,j,k)-zs(i,j)
            enddo
            enddo
            enddo
          else
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
          do j=1,nj
          do i=1,ni
            ds2(i,j)=dum1(i,j,1)
          enddo
          enddo

          varname='uh'
          description='integrated updraft helicity'
          units='m^2/s^2'
          call h5_write_2d_float(twod_group_id,varname,description,units,ds2,ni,nj)

      endif!output_uh.eq.1

      if(output_vort.eq.1)then
          call calcvort(xh,xf,uf,vf,zh,mf,zf,ua,va,wa,dum1,dum2,dum3,dum4)

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=dum1(i,j,k)
          enddo
          enddo
          enddo

          varname='xvort'
          description='horizontal vorticity (x), usually called xi'
          units='s^{-1}'
          call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=dum2(i,j,k)
          enddo
          enddo
          enddo

          varname='yvort'
          description='horizontal vorticity (y), usually called eta'
          units='s^{-1}'
          call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=dum3(i,j,k)
          enddo
          enddo
          enddo

          varname='zvort'
          description='vertical vorticity (zeta)'
          units='s^{-1}'
          call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)

      endif !output_vort.eq.1

      if(output_dissten.eq.1) then

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ds(i,j,k)=dissten(i,j,k)
          enddo
          enddo
          enddo

          varname='dissten'
          description='dissipative heating tendency'
          units='K/s'
          call h5_write_3d_float(threed_group_id,varname,description,units,6,ds,ni,nj,nk)

       endif


!--------------------------------------------------------------
!  writeout data on u-points

! Note velocity data is staggered and requires 1 extra point

       if(output_u    .eq.1) then
           do k=1,nk
           do j=1,nj
           do i=1,ni+1
             dumu(i,j,k)=ua(i,j,k)
           enddo
           enddo
           enddo
           varname='u'
           description = 'velocity in x direction on native Arakawa C grid'
           units = 'm/s'
           call h5_write_3d_float(threed_group_id,varname,description,units,3,dumu(1:ni+1,1:nj,1:nk),ni+1,nj,nk)
       endif


      if(output_upert.eq.1)then
          do k=1,nk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k)=ua(i,j,k)-u0(i,j,k)
          enddo
          enddo
          enddo
          varname='upert'
          description = 'perturbation velocity in x direction on native Arakawa C grid'
          units = 'm/s'
          call h5_write_3d_float(threed_group_id,varname,description,units,3,dumu(1:ni+1,1:nj,1:nk),ni+1,nj,nk)
      endif

!--------------------------------------------------------------
!  writeout data on v-points

      if(output_v    .eq.1) then
          do k=1,nk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k)=va(i,j,k)
          enddo
          enddo
          enddo
          varname='v'
          description = 'velocity in y direction on native Arakawa C grid'
          units = 'm/s'
          call h5_write_3d_float(threed_group_id,varname,description,units,3,dumv(1:ni,1:nj+1,1:nk),ni,nj+1,nk)
      endif


      if(output_vpert.eq.1)then
          do k=1,nk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k)=va(i,j,k)-v0(i,j,k)
          enddo
          enddo
          enddo
          varname='vpert'
          description = 'perturbation velocity in y direction on native Arakawa C grid'
          units = 'm/s'
          call h5_write_3d_float(threed_group_id,varname,description,units,3,dumv(1:ni,1:nj+1,1:nk),ni,nj+1,nk)
      endif

!--------------------------------------------------------------
!  writeout data on w-points

      if(output_w    .eq.1) then
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k)=wa(i,j,k)
          enddo
          enddo
          enddo
          varname='w'
          description = 'velocity in z direction on native Arakawa C grid'
          units = 'm/s'
          call h5_write_3d_float(threed_group_id,varname,description,units,3,dumw(1:ni,1:nj,1:nk+1),ni,nj,nk+1)
      endif

!--------------------------------------------------------------

      write(outfile,*) 'Done Writing Data to File: ',trim(filename)
! NOTE if compiler says undefined symbol here, try flush_ (IBM)
!     call flush(outfile)

! close everything
      call h5gclose_f(twod_group_id,ierror)
      call h5gclose_f(threed_group_id,ierror)
      call h5fclose_f(file_id,ierror)
      call h5close_f(ierror)

#ifdef MPI
!     call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

      return
      end

!--------------------------------------------------------------------------------------------MARK

! Have messed around with chunking. For 60x60x80 (per hdf file),
! 60x60x20 chunking works well. I think with small datasets you really
! don't gain much with chunking. Chunking horizontally (e.g., 30x30x20)
! leads to a penalty with small files. You may wish to experiment
! with chunking if you are interested in saving space and are using
! compression, or wish to increase file access and write performance
! (see hdf5 docs).

! Also note that scale_factor is passed to this routine, you may wish
! to change this value if you find that you are seeing compression
! artifacts or you are getting too much (or too little) compression. The
! value I am passing is the one I chose based upon the expected range of
! the data and the level of nrecision I wanted.

      subroutine h5_write_3d_float(file_id,varname,description,units,scale_factor,fl_array_3d,numi,numj,numk)
      use HDF5

      implicit none
      include 'input.incl'
      character*120  description
      character*20  units
      integer(HSIZE_T), DIMENSION(3) :: dims,chunkdims
      integer(HSIZE_T), DIMENSION(1) :: strdims
      integer(HSIZE_T) level
      integer :: arank = 1
      integer(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      integer(SIZE_T) :: stringlen
      integer(HID_T) :: file_id,dset_id,dspace_id,chunk_id,datatype_id,desctype_id,unitstype_id,desc_id,units_id,aspace_id
      integer :: scale_type, scale_factor
      integer :: options_mask, pix_per_block
      integer :: rank,numi,numj,numk,ierror
      character(LEN=20) :: varname
      real, dimension(numi,numj,numk) :: fl_array_3d
      integer, dimension(2) :: szip_filter_values
      integer :: gzip_compression_level
      integer(size_t) :: n0,n1,n2

      gzip_compression_level = 1 ! Bump this up if you want more compression

      rank=3
      dims(1)=numi
      dims(2)=numj
      dims(3)=numk
      strdims(1)=1
      level=20
      chunkdims(1) = dims(1) ! do not make this smaller or file gets much(?) larger
      chunkdims(2) = dims(2) ! do not make this smaller or file gets much(?) larger
      chunkdims(3) = min(level,dims(3)) ! Tried a few values, 20 worked nicely

! ORF 1/5/11 add description and units attributes

      call h5screate_simple_f(arank,adims,aspace_id,ierror)
      call h5tcopy_f(H5T_NATIVE_CHARACTER,desctype_id,ierror)
      call h5tcopy_f(H5T_NATIVE_CHARACTER,unitstype_id,ierror)
      stringlen = len_trim(description)
      call h5tset_size_f(desctype_id,stringlen,ierror)
      stringlen = len_trim(units)
      call h5tset_size_f(unitstype_id,stringlen,ierror)

!     call h5tcopy_f(H5T_IEEE_F32LE,datatype_id,ierror)
      call h5screate_simple_f(rank,dims,dspace_id,ierror)
      call h5pcreate_f(H5P_DATASET_CREATE_F,chunk_id,ierror)
      call h5pset_chunk_f(chunk_id,rank,chunkdims,ierror)

! sometimes shuffling first makes compression more efficient
!     call h5pset_shuffle_f (chunk_id,ierror)

      if (output_format.eq.3) call h5pset_scaleoffset_f (chunk_id,0,scale_factor,ierror)
      if (output_format.eq.3.or.output_format.eq.4)call h5pset_deflate_f (chunk_id,gzip_compression_level,ierror)

!     call h5dcreate_f(file_id,trim(varname),datatype_id,dspace_id,dset_id,ierror,chunk_id)
      call h5dcreate_f(file_id,trim(varname),H5T_NATIVE_REAL,dspace_id,dset_id,ierror,chunk_id)
      call h5acreate_f(dset_id,"description",desctype_id,aspace_id,desc_id,ierror,H5P_DEFAULT_F,H5P_DEFAULT_F)
      call h5awrite_f(desc_id,desctype_id,description,strdims,ierror)
      call h5acreate_f(dset_id,"units",unitstype_id,aspace_id,units_id,ierror,H5P_DEFAULT_F,H5P_DEFAULT_F)
      call h5awrite_f(units_id,unitstype_id,units,strdims,ierror)

!     call h5dwrite_f(dset_id,datatype_id,fl_array_3d,dims,ierror)
      call h5dwrite_f(dset_id,H5T_NATIVE_REAL,fl_array_3d,dims,ierror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)

      call H5aclose_f (desc_id,ierror)
      call H5aclose_f (units_id,ierror)
      call H5sclose_f (aspace_id,ierror)
      call H5tclose_f (desctype_id,ierror)
      call H5tclose_f (unitstype_id,ierror)
!     call H5tclose_f (datatype_id,ierror)
      call h5dclose_f (dset_id,ierror)
      call h5sclose_f (dspace_id,ierror)
      call H5pclose_f (chunk_id,ierror)

      return
      end

!--------------------------------------------------------------------------------------------MARK
      subroutine h5_write_float (file_id,varname,description,units,value)
      use HDF5
      implicit none
      character(LEN=20) :: varname
      character*120  description
      character*20  units

      integer(HSIZE_T), DIMENSION(1) :: strdims
      integer :: arank = 1
      integer(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      integer(SIZE_T) :: stringlen

      real value
      integer rank,ierror
      integer(HSIZE_T), DIMENSION(1) :: dims
      integer(HID_T) :: file_id,dset_id,dspace_id,aspace_id,desc_id,desctype_id,unitstype_id,units_id

      rank=1
      strdims(1)=1
      dims(1)=1

      call h5screate_simple_f(arank,adims,aspace_id,ierror)
      call h5tcopy_f(H5T_NATIVE_CHARACTER,desctype_id,ierror)
      call h5tcopy_f(H5T_NATIVE_CHARACTER,unitstype_id,ierror)
      stringlen = len_trim(description)
      call h5tset_size_f(desctype_id,stringlen,ierror)
      stringlen = len_trim(units)
      call h5tset_size_f(unitstype_id,stringlen,ierror)

      call h5screate_simple_f(rank,dims,dspace_id,ierror)
      call h5dcreate_f(file_id,trim(varname),H5T_NATIVE_REAL,dspace_id,dset_id,ierror)
      call h5dwrite_f(dset_id,H5T_NATIVE_REAL,value,dims,ierror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)

      call h5acreate_f(dset_id,"description",desctype_id,aspace_id,desc_id,ierror,H5P_DEFAULT_F,H5P_DEFAULT_F)
      call h5awrite_f(desc_id,desctype_id,trim(description),strdims,ierror)
      call h5acreate_f(dset_id,"units",unitstype_id,aspace_id,units_id,ierror,H5P_DEFAULT_F,H5P_DEFAULT_F)
      call h5awrite_f(units_id,unitstype_id,trim(units),strdims,ierror)


      call H5aclose_f (desc_id,ierror)
      call H5aclose_f (units_id,ierror)
      call H5sclose_f (aspace_id,ierror)
      call H5tclose_f (desctype_id,ierror)
      call H5tclose_f (unitstype_id,ierror)
      call h5dclose_f (dset_id,ierror)
      call h5sclose_f (dspace_id,ierror)
      return
      end
!--------------------------------------------------------------------------------------------MARK
      subroutine h5_write_int (file_id,varname,description,units,value)
      use HDF5
      implicit none
      character(LEN=20) :: varname
      character*120  description
      character*20  units

      integer(HSIZE_T), DIMENSION(1) :: strdims
      integer :: arank = 1
      integer(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      integer(SIZE_T) :: stringlen

      integer value
      integer rank,ierror
      integer(HSIZE_T), DIMENSION(1) :: dims
      integer(HID_T) :: file_id,dset_id,dspace_id,aspace_id,desc_id,desctype_id,unitstype_id,units_id

      rank=1
      strdims(1)=1
      dims(1)=1

      call h5screate_simple_f(arank,adims,aspace_id,ierror)
      call h5tcopy_f(H5T_NATIVE_CHARACTER,desctype_id,ierror)
      call h5tcopy_f(H5T_NATIVE_CHARACTER,unitstype_id,ierror)
      stringlen = len_trim(description)
      call h5tset_size_f(desctype_id,stringlen,ierror)
      stringlen = len_trim(units)
      call h5tset_size_f(unitstype_id,stringlen,ierror)

      call h5screate_simple_f(rank,dims,dspace_id,ierror)
      call h5dcreate_f(file_id,trim(varname),H5T_NATIVE_INTEGER,dspace_id,dset_id,ierror)
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,value,dims,ierror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)

      call h5acreate_f(dset_id,"description",desctype_id,aspace_id,desc_id,ierror,H5P_DEFAULT_F,H5P_DEFAULT_F)
      call h5awrite_f(desc_id,desctype_id,trim(description),strdims,ierror)
      call h5acreate_f(dset_id,"units",unitstype_id,aspace_id,units_id,ierror,H5P_DEFAULT_F,H5P_DEFAULT_F)
      call h5awrite_f(units_id,unitstype_id,trim(units),strdims,ierror)


      call H5aclose_f (desc_id,ierror)
      call H5aclose_f (units_id,ierror)
      call H5sclose_f (aspace_id,ierror)
      call H5tclose_f (desctype_id,ierror)
      call H5tclose_f (unitstype_id,ierror)
      call h5dclose_f (dset_id,ierror)
      call h5sclose_f (dspace_id,ierror)
      return
      end
!--------------------------------------------------------------------------------------------MARK
      subroutine h5_write_1d_float (file_id,varname,description,units,fl_array_1d,numk)
      use HDF5
      implicit none
      character(LEN=20) :: varname
      character*120  description
      character*20  units

      integer(HSIZE_T), DIMENSION(1) :: strdims
      integer :: arank = 1
      integer(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      integer(SIZE_T) :: stringlen

      real value
      integer rank,numk,ierror
      real, dimension(1:numk) :: fl_array_1d
      integer(HSIZE_T), DIMENSION(1) :: dims
      integer(HID_T) :: file_id,dset_id,dspace_id,aspace_id,desc_id,desctype_id,unitstype_id,units_id

      rank=1
      strdims(1)=1
      dims(1)=numk

      call h5screate_simple_f(arank,adims,aspace_id,ierror)
      call h5tcopy_f(H5T_NATIVE_CHARACTER,desctype_id,ierror)
      call h5tcopy_f(H5T_NATIVE_CHARACTER,unitstype_id,ierror)
      stringlen = len_trim(description)
      call h5tset_size_f(desctype_id,stringlen,ierror)
      stringlen = len_trim(units)
      call h5tset_size_f(unitstype_id,stringlen,ierror)

      call h5screate_simple_f(rank,dims,dspace_id,ierror)
      call h5dcreate_f(file_id,trim(varname),H5T_NATIVE_REAL,dspace_id,dset_id,ierror)
      call h5dwrite_f(dset_id,H5T_NATIVE_REAL,fl_array_1d,dims,ierror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)

      call h5acreate_f(dset_id,"description",desctype_id,aspace_id,desc_id,ierror,H5P_DEFAULT_F,H5P_DEFAULT_F)
      call h5awrite_f(desc_id,desctype_id,trim(description),strdims,ierror)
      call h5acreate_f(dset_id,"units",unitstype_id,aspace_id,units_id,ierror,H5P_DEFAULT_F,H5P_DEFAULT_F)
      call h5awrite_f(units_id,unitstype_id,trim(units),strdims,ierror)


      call H5aclose_f (desc_id,ierror)
      call H5aclose_f (units_id,ierror)
      call H5sclose_f (aspace_id,ierror)
      call H5tclose_f (desctype_id,ierror)
      call H5tclose_f (unitstype_id,ierror)
      call h5dclose_f (dset_id,ierror)
      call h5sclose_f (dspace_id,ierror)
      return
      end
!--------------------------------------------------------------------------------------------MARK

      subroutine h5_write_2d_float(file_id,varname,description,units,fl_array_2d,numi,numj)
      use HDF5

      implicit none
      character(LEN=20) :: varname
      character*20  units
      character*120  description
      integer(HSIZE_T), DIMENSION(2) :: dims
      integer(HSIZE_T), DIMENSION(1) :: strdims
      integer :: arank = 1
      integer(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      integer(SIZE_T) :: stringlen
      integer(HID_T) :: file_id,dset_id,dspace_id,chunk_id,datatype_id,desctype_id,unitstype_id,desc_id,units_id,aspace_id
      integer :: rank,numi,numj,ierror
      real, dimension(numi,numj) :: fl_array_2d

      rank=2
      dims(1)=numi
      dims(2)=numj
      strdims(1)=1

! ORF TODO: Add chunking and gzip as we do in the other 2d write code with
! unlimited

      call h5screate_simple_f(arank,adims,aspace_id,ierror)
      call h5tcopy_f(H5T_NATIVE_CHARACTER,desctype_id,ierror)
      call h5tcopy_f(H5T_NATIVE_CHARACTER,unitstype_id,ierror)
      stringlen = len_trim(description)
      call h5tset_size_f(desctype_id,stringlen,ierror)
      stringlen = len_trim(units)
      call h5tset_size_f(unitstype_id,stringlen,ierror)

!     call h5tcopy_f(H5T_IEEE_F32LE,datatype_id,ierror)
      call h5screate_simple_f(rank,dims,dspace_id,ierror)
      call h5pcreate_f(H5P_DATASET_CREATE_F,chunk_id,ierror)

!     call h5dcreate_f(file_id,trim(varname),datatype_id,dspace_id,dset_id,ierror,chunk_id)
      call h5dcreate_f(file_id,trim(varname),H5T_NATIVE_REAL,dspace_id,dset_id,ierror,chunk_id)
      call h5acreate_f(dset_id,"description",desctype_id,aspace_id,desc_id,ierror,H5P_DEFAULT_F,H5P_DEFAULT_F)
      call h5awrite_f(desc_id,desctype_id,description,strdims,ierror)
      call h5acreate_f(dset_id,"units",unitstype_id,aspace_id,units_id,ierror,H5P_DEFAULT_F,H5P_DEFAULT_F)
      call h5awrite_f(units_id,unitstype_id,units,strdims,ierror)

!     call h5dwrite_f(dset_id,datatype_id,fl_array_2d,dims,ierror)
      call h5dwrite_f(dset_id,H5T_NATIVE_REAL,fl_array_2d,dims,ierror,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)

      call H5aclose_f (desc_id,ierror)
      call H5aclose_f (units_id,ierror)
      call H5sclose_f (aspace_id,ierror)
      call H5tclose_f (desctype_id,ierror)
      call H5tclose_f (unitstype_id,ierror)
!     call H5tclose_f (datatype_id,ierror)
      call h5dclose_f (dset_id,ierror)
      call h5sclose_f (dspace_id,ierror)
      call H5pclose_f (chunk_id,ierror)

      return
      end

!--------------------------------------------------------------------------------------------MARK

! Fields along parcels are hardwired in parcel.F. If you add more fields, change
! here accordingly.

subroutine writepdata_hdf5(nrec,rtime,pdata)

      use HDF5
      implicit none
      include 'input.incl'

      character*160 filename
      integer :: nrec
      real :: rtime
      real, dimension(npvals,nparcels) :: pdata
      integer(SIZE_T) :: stringlen
      integer(HID_T) :: file_id,dset_id,dspace_id,chunk_id,datatype_id
      integer(HID_T) :: desctype_id,unitstype_id,desc_id,units_id,aspace_id
      integer(HID_T) :: parcel_group_id,memoryspace,filespace,time_id
      integer(HSIZE_T), dimension(2) :: maxdims,chunkdims,currentdims,singledims,offset
      integer(HSIZE_T), dimension(1) :: maxtimedims,timechunkdims,timecurrentdims,timesingledims,timeoffset
      integer(HID_T) :: timechunk_id,timememoryspace,timefilespace,timedset_id
      character(len=6), parameter :: parcel_group_name="parcel"
      character(len=120), parameter :: timedesc="model time"
      character(len=20), parameter :: timevar="time",timeunits="s"
      integer :: i,ierror
      integer :: rank,nvals

      character(len=20), parameter :: varname(14) = (/character(len=20) :: &
      'x', &
      'y', &
      'z', &
      'qv', &
      'qc', &
      'qr', &
      'nsq', &
      'u', &
      'v', &
      'w', &
      'kh', &
      'the', &
      'b', &
      'dpdz'/)

      character(len=120), parameter :: vardesc(14) = (/character(len=120) :: &
      'x position of parcel', &
      'y position of parcel', &
      'z position of parcel', &
      'water vapor mixing ratio', &
      'cloud water mixing ratio', &
      'rain water mixing ratio', &
      'Brunt Vaisala frequency squared', &
      'u component of wind', &
      'v component of wind', &
      'w component of wind', &
      'turbulence coefficient for scalars', &
      'theta-e', &
      'buoyancy', &
      'dpdz tendency'/)

      character(len=20), parameter :: varunits(14) = (/character(len=20) :: &
      'm','m','m', &
      'kg/kg','kg/kg','kg/kg', &
      'dimensionless', &
      'm/s','m/s','m/s', &
      'm^2/s','K', &
      'm/s^2', &
      'm/s^2'/)

      call h5open_f (ierror)

      maxdims(1) = nparcels
      maxdims(2) = H5S_UNLIMITED_F
      currentdims(1)=nparcels
      currentdims(2)=1
      singledims(1)=nparcels
      singledims(2)=1
      chunkdims(1) =nparcels
      chunkdims(2) = 1
      offset(1) = 0
      offset(2) = 0

      maxtimedims(1) = H5S_UNLIMITED_F
      timechunkdims(1) = 100
      timecurrentdims(1)=1
      timesingledims(1)=1
      timeoffset(1)=0

      timeoffset(1)=nrec-1
      timecurrentdims(1)=nrec
      offset(2) = nrec-1
      currentdims(2) = nrec

      filename = trim(output_basename)//'_pdata.hdf5'
      if(nrec.eq.1) then 
             call h5fcreate_f(trim(filename),H5F_ACC_TRUNC_F,file_id,ierror)
             call h5gcreate_f(file_id,parcel_group_name,parcel_group_id,ierror)      
      else
!            print *,'myid = ',myid,' nrec= ',nrec,' filename = ',trim(filename)
             call h5fopen_f(trim(filename),H5F_ACC_RDWR_F, file_id,ierror)
             call h5gopen_f(file_id,parcel_group_name,parcel_group_id,ierror)
      endif

      rank=1
      nvals=1
      call writefloat_hdf_unlimited(rank,nvals,nrec,parcel_group_id,timevar,timedesc,timeunits,time_id,timefilespace, &
            timememoryspace,timedset_id,timechunk_id,timechunkdims,maxtimedims,timecurrentdims,timesingledims,timeoffset,rtime)

      rank=2
      nvals=nparcels
      do i=1,14
            call writefloat_hdf_unlimited(rank,nvals,nrec,parcel_group_id,varname(i),vardesc(i), &
                  varunits(i),dspace_id,filespace,memoryspace,dset_id,chunk_id, &
                  chunkdims,maxdims,currentdims,singledims,offset,pdata(i,:))
            enddo

      nrec = nrec + 1
      call h5close_f (ierror)

      return
      end 


!--------------------------------------------------------------------------------------------MARK
! ORF 1/11/11
! Following routine works on data with an unlimited dimension,
! i.e., for creating and appending variables in parcel or stats hdf5 files

subroutine writefloat_hdf_unlimited(rank,nvals,irec,file_id,varname,vardesc,varunits,dspace_id, &
      filespace,memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,floatval)

      use HDF5
      implicit none
      include 'input.incl'

      integer(HID_T) :: file_id,dset_id,dspace_id,chunk_id,datatype_id
      integer(HID_T) :: desctype_id,unitstype_id,desc_id,units_id,aspace_id
      integer(HID_T) :: parcel_group_id,memoryspace,filespace
      integer rank,nvals
      integer(HSIZE_T),dimension(rank) :: chunkdims,maxdims,currentdims,singledims,offset
      integer(HSIZE_T), DIMENSION(1) :: adims = (/1/)
      integer(HSIZE_T), DIMENSION(1) :: strdims = (/1/)
      real, dimension(nvals) :: floatval
      character(len=20) ::  varname
      character(len=120) :: vardesc
      character(len=20) :: varunits
      integer irec,ierror
      integer(SIZE_T) :: stringlen

!     print *, 'myid: ',myid,'irec: ',irec, varname
      if (irec.eq.1) then

            call h5screate_simple_f(rank, singledims, dspace_id, ierror, maxdims)
            call h5screate_simple_f(rank, singledims, memoryspace, ierror, singledims)
            call h5pcreate_f(H5P_DATASET_CREATE_F,chunk_id,ierror)
            call h5pset_chunk_f(chunk_id,rank,chunkdims,ierror)
! This is an example of how you would add gzip compression to unlimited dimension arrays for parcel and stats hdf5 files
!           if (output_format.eq.3.or.output_format.eq.4)call h5pset_deflate_f (chunk_id,1,ierror)

            call h5dcreate_f(file_id,trim(varname),H5T_NATIVE_REAL,dspace_id,dset_id,ierror,chunk_id)
            call h5dset_extent_f(dset_id,currentdims,ierror)
            call h5dget_space_f(dset_id,filespace,ierror)
            call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,singledims,ierror)
            call h5dwrite_f(dset_id,H5T_NATIVE_REAL,floatval,singledims,ierror,memoryspace,filespace,H5P_DEFAULT_F)
            call h5sclose_f(memoryspace,ierror)
            call h5sclose_f(filespace,ierror)
            call h5sclose_f(dspace_id,ierror)
            call H5pclose_f(chunk_id,ierror)

            call h5tcopy_f(H5T_NATIVE_CHARACTER,desctype_id,ierror)
            stringlen = len_trim(vardesc)
            call h5tset_size_f(desctype_id,stringlen,ierror)
            call h5screate_simple_f(1,adims,aspace_id,ierror)
            call h5acreate_f(dset_id,"description",desctype_id,aspace_id,desc_id,ierror,H5P_DEFAULT_F,H5P_DEFAULT_F)
            call h5awrite_f(desc_id,desctype_id,vardesc,strdims,ierror)
            call H5aclose_f (desc_id,ierror)
            call H5sclose_f (aspace_id,ierror)
            call H5tclose_f (desctype_id,ierror)

            call h5tcopy_f(H5T_NATIVE_CHARACTER,desctype_id,ierror)
            stringlen = len_trim(varunits)
            call h5tset_size_f(desctype_id,stringlen,ierror)
            call h5screate_simple_f(1,adims,aspace_id,ierror)
            call h5acreate_f(dset_id,"units",desctype_id,aspace_id,desc_id,ierror,H5P_DEFAULT_F,H5P_DEFAULT_F)
            call h5awrite_f(desc_id,desctype_id,varunits,strdims,ierror)
            call H5aclose_f (desc_id,ierror)
            call H5sclose_f (aspace_id,ierror)
            call H5tclose_f (desctype_id,ierror)

            call h5dclose_f(dset_id,ierror)
      else
!           print *,'myid = ',myid,' irec= ',irec,' varname = ',varname
            call h5dopen_f(file_id,trim(varname),dset_id,ierror)
            call h5dset_extent_f(dset_id,currentdims,ierror)
            call h5dget_space_f(dset_id,filespace,ierror)
            call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,singledims,ierror)
            call h5screate_simple_f(rank, singledims, memoryspace, ierror, singledims)
            call h5dwrite_f(dset_id,H5T_NATIVE_REAL,floatval,singledims,ierror,memoryspace,filespace,H5P_DEFAULT_F)
            call h5sclose_f(memoryspace,ierror)
            call h5sclose_f(filespace,ierror)
            call h5dclose_f(dset_id,ierror)
      endif
      return
end 


!--------------------------------------------------------------------------------------------MARK
subroutine writestat_hdf5(nrec,rtime,nstat,rstat,qname,budname)

      use HDF5
      implicit none
      include 'input.incl'

      integer, intent(inout) :: nrec
      real,    intent(in)    :: rtime
      integer, intent(in)    :: nstat
      real, dimension(stat_out), intent(in) :: rstat
      character*3, dimension(maxq), intent(in) :: qname
      character*6, dimension(maxq), intent(in) :: budname

      character*160 filename
      integer(SIZE_T) :: stringlen
      integer(HID_T) :: file_id,dset_id,dspace_id,chunk_id,datatype_id
      integer(HID_T) :: desctype_id,unitstype_id,desc_id,units_id,aspace_id
      integer(HID_T) :: stats_group_id,memoryspace,filespace,var_id
      integer(HSIZE_T), dimension(1) :: maxdims,chunkdims,currentdims,singledims,offset
      character(len=6), parameter :: stats_group_name="stats"
      integer :: i,n,ierror
      integer :: rank,nvals
      character(len=20) :: varname
      character(len=20) :: units
      character(len=120)::  description
      character(len=8) :: text1
      character(len=30) :: text2

      call h5open_f (ierror)

      maxdims(1) = H5S_UNLIMITED_F
      chunkdims(1) = 100
      currentdims(1)=1
      singledims(1)=1
      offset(1)=0

      rank=1
      nvals=1

!     print *,'nrec = ',nrec
!     print *,'rtime = ',rtime

      filename = trim(output_basename)//'_stats.hdf5'
      if(nrec.eq.1) then 
             call h5fcreate_f(trim(filename),H5F_ACC_TRUNC_F,file_id,ierror)
             call h5gcreate_f(file_id,stats_group_name,stats_group_id,ierror)      
      else
!            print *,'myid = ',myid,' nrec= ',nrec,' filename = ',trim(filename)
             call h5fopen_f(trim(filename),H5F_ACC_RDWR_F, file_id,ierror)
             call h5gopen_f(file_id,stats_group_name,stats_group_id,ierror)
      endif

      offset(1)=nrec-1
      currentdims(1)=nrec

      varname='time'
      description='model time'
      units='s'
      call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
            memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rtime)
 
      i=0
      if(adapt_dt.eq.1) then
            i=i+1
            varname='dt'
            description='average timestep (this simulation used adaptive timestepping)'
            units='s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_w.eq.1) then
            i=i+1
            varname='wmax'
            description='maximum vertical velocity'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='wmin'
            description='minimum vertical velocity'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_u.eq.1) then
            i=i+1
            varname='umax'
            description='maximum E/W wind speed'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='umin'
            description='minimum E/W wind speed'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='sumax'
            description='maximum surface E/W wind speed'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='sumin'
            description='minimum surface E/W wind speed'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_v.eq.1) then
            i=i+1
            varname='vmax'
            description='maximum N/S wind speed'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='vmin'
            description='minimum N/S wind speed'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='svmax'
            description='maximum surface N/S wind speed'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='svmin'
            description='minimum surface N/S wind speed'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_rmw.eq.1) then
            i=i+1
            varname='rmw'
            description='radius of maximum V'
            units='m'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_pipert.eq.1) then
            i=i+1
            varname='ppimax'
            description='maximum pi perturbation'
            units='dimensionless'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='ppimin'
            description='minimum pi perturbation'
            units='dimensionless'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif

      if(stat_prspert.eq.1) then
            i=i+1
            varname='ppmax'
            description='maximum pressure perturbation'
            units='Pa'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='ppmin'
            description='minimum pressure perturbation'
            units='Pa'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_thpert.eq.1) then
            i=i+1
            varname='thpmax'
            description='maximum potential temperature perturbation'
            units='K'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='thpmin'
            description='minimum potential temperature perturbation'
            units='K'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='sthpmax'
            description='maximum surface potential temperature perturbation'
            units='K'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='sthpmin'
            description='minimum surface potential temperature perturbation'
            units='K'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_q.eq.1) then
            do n=1,numq
                  text1='max     '
                  text2='max                           '
                  write(text1(4:6),156) qname(n)
                  write(text2(5:7),156) qname(n)
156 format(a3)
                  i=i+1
                  varname=text1
                  description=text2
                  units='kg/kg'
                  call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                        memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
 
                  text1='min     '
                  text2='min                           '
                  write(text1(4:6),156) qname(n)
                  write(text2(5:7),156) qname(n)
                  i=i+1
                  varname=text1
                  description=text2
                  units='kg/kg'
                  call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                        memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            enddo
      endif

      if(stat_tke.eq.1) then
            i=i+1
            varname='tkemax'
            description='maximum turbulent kinetic energy'
            units='m^2/s^2'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='tkemin'
            description='minimum turbulent kinetic energy'
            units='m^2/s^2'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_km.eq.1) then
            i=i+1
            varname='kmhmax'
            description='maximum horizontal eddy mixing coefficient for momentum'
            units='m^2/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='kmhmin'
            description='minimum horizontal eddy mixing coefficient for momentum'
            units='m^2/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='kmvmax'
            description='maximum vertical eddy mixing coefficient for momentum'
            units='m^2/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='kmvmin'
            description='minimum vertical eddy mixing coefficient for momentum'
            units='m^2/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_kh.eq.1) then
            i=i+1
            varname='khhmax'
            description='maximum horizontal eddy mixing coefficient for scalars'
            units='m^2/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='khhmin'
            description='minimum horizontal eddy mixing coefficient for scalars'
            units='m^2/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='khvmax'
            description='maximum vertical eddy mixing coefficient for scalars'
            units='m^2/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='khvmin'
            description='minimum vertical eddy mixing coefficient for scalars'
            units='m^2/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_div.eq.1) then
            i=i+1
            varname='divmax'
            description='maximum 3D divergence'
            units='s^{-1}'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='divmin'
            description='minimum 3D divergence'
            units='s^{-1}'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_rh.eq.1) then
            i=i+1
            varname='rhmax'
            description='maximum relative humidity'
            units='dimensionless'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='rhmin'
            description='minimum relative humidity'
            units='dimensionless'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_rhi.eq.1) then
            i=i+1
            varname='rhimax'
            description='maximum relative humidity wrt ice'
            units='dimensionless'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='rhimin'
            description='minimum relative humidity wrt ice'
            units='dimensionless'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if (iptra.eq.1) then
            do n=1,npt
                  text1='maxpt   '
                  text2='max pt                        '
                  write(text1(6:6),157) n
                  write(text2(7:7),157) n
157 format(i1)
                  i=i+1
                  varname=text1
                  description=text2
                  units='kg/kg'
                  call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                        memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
 
                  text1='minpt   '
                  text2='min pt                        '
                  write(text1(6:6),157) n
                  write(text2(7:7),157) n
                  i=i+1
                  varname=text1
                  description=text2
                  units='kg/kg'
                  call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                        memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            enddo
      endif
 
      if(stat_the.eq.1) then
            i=i+1
            varname='themax'
            description='maximum theta-e below 10 km'
            units='K'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='themin'
            description='minimum theta-e below 10 km'
            units='K'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='sthemax'
            description='maximum surface theta-e'
            units='K'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='sthemin'
            description='minimum surface theta-e'
            units='K'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_cloud.eq.1) then
            i=i+1
            varname='qctop'
            description='maximum cloud top height'
            units='m'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='qcbot'
            description='minimum cloud base height'
            units='m'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_sfcprs.eq.1) then
            i=i+1
            varname='sprsmax'
            description='maximum surface pressure'
            units='Pa'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='sprsmin'
            description='minumum surface pressure'
            units='Pa'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_wsp.eq.1) then
            i=i+1
            varname='wspmax'
            description='maximum wind speed'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='wspmin'
            description='minumum wind speed'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='swspmax'
            description='maximum surface wind speed'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='swspmin'
            description='minumum surface wind speed'
            units='m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            if(idrag.eq.1)then
                  i=i+1
                  varname='wsp10max'
                  description='maximum 10 m wind speed'
                  units='m/s'
                  call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                        memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
                  i=i+1
                  varname='wsp10min'
                  description='minimum 10 m wind speed'
                  units='m/s'
                  call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                        memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            endif
      endif
 
      if(stat_cfl.eq.1) then
            if(adapt_dt.eq.1) then
                  i=i+1
                  varname='cflmax'
                  description='maximum Courant number (average - model ran with adaptive time stepping)'
                  units='dimensionless'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            else
                  i=i+1
                  varname='cflmax'
                  description='maximum Courant number'
                  units='dimensionless'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            endif
 
            i=i+1
            varname='kshmax'
            description='maximum horizontal K stability factor'
            units='dimensionless'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
 
            i=i+1
            varname='ksvmax'
            description='maximum vertical K stability factor'
            units='dimensionless'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if (stat_vort.eq.1) then
            i=i+1
            varname='vortsfc'
            description='maximum surface vertical vorticity'
            units='s^{-1}'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
 
            i=i+1
            varname='vort1km'
            description='maximum 1 km vertical vorticity'
            units='s^{-1}'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
 
            i=i+1
            varname='vort2km'
            description='maximum 2 km vertical vorticity'
            units='s^{-1}'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
 
            i=i+1
            varname='vort3km'
            description='maximum 3 km vertical vorticity'
            units='s^{-1}'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
 
            i=i+1
            varname='vort4km'
            description='maximum 4 km vertical vorticity'
            units='s^{-1}'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
 
            i=i+1
            varname='vort5km'
            description='maximum 5 km vertical vorticity'
            units='s^{-1}'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_tmass.eq.1) then
            i=i+1
            varname='tmass'
            description='total mass of dry air'
            units='kg'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
      if(stat_tmois.eq.1) then
            i=i+1
            varname='tmois'
            description='total mass of water vapor in air (doublecheck)'
            units='kg'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
! Note, tmois seems to be same as massqv 
      if(stat_qmass.eq.1)then
            do n=1,numq
                  if( (n.eq.nqv) .or.                                 &
                      (n.ge.nql1.and.n.le.nql2) .or.                  &
                      (n.ge.nqs1.and.n.le.nqs2.and.iice.eq.1) )then
                            text1='mass    '
                            text2='total mass of                 '
                            write(text1( 5: 7),156) qname(n)
                            write(text2(15:17),156) qname(n)
 
                            i=i+1
                            varname=text1
                            description=text2
                            units='kg'
                            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
                   endif
             enddo
       endif
 
      if(stat_tenerg.eq.1) then
            i=i+1
            varname='ek'
            description='total kinetic energy'
            units='J'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='ei'
            description='total internal energy'
            units='J'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='ep'
            description='total potential energy'
            units='J'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='le'
            description='total latent energy'
            units='J'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='et'
            description='total energy'
            units='J'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_mo.eq.1) then
            i=i+1
            varname='tmu'
            description='total E-W momentum '
            units='kg m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='tmv'
            description='total N-S momentum '
            units='kg m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='tmw'
            description='total vertical momentum '
            units='kg m/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
 
      if(stat_tmf.eq.1) then
            i=i+1
            varname='tmfu'
            description='total upward mass flux'
            units='kg/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            i=i+1
            varname='tmfd'
            description='total downward mass flux'
            units='kg/s'
            call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                  memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
      endif
      if(stat_pcn.eq.1)then
            do n=1,nbudget
                  text1='        '
                  text2='                              '
                  write(text1(1:6),158) budname(n)
                  write(text2(1:6),158) budname(n)
158 format(a6)
                  i=i+1
                  varname=text1
                  description=trim(text2)//' (microphysics water budget in qbudget array - try `grep qbudget solve.F`)'
                  units='kg'
                  call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                        memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            enddo
      endif
 
      IF(stat_qsrc.eq.1)THEN
            do n=1,numq
                  text1='as      '
                  text2='artificial source of          '
                  write(text1( 3: 5),156) qname(n)
                  write(text2(22:24),156) qname(n)
                  i=i+1
                  varname=text1
                  description=text2
                  units='kg'
                  call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                        memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            enddo
            do n=1,numq
                  text1='bs      '
                  text2='bndry source/sink of          '
                  write(text1( 3: 5),156) qname(n)
                  write(text2(22:24),156) qname(n)
                  i=i+1
                  varname=text1
                  description=text2
                  units='kg'
                  call writefloat_hdf_unlimited(rank,nvals,nrec,stats_group_id,varname,description,units,var_id,filespace, &
                        memoryspace,dset_id,chunk_id,chunkdims,maxdims,currentdims,singledims,offset,rstat(i))
            enddo
      endif

! Note: i should now equal stat_out (see procXXXXX.out file) or something is wrong
      if (i .ne. stat_out) then
            print *, 'Number of stats written in writestat_hdf5: ',i 
            print *, 'Number of stats caculated in param: ',stat_out
            print *, 'These two must be equal, something is wrong, likely in writestat_hdf5'
            call h5close_f (ierror)
            call stopcm1
      endif

      nrec = nrec + 1
      call h5close_f (ierror)
#endif
! above endif matches #ifdef HDFOUT
      return
      end 
