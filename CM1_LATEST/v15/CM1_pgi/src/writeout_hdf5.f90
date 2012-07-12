

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



      implicit none

      include 'input.incl'
      include 'constants.incl'




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


! above endif matches #ifdef HDFOUT
      return
      end 
