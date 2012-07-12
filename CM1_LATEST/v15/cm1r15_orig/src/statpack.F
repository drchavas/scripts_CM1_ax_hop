

      subroutine statpack(nrec,ndt,dt,rtime,adt,acfl,cloudvar,qname,budname,qbudget,asq,bsq, &
                          xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,vf,zh,mh,rmh,mf,    &
                          rstat,pi0,rho0,thv0,th0,qv0,u0,v0,                &
                          dum1,dum2,dum3,dum4,divx,rho,prs,                 &
                          ua,va,wa,ppi,tha,qa,vq,kmh,kmv,khh,khv,tkea,pta,u10,v10)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      integer :: nrec,ndt
      real :: dt,rtime
      real*8 :: adt,acfl
      logical, dimension(maxq) :: cloudvar
      character*3, dimension(maxq) :: qname
      character*6, dimension(maxq) :: budname
      real*8, dimension(nbudget) :: qbudget
      real*8, dimension(numq) :: asq,bsq
      real, dimension(ib:ie) :: xh,rxh,uh,ruh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je) :: yh,vh,rvh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(stat_out) :: rstat
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,thv0,th0,qv0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,divx,rho,prs
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u0,ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v0,va
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa,vq
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      real, intent(in), dimension(ibl:iel,jbl:jel) :: u10,v10

!-----------------------------------------------------------------------

      integer i,j,k,n,nstat
      character*6 :: text1,text2
      real qvs
      real rslf,rsif
#ifdef MPI
      real*8 var
#endif

!-----------------------------------------------------------------------

      nstat = 0

    IF( adapt_dt.eq.1 )THEN
      nstat = 1
      rstat(nstat) = sngl(  adt/float(max(1,ndt)) )
      acfl         = sngl( acfl/float(max(1,ndt)) )
      ndt  = 0
      adt  = 0.0
    ENDIF

      if(stat_w.eq.1) call maxmin(ni,nj,nk+1,wa,nstat,rstat,'WMAX  ','WMIN  ')
      if(stat_u.eq.1)then
        call maxmin(ni+1,nj,nk,ua,nstat,rstat,'UMAX  ','UMIN  ')
        call maxmin2d(ni+1,nj,ua(ib,jb,1),nstat,rstat,'SUMAX ','SUMIN ')
      endif
      if(stat_v.eq.1)then
        call maxmin(ni,nj+1,nk,va,nstat,rstat,'VMAX  ','VMIN  ')
!!!      print *,'  umax:',rstat(nstat)+rstat(nstat-1),rstat(nstat-4)+rstat(nstat-5),rstat(nstat-1)-rstat(nstat-5)
        call maxmin2d(ni,nj+1,va(ib,jb,1),nstat,rstat,'SVMAX ','SVMIN ')
      endif
      if(stat_rmw.eq.1)then
        call getrmw(nstat,rstat,xh,va)
      endif
 
      if(stat_pipert.eq.1) call maxmin(ni,nj,nk,ppi,nstat,rstat,'PPIMAX','PPIMIN')

      if(stat_prspert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k)=prs(i,j,k)-p00*(pi0(i,j,k)**cpdrd)
        enddo
        enddo
        enddo
        call maxmin(ni,nj,nk,dum2,nstat,rstat,'PPMAX ','PPMIN ')
      endif

      if(stat_thpert.eq.1)then
        call maxmin(ni,nj,nk,tha,nstat,rstat,'THPMAX','THPMIN')
        call maxmin2d(ni,nj,tha(ib,jb,1),nstat,rstat,'STHPMX','STHPMN')
      endif

      if(imoist.eq.1.and.stat_q.eq.1)then
        do n=1,numq
          text1='MAX   '
          text2='MIN   '
          write(text1(4:6),121) qname(n)
          write(text2(4:6),121) qname(n)
121       format(a3)
          call maxmin(ni,nj,nk,qa(ib,jb,kb,n),nstat,rstat,text1,text2)
        enddo
      endif

      if(iturb.eq.1)then
        if(stat_tke.eq.1) call maxmin(ni,nj,nkt,tkea,nstat,rstat,'TKEMAX','TKEMIN')
      endif

      if(iturb.ge.1)then
        if(stat_km.eq.1) call maxmin(ni,nj,nkt,kmh,nstat,rstat,'KMHMAX','KMHMIN')
        if(stat_km.eq.1) call maxmin(ni,nj,nkt,kmv,nstat,rstat,'KMVMAX','KMVMIN')
        if(stat_kh.eq.1) call maxmin(ni,nj,nkt,khh,nstat,rstat,'KHHMAX','KHHMIN')
        if(stat_kh.eq.1) call maxmin(ni,nj,nkt,khv,nstat,rstat,'KHVMAX','KHVMIN')
      endif

      if(stat_div.eq.1)then
      IF(axisymm.eq.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          divx(i,j,k)=                                                     &
              0.5*( (rho0(i,j,k)+rho0(i+1,j,k))*ua(i+1,j,k)                &
                   -(rho0(i,j,k)+rho0(i-1,j,k))*ua(i  ,j,k) )*rdx*uh(i)    &
             +0.5*( (rho0(i,j,k)+rho0(i,j+1,k))*va(i,j+1,k)                &
                   -(rho0(i,j,k)+rho0(i,j-1,k))*va(i,j  ,k) )*rdy*vh(j)    &
             +0.5*( (rho0(i,j,k)+rho0(i,j,k+1))*wa(i,j,k+1)                &
                   -(rho0(i,j,k)+rho0(i,j,k-1))*wa(i,j,k  ) )*rdz*mh(i,j,k)
        enddo
        enddo
        enddo
      ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          divx(i,j,k)=                                                     &
              rho0(1,1,k)*( xf(i+1)*ua(i+1,j,k)                            &
                           -xf(i  )*ua(i  ,j,k) )*rdx*uh(i)*rxh(i)         &
             +0.5*( (rho0(i,j,k)+rho0(i,j,k+1))*wa(i,j,k+1)                &
                   -(rho0(i,j,k)+rho0(i,j,k-1))*wa(i,j,k  ) )*rdz*mh(i,j,k)
        enddo
        enddo
        enddo
      ENDIF
        call maxmin(ni,nj,nk,divx,nstat,rstat,'DIVMAX','DIVMIN')
      endif

      IF(imoist.eq.1)THEN

        if(stat_rh.eq.1 .or. stat_the.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k,qvs)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qvs=rslf( prs(i,j,k) , (th0(i,j,k)+tha(i,j,k))*(pi0(i,j,k)+ppi(i,j,k)) )
            dum2(i,j,k)=qa(i,j,k,nqv)*(1.0+qvs*reps)    &
                       /(qvs*(1.0+qa(i,j,k,nqv)*reps))
          enddo
          enddo
          enddo
        endif

        if(stat_rh.eq.1)then
          call maxmin(ni,nj,nk,dum2,nstat,rstat,'RHMAX ','RHMIN ')
        endif

        if(iice.eq.1 .and. stat_rhi.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k,qvs)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qvs=rsif( prs(i,j,k) , (th0(i,j,k)+tha(i,j,k))*(pi0(i,j,k)+ppi(i,j,k)) )
            dum3(i,j,k)=qa(i,j,k,nqv)*(1.0+qvs*reps)    &
                       /(qvs*(1.0+qa(i,j,k,nqv)*reps))
          enddo
          enddo
          enddo
          call maxmin(ni,nj,nk,dum3,nstat,rstat,'RHIMAX','RHIMIN')
        endif

      ENDIF

        if(iptra.eq.1)then
          do n=1,npt
            text1='MAXPT '
            text2='MINPT '
            write(text1(6:6),122) n
            write(text2(6:6),122) n
122         format(i1)
            call maxmin(ni,nj,nk,pta(ib,jb,kb,n),nstat,rstat,text1,text2)
          enddo
        endif

      IF(imoist.eq.1)THEN

        if(stat_the.eq.1)then
          call calcthe(zh,pi0,th0,dum4,dum2,prs,ppi,tha,qa)
          call maxmin(ni,nj,nk,dum4,nstat,rstat,'THEMAX','THEMIN')
          call maxmin2d(ni,nj,dum4(ib,jb,1),nstat,rstat,'STHEMX','STHEMN')
        endif

        if(stat_cloud.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k)=0.0
          enddo
          enddo
          enddo
          do n=1,numq
            if(cloudvar(n))then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
              do k=1,nk
              do j=1,nj
              do i=1,ni
                dum1(i,j,k)=dum1(i,j,k)+qa(i,j,k,n)
              enddo
              enddo
              enddo
            endif
          enddo
          call cloud(nstat,rstat,zh,dum1)
        endif
      ENDIF

      if(stat_sfcprs.eq.1)then
        call maxmin2d(ni,nj,prs(ib,jb,1),nstat,rstat,'SFPMAX','SFPMIN')
      endif

      if(stat_wsp.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=sqrt( (umove+0.5*(ua(i,j,k)+ua(i+1,j,k)))**2     &
                           +(vmove+0.5*(va(i,j,k)+va(i,j+1,k)))**2 )
        enddo
        enddo
        enddo
        call maxmin(ni,nj,nk,dum1,nstat,rstat,'WSPMAX','WSPMIN')
        call maxmin2d(ni,nj,dum1(ib,jb,1),nstat,rstat,'SWSPMX','SWSPMN')
      IF(idrag.eq.1)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1)=sqrt( u10(i,j)**2 + v10(i,j)**2 )
        enddo
        enddo
        call maxmin2d(ni,nj,dum1(ib,jb,1),nstat,rstat,'10MWMX','10MWMN')
      ENDIF
      endif

      if(stat_cfl.eq.1) call calccfl(nstat,rstat,dt,acfl,uf,vf,mf,ua,va,wa,1)

      if(stat_cfl.eq.1) call calcksmax(nstat,rstat,dt,uh,vh,mf,kmh,kmv,khh,khv)

      if(stat_vort.eq.1) call vertvort(nstat,rstat,xh,xf,uf,vf,zh,ua,va)

      if(stat_tmass.eq.1) call calcmass(nstat,rstat,ruh,rvh,rmh,rho)

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum1(i,j,k)=0.0
        dum2(i,j,k)=0.0
        dum3(i,j,k)=0.0
      enddo
      enddo
      enddo
 
      IF(imoist.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=qa(i,j,k,nqv)
        enddo
        enddo
        enddo

        call getqli(0,qa,dum2,dum3)

        if(stat_tmois.eq.1)then
          call totmois(nstat,rstat,qbudget(budrain),ruh,rvh,rmh,dum1,dum2,dum3,rho)
        endif

        if(stat_qmass.eq.1)then
          do n=1,numq
            IF( (n.eq.nqv) .or.                                 &
                (n.ge.nql1.and.n.le.nql2) .or.                  &
                (n.ge.nqs1.and.n.le.nqs2.and.iice.eq.1) )THEN
              text1='   MAS'
              write(text1(1:3),121) qname(n)
              call totq(nstat,rstat,ruh,rvh,rmh,qa(ib,jb,kb,n),rho,text1)
            ENDIF
          enddo
        endif

      ENDIF

        if(imoist.eq.1)then
          if(ptype.eq.1.or.ptype.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum4(i,j,k)=vq(i,j,k,3)
            enddo
            enddo
            enddo
          elseif(ptype.eq.6)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum4(i,j,k)=vq(i,j,k,2)
            enddo
            enddo
            enddo
          else
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum4(i,j,k)=0.0
            enddo
            enddo
            enddo
          endif
        else
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum4(i,j,k)=0.0
          enddo
          enddo
          enddo
        endif
 
      if(stat_tenerg.eq.1)then
        call calcener(nstat,rstat,ruh,rvh,zh,rmh,pi0,th0,rho,ua,va,wa,ppi,tha,    &
                      dum1,dum2,dum3,dum4)
      endif

      if(stat_mo.eq.1)then
        call calcmoe(nstat,rstat,ruh,rvh,rmh,rho,ua,va,wa,dum1,dum2,dum3,dum4)
      endif

      if(stat_tmf.eq.1) call tmf(nstat,rstat,ruh,rvh,rho,wa)

!----------

      IF(imoist.eq.1 .and. stat_pcn.eq.1)THEN
#ifdef MPI
        do n=1,nbudget
          var=0.0d0
          call MPI_REDUCE(qbudget(n),var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,   &
                          MPI_COMM_WORLD,ierr)
          qbudget(n)=var
        enddo
        if(myid.eq.0)then
#endif
100     format(2x,a6,':',1x,e13.6)
 
        do n=1,nbudget
          write(6,100) budname(n),qbudget(n)
          nstat = nstat + 1
          rstat(nstat) = qbudget(n)
        enddo
#ifdef MPI
        endif
#endif
      ENDIF

      IF(imoist.eq.1 .and. stat_qsrc.eq.1)THEN
#ifdef MPI
        do n=1,numq
          var=0.0d0
          call MPI_REDUCE(asq(n),var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,   &
                          MPI_COMM_WORLD,ierr)
          asq(n)=var
        enddo
        do n=1,numq
          var=0.0d0
          call MPI_REDUCE(bsq(n),var,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,   &
                          MPI_COMM_WORLD,ierr)
          bsq(n)=var
        enddo
        if(myid.eq.0)then
#endif
        do n=1,numq
          text1='as    '
          write(text1(3:5),121) qname(n)
          write(6,100) text1,asq(n)
          nstat = nstat + 1
          rstat(nstat) = asq(n)
        enddo
        do n=1,numq
          text1='bs    '
          write(text1(3:5),121) qname(n)
          write(6,100) text1,bsq(n)
          nstat = nstat + 1
          rstat(nstat) = bsq(n)
        enddo
#ifdef MPI
        endif
#endif
      ENDIF

  IF(myid.eq.0)THEN

!-----------------------------------------------------------------------
!  writeitout:  GrADS format

    IF(output_format.eq.1)THEN

      open(unit=60,file=statfile,form='unformatted',access='direct',   &
           recl=4,status='unknown')
      if( nstat.ne.stat_out )   stop 12998
      do n=1,nstat
#ifdef DP
        write(60,rec=nrec) sngl(rstat(n))
#else
        write(60,rec=nrec) rstat(n)
#endif
        nrec = nrec + 1
      enddo
      close(unit=60)

!-----------------------------------------------------------------------
!  writeitout:  netcdf format

#ifdef NETCDF
    ELSEIF(output_format.eq.2)THEN

      call writestat_nc(nrec,rtime,nstat,rstat,qname,budname)

#endif

!-----------------------------------------------------------------------
!  writeout:  hdf5 format

#ifdef HDFOUT
    ELSEIF(output_format.ge.3)THEN

      call writestat_hdf5(nrec,rtime,nstat,rstat,qname,budname)

#endif

!-----------------------------------------------------------------------

    ENDIF

  ENDIF

      if(timestats.ge.1) time_stat=time_stat+mytime()

      return
      end


