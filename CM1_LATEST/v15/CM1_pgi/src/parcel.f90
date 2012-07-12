



      subroutine parcel_driver(prec,dt,xh,uh,ruh,yh,vh,rvh,zh,mh,rmh,mf,   &
                               pi0,thv0,th0,the,b,dpdz,thv,qt,prs,      &
                               ua,va,wa,ppi,nm,tha,qa,kh,pdata,rtime,   &
                               ploc,packet,reqs_s,                      &
                               sw1,sw2,se1,se2,ss1,ss2,sn1,sn2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      integer :: prec
      real :: dt
      real, dimension(ib:ie) :: xh,uh,ruh
      real, dimension(jb:je) :: yh,vh,rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,th0
      real, dimension(ib:ie,jb:je,kb:ke) :: the,b,dpdz,thv,qt,prs
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,nm,tha
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kh
      real, dimension(npvals,nparcels) :: pdata
      real :: rtime
      real, dimension(3,nparcels) :: ploc
      real, dimension(npvals+1,nparcels) :: packet
      integer, dimension(rmp) :: reqs_s
      real, dimension(jmp,kmp) :: sw1,sw2,se1,se2
      real, dimension(imp,kmp) :: ss1,ss2,sn1,sn2

      integer :: n,np,i,j,k,iflag,jflag,kflag
      real :: tx,cpm,qvs,tem
      real :: uval,vval,wval,rx,ry,rz,w1,w2,w3,w4,w5,w6,w7,w8
      real :: rslf







!----------------------------------------------------------------------
!  Calculate derived variables
!
!  Note:  n-squared should be on w-points

    call bcs(prs)




    call bcs(nm)







    IF(imoist.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        qt(i,j,k)=0.0
      enddo
      enddo
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        qt(i,j,k)=qt(i,j,k)+qa(i,j,k,nqv)
      enddo
      enddo
      enddo

      do n=nql1,nql2
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
          qt(i,j,k)=qt(i,j,k)+qa(i,j,k,n)
        enddo
        enddo
        enddo
      enddo

      do n=nqs1,nqs2
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
          qt(i,j,k)=qt(i,j,k)+qa(i,j,k,n)
        enddo
        enddo
        enddo
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,tx,cpm,qvs)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        tx=(th0(i,j,k)+tha(i,j,k))*(pi0(i,j,k)+ppi(i,j,k))
        cpm=cp+cpl*qa(i,j,k,nqv)
        do n=nql1,nql2
          cpm=cpm+cpl*qa(i,j,k,n)
        enddo
        qvs=rslf( prs(i,j,k) , tx )
        the(i,j,k)=tx*((p00*(1.0+qa(i,j,k,nqv)*reps)      &
                        /prs(i,j,k))**(rd/cpm))       &
           *((qa(i,j,k,nqv)/qvs)**(-qa(i,j,k,nqv)*rv/cpm))    &
           *exp((lv1-lv2*tx)*qa(i,j,k,nqv)/(cpm*tx))
        thv(i,j,k)=(th0(i,j,k)+tha(i,j,k))*(1.0+reps*qa(i,j,k,nqv))/(1.0+qt(i,j,k))
        b(i,j,k)=g*( thv(i,j,k)/thv0(i,j,k) - 1.0 )
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        the(i,j,k)=th0(i,j,k)+tha(i,j,k)
        thv(i,j,k)=th0(i,j,k)+tha(i,j,k)
        b(i,j,k)=g*( thv(i,j,k)/thv0(i,j,k) - 1.0 )
      enddo
      enddo
      enddo

    ENDIF

    tem=rdz*cp*0.5
!$omp parallel do default(shared)  &
!$omp private(i,j,k,tem)
    do k=2,nk
      do j=0,nj+1
      do i=0,ni+1
        dpdz(i,j,k)=-( tem*(ppi(i,j,k)-ppi(i,j,k-1))*mf(i,j,k)       &
                          *(thv(i,j,k)+thv(i,j,k-1)) )
      enddo
      enddo
    enddo

!$omp parallel do default(shared)  &
!$omp private(i,j)
    do j=0,nj+1
    do i=0,ni+1
      dpdz(i,j,1   )=0.0
      dpdz(i,j,nk+1)=0.0
    enddo
    enddo

    if(timestats.ge.1) time_parcels=time_parcels+mytime()






!----------------------------------------------------------------------
!  get corner info for MPI runs
!  (may not parallelize correctly if this is not done)


      do j=0,nj+1
      do i=0,ni+2
        ua(i,j,0)    = ua(i,j,1)
        ua(i,j,nk+1) = ua(i,j,nk)
      enddo
      enddo

      do j=0,nj+2
      do i=0,ni+1
        va(i,j,0)    = va(i,j,1)
        va(i,j,nk+1) = va(i,j,nk)
      enddo
      enddo

      do j=0,nj+1
      do i=0,ni+1
        wa(i,j,0)    = -wa(i,j,2)
        wa(i,j,nk+2) = -wa(i,j,nk)
      enddo
      enddo

      if(timestats.ge.1) time_parcels=time_parcels+mytime()
      call prepcorners(nm)
      call prepcorners(dpdz)
      if(imoist.eq.1)then
        if(nqv.ne.0) call prepcorners(qa(ib,jb,kb,nqv))
        if(nqc.ne.0) call prepcorners(qa(ib,jb,kb,nqc))
        if(nqr.ne.0) call prepcorners(qa(ib,jb,kb,nqr))
      endif
      call prepcorners(the)
      call prepcorners(b)

!----------------------------------------------------------------------
!
!  Currently, pdata( 1) = x
!                  ( 2) = y
!                  ( 3) = z
!                  ( 4) = qv
!                  ( 5) = qc  ( = qa(2) for kessler and goddard)
!                  ( 6) = qr  ( = qa(3) for kessler and goddard)
!                       (may have to be changed for other microphysics schemes)
!                  ( 7) = n-squared
!                  ( 8) = u
!                  ( 9) = v
!                  (10) = w
!                  (11) = kh
!                  (12) = theta-e
!                  (13) = b
!                  (14) = dpdz
!
!    npvals should be 17 (14 + 3) in param.F
!
!----------------------------------------------------------------------

    DO np=1,nparcels

      pdata(1,np)=pdata(ifx,np)
      pdata(2,np)=pdata(ify,np)
      pdata(3,np)=pdata(ifz,np)

      iflag=0
      jflag=0

      do i=1,ni
        if( abs(xh(i)-pdata(1,np)).le.0.5*dx*ruh(i) ) iflag=i
      enddo

    IF(axisymm.eq.1)THEN
      jflag = 1
    ELSE
      do j=1,nj
        if( abs(yh(j)-pdata(2,np)).le.0.5*dy*rvh(j) ) jflag=j
      enddo
    ENDIF

      IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
          (jflag.ge.1.and.jflag.le.nj) )THEN

        i=iflag
        j=jflag

        do k=1,nk
          if( abs(zh(i,j,k)-pdata(3,np)).le.0.5*dz*rmh(i,j,k) ) kflag=k
        enddo

!----------------------------------------------------------------------
!  Data on u points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(2,np).lt.yh(j) )then
          j=j-1
        endif
        if( pdata(3,np).lt.zh(i,j,k) )then
          k=k-1
        endif

        rx = ( pdata(1,np)-xh(i)+0.5*dx*ruh(i) )*rdx*uh(i)
        ry = ( pdata(2,np)-yh(j) )*rdy*vh(j)
        rz = ( pdata(3,np)-zh(iflag,jflag,k) )*rdz*mh(i,j,k)

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        call tri_interp(ni+1,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,ua,uval)

!----------------------------------------------------------------------
!  Data on v points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(1,np).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(3,np).lt.zh(i,j,k) )then
          k=k-1
        endif

        rx = ( pdata(1,np)-xh(i) )*rdx*uh(i)
        ry = ( pdata(2,np)-yh(j)+0.5*dy*rvh(j) )*rdy*vh(j)
        rz = ( pdata(3,np)-zh(iflag,jflag,k) )*rdz*mh(i,j,k)

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        call tri_interp(ni,nj+1,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,va,vval)

!----------------------------------------------------------------------
!  Data on w points


        i=iflag
        j=jflag
        k=kflag

        if( pdata(1,np).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(2,np).lt.yh(j) )then
          j=j-1
        endif

        rx = ( pdata(1,np)-xh(i) )*rdx*uh(i)
        ry = ( pdata(2,np)-yh(j) )*rdy*vh(j)
        rz = ( pdata(3,np)-(zh(iflag,jflag,k)-0.5*dz*rmh(i,j,k)) )*rdz*mh(i,j,k)

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,wa  ,wval)
        call tri_interp(ni,nj,nk  ,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,nm  ,pdata( 7,np))
        call tri_interp(ni,nj,nk  ,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dpdz,pdata(14,np))

!----------------------------------------------------------------------
!  Data on scalar points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(1,np).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(2,np).lt.yh(j) )then
          j=j-1
        endif
        if( pdata(3,np).lt.zh(i,j,k) )then
          k=k-1
        endif

        rx = ( pdata(1,np)-xh(i) )*rdx*uh(i)
        ry = ( pdata(2,np)-yh(j) )*rdy*vh(j)
        rz = ( pdata(3,np)-zh(iflag,jflag,k) )*rdz*mh(i,j,k)

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

      if(imoist.eq.1)then
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qa(ib,jb,kb,nqv),pdata( 4,np))
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qa(ib,jb,kb,nqc),pdata( 5,np))
        if(ptype.ne.6)then
          call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qa(ib,jb,kb,nqr),pdata( 6,np))
        endif
      endif
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,the,pdata(12,np))
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,b  ,pdata(13,np))

!----------------------------------------------------------------------

        pdata( 8,np)=uval
        pdata( 9,np)=vval
        pdata(10,np)=wval

        pdata(ifx,np)=pdata(1,np)+dt*uval
      IF(axisymm.eq.1)THEN
        pdata(ify,np)=0.0
      ELSE
        pdata(ify,np)=pdata(2,np)+dt*vval
      ENDIF
        pdata(ifz,np)=pdata(3,np)+dt*wval

        if(pdata(ifx,np).lt. 0.0)then
          if(wbc.eq.1)then
            pdata(ifx,np)=pdata(ifx,np)+(maxx-minx)
          else
            pdata(ifx,np)=minx
          endif
        endif
        if(pdata(ifx,np).ge.maxx)then
           if(ebc.eq.1)then
             pdata(ifx,np)=pdata(ifx,np)-(maxx-minx)
           else
             pdata(ifx,np)=maxx
           endif
        endif

        if((pdata(ify,np).ge.maxy).and.(axisymm.ne.1))then
          if(nbc.eq.1)then
            pdata(ify,np)=pdata(ify,np)-(maxy-miny)
          else
            pdata(ify,np)=maxy
          endif
        endif
        if((pdata(ify,np).lt. 0.0).and.(axisymm.ne.1))then
          if(sbc.eq.1)then
            pdata(ify,np)=pdata(ify,np)+(maxy-miny)
          else
            pdata(ify,np)=miny
          endif
        endif

        pdata(ifz,np)=max(pdata(ifz,np),0.0)
        pdata(ifz,np)=min(pdata(ifz,np),maxz)


      ENDIF

    ENDDO

!----------------------------------------------------------------------
!  communicate data


!----------------------------------------------------------------------
!  write out data

    IF( rtime.ge.prcltim .or. prclfrq.lt.0.0 )THEN

      IF(myid.eq.0)THEN

      IF(output_format.eq.1)THEN
        ! GrADS format:

        string(totlen+1:totlen+1+12) = '_pdata.dat  '
        write(outfile,*)
        write(outfile,*) string
        write(outfile,*)
        open(unit=61,file=string,form='unformatted',access='direct',   &
             recl=4,status='unknown')

          do n=1,npvals-3
          do np=1,nparcels
            write(61,rec=prec) pdata(n,np)
            prec=prec+1
          enddo
          enddo

        close(unit=61)

      ELSEIF(output_format.eq.2)THEN

        call writepdata_nc(prec,rtime,pdata)

      ENDIF

      ENDIF

      prcltim = prcltim + prclfrq

    ENDIF

!----------------------------------------------------------------------

      if(timestats.ge.1) time_parcels=time_parcels+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine tri_interp(iz,jz,kz,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,s,pdata)
      implicit none

      include 'input.incl'

      integer :: iz,jz,kz,i,j,k
      real :: w1,w2,w3,w4,w5,w6,w7,w8
      real, dimension(1-ngxy:iz+ngxy,1-ngxy:jz+ngxy,1-ngz:kz+ngz) :: s
      real :: pdata

      pdata=s(i  ,j  ,k  )*w1    &
           +s(i+1,j  ,k  )*w2    &
           +s(i  ,j+1,k  )*w3    &
           +s(i  ,j  ,k+1)*w4    &
           +s(i+1,j  ,k+1)*w5    &
           +s(i  ,j+1,k+1)*w6    &
           +s(i+1,j+1,k  )*w7    &
           +s(i+1,j+1,k+1)*w8

      return
      end


