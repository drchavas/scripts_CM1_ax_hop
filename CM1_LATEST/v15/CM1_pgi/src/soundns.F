

      subroutine soundns(xh,rxh,uh,xf,uf,yh,vh,yf,vf,zh,mh,mf,pi0,thv0,       &
                         radbcw,radbce,radbcs,radbcn,                &
                         divx,u0,ua,u3d,uten,v0,va,v3d,vten,wa,w3d,wten,   &
                         ppi,pp3d,ppten,thv,ppterm,dttmp,nrk,rtime,  &
                         reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,         &
                         uw31,uw32,ue31,ue32,us31,us32,un31,un32,    &
                         vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,    &
                         ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,    &
                         sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,    &
                         pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie) :: xh,rxh,uh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je) :: yh,vh
      real, dimension(jb:je+1) :: yf,vf
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0
      real, dimension(jb:je,kb:ke) :: radbcw,radbce
      real, dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, dimension(ib:ie,jb:je,kb:ke) :: divx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u0,ua,u3d,uten
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v0,va,v3d,vten
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa,w3d,wten
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppten
      real, dimension(ib:ie,jb:je,kb:ke) :: thv,ppterm
      real dttmp
      integer nrk
      real, intent(in) :: rtime
      integer, dimension(rmp) :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p
      real, dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, dimension(imp,cmp,kmp)   :: vs31,vs32,vn31,vn32
      real, dimension(cmp,jmp,kmp-1) :: ww31,ww32,we31,we32
      real, dimension(imp,cmp,kmp-1) :: ws31,ws32,wn31,wn32
      real, dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, dimension(imp,kmp) :: ps1,ps2,pn1,pn2

!-----

      integer :: i,j,k
      real :: tem,tem2

#ifdef MPI
      include 'mpif.h'
#endif

!---------------------------------------------------------------------

        if(irbc.eq.2)then
 
          if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,ua)
 
          if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,va)
 
        endif

!-----

        if(wbc.eq.2.and.ibw.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            u3d(1,j,k)=ua(1,j,k)-dttmp*radbcw(j,k)       &
                      *(ua(2,j,k)-ua(1,j,k))*rdx*uh(1)   &
                         +dttmp*uten(1,j,k)
          enddo
          enddo
        endif

        if(ebc.eq.2.and.ibe.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            u3d(ni+1,j,k)=ua(ni+1,j,k)-dttmp*radbce(j,k)           &
                         *(ua(ni+1,j,k)-ua(ni  ,j,k))*rdx*uh(ni)   &
                         +dttmp*uten(ni+1,j,k)
          enddo
          enddo
        endif

        tem=dttmp*rdx*cp*0.5

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1+ibw,ni+1-ibe
          u3d(i,j,k)=ua(i,j,k)+dttmp*uten(i,j,k)             &
                  -(tem*(ppi(i,j,k)-ppi(i-1,j,k))*uf(i)*     &
                        (thv(i,j,k)+thv(i-1,j,k)))
        enddo
        enddo
        enddo

!----------------------------------------------
!  convergence forcing:

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. nx.gt.1 )THEN
            call convinitu(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibw,ibe,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xf,yh,zh,u0,u3d)
          ENDIF
        ENDIF

!----------------------------------------------

        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcu(u3d)

#ifdef MPI
        call comm_3u_start(u3d,uw31,uw32,ue31,ue32,   &
                               us31,us32,un31,un32,reqs_u)
#endif

!-----

      IF(axisymm.eq.0)THEN

        if(sbc.eq.2.and.ibs.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            v3d(i,1,k)=va(i,1,k)-dttmp*radbcs(i,k)       &
                      *(va(i,2,k)-va(i,1,k))*rdy*vh(1)   &
                      +dttmp*vten(i,1,k)
          enddo
          enddo
        endif
 
        if(nbc.eq.2.and.ibn.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            v3d(i,nj+1,k)=va(i,nj+1,k)-dttmp*radbcn(i,k)           &
                         *(va(i,nj+1,k)-va(i,nj  ,k))*rdy*vh(nj)   &
                         +dttmp*vten(i,nj+1,k)
          enddo
          enddo
        endif

        tem=dttmp*rdy*cp*0.5

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1+ibs,nj+1-ibn
        do i=1,ni
          v3d(i,j,k)=va(i,j,k)+dttmp*vten(i,j,k)             &
                  -(tem*(ppi(i,j,k)-ppi(i,j-1,k))*vf(j)*     &
                        (thv(i,j,k)+thv(i,j-1,k)))
        enddo
        enddo
        enddo

!----------------------------------------------
!  convergence forcing:

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. ny.gt.1 )THEN
            call convinitv(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibs,ibn,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xh,yf,zh,v0,v3d)
          ENDIF
        ENDIF

!----------------------------------------------

        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcv(v3d)

#ifdef MPI
        call comm_3v_start(v3d,vw31,vw32,ve31,ve32,   &
                               vs31,vs32,vn31,vn32,reqs_v)
#endif

      ENDIF

!-----

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          tem=dttmp*rdz*mf(i,j,k)*cp*0.5
          w3d(i,j,k)=wa(i,j,k)+dttmp*wten(i,j,k)             &
                  -(tem*(ppi(i,j,k)-ppi(i,j,k-1))*    &
                        (thv(i,j,k)+thv(i,j,k-1)))
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcw(w3d,1)

#ifdef MPI
        call comm_3w_start(w3d,ww31,ww32,we31,we32,   &
                               ws31,ws32,wn31,wn32,reqs_w)
#endif

!-----

    IF(axisymm.eq.0)THEN
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        divx(i,j,k)=(u3d(i+1,j,k)-u3d(i,j,k))*rdx*uh(i)    &
                   +(v3d(i,j+1,k)-v3d(i,j,k))*rdy*vh(j)    &
                   +(w3d(i,j,k+1)-w3d(i,j,k))*rdz*mh(i,j,k)
      enddo
      enddo
      enddo

    ELSE
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        divx(i,j,k)=(xf(i+1)*u3d(i+1,j,k)-xf(i)*u3d(i,j,k))*rdx*uh(i)*rxh(i)   &
                   +(w3d(i,j,k+1)-w3d(i,j,k))*rdz*mh(i,j,k)
      enddo
      enddo
      enddo

    ENDIF

        tem=dttmp*g*0.5/cp
        tem2=dttmp*rddcv

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          pp3d(i,j,k)=ppi(i,j,k)+dttmp*ppten(i,j,k)              &
              +( tem*(w3d(i,j,k)+w3d(i,j,k+1))/thv0(i,j,k) )         &
              -( (pi0(i,j,k)+ppi(i,j,k))*dttmp*ppterm(i,j,k)*divx(i,j,k) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
        enddo
        enddo
        enddo

        if(timestats.ge.1) time_sound=time_sound+mytime()
 
        call bcs(pp3d)

#ifdef MPI
        if(nrk.lt.3.or.imoist.eq.0)then
          call comm_3s_start(pp3d,sw31,sw32,se31,se32,   &
                                  ss31,ss32,sn31,sn32,reqs_s)
        endif

!--------------------------------------------------------------------
!  Finish comm

        call comm_3u_end(u3d,uw31,uw32,ue31,ue32,   &
                             us31,us32,un31,un32,reqs_u)

!-----

        call comm_3v_end(v3d,vw31,vw32,ve31,ve32,   &
                             vs31,vs32,vn31,vn32,reqs_v)

!-----

        call comm_3w_end(w3d,ww31,ww32,we31,we32,   &
                             ws31,ws32,wn31,wn32,reqs_w)

!-----

        if(nrk.lt.3.or.imoist.eq.0)then
          call comm_3s_end(pp3d,sw31,sw32,se31,se32,   &
                                ss31,ss32,sn31,sn32,reqs_s)
        endif
#endif
!-------------------------------------------------------------------- 

      return
      end


