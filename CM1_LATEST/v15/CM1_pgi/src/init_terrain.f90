

      subroutine init_terrain(xh,xf,yh,yf,sigma,sigmaf,   &
                              zh,mh,rmh,zf,mf,rmf,zs,gz,dzdx,dzdy,gx,gy)
      implicit none

      include 'input.incl'
      include 'constants.incl'




      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(inout), dimension(kb:ke) :: sigma
      real, intent(inout), dimension(kb:ke+1) :: sigmaf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf,rmf
      real, intent(inout), dimension(itb:ite,jtb:jte) :: zs,gz,dzdx,dzdy
      real, intent(inout), dimension(itb:ite+1,jtb:jte,ktb:kte) :: gx
      real, intent(inout), dimension(itb:ite,jtb:jte+1,ktb:kte) :: gy

      integer :: i,j,k,irec
      real :: hh,aa,xval,xc


!-----------------------------------------------------------------------
!     SPECIFY TERRAIN HERE
!-----------------------------------------------------------------------

!----------------------------------------------------------
!  itern = 1
!  bell-shaped

        IF(itern.eq.1)THEN

          hh =      400.0              ! max. height (m)
          aa =     1000.0              ! half width (m)
          xc =        0.0 + 0.5*dx     ! x-location (m)

          do j=jb,je
          do i=ib,ie
            zs(i,j)=hh/( 1.0+( (xh(i)-xc)/aa )**2 )
          enddo
          enddo

!---------------
!  itern = 2
!  Schaer case

        ELSEIF(itern.eq.2)THEN

          do j=jb,je
          do i=ib,ie
            xval=dx*(i-ni/2)
            zs(i,j)=250.0*exp(-(xval/5000.0)**2)*(cos(pi*xval/4000.0)**2)
          enddo
          enddo

!---------------

        ELSEIF(itern.eq.3)THEN

          hh =      500.0     ! max. height (m)
          aa =    20000.0     ! half width (m)

          do j=jb,je
          do i=ib,ie
            xval = sqrt( (xh(i)-129000.0)**2   &
                        +(yh(j)-129000.0)**2   &
                                             )
            zs(i,j)=hh*( (1.0+(xval/aa)**2 )**(-1.5) )
          enddo
          enddo

!----------------------------------------------------------
!  itern = 4
!  read from GrADS file "perts.dat"

        ELSEIF(itern.eq.4)THEN

          open(unit=73,file='perts.dat',status='old',   &
               form='unformatted',access='direct',recl=4)

          do j=1,nj
          do i=1,ni
            irec=(myj-1)*nx*nj   &
                +(j-1)*nx        &
                +(myi-1)*ni      &
                +i
            read(73,rec=irec) zs(i,j)
          enddo
          enddo

          close(unit=73)

!----------------------------------------------------------

        ENDIF

!--------------------------------------------------------------
!  DO NOT CHANGE ANYTHING BELOW HERE !
!--------------------------------------------------------------

        call bc2d(zs)

        zt = maxz

        if(stretch_z.ne.1)then

          do k=kb,ke+1
            sigmaf(k)=dz*(k-1)
          enddo

        endif

        write(outfile,*)
        do k=1,nk+1
          write(outfile,*) '  sigmaf:',k,sigmaf(k)
        enddo
        write(outfile,*)

        do k=kb,ke
          sigma(k)=0.5*(sigmaf(k)+sigmaf(k+1))
        enddo

        do k=1,nk
        do j=1,nj
        do i=1,ni
          zh(i,j,k)=zs(i,j)+sigma(k)*(zt-zs(i,j))/zt
        enddo
        enddo
        enddo

        do k=kb,ke+1
        do j=jb,je
        do i=ib,ie
          zf(i,j,k)=zs(i,j)+sigmaf(k)*(zt-zs(i,j))/zt
        enddo
        enddo
        enddo

        do j=1,nj
        do i=1,ni
          gz(i,j)=zt/(zt-zs(i,j))
        enddo
        enddo

        do j=1,nj
        do i=1,ni
          dzdx(i,j)=( 45.0*( zs(i+1,j)-zs(i-1,j) )                &
                      -9.0*( zs(i+2,j)-zs(i-2,j) )                &
                          +( zs(i+3,j)-zs(i-3,j) ) )/(60.0*dx)
          dzdy(i,j)=( 45.0*( zs(i,j+1)-zs(i,j-1) )                &
                      -9.0*( zs(i,j+2)-zs(i,j-2) )                &
                          +( zs(i,j+3)-zs(i,j-3) ) )/(60.0*dx)
        enddo
        enddo

!--------------------------------
!  set boundary points

        call bc2d(gz)
        call bc2d(dzdx)
        call bc2d(dzdy)
        call bcs(zh)


!--------------------------------

        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          gx(i,j,k)=(zs(i,j)-zs(i-1,j))*rdx*(sigma(k)-zt)    &
                   /(zt-0.5*(zs(i-1,j)+zs(i,j)))
        enddo
        enddo
        enddo

        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          gy(i,j,k)=(zs(i,j)-zs(i,j-1))*rdx*(sigma(k)-zt)    &
                   /(zt-0.5*(zs(i,j-1)+zs(i,j)))
        enddo
        enddo
        enddo

        call bcu(gx)
        call bcv(gy)


        do j=jb,je
        do i=ib,ie+1
          gx(i,j, 0)=gx(i,j,1)
          gx(i,j,nk+1)=gx(i,j,nk  )
        enddo
        enddo

        do j=jb,je+1
        do i=ib,ie
          gy(i,j, 0)=gy(i,j,1)
          gy(i,j,nk+1)=gy(i,j,nk  )
        enddo
        enddo

!--------------------------------

        do j=jb,je
        do i=ib,ie
          zf(i,j,0)=zf(i,j,1)-(zf(i,j,2)-zf(i,j,1))
          zf(i,j,nk+2)=zf(i,j,nk+1)+(zf(i,j,nk+1)-zf(i,j,nk))
          zh(i,j,0)=0.5*(zf(i,j,0)+zf(i,j,1))
          zh(i,j,nk+1)=0.5*(zf(i,j,nk+1)+zf(i,j,nk+2))
        enddo
        enddo

        write(outfile,*)
        do i=ib,ie
          write(outfile,*) '  zs at nj/2:',i,zs(i,nj/2)
        enddo
        write(outfile,*)

        write(outfile,*)
        do j=jb,je
          write(outfile,*) '  zs at ni/2:',j,zs(ni/2,j)
        enddo
        write(outfile,*)

!---------------------------------------

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          mh(i,j,k)=mh(i,j,k)*gz(i,j)
          rmh(i,j,k)=1.0/mh(i,j,k)
        enddo
        enddo
        enddo

        do k=kb,ke+1
        do j=jb,je
        do i=ib,ie
          mf(i,j,k)=mf(i,j,k)*gz(i,j)
          rmf(i,j,k)=1.0/mf(i,j,k)
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------

      end subroutine init_terrain
