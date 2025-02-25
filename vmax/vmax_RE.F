	!----------------------------------------------------------------------------------------
	! This program calculates the vmax and rmax needed for RE87 model to give desired actual vmax and rmax
 	! and outputs it to the terminal: first line rmax [m], second line vmax[ms-1]
 	!
	! Created: 22 Mar 2011, Dan Chavas, MIT
	! 
	! Last modified:  23 Mar 2011, Dan Chavas, MIT
	! 
	!----------------------------------------------------------------------------------------

	 PROGRAM vmax_RE

	  implicit none

	!variables
	integer :: nr,nr2
	integer :: i,j,k,i_good
	real :: fcor
	real :: vmax,r0,rmax,rm_min,dr,dd1,dd2,r0drmax
	real :: rref,vr,vr_old
	real :: rm_RE,vm_RE,vmax_real,rmax_real
	real, dimension(:), allocatable :: rm_REs,vm_REs,dvdr


!!!!!!! user input !!!!!!!!!!!!!!!!!!!!!!!!!!!
!	  write(*,"(/,'Desired vmax for initial vortex [ms-1]')")
	  read(*,*) vmax

!	  write(*,"(/,'Desired outer for initial vortex [m]')")
	  read(*,*) r0

!	  write(*,"(/,'Ratio of r0:rmax for initial vortex [-]')")
	  read(*,*) r0drmax

!	  write(*,"(/,'Coriolis parameter [s-1]')")
	  read(*,*) fcor

!        vmax=12.5	!ms-1
!        r0=100000	!m
!	fcor=.00005	!s-1
!	r0drmax=5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        rmax=r0/r0drmax !assume rmax scales with r0

	dr=100	!m
	rm_min=5000	!m
	nr=(r0/2-rm_min)/dr+1

	allocate(rm_REs(nr))
	rm_REs(:)=0

	DO i=1,nr
		rm_REs(i)=rm_min+dr*(i-1)
	ENDDO

    !!! Solve for RE87 input vm_RE and rm_RE that will give desired true vmax and rmax

     !! Step 1: vm_RE as a function of rm_REs; from RE87: v(rmax)=vmax
	allocate(vm_REs(nr))
	vm_REs(:)=0

	DO i=1,nr        
		dd1=2*rm_REs(i)/(rmax+rm_REs(i))
		dd2=2*rm_REs(i)/(r0+rm_REs(i))
	        vm_REs(i)=sqrt((((vmax+.5*fcor*rmax)**2)-.25*(fcor**2)*(rmax**2))/(((rmax/rm_REs(i))**2)*((dd1**3)-(dd2**3))))
	ENDDO

     !! Step 2: dv/dr = 0 at r=rmax; analytic dvdr eqn found in matlab
	allocate(dvdr(nr))
	dvdr(:)=0

        rref=rmax
	DO i=1,nr
		dvdr(i) = - fcor/2 - ((24*rm_REs(i)*(rref**2)*(vm_REs(i)**2))/((rm_REs(i) + rref)**4) &
			- ((fcor**2)*rref)/2 + (2*rref*(vm_REs(i)**2)*((8*(rm_REs(i)**3))/((r0 + rm_REs(i))**3) &
			- (8*(rm_REs(i)**3))/((rm_REs(i) + rref)**3)))/(rm_REs(i)**2))/(2*(((fcor**2)*(rref**2))/4 &
			- ((rref**2)*(vm_REs(i)**2)*((8*(rm_REs(i)**3))/((r0 + rm_REs(i))**3) &
			- (8*(rm_REs(i)**3))/((rm_REs(i) + rref)**3)))/(rm_REs(i)**2))**(.5))
	ENDDO

     !! Step 3: desired value occurs when dv/dr = 0

	i=1
	DO WHILE(dvdr(i)<0)
		i=i+1
	ENDDO

        i_good=i
        rm_RE=rm_REs(i_good)
        vm_RE=vm_REs(i_good)

     !! Step 4: double check that it worked by plugging back into original RE87 equation
        dd2=2*rm_RE/(r0+rm_RE)
	
	vr_old=-1	!less than vr and 0
	vr=0
	i=0
	rref=0
	DO WHILE(vr>vr_old)	!v increases moving out from center til reaches rmax
		i=i+1
		vr_old=vr
		rref=dr*i
	        dd1=2*rm_RE/(rref+rm_RE)
	        vr=sqrt((vm_RE**2)*((rref/rm_RE)**2)*((dd1**3)-(dd2**3))+.25*(fcor**2)*(rref**2))-.5*fcor*rref
	ENDDO

        vmax_real=vr_old
        rmax_real=rref

!	write(*,*)rmax,vmax,rm_RE,vm_RE,rmax_real,vmax_real
	write(*,*)rm_RE
	write(*,*)vm_RE

	END PROGRAM
