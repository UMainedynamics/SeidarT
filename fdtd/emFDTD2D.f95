! The electromagFDTD2D module is designed to be integrated with python via f2py. 
!
! Compile using
!     f2py3 -c --fcompiler=gnu95 -m emfdtd2d emFDTD2d.f95
!
! Created by Steven Bernsen with T-minus one week to AGU
! University of Maine
! Department of Earth and Environmental Sciences 
! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


module electromagFDTD2d

implicit none
  
contains

  
  !==============================================================================
  subroutine permittivity_write(im, mlist, npoints_pml, nx, nz) 
  ! STIFFNESS_ARRAYS takes a matrix containing the material integer identifiers 
  ! and creates the same size array for each independent coefficient of the 
  ! stiffness matrix along with a density matrix. Since we ae using PML
  ! boundaries, we will extend the the boundary values through the PML region.
  ! 
  ! INPUT 
  !   im (INTEGER)  
  !   mlist (REAL)
  !   epsilonx(i,j), sigmax(i,j), epsilony(i,j), sigmay, (REAL) -
  !   npoints_pml (INTEGER) - the 
  !   
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  implicit none 
  
  integer :: nx,nz
  integer,parameter :: dp = kind(0.d0)
  integer,dimension(nx,nz) :: im
  integer :: i, j, npoints_pml
  real(kind=dp), dimension(:,:) :: mlist
  real(kind=dp), dimension(2*npoints_pml+nx,2*npoints_pml+nz) :: &
            epsilonx, epsilony, epsilonz, sigmax, sigmay, sigmaz
  
  !f2py3 intent(in):: im, mlist, npoints_pml, nx, nz
  
  ! Allocate space for permittivity and conductivity values
  epsilonx(:,:) = 0.d0
  epsilony(:,:) = 0.d0
  epsilonz(:,:) = 0.d0
  sigmax(:,:) = 0.d0
  sigmay(:,:) = 0.d0
  sigmaz(:,:) = 0.d0
  
  do i=npoints_pml+1,nx + npoints_pml
    do j=npoints_pml+1,nz + npoints_pml
      epsilonx(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 2)
      epsilony(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 3)
      epsilonz(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 4)
      
      sigmax(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 5) 
      sigmay(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 6)
      sigmaz(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 7)
    end do
  end do
  
  ! Extend the boundary values of the stiffnesses into the PML region
  do i = 1,npoints_pml+1
    ! top and bottom
    epsilonx( i, : ) = epsilonx(npoints_pml+1,:)
    epsilony( i, : ) = epsilony(npoints_pml+1,:)
    epsilonz( i, : ) = epsilonz(npoints_pml+1,:) 
    sigmax( i, : ) = sigmax(npoints_pml+1,:)
    sigmay( i, : ) = sigmay(npoints_pml+1,:)
    sigmaz( i, : ) = sigmaz(npoints_pml+1,:)
    
    epsilonx( nx+npoints_pml-1+i, : ) = epsilonx(nx+npoints_pml-1,:)
    epsilony( nx+npoints_pml-1+i, : ) = epsilony(nx+npoints_pml-1,:)
    epsilonz( nx+npoints_pml-1+i, : ) = epsilonz(nx+npoints_pml-1,:)
    sigmax( nx+npoints_pml-1+i, : ) = sigmax(nx+npoints_pml-1,:)
    sigmay( nx+npoints_pml-1+i, : ) = sigmay(nx+npoints_pml-1,:)
    sigmaz( nx+npoints_pml-1+i, : ) = sigmaz(nx+npoints_pml-1,:)
  
    ! left and right
    epsilonx( :, i ) = epsilonx(:, npoints_pml+1)
    epsilony( :, i ) = epsilony(:, npoints_pml+1)
    epsilonz( :, i ) = epsilonz(:, npoints_pml+1)
    sigmax( :, i ) = sigmax(:, npoints_pml+1)
    sigmay( :, i ) = sigmay(:, npoints_pml+1)
    sigmaz( :, i ) = sigmaz(:, npoints_pml+1)
  
    epsilonx( :, nz+npoints_pml-1+i ) = epsilonx(:,nz+npoints_pml-1)    
    epsilony( :, nz+npoints_pml-1+i ) = epsilony(:,nz+npoints_pml-1)
    epsilonz( :, nz+npoints_pml-1+i ) = epsilonz(:,nz+npoints_pml-1)
    sigmax( :, nz+npoints_pml-1+i ) = sigmax(:,nz+npoints_pml-1)    
    sigmay( :, nz+npoints_pml-1+i ) = sigmay(:,nz+npoints_pml-1)
    sigmaz( :, nz+npoints_pml-1+i ) = sigmaz(:,nz+npoints_pml-1)
  end do 
  
  ! Write each of the matrices to file
  call material_rw('epsx.dat', epsilonx, .FALSE.)
  call material_rw('epsy.dat', epsilony, .FALSE.)
  call material_rw('epsz.dat', epsilonz, .FALSE.)
  call material_rw('sigx.dat', sigmax, .FALSE.)
  call material_rw('sigy.dat', sigmay, .FALSE.)
  call material_rw('sigz.dat', sigmaz, .FALSE.)
  
end subroutine permittivity_write


! ---------------------------------------------------------------------------
  subroutine material_rw(filename, image_data, readfile)

    implicit none
    
    integer,parameter :: dp = kind(0.d0)
    character(len=8) :: filename
    real(kind=dp),dimension(:,:) :: image_data
    logical :: readfile
    
    
    open(unit = 13, form="unformatted", file = trim(filename))
    
    if ( readfile ) then
      read(13) image_data
    else
      write(13) image_data
    endif
    
    close(unit = 13)
  
  end subroutine material_rw
  


!==============================================================================
  subroutine loadsource(filename, N, srcfn)
    
    implicit none

    integer,parameter :: dp = kind(0.d0)
    character(len=26) :: filename
    integer :: N
    real(kind=dp),dimension(N) :: srcfn
    
    open(unit = 13, form="unformatted", file = trim(filename))
    read(13) srcfn
    
    close(unit = 13)

  end subroutine loadsource

! -----------------------------------------------------------------------------

  subroutine cpml_coeffs(nx, dx, dt, npml, k_max, alpha_max, &
              kappa, alpha, acoeff, bcoeff, HALF)

    implicit none

    integer,parameter :: dp=kind(0.d0)
    integer :: i

    ! Define real inputs 
    real(kind=dp) :: dx, dt, sig_max, k_max, alpha_max 
    integer :: nx, npml
    logical :: HALF

    ! define the output arrays
    real(kind=dp),dimension(nx) :: kappa, alpha, acoeff, bcoeff

    ! Define all other variables needed in the program
    real(kind=dp) :: xoriginleft, xoriginright
    real(kind=dp),dimension(nx) :: xval, sigma
    real(kind=dp),parameter :: PI = 3.141592653589793238462643d0
    real(kind=dp),parameter :: eps0 = 8.85418782d-12, mu0 = 4.0d0*pi*1.0d-7
    integer,parameter :: NP = 2, NPA = 2

    real(kind=dp) :: abscissa_in_pml, abscissa_normalized

    ! ===========================================================
    sig_max = ( 0.8d0 * ( dble(NP+1) ) / ( dx * ( mu0 / eps0 )**0.5d0 ) )

    sigma(:) = 0.d0

    do i=1,nx 
      xval(i) = dx * dble(i - 1)
    enddo

    if (HALF) then 
        xval = xval + dx/2.0
    endif

    xoriginleft = dx * dble( npml )
    xoriginright = dx * dble( (NX-1) - npml )

    do i=1,nx
        !---------- left edge
        abscissa_in_PML = xoriginleft - xval(i)
        if (abscissa_in_PML >= 0.d0) then
            abscissa_normalized = abscissa_in_PML / dble(dx * npml)
            sigma(i) = sig_max * abscissa_normalized**NP
            ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
            kappa(i) = 1.d0 + (K_MAX - 1.d0) * abscissa_normalized**NP
            alpha(i) = ALPHA_MAX * (1.d0 - abscissa_normalized)**NPA
        endif

        !---------- right edge
        ! define damping profile at the grid points
        abscissa_in_PML = xval(i) - xoriginright
        if (abscissa_in_PML >= 0.d0) then
          abscissa_normalized = abscissa_in_PML / dble(dx * npml)
          sigma(i) = sig_max * abscissa_normalized**NP
          kappa(i) = 1.d0 + (k_max - 1.d0) * abscissa_normalized**NP
          alpha(i) = alpha_max * (1.d0 - abscissa_normalized)**NPA
        endif

        ! just in case, for -5 at the end
        if (alpha(i) < 0.d0) alpha(i) = 0.d0
        ! Compute the b_i coefficents
        bcoeff(i) = exp( - (sigma(i) / kappa(i) + alpha(i)) * DT/eps0 )
        
        ! Compute the a_i coefficients
        ! this to avoid division by zero outside the PML
        if (abs(sigma(i)) > 1.d-6) then 
          acoeff(i) = sigma(i) * (bcoeff(i) - 1.d0) / ( (sigma(i) + kappa(i) * alpha(i)) ) / kappa(i)
        endif

    enddo 

  end subroutine cpml_coeffs

!==============================================================================


subroutine electromag_cpml_2d(nx, ny, dx, dy, npoints_pml, src, f0, nstep)

! 2D elastic finite-difference code in velocity and stress formulation
! with Convolutional-PML (C-PML) absorbing conditions for an anisotropic medium

! Dimitri Komatitsch, University of Pau, France, April 2007.
! Anisotropic implementation by Roland Martin and Dimitri Komatitsch, University of Pau, France, April 2007.

! The second-order staggered-grid formulation of Madariaga (1976) and Virieux (1986) is used:
!
!            ^ y
!            |
!            |
!
!            +-------------------+
!            |                   |
!            |                   |
!            |                   |
!            |                   |
!            |        v_y        |
!   sigma_xy +---------+         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            +---------+---------+  ---> x
!           v_x    sigma_xx
!                  sigma_yy
!
! IMPORTANT : all our CPML codes work fine in single precision as well (which is significantly faster).
!             If you want you can thus force automatic conversion to single precision at compile time
!             or change all the declarations and constants in the code from real(kind=dp) to single.
!
! INPUT
!   im (INTEGER)
!   nx, ny (INTEGER)
!   epsilonx, sigmax, epsilony, sigmay (REAL)
!   dx, dy (REAL)
!   npoints_pml (INTEGER) - the thickness of the pml
!   rcx (INTEGER) - the x and y indices for an array of recievers
!


implicit none

integer,parameter :: dp=kind(0.d0)

! total number of grid points in each direction of the grid
integer :: nx
integer :: ny
! integer, dimension(:,:) :: rcx 

! thickness of the PML layer in grid points
integer :: npoints_pml
! integer, dimension(nx,ny)
real(kind=dp), dimension(nx,ny) :: epsilonx, sigmax, epsilony, sigmay

! total number of time steps
integer :: nstep

! time step in seconds. decreasing the time step improves the pml attenuation
! but it should be inversely proportional to the center frequency of the 
! source frequency 
real(kind=dp) :: DT
real(kind=dp) :: dx, dy 

! source
integer,dimension(:) :: src
real(kind=dp),dimension(nstep) :: srcx, srcy
integer :: isource, jsource

! value of PI
real(kind=dp), parameter :: PI = 3.141592653589793238462643d0

! speed of mother fuckin' light 
real(kind=dp), parameter :: Clight = 2.9979458d+8

! permability and permittivity of free space 
real(kind=dp), parameter :: mu0 = 4.0d0*pi*1.0d-7, eps0 = 8.85418782d-12, epsR = 1.0d0

! typical relative permeability of common materials is close to unity but for
! the more specific case we can edit the following line to input permeability 
! as a 2D array 
real(kind=dp), parameter :: mu = 1.d0

! conversion from degrees to radians
real(kind=dp), parameter :: DEGREES_TO_RADIANS = PI / 180.d0

! E-field threshold above which we consider that the code became unstable
real(kind=dp), parameter :: STABILITY_THRESHOLD = 1.d+25

! main arrays
real(kind=dp), dimension(nx-1,ny) :: Ex
real(kind=dp), dimension(nx,ny-1) :: Ey
real(kind=dp), dimension(nx-1,ny-1) :: Hz

! we will compute the coefficients for the finite difference scheme 
real(kind=dp), dimension(nx, ny) :: caEx, cbEx
real(kind=dp), dimension(nx, ny) :: caEy, cbEy
real(kind=dp) :: daHz, dbHz

real(kind=dp) :: value_dEx_dy, value_dEy_dx, value_dHz_dy, value_dHz_dx

! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
real(kind=dp), dimension(nx,ny) ::  memory_dEy_dx,  memory_dEx_dy 
real(kind=dp), dimension(nx,ny) :: memory_dHz_dx
real(kind=dp), dimension(nx,ny) ::  memory_dHz_dy

! parameters for the source
! angle of source force clockwise with respect to vertical (Y) axis
! this will later be treated as an input
real(kind=dp) :: f0
real(kind=dp), parameter :: factor = 1.d0
! character(len=6) :: src_type
integer :: i,j, it

real(kind=dp) :: velocnorm

! -------------------------------- PML parameters 
! power to compute d0 profile. Increasing this value allows for a larger dampening gradient in the PML
real(kind=dp), parameter :: NP = 2.d0, NPA = 1.d0

! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
real(kind=dp) :: k_max!, parameter :: k_max = 1.2d5
real(kind=dp) :: alpha_max !, parameter :: alpha_max = pi*f0 !7.0d-1 ! From Taflove ! PI*f0 ! from Festa and Vilotte 


! 1D arrays for the damping profiles
real(kind=dp), dimension(nx) :: sigh_x,K_x,alpha_x,a_x,b_x
real(kind=dp), dimension(nx) :: sige_x_half,K_x_half, alpha_x_half,a_x_half,b_x_half
real(kind=dp), dimension(ny) :: sigh_y,K_y,alpha_y,a_y,b_y
real(kind=dp), dimension(ny) :: sige_y_half,K_y_half, alpha_y_half,a_y_half,b_y_half

integer :: thickness_PML_x,thickness_PML_y


! ------------------- Name the f2py inputs 
!f2py3 intent(in) :: nx, ny, dx, dy,
!f2py3 intent(in) :: noints_pml, src, f0, nstep

! =============================================================================
! ----------------------- Load Permittivity Coefficients ----------------------

call material_rw('epsx.dat', epsilonx, .TRUE.)
call material_rw('epsz.dat', epsilony, .TRUE.) ! We will change y to z soon
call material_rw('sigx.dat', sigmax, .TRUE.)
call material_rw('sigz.dat', sigmay, .TRUE.)

! ------------------------ Assign some constants -----------------------

! Assign the source location indices
isource = int(src(1)) + npoints_pml
jsource = int(src(2)) + npoints_pml

! nrec = size(rcx,1) 
! allocate( sisEx(nstep, nrec ), &
!   sisEy(nstep, nrec ) )

! Define the 
DT = minval( (/dx, dy/) )/ ( 2.d0 * Clight/sqrt( minval( (/ epsilonx, epsilony /) ) ) ) 

alpha_max = 2*pi*eps0*f0/10
k_max = 1.5d1

! Compute the coefficients of the FD scheme. First scale the relative 
! permittivity and permeabilities to get the absolute values 
epsilonx(:,:) = epsilonx*eps0
epsilony(:,:) = epsilony*eps0


! We need to change sigma to dsigma, same for epsilon

caEx(:,:) = ( 1.0d0 - sigmax * dt / ( 2.0d0 * epsilonx ) ) / &
            ( 1.0d0 + sigmax * dt / (2.0d0 * epsilonx ) )
cbEx(:,:) = (dt / epsilonx ) / ( 1.0d0 + sigmax * dt / ( 2.0d0 * epsilonx ) )

caEy(:,:) = ( 1.0d0 - sigmay * dt / ( 2.0d0 * epsilony ) ) / &
            ( 1.0d0 + sigmay * dt / (2.0d0 * epsilony ) )
cbEy(:,:) = (dt / epsilony ) / ( 1.0d0 + sigmay * dt / ( 2.0d0 * epsilony ) )

daHz = dt/(4.0d0*mu0*mu)
dbHz = dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
daHz = 1.0d0 ! (1-daHz)/(1+daHz) ! 


! ----------------------------------------------------------------------
!---
!--- program starts here
!---

! ================================ LOAD SOURCE ================================

call loadsource('electromagneticsourcex.dat', nstep, srcx)
call loadsource('electromagneticsourcez.dat', nstep, srcy)


! Initialize CPML damping variables
  sigh_x(:) = 0.d0
  sige_x_half(:) = 0.d0
  K_x(:) = 1.d0
  K_x_half(:) = 1.d0
  alpha_x(:) = 0.d0
  alpha_x_half(:) = 0.d0
  a_x(:) = 0.d0
  a_x_half(:) = 0.d0

  sigh_y(:) = 0.d0
  sige_y_half(:) = 0.d0
  K_y(:) = 1.d0
  K_y_half(:) = 1.d0
  alpha_y(:) = 0.d0
  alpha_y_half(:) = 0.d0
  a_y(:) = 0.d0
  a_y_half(:) = 0.d0


! thickness of the PML layer in meters. 
! This can be varied for either axis but for now it constant. 
thickness_PML_x = NPOINTS_PML
thickness_PML_y = NPOINTS_PML

! The source and reciever values are given relative to the image so we
! need to add the PML thickness to the coordinates
src(1) = src(1) + thickness_PML_x 
src(2) = src(2) + thickness_PML_y



call cpml_coeffs(nx, dx, dt, thickness_PML_x, k_max, alpha_max, &
            K_x, alpha_x, a_x, b_x, .FALSE.)

call cpml_coeffs(nx, dx, dt, thickness_PML_x, k_max, alpha_max, &
            K_x_half, alpha_x_half, a_x_half, b_x_half, .TRUE.)

call cpml_coeffs(ny, dy, dt, thickness_PML_y, k_max, alpha_max, &
            K_y, alpha_y, a_y, b_y, .FALSE.)

call cpml_coeffs(ny, dy, dt, thickness_PML_y, k_max, alpha_max, &
            K_y_half, alpha_y_half, a_y_half, b_y_half, .TRUE.)



! initialize arrays
Ex(:,:) = 0.0d0
Ey(:,:) = 0.0d0
Hz(:,:) = 0.0d0

! PML
memory_dEx_dy(:,:) = 0.0d0
memory_dEy_dx(:,:) = 0.0d0

memory_dHz_dx(:,:) = 0.0d0
memory_dHz_dy(:,:) = 0.0d0

!---
!---  beginning of time loop
!---

do it = 1,NSTEP
  
  !--------------------------------------------------------
  ! compute magnetic field and update memory variables for C-PML
  !--------------------------------------------------------
  do i = 1,nx-1  
    do j = 1,ny-1
      
      ! Values needed for the magnetic field updates
      value_dEx_dy = ( Ex(i,j+1) - Ex(i,j) )/dy
      memory_dEx_dy(i,j) = b_y(j) * memory_dEx_dy(i,j) + a_y(j) * value_dEx_dy
      value_dEx_dy = value_dEx_dy/ K_y(j) + memory_dEx_dy(i,j)

      ! The rest of the equation needed for agnetic field updates
      value_dEy_dx = ( Ey(i+1,j) - Ey(i,j) )/dx
      memory_dEy_dx(i,j) = b_x(i) * memory_dEy_dx(i,j) + a_x(i) * value_dEy_dx
      value_dEy_dx = value_dEy_dx/ K_x(i) + memory_dEy_dx(i,j)

      ! Now update the Magnetic field
      Hz(i,j) = daHz*Hz(i,j) + dbHz*( value_dEy_dx + value_dEx_dy )

    enddo  
  enddo

  !--------------------------------------------------------
  ! compute electric field and update memory variables for C-PML
  !--------------------------------------------------------
  
  ! Compute the differences in the y-direction
  do j = 2,ny-1
    do i = 1,nx-1
      ! Update the Ex field
      value_dHz_dy = ( Hz(i,j) - Hz(i,j-1) )/dy ! this is ny-1 length vector
      memory_dHz_dy(i,j) = b_y_half(j) * memory_dHz_dy(i,j) + a_y_half(j) * value_dHz_dy
      value_dHz_dy = value_dHz_dy/K_y_half(j) + memory_dHz_dy(i,j)

      Ex(i,j) = (( caEx(i,j) + caEx(i,j-1) )/2) * Ex(i,j) + &
        (( cbEx(i,j) + cbEx(i,j-1) )/2 ) * value_dHz_dy
    enddo
  enddo

  do j = 1,ny-1
    do i = 2,nx-1
      ! Update the Ey field
      value_dHz_dx = ( Hz(i,j) - Hz(i-1,j) )/dx
      memory_dHz_dx(i,j) = b_x_half(i) * memory_dHz_dx(i,j) + a_x_half(i) * value_dHz_dx
      value_dHz_dx = value_dHz_dx/K_x_half(i) + memory_dHz_dx(i,j)
      
      Ey(i,j) = (( caEy(i,j) + caEy(i-1,j) )/2) * Ey(i,j) + &
        (( cbEy(i,j) + cbEy(i-1,j) )/2) * value_dHz_dx 
    enddo
  enddo


  !----------------------------------------------------------------------------

  Ex(isource,jsource) = Ex(isource,jsource) + srcx(it) * DT / epsilonx(isource,jsource)
  Ey(isource,jsource) = Ey(isource,jsource) + srcy(it) * DT / epsilony(isource,jsource) 
  
  ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
  Ex(1,:) = 0.d0
  Ex(nx-1,:) = 0.d0
  Ex(:,1) = 0.d0
  Ex(:,ny) = 0.d0

  Ey(1,:) = 0.d0
  Ey(nx,:) = 0.d0
  Ey(:,1) = 0.d0
  Ey(:,ny-1) = 0.d0

  Hz(1,:) = 0.d0
  Hz(nx-1,:) = 0.d0
  Hz(:,1) = 0.d0
  Hz(:,ny-1) = 0.d0

  ! print maximum of norm of velocity
  velocnorm = maxval(sqrt(Ex**2 + Ey**2))
  if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'

  call write_image(Ex, nx-1, ny, it, 'Ex')
  call write_image(Ey, nx, ny-1, it, 'Ez')

enddo   ! end of time loop


end subroutine electromag_cpml_2d


!==============================================================================
!==============================================================================
!==============================================================================
subroutine write_image(image_data, nx, ny, it, channel)

implicit none

integer, parameter :: dp = kind(0.d0)
integer :: nx, ny, it
real(kind=dp) :: image_data(nx, ny)
character(len=2) :: channel
character(len=100) :: filename

WRITE (filename, "(a2, i6.6, '.dat')" ) channel, it

open(unit = 10, form = 'unformatted', file = trim(filename) )
write(10) sngl(image_data)

close(unit = 10)

end subroutine write_image


end module electromagFDTD2d