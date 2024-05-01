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


! module electromagFDTD2d

!   implicit none
    
!   contains

  
  !==============================================================================
  !   subroutine permittivity_write(im, mlist, npoints_pml, nx, nz) 
  !     ! STIFFNESS_ARRAYS takes a matrix containing the material integer identifiers 
  !     ! and creates the same size array for each independent coefficient of the 
  !     ! stiffness matrix along with a density matrix. Since we ae using PML
  !     ! boundaries, we will extend the the boundary values through the PML region.
  !     ! 
  !     ! INPUT 
  !     !   im (INTEGER)  
  !     !   mlist (REAL)
  !     !   eps11(i,j), sig11(i,j), eps22(i,j), sig22, (REAL) -
  !     !   npoints_pml (INTEGER) - the 
  !     !   
  !     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
  !     implicit none 
      
  !     integer :: nx,nz
  !     integer,parameter :: dp = kind(0.d0)
  !     integer,dimension(nx,nz) :: im
  !     integer :: i, j, npoints_pml
  !     real(kind=dp), dimension(:,:) :: mlist
  !     real(kind=dp), dimension(2*npoints_pml+nx,2*npoints_pml+nz) :: &
  !               eps11, eps22, eps33, &
  !               eps12, eps13, eps23, &
  !               sig11, sig22, sig33, &
  !               sig12, sig13, sig23
      
  !     !f2py3 intent(in):: im, mlist, npoints_pml, nx, nz
      
  !     ! Allocate space for permittivity and conductivity values
  !     eps11(:,:) = 0.d0
  !     eps12(:,:) = 0.d0
  !     eps13(:,:) = 0.d0
  !     eps22(:,:) = 0.d0
  !     eps23(:,:) = 0.d0
  !     eps33(:,:) = 0.d0
  !     sig11(:,:) = 0.d0
  !     sig12(:,:) = 0.d0
  !     sig13(:,:) = 0.d0
  !     sig22(:,:) = 0.d0
  !     sig23(:,:) = 0.d0
  !     sig33(:,:) = 0.d0
      
  !     do i=npoints_pml+1,nx + npoints_pml
  !       do j=npoints_pml+1,nz + npoints_pml
  !         eps11(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 2)
  !         eps12(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),3)
  !         eps13(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),4)
  !         eps22(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 5)
  !         eps23(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),6)
  !         eps33(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 7)
          
  !         sig11(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 8) 
  !         sig12(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),9)
  !         sig13(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),10)
  !         sig22(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 11)
  !         sig23(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),12)
  !         sig33(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 13)
  !       end do
  !     end do
      
  !     ! Extend the boundary values of the stiffnesses into the PML region
  !     do i = 1,npoints_pml+1
  !       ! top and bottom
  !       eps11( i, : ) = eps11(npoints_pml+1,:)
  !       eps22( i, : ) = eps22(npoints_pml+1,:)
  !       eps33( i, : ) = eps33(npoints_pml+1,:)
  !       eps12( i, : ) = eps12(npoints_pml+1,:)
  !       eps13( i, : ) = eps13(npoints_pml+1,:)
  !       eps23( i, : ) = eps23(npoints_pml+1,:)

  !       eps11( nx+npoints_pml-1+i, : ) = eps11(nx+npoints_pml-1,:)
  !       eps22( nx+npoints_pml-1+i, : ) = eps22(nx+npoints_pml-1,:)
  !       eps33( nx+npoints_pml-1+i, : ) = eps33(nx+npoints_pml-1,:)
  !       eps12( nx+npoints_pml-1+i, : ) = eps12(nx+npoints_pml-1,:)
  !       eps13( nx+npoints_pml-1+i, : ) = eps13(nx+npoints_pml-1,:)
  !       eps23( nx+npoints_pml-1+i, : ) = eps23(nx+npoints_pml-1,:)
        
  !       sig11( i, : ) = sig11(npoints_pml+1,:)
  !       sig22( i, : ) = sig22(npoints_pml+1,:)
  !       sig33( i, : ) = sig33(npoints_pml+1,:)
  !       sig12( i, : ) = sig12(npoints_pml+1,:)
  !       sig13( i, : ) = sig13(npoints_pml+1,:)
  !       sig23( i, : ) = sig23(npoints_pml+1,:)

  !       sig11( nx+npoints_pml-1+i, : ) = sig11(nx+npoints_pml-1,:)
  !       sig22( nx+npoints_pml-1+i, : ) = sig22(nx+npoints_pml-1,:)
  !       sig33( nx+npoints_pml-1+i, : ) = sig33(nx+npoints_pml-1,:)
  !       sig13( nx+npoints_pml-1+i, : ) = sig12(nx+npoints_pml-1,:)
  !       sig13( nx+npoints_pml-1+i, : ) = sig13(nx+npoints_pml-1,:)
  !       sig23( nx+npoints_pml-1+i, : ) = sig23(nx+npoints_pml-1,:)
      
  !     !!!!!  ! left and right
  !       eps11( :, i ) = eps11(:, npoints_pml+1)
  !       eps22( :, i ) = eps22(:, npoints_pml+1)
  !       eps33( :, i ) = eps33(:, npoints_pml+1)
  !       eps12( :, i ) = eps12(:, npoints_pml+1)
  !       eps13( :, i ) = eps13(:, npoints_pml+1)
  !       eps23( :, i ) = eps23(:, npoints_pml+1)

  !       eps11( :, nz+npoints_pml-1+i ) = eps11(:,nz+npoints_pml-1)    
  !       eps22( :, nz+npoints_pml-1+i ) = eps22(:,nz+npoints_pml-1)
  !       eps33( :, nz+npoints_pml-1+i ) = eps33(:,nz+npoints_pml-1)
  !       eps12( :, nz+npoints_pml-1+i ) = eps12(:,nz+npoints_pml-1)    
  !       eps13( :, nz+npoints_pml-1+i ) = eps13(:,nz+npoints_pml-1)
  !       eps23( :, nz+npoints_pml-1+i ) = eps23(:,nz+npoints_pml-1)
        
  !       sig11( :, i ) = sig11(:, npoints_pml+1)
  !       sig22( :, i ) = sig22(:, npoints_pml+1)
  !       sig33( :, i ) = sig33(:, npoints_pml+1)
  !       sig12( :, i ) = sig11(:, npoints_pml+1)
  !       sig13( :, i ) = sig13(:, npoints_pml+1)
  !       sig23( :, i ) = sig33(:, npoints_pml+1)
      
  !       sig11( :, nz+npoints_pml-1+i ) = sig11(:,nz+npoints_pml-1)    
  !       sig22( :, nz+npoints_pml-1+i ) = sig22(:,nz+npoints_pml-1)
  !       sig33( :, nz+npoints_pml-1+i ) = sig33(:,nz+npoints_pml-1)
  !       sig12( :, nz+npoints_pml-1+i ) = sig12(:,nz+npoints_pml-1)    
  !       sig13( :, nz+npoints_pml-1+i ) = sig13(:,nz+npoints_pml-1)
  !       sig23( :, nz+npoints_pml-1+i ) = sig23(:,nz+npoints_pml-1)
  !     end do 

  !     ! Write each of the matrices to file
  !     call material_rw('eps11.dat', eps11, .FALSE.)
  !     call material_rw('eps12.dat', eps12, .FALSE.)
  !     call material_rw('eps13.dat', eps13, .FALSE.)
  !     call material_rw('eps22.dat', eps22, .FALSE.)
  !     call material_rw('eps23.dat', eps23, .FALSE.)
  !     call material_rw('eps33.dat', eps33, .FALSE.)
  !     call material_rw('sig11.dat', sig11, .FALSE.)
  !     call material_rw('sig12.dat', sig12, .FALSE.)
  !     call material_rw('sig13.dat', sig13, .FALSE.)
  !     call material_rw('sig22.dat', sig22, .FALSE.)
  !     call material_rw('sig23.dat', sig23, .FALSE.)
  !     call material_rw('sig33.dat', sig33, .FALSE.)

  !   end subroutine permittivity_write


  ! ! ---------------------------------------------------------------------------
  !   subroutine material_rw(filename, image_data, readfile)

  !     implicit none
      
  !     integer,parameter :: dp = kind(0.d0)
  !     character(len=9) :: filename
  !     real(kind=dp),dimension(:,:) :: image_data
  !     logical :: readfile
      
      
  !     open(unit = 13, form="unformatted", file = trim(filename))
      
  !     if ( readfile ) then
  !       read(13) image_data
  !     else
  !       write(13) image_data
  !     endif
      
  !     close(unit = 13)
    
  !   end subroutine material_rw
    


  ! !==============================================================================
  !   subroutine loadsource(filename, N, srcfn)
      
  !     implicit none

  !     integer,parameter :: dp = kind(0.d0)
  !     character(len=26) :: filename
  !     integer :: N
  !     real(kind=dp),dimension(N) :: srcfn
      
  !     open(unit = 13, form="unformatted", file = trim(filename))
  !     read(13) srcfn
      
  !     close(unit = 13)

  !   end subroutine loadsource

  ! ! -----------------------------------------------------------------------------

  !   subroutine cpml_coeffs(nx, dx, dt, npml, k_max, alpha_max, &
  !               kappa, alpha, acoeff, bcoeff, HALF)

  !     implicit none

  !     integer,parameter :: dp=kind(0.d0)
  !     integer :: i

  !     ! Define real inputs 
  !     real(kind=dp) :: dx, dt, sig_max, k_max, alpha_max 
  !     integer :: nx, npml
  !     logical :: HALF

  !     ! define the output arrays
  !     real(kind=dp),dimension(nx) :: kappa, alpha, acoeff, bcoeff

  !     ! Define all other variables needed in the program
  !     real(kind=dp) :: xoriginleft, xoriginright
  !     real(kind=dp),dimension(nx) :: xval, sigma
  !     real(kind=dp),parameter :: PI = 3.141592653589793238462643d0
  !     real(kind=dp),parameter :: eps0 = 8.85418782d-12, mu0 = 4.0d0*pi*1.0d-7
  !     integer,parameter :: NP = 2, NPA = 2

  !     real(kind=dp) :: abscissa_in_pml, abscissa_normalized

  !     ! ===========================================================
  !     sig_max = ( 0.8d0 * ( dble(NP+1) ) / ( dx * ( mu0 / eps0 )**0.5d0 ) )

  !     sigma(:) = 0.d0

  !     do i=1,nx 
  !       xval(i) = dx * dble(i - 1)
  !     enddo

  !     if (HALF) then 
  !         xval = xval + dx/2.0
  !     endif

  !     xoriginleft = dx * dble( npml )
  !     xoriginright = dx * dble( (NX-1) - npml )

  !     do i=1,nx
  !         !---------- left edge
  !         abscissa_in_PML = xoriginleft - xval(i)
  !         if (abscissa_in_PML >= 0.d0) then
  !             abscissa_normalized = abscissa_in_PML / dble(dx * npml)
  !             sigma(i) = sig_max * abscissa_normalized**NP
  !             ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
  !             kappa(i) = 1.d0 + (K_MAX - 1.d0) * abscissa_normalized**NP
  !             alpha(i) = ALPHA_MAX * (1.d0 - abscissa_normalized)**NPA
  !         endif

  !         !---------- right edge
  !         ! define damping profile at the grid points
  !         abscissa_in_PML = xval(i) - xoriginright
  !         if (abscissa_in_PML >= 0.d0) then
  !           abscissa_normalized = abscissa_in_PML / dble(dx * npml)
  !           sigma(i) = sig_max * abscissa_normalized**NP
  !           kappa(i) = 1.d0 + (k_max - 1.d0) * abscissa_normalized**NP
  !           alpha(i) = alpha_max * (1.d0 - abscissa_normalized)**NPA
  !         endif

  !         ! just in case, for -5 at the end
  !         if (alpha(i) < 0.d0) alpha(i) = 0.d0
  !         ! Compute the b_i coefficents
  !         bcoeff(i) = exp( - (sigma(i) / kappa(i) + alpha(i)) * DT/eps0 )
          
  !         ! Compute the a_i coefficients
  !         ! this to avoid division by zero outside the PML
  !         if (abs(sigma(i)) > 1.d-6) then 
  !           acoeff(i) = sigma(i) * (bcoeff(i) - 1.d0) / ( (sigma(i) + kappa(i) * alpha(i)) ) / kappa(i)
  !         endif

  !     enddo 

  !   end subroutine cpml_coeffs


  !   !==============================================================================

  !   subroutine loadcpml(filename, image_data)

  !     implicit none
      
  !     integer,parameter :: dp = kind(0.d0)
  !     character(len=*) :: filename
  !     real(kind=dp),dimension(:) :: image_data
      
  !     open(unit = 13, form="unformatted", file = trim(filename), access='stream')
  !     read(13) image_data
  !     close(unit = 13)
  !   end subroutine loadcpml
  !   !==============================================================================


  subroutine electromag_cpml_2d(nx, nz, dx, dz, npoints_pml, src, nstep)

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
    !   nx, nz (INTEGER)
    !   eps11, sig11, eps22, sig22 (REAL)
    !   dx, dz (REAL)
    !   npoints_pml (INTEGER) - the thickness of the pml
    !   rcx (INTEGER) - the x and y indices for an array of recievers
    !


    implicit none

    integer,parameter :: dp=kind(0.d0)

    ! total number of grid points in each direction of the grid
    integer :: nx, nz
    ! integer, dimension(:,:) :: rcx 

    ! thickness of the PML layer in grid points
    integer :: npoints_pml, nstep

    ! integer, dimension(nx,nz)
    real(kind=dp), dimension(nx,nz) :: eps11, eps13, eps33, &
                                      sig11, sig13, sig33, &
                                      epsilonx, epsilonz, &
                                      sigmax, sigmaz


    ! time step in seconds. decreasing the time step improves the pml attenuation
    ! but it should be inversely proportional to the center frequency of the 
    ! source frequency 
    real(kind=dp) :: DT, dx, dz 

    ! source
    integer,dimension(:) :: src
    real(kind=dp),dimension(nstep) :: srcx, srcz
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
    real(kind=dp), dimension(nx-1,nz) :: Ex
    real(kind=dp), dimension(nx,nz-1) :: Ez
    real(kind=dp), dimension(nx-1,nz-1) :: Hy

    ! we will compute the coefficients for the finite difference scheme 
    real(kind=dp), dimension(nx, nz) :: caEx, cbEx
    real(kind=dp), dimension(nx, nz) :: caEz, cbEz
    real(kind=dp) :: daHy, dbHy

    real(kind=dp) :: value_dEx_dz, value_dEz_dx, value_dHy_dz, value_dHy_dx

    ! arrays for the memory variables
    ! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
    real(kind=dp), dimension(nx,nz) ::  memory_dEz_dx,  memory_dEx_dz 
    real(kind=dp), dimension(nx,nz) :: memory_dHy_dx
    real(kind=dp), dimension(nx,nz) ::  memory_dHy_dz

    ! parameters for the source
    ! angle of source force clockwise with respect to vertical (Y) axis
    ! this will later be treated as an input
    real(kind=dp), parameter :: factor = 1.d0
    ! character(len=6) :: src_type
    integer :: i,j, it

    real(kind=dp) :: velocnorm

    ! -------------------------------- PML parameters 
    ! power to compute d0 profile. Increasing this value allows for a larger dampening gradient in the PML
    real(kind=dp), parameter :: NP = 2.d0, NPA = 1.d0

    ! 1D arrays for the damping profiles
    real(kind=dp), dimension(nx) :: K_x,alpha_x,a_x,b_x, &
                                    K_x_half, alpha_x_half,a_x_half,b_x_half
    real(kind=dp), dimension(nz) :: K_z,alpha_z,a_z,b_z, &
                                    K_z_half, alpha_z_half,a_z_half,b_z_half

    ! ------------------- Name the f2py inputs 
    !f2py3 intent(in) :: nx, nz, dx, dz,
    !f2py3 intent(in) :: noints_pml, src, nstep

    ! =============================================================================
    ! ----------------------- Load Permittivity Coefficients ----------------------

    call material_rw('eps11.dat', eps11, .TRUE.)
    call material_rw('eps13.dat', eps13, .TRUE.)
    call material_rw('eps33.dat', eps33, .TRUE.) ! We will change y to z soon
    call material_rw('sig11.dat', sig11, .TRUE.)
    call material_rw('sig13.dat', sig13, .TRUE.)
    call material_rw('sig33.dat', sig33, .TRUE.)

    ! ------------------------ Assign some constants -----------------------

    ! Assign the source location indices
    isource = int(src(1)) + npoints_pml
    jsource = int(src(2)) + npoints_pml

    ! Define the 
    DT = minval( (/dx, dz/) )/ ( 2.d0 * Clight/sqrt( minval( (/ eps11, eps33 /) ) ) ) 

    ! Compute the coefficients of the FD scheme. First scale the relative 
    ! permittivity and permeabilities to get the absolute values 
    epsilonx(:,:) = (eps11 + eps13)*eps0
    epsilonz(:,:) = (eps13 + eps33)*eps0
    sigmax(:,:) = sig11 + sig13 
    sigmaz(:,:) = sig13 + sig33

    ! We need to change sigma to dsigma, same for epsilon

    caEx(:,:) = ( 1.0d0 - sigmax * dt / ( 2.0d0 * epsilonx ) ) / &
                ( 1.0d0 + sigmax * dt / (2.0d0 * epsilonx ) )
    cbEx(:,:) = (dt / epsilonx ) / ( 1.0d0 + sigmax * dt / ( 2.0d0 * epsilonx ) )

    caEz(:,:) = ( 1.0d0 - sigmaz * dt / ( 2.0d0 * epsilonz ) ) / &
                ( 1.0d0 + sigmaz * dt / (2.0d0 * epsilonz ) )
    cbEz(:,:) = (dt / epsilonz ) / ( 1.0d0 + sigmaz * dt / ( 2.0d0 * epsilonz ) )

    daHy = dt/(4.0d0*mu0*mu)
    dbHy = dt/mu0 !dt/(mu*mu*dx*(1+daHy) ) 
    daHy = 1.0d0 ! (1-daHy)/(1+daHy) ! 


    ! ----------------------------------------------------------------------

    ! ================================ LOAD SOURCE ================================

    call loadsource('electromagneticsourcex.dat', nstep, srcx)
    call loadsource('electromagneticsourcez.dat', nstep, srcz)



    ! ----------------------------------------------------------------------

    ! Initialize CPML damping variables
    K_x(:) = 1.0d0
    K_x_half(:) = 1.0d0
    alpha_x(:) = 0.0d0
    alpha_x_half(:) = 0.0d0
    a_x(:) = 0.0d0
    a_x_half(:) = 0.0d0
    b_x(:) = 0.0d0 
    b_x_half(:) = 0.0d0 

    K_z(:) = 1.0d0
    K_z_half(:) = 1.0d0
    alpha_z(:) = 0.0d0
    alpha_z_half(:) = 0.0d0
    a_z(:) = 0.0d0
    a_z_half(:) = 0.0d0


    call loadcpml('kappax_cpml.dat', K_x)
    call loadcpml('alphax_cpml.dat', alpha_x)
    call loadcpml('acoefx_cpml.dat', a_x)
    call loadcpml('bcoefx_cpml.dat', b_x)

    call loadcpml('kappaz_cpml.dat', K_z)
    call loadcpml('alphaz_cpml.dat', alpha_z)
    call loadcpml('acoefz_cpml.dat', a_z)
    call loadcpml('bcoefz_cpml.dat', b_z)

    call loadcpml('kappax_half_cpml.dat', K_x_half)
    call loadcpml('alphax_half_cpml.dat', alpha_x_half)
    call loadcpml('acoefx_half_cpml.dat', a_x_half)
    call loadcpml('bcoefx_half_cpml.dat', b_x_half)

    call loadcpml('kappaz_half_cpml.dat', K_z_half)
    call loadcpml('alphaz_half_cpml.dat', alpha_z_half)
    call loadcpml('acoefz_half_cpml.dat', a_z_half)
    call loadcpml('bcoefz_half_cpml.dat', b_z_half)


    ! initialize arrays
    Ex(:,:) = 0.0d0
    Ez(:,:) = 0.0d0
    Hy(:,:) = 0.0d0

    ! PML
    memory_dEx_dz(:,:) = 0.0d0
    memory_dEz_dx(:,:) = 0.0d0

    memory_dHy_dx(:,:) = 0.0d0
    memory_dHy_dz(:,:) = 0.0d0

    !---
    !---  beginning of time loop
    !---

    do it = 1,NSTEP
      
      !--------------------------------------------------------
      ! compute magnetic field and update memory variables for C-PML
      !--------------------------------------------------------
      do i = 1,nx-1  
        do j = 1,nz-1
          
          ! Values needed for the magnetic field updates
          value_dEx_dz = ( Ex(i,j+1) - Ex(i,j) )/dz
          memory_dEx_dz(i,j) = b_z(j) * memory_dEx_dz(i,j) + a_z(j) * value_dEx_dz
          value_dEx_dz = value_dEx_dz/ K_z(j) + memory_dEx_dz(i,j)

          ! The rest of the equation needed for agnetic field updates
          value_dEz_dx = ( Ez(i+1,j) - Ez(i,j) )/dx
          memory_dEz_dx(i,j) = b_x(i) * memory_dEz_dx(i,j) + a_x(i) * value_dEz_dx
          value_dEz_dx = value_dEz_dx/ K_x(i) + memory_dEz_dx(i,j)

          ! Now update the Magnetic field
          Hy(i,j) = daHy*Hy(i,j) + dbHy*( value_dEz_dx + value_dEx_dz )

        enddo  
      enddo

      !--------------------------------------------------------
      ! compute electric field and update memory variables for C-PML
      !--------------------------------------------------------
      
      ! Compute the differences in the y-direction
      do j = 2,nz-1
        do i = 1,nx-1
          ! Update the Ex field
          value_dHy_dz = ( Hy(i,j) - Hy(i,j-1) )/dz ! this is nz-1 length vector
          memory_dHy_dz(i,j) = b_z_half(j) * memory_dHy_dz(i,j) + a_z_half(j) * value_dHy_dz
          value_dHy_dz = value_dHy_dz/K_z_half(j) + memory_dHy_dz(i,j)

          Ex(i,j) = (( caEx(i,j) + caEx(i,j-1) )/2) * Ex(i,j) + &
            (( cbEx(i,j) + cbEx(i,j-1) )/2 ) * value_dHy_dz
        enddo
      enddo

      do j = 1,nz-1
        do i = 2,nx-1
          ! Update the Ez field
          value_dHy_dx = ( Hy(i,j) - Hy(i-1,j) )/dx
          memory_dHy_dx(i,j) = b_x_half(i) * memory_dHy_dx(i,j) + a_x_half(i) * value_dHy_dx
          value_dHy_dx = value_dHy_dx/K_x_half(i) + memory_dHy_dx(i,j)
          
          Ez(i,j) = (( caEz(i,j) + caEz(i-1,j) )/2) * Ez(i,j) + &
            (( cbEz(i,j) + cbEz(i-1,j) )/2) * value_dHy_dx 
        enddo
      enddo


      !----------------------------------------------------------------------------

      Ex(isource,jsource) = Ex(isource,jsource) + srcx(it) * DT / eps11(isource,jsource)
      Ez(isource,jsource) = Ez(isource,jsource) + srcz(it) * DT / eps33(isource,jsource) 
      
      ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
      Ex(1,:) = 0.d0
      Ex(nx-1,:) = 0.d0
      Ex(:,1) = 0.d0
      Ex(:,nz) = 0.d0

      Ez(1,:) = 0.d0
      Ez(nx,:) = 0.d0
      Ez(:,1) = 0.d0
      Ez(:,nz-1) = 0.d0

      Hy(1,:) = 0.d0
      Hy(nx-1,:) = 0.d0
      Hy(:,1) = 0.d0
      Hy(:,nz-1) = 0.d0

      ! print maximum of norm of velocity
      velocnorm = maxval(sqrt(Ex**2 + Ez**2))
      if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'

      call write_image2(Ex, nx-1, nz, it, 'Ex')
      call write_image2(Ez, nx, nz-1, it, 'Ez')

    enddo   ! end of time loop


  end subroutine electromag_cpml_2d

  !==============================================================================
  ! subroutine write_image2(image_data, nx, nz, src, it, channel)

  !   implicit none

  !   integer, parameter :: dp = kind(0.d0)
  !   integer :: nx, nz, it
  !   real(kind=dp) :: image_data(nx, nz)
  !   character(len=2) :: channel
  !   character(len=100) :: filename

  !   WRITE (filename, "(a2, i6.6, '.dat')" ) channel, it

  !   open(unit = 10, form = 'unformatted', file = trim(filename) )
  !   write(10) sngl(image_data)

  !   close(unit = 10)

  ! end subroutine write_image2


! end module electromagFDTD2d