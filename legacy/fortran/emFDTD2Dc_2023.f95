
subroutine electromag_cpml_2d(nx, nz, dx, dz, npoints_pml, src, nstep)

  implicit none

  integer,parameter :: dp=kind(0.d0)

  ! total number of grid points in each direction of the grid
  integer :: nx, nz
  ! integer, dimension(:,:) :: rcx 

  ! thickness of the PML layer in grid points
  integer :: npoints_pml, nstep

  ! integer, dimension(nx,nz)
  real(kind=dp), dimension(nx,nz) :: eps11, eps13, eps33, epsilonx, epsilonz
  real(kind=dp), dimension(nx,nz) :: sig11, sig13, sig33, sigmax, sigmaz

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
  complex(kind=dp), dimension(nx-1,nz) :: Ex
  complex(kind=dp), dimension(nx,nz-1) :: Ez
  complex(kind=dp), dimension(nx-1,nz-1) :: Hy

  ! we will compute the coefficients for the finite difference scheme 
  complex`(kind=dp), dimension(nx, nz) :: caEx, cbEx
  complex`(kind=dp), dimension(nx, nz) :: caEz, cbEz
  real(kind=dp) :: daHy, dbHy

  complex(kind=dp) :: value_dEx_dz, value_dEz_dx, value_dHy_dz, value_dHy_dx

  ! arrays for the memory variables
  ! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
  complex(kind=dp), dimension(nx,nz) ::  memory_dEz_dx,  memory_dEx_dz 
  complex(kind=dp), dimension(nx,nz) :: memory_dHy_dx
  complex(kind=dp), dimension(nx,nz) ::  memory_dHy_dz

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
  DT = minval( (/dx, dz/) )/ ( 2.d0 * Clight/sqrt( minval( (/ REAL(eps11), REAL(eps33) /) ) ) ) 

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
  Ex(:,:) = complex(0.d0, 0.0d0)
  Ez(:,:) = complex(0.d0, 0.0d0)
  Hy(:,:) = complex(0.d0, 0.0d0)

  ! PML
  memory_dEx_dz(:,:) = complex(0.d0, 0.0d0)
  memory_dEz_dx(:,:) = complex(0.d0, 0.0d0)

  memory_dHy_dx(:,:) = complex(0.d0, 0.0d0)
  memory_dHy_dz(:,:) = complex(0.d0, 0.0d0)

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

  !   call write_image2(Ex, nx-1, nz, it, 'Ex')
  !   call write_image2(Ez, nx, nz-1, it, 'Ez')
    call write(Ex, nx-1, nz, it, 'Ex')
    call write(Ez, nx, nz-1, it, 'Ez')

  enddo   ! end of time loop


end subroutine electromag_cpml_2d

