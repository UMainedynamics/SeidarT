! The seismicFDTD2D module is designed to be integrated with python via f2py. 
!
! Compile using
!     f2py3 -c --fcompiler=gnu95 -m seismicfdtd2d_dp seismicFDTD2D_dp.f95
!
! Created by Steven Bernsen with T-minus one week to AGU
! University of Maine
! Department of Earth and Environmental Sciences 
! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module seismicFDTD25d

implicit none

contains

  !==============================================================================
subroutine stiffness_write(im, mlist, npoints_pml, nx, nz, gradient) 
  ! STIFFNESS_ARRAYS takes a matrix containing the material integer identifiers 
  ! and creates the same size array for each independent coefficient of the 
  ! stiffness matrix along with a density matrix. Since we ae using PML
  ! boundaries, we will extend the the boundary values through the PML region.
  ! 
  ! INPUT 
  !   im (INTEGER)  
  !   mlist (REAL)
  !   c11(i,j), c12(i,j), c22(i,j), c66, rho(i,j) (REAL) -
  !   npoints_pml (INTEGER) - the 
  !   
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  implicit none 

  integer,parameter :: dp = kind(0.d0)
  integer :: nx, nz
  integer,dimension(nx,nz) :: im
  integer :: i, j, npoints_pml
  real(kind=dp),dimension(:,:) :: mlist
  real(kind=dp),dimension(nx,nz) :: gradient
  real(kind=dp),dimension(2*npoints_pml+nx,2*npoints_pml+nz) :: c11,c12,c13,&
                                                              c14,c15,c16, &
                                                              c22,c23,c24,c25,c26,&
                                                              c33,c34,c35,c36, &
                                                              c44,c45,c46, &
                                                              c55,c56,c66, &
                                                              rho
  

  !f2py3 intent(in) :: im, mlist, npoints_pml, nx, nz, gradient

  c11(:,:) = 0.d0 
  c12(:,:) = 0.d0 
  c13(:,:) = 0.d0
  c14(:,:) = 0.d0 
  c15(:,:) = 0.d0 
  c16(:,:) = 0.d0 
  c22(:,:) = 0.d0 
  c23(:,:) = 0.d0 
  c24(:,:) = 0.d0 
  c25(:,:) = 0.d0 
  c26(:,:) = 0.d0 
  c33(:,:) = 0.d0 
  c34(:,:) = 0.d0 
  c35(:,:) = 0.d0 
  c36(:,:) = 0.d0 
  c44(:,:) = 0.d0 
  c45(:,:) = 0.d0 
  c46(:,:) = 0.d0 
  c55(:,:) = 0.d0 
  c56(:,:) = 0.d0 
  c66(:,:) = 0.d0 
  rho(:,:) = 0.d0 

  !Assign between the PML regions
  do i = npoints_pml+1, nx+npoints_pml
    do j = npoints_pml+1, nz+npoints_pml
      c11(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 2)
      c12(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 3)
      c13(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 4)
      c14(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 5)
      c15(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 6)
      c16(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 7)
      c22(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 8)
      c23(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 9)
      c24(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 10)
      c25(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 11)
      c26(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 12)
      c33(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 13)
      c34(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 14)
      c35(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 15)
      c36(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 16)
      c44(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 17)
      c45(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 18)
      c46(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 19)
      c55(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 20)
      c56(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 21)
      c66(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 22)
      rho(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 23) 
    enddo
  enddo

  rho(npoints_pml+1:nx+npoints_pml,npoints_pml+1:nz+npoints_pml) = &
      rho(npoints_pml+1:nx+npoints_pml,npoints_pml+1:nz+npoints_pml)*gradient

  ! Extend the boundary values of the stiffnesses into the PML region
  do i = 1,npoints_pml+1
    ! top 
    c11( i, :) = c11(npoints_pml+1,:)
    c12( i, :) = c12(npoints_pml+1,:)
    c13( i, :) = c13(npoints_pml+1,:)
    c14( i, :) = c14(npoints_pml+1,:)
    c15( i, :) = c15(npoints_pml+1,:)
    c16( i, :) = c16(npoints_pml+1,:)
    c22( i, :) = c22(npoints_pml+1,:)
    c23( i, :) = c23(npoints_pml+1,:)
    c24( i, :) = c24(npoints_pml+1,:)
    c25( i, :) = c25(npoints_pml+1,:)
    c26( i, :) = c26(npoints_pml+1,:)
    c33( i, :) = c33(npoints_pml+1,:)
    c34( i, :) = c34(npoints_pml+1,:)
    c35( i, :) = c35(npoints_pml+1,:)
    c36( i, :) = c36(npoints_pml+1,:)
    c44( i, :) = c44(npoints_pml+1,:)
    c45( i, :) = c45(npoints_pml+1,:)
    c46( i, :) = c46(npoints_pml+1,:)
    c55( i, :) = c55(npoints_pml+1,:)
    c56( i, :) = c56(npoints_pml+1,:)
    c66( i, :) = c66(npoints_pml+1,:)
    rho( i, :) = rho(npoints_pml+1,:)

    ! bottom
    c11( nx+npoints_pml-1+i, :) = c11(nx+npoints_pml-1,:)
    c12( nx+npoints_pml-1+i, :) = c12(nx+npoints_pml-1,:)
    c13( nx+npoints_pml-1+i, :) = c13(nx+npoints_pml-1,:)
    c14( nx+npoints_pml-1+i, :) = c14(nx+npoints_pml-1,:)
    c15( nx+npoints_pml-1+i, :) = c15(nx+npoints_pml-1,:)
    c16( nx+npoints_pml-1+i, :) = c16(nx+npoints_pml-1,:)
    c22( nx+npoints_pml-1+i, :) = c22(nx+npoints_pml-1,:)
    c23( nx+npoints_pml-1+i, :) = c23(nx+npoints_pml-1,:)
    c24( nx+npoints_pml-1+i, :) = c24(nx+npoints_pml-1,:)
    c25( nx+npoints_pml-1+i, :) = c25(nx+npoints_pml-1,:)
    c26( nx+npoints_pml-1+i, :) = c26(nx+npoints_pml-1,:)
    c33( nx+npoints_pml-1+i, :) = c33(nx+npoints_pml-1,:)
    c34( nx+npoints_pml-1+i, :) = c34(nx+npoints_pml-1,:)
    c35( nx+npoints_pml-1+i, :) = c35(nx+npoints_pml-1,:)
    c36( nx+npoints_pml-1+i, :) = c36(nx+npoints_pml-1,:)
    c44( nx+npoints_pml-1+i, :) = c44(nx+npoints_pml-1,:)
    c45( nx+npoints_pml-1+i, :) = c45(nx+npoints_pml-1,:)
    c46( nx+npoints_pml-1+i, :) = c46(nx+npoints_pml-1,:)
    c55( nx+npoints_pml-1+i, :) = c55(nx+npoints_pml-1,:)
    c56( nx+npoints_pml-1+i, :) = c56(nx+npoints_pml-1,:)
    c66( nx+npoints_pml-1+i, :) = c66(nx+npoints_pml-1,:)
    rho( nx+npoints_pml-1+i, :) = rho(nx+npoints_pml-1,:)

    ! left 
    c11( :, i) = c11(:, npoints_pml+1)
    c12( :, i) = c12(:, npoints_pml+1)
    c13( :, i) = c13(:, npoints_pml+1)
    c14( :, i) = c14(:, npoints_pml+1)
    c15( :, i) = c15(:, npoints_pml+1)
    c16( :, i) = c16(:, npoints_pml+1)
    c22( :, i) = c22(:, npoints_pml+1)
    c23( :, i) = c23(:, npoints_pml+1)
    c24( :, i) = c24(:, npoints_pml+1)
    c25( :, i) = c25(:, npoints_pml+1)
    c26( :, i) = c26(:, npoints_pml+1)
    c33( :, i) = c33(:, npoints_pml+1)
    c34( :, i) = c34(:, npoints_pml+1)
    c35( :, i) = c35(:, npoints_pml+1)
    c36( :, i) = c36(:, npoints_pml+1)
    c44( :, i) = c44(:, npoints_pml+1)
    c45( :, i) = c45(:, npoints_pml+1)
    c46( :, i) = c46(:, npoints_pml+1)
    c55( :, i) = c55(:, npoints_pml+1)
    c56( :, i) = c56(:, npoints_pml+1)
    c66( :, i) = c66(:, npoints_pml+1)
    rho( :, i) = rho(:, npoints_pml+1)

    ! right
    c11( :, nz+npoints_pml-1+i) = c11(:,nz+npoints_pml-1)
    c12( :, nz+npoints_pml-1+i) = c12(:,nz+npoints_pml-1)
    c13( :, nz+npoints_pml-1+i) = c13(:,nz+npoints_pml-1)      
    c14( :, nz+npoints_pml-1+i) = c14(:,nz+npoints_pml-1)      
    c15( :, nz+npoints_pml-1+i) = c15(:,nz+npoints_pml-1)      
    c16( :, nz+npoints_pml-1+i) = c16(:,nz+npoints_pml-1)      
    c22( :, nz+npoints_pml-1+i) = c22(:,nz+npoints_pml-1)
    c23( :, nz+npoints_pml-1+i) = c23(:,nz+npoints_pml-1)
    c24( :, nz+npoints_pml-1+i) = c24(:,nz+npoints_pml-1)
    c25( :, nz+npoints_pml-1+i) = c25(:,nz+npoints_pml-1)
    c26( :, nz+npoints_pml-1+i) = c26(:,nz+npoints_pml-1)
    c33( :, nz+npoints_pml-1+i) = c33(:,nz+npoints_pml-1)
    c34( :, nz+npoints_pml-1+i) = c34(:,nz+npoints_pml-1)
    c35( :, nz+npoints_pml-1+i) = c35(:,nz+npoints_pml-1)
    c36( :, nz+npoints_pml-1+i) = c36(:,nz+npoints_pml-1)
    c44( :, nz+npoints_pml-1+i) = c44(:,nz+npoints_pml-1)
    c45( :, nz+npoints_pml-1+i) = c45(:,nz+npoints_pml-1)
    c46( :, nz+npoints_pml-1+i) = c46(:,nz+npoints_pml-1)
    c55( :, nz+npoints_pml-1+i) = c55(:,nz+npoints_pml-1)      
    c56( :, nz+npoints_pml-1+i) = c56(:,nz+npoints_pml-1)      
    c66( :, nz+npoints_pml-1+i) = c66(:,nz+npoints_pml-1)
    rho( :, nz+npoints_pml-1+i) = rho(:,nz+npoints_pml-1)

  end do 

  ! Write each of the matrices to file
  call material_rw('c11.dat', c11, .FALSE.)
  call material_rw('c12.dat', c12, .FALSE.)
  call material_rw('c13.dat', c13, .FALSE.)
  call material_rw('c14.dat', c14, .FALSE.)
  call material_rw('c15.dat', c15, .FALSE.)
  call material_rw('c16.dat', c16, .FALSE.)
  call material_rw('c22.dat', c22, .FALSE.)
  call material_rw('c23.dat', c23, .FALSE.)
  call material_rw('c24.dat', c24, .FALSE.)
  call material_rw('c25.dat', c25, .FALSE.)
  call material_rw('c26.dat', c26, .FALSE.)
  call material_rw('c33.dat', c33, .FALSE.)
  call material_rw('c34.dat', c34, .FALSE.)
  call material_rw('c35.dat', c35, .FALSE.)
  call material_rw('c36.dat', c36, .FALSE.)
  call material_rw('c44.dat', c44, .FALSE.)
  call material_rw('c45.dat', c45, .FALSE.)
  call material_rw('c46.dat', c46, .FALSE.)
  call material_rw('c55.dat', c55, .FALSE.)
  call material_rw('c56.dat', c56, .FALSE.)
  call material_rw('c66.dat', c66, .FALSE.)
  call material_rw('rho.dat', rho, .FALSE. )

end subroutine stiffness_write

  ! ---------------------------------------------------------------------------
  subroutine material_rw(filename, image_data, readfile)

  implicit none

  integer,parameter :: dp = kind(0.d0)
  character(len=7) :: filename
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
    character(len=18) :: filename
    integer :: N
    real(kind=dp),dimension(N) :: srcfn
    
    open(unit = 13, form="unformatted", file = trim(filename))
    read(13) srcfn
    
    close(unit = 13)
  
  end subroutine loadsource

! -----------------------------------------------------------------------------

  subroutine cpml_coeffs(nx, dx, dt, npml, sig_max, k_max, alpha_max, &
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
  integer,parameter :: NP = 2, NPA = 2

  real(kind=dp) :: abscissa_in_pml, abscissa_normalized

  ! ---------------- End of memory definitions ----------------
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
          ! from Stephen Gedney's unpublished class notes for class EE699, 
          ! lecture 8, slide 8-2
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
      bcoeff(i) = exp( - (sigma(i) / kappa(i) + alpha(i)) * DT )
      
      ! Compute the a_i coefficients
      ! this to avoid division by zero outside the PML
      if (abs(sigma(i)) > 1.d-6) then 
        acoeff(i) = sigma(i) * (bcoeff(i) - 1.d0) / &
        (kappa(i) * sigma(i) + kappa(i) * alpha(i) ) !( (sigma(i) + kappa(i) * alpha(i)) ) / kappa(i)
      endif

  enddo 

  end subroutine cpml_coeffs
! -----------------------------------------------------------------------------

  subroutine loadcpml(filename, image_data)

    implicit none
    
    integer,parameter :: dp = kind(0.d0)
    character(len=*) :: filename
    real(kind=dp),dimension(:) :: image_data
    
    open(unit = 13, form="unformatted", file = trim(filename), access='stream')
    read(13) image_data
    close(unit = 13)
  end subroutine loadcpml


!==============================================================================
subroutine seismic_cpml_25d(nx, ny, nz, dx, dy, dz, &
                      npoints_pml, src, f0, nstep)

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
!   c11, c12, c22, c66, rho (REAL)
!   dx, dy (REAL)
!   npoints_pml (INTEGER) - the thickness of the pml


implicit none

integer, parameter :: dp=kind(0.d0)

! total number of grid points in each direction of the grid
integer :: nx, ny, nz

! thickness of the PML layer in grid points
integer :: npoints_pml
real(kind=dp), dimension(nx,nz) ::  c11, c12, c13, c14, c15, c16, &
                                    c22, c23, c24, c25, c26, &
                                    c33, c34, c35, c36, &
                                    c44, c45, c46, &
                                    c55, c56, &
                                    c66, &
                                    rho
real(kind=dp) :: f0, deltarho

! total number of time steps
integer :: nstep

! time step in seconds. decreasing the time step improves the pml attenuation
! but it should be inversely proportional to the center frequency of the 
! source frequency 
real(kind=dp) :: DT
real(kind=dp) :: dx, dy, dz
! parameters for the source
real(kind=dp), parameter :: factor = 1.d07

! source
integer,dimension(:) :: src
integer :: isource, jsource, ksource

! value of PI
real(kind=dp), parameter :: PI = 3.141592653589793238462643d0

! conversion from degrees to radians
real(kind=dp), parameter :: DEGREES_TO_RADIANS = PI / 180.d0

! large value for maximum
real(kind=dp), parameter :: HUGEVAL = 1.d+30

! velocity threshold above which we consider that the code became unstable
real(kind=dp), parameter :: STABILITY_THRESHOLD = 1.d+25

! main arrays
real(kind=dp), dimension(nx,ny,nz) :: vx,vy,vz, &
                                      sigmaxx,sigmayy,sigmazz, &
                                      sigmaxy,sigmaxz,sigmayz

! power to compute d0 profile. Increasing this value allows for a larger dampening gradient in the PML
real(kind=dp), parameter :: NPOWER = 2.d0

! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
real(kind=dp), parameter :: K_MAX = 1.d0
real(kind=dp) :: ALPHA_MAX

! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
real(kind=dp), dimension(nx,ny,nz) :: &
    memory_dvx_dx, memory_dvx_dy, memory_dvx_dz, &
    memory_dvy_dx, memory_dvy_dy, memory_dvy_dz, &
    memory_dvz_dx, memory_dvz_dy, memory_dvz_dz, &
    memory_dsigmaxx_dx, memory_dsigmayy_dy, memory_dsigmazz_dz, &
    memory_dsigmaxy_dx, memory_dsigmaxy_dy, &
    memory_dsigmaxz_dx, memory_dsigmaxz_dz, &
    memory_dsigmayz_dy, memory_dsigmayz_dz

real(kind=dp) :: &
    dvx_dx, dvx_dy, dvx_dz, &
    dvy_dx, dvy_dy, dvy_dz, &
    dvz_dx, dvz_dy, dvz_dz, &
    dsigmaxx_dx, dsigmayy_dy, dsigmazz_dz, &
    dsigmaxy_dx, dsigmaxy_dy, &
    dsigmaxz_dx, dsigmaxz_dz, &
    dsigmayz_dy, dsigmayz_dz

! 1D arrays for the damping profiles
real(kind=dp), dimension(nx) :: d_x,K_x,alpha_x,a_x,b_x,d_x_half, &
                                K_x_half,alpha_x_half,a_x_half,b_x_half
real(kind=dp), dimension(ny) :: d_y,K_y,alpha_y,a_y,b_y,d_y_half, &
                                K_y_half,alpha_y_half,a_y_half,b_y_half
real(kind=dp), dimension(nz) :: d_z,K_z,alpha_z,a_z,b_z,d_z_half, &
                                K_z_half,alpha_z_half,a_z_half,b_z_half

real(kind=dp) :: thickness_PML_x,thickness_PML_y,thickness_PML_z
real(kind=dp) :: Rcoef,d0_x,d0_y,d0_z

! for the source
real(kind=dp) :: a!, t
real(kind=dp),dimension(nstep) :: srcx, srcy, srcz

integer :: i,j,k,it

real(kind=dp) :: velocnorm

! for stability estimate
real(kind=dp) :: quasi_cp_max

! Name the f2py inputs 
!f2py3 intent(in) :: nx, ny, nz, dx, dy, dz,
!f2py3 intent(in) :: noints_pml, src, f0, nstep

! ------------------------ Load Stiffness Coefficients ------------------------

call material_rw('c11.dat', c11, .TRUE.)
call material_rw('c12.dat', c12, .TRUE.)
call material_rw('c13.dat', c13, .TRUE.)
call material_rw('c14.dat', c14, .TRUE.)
call material_rw('c15.dat', c15, .TRUE.)
call material_rw('c16.dat', c16, .TRUE.)
call material_rw('c22.dat', c22, .TRUE.)
call material_rw('c23.dat', c23, .TRUE.)
call material_rw('c24.dat', c24, .TRUE.)
call material_rw('c25.dat', c25, .TRUE.)
call material_rw('c26.dat', c26, .TRUE.)
call material_rw('c33.dat', c33, .TRUE.)
call material_rw('c34.dat', c34, .TRUE.)
call material_rw('c35.dat', c35, .TRUE.)
call material_rw('c36.dat', c36, .TRUE.)
call material_rw('c44.dat', c44, .TRUE.)
call material_rw('c45.dat', c45, .TRUE.)
call material_rw('c46.dat', c46, .TRUE.)
call material_rw('c55.dat', c55, .TRUE.)
call material_rw('c56.dat', c56, .TRUE.)
call material_rw('c66.dat', c66, .TRUE.)
call material_rw('rho.dat', rho, .TRUE.)

! ------------------------ Assign some constants -----------------------
isource = src(1)+npoints_pml
jsource = src(2)+npoints_pml
ksource = src(3)+npoints_pml

! To ensure a courant number <= 1.0, we can calculate the time step from
! the velocity
DT = 0.7 * minval( (/dx,dy,dz/) )/ &
    ( sqrt( 3.d0 * ( maxval( (/ c11/rho, c22/rho, c33/rho /) ) ) ) )

ALPHA_MAX = pi*f0  ! from Festa and Vilotte
a = pi*pi*f0*f0

! ----------------------------------------------------------------------
!---
!--- program starts here
!---

! ================================ LOAD SOURCE ================================

call loadsource('seismicsourcex.dat', nstep, srcx)
call loadsource('seismicsourcey.dat', nstep, srcy)
call loadsource('seismicsourcez.dat', nstep, srcz)

! =============================================================================

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
Rcoef = 0.001d0

! check that NPOWER is okay
  if (NPOWER < 1) stop 'NPOWER must be greater than 1'
  
! ==================================== PML ====================================
! ---------------- define profile of absorption in PML region -----------------

! thickness of the PML layer in meters
! thickness_PML_x = npoints_pml * dx
! thickness_PML_y = npoints_pml * dy
! thickness_PML_z = npoints_pml * dz

! Initialize PML 
  d_x(:) = 0.d0
  d_x_half(:) = 0.d0
  K_x(:) = 1.d0
  K_x_half(:) = 1.d0
  alpha_x(:) = 0.d0
  alpha_x_half(:) = 0.d0
  a_x(:) = 0.d0
  a_x_half(:) = 0.d0

  d_y(:) = 0.d0
  d_y_half(:) = 0.d0
  K_y(:) = 1.d0
  K_y_half(:) = 1.d0
  alpha_y(:) = 0.d0
  alpha_y_half(:) = 0.d0
  a_y(:) = 0.d0
  a_y_half(:) = 0.d0

  d_z(:) = 0.d0
  d_z_half(:) = 0.d0
  K_z(:) = 1.d0
  K_z_half(:) = 1.d0 
  alpha_z(:) = 0.d0
  alpha_z_half(:) = 0.d0
  a_z(:) = 0.d0
  a_z_half(:) = 0.d0


  ! to compute d0 below, and for stability estimate
  ! quasi_cp_max = sqrt( maxval( (/ c11/rho, c22/rho, c33/rho /) ) )
  quasi_cp_max = ( minval( (/dx,dy,dz/) )/ ( 2.0 * dt) )

  ! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_x)
  d0_y = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_y)
  d0_z = - dble(NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_z)


! ------------------------- damping in the X direction-------------------------

! ! Full nodes
! call cpml_coeffs(nx, dx, dt, npoints_pml, d0_x, k_max, alpha_max, &
!             K_x, alpha_x, a_x, b_x, .FALSE.)

! ! Half nodes
! call cpml_coeffs(nx, dx, dt, npoints_pml, d0_x, k_max, alpha_max, &
!             K_x_half, alpha_x_half, a_x_half, b_x_half, .TRUE.)

! ! ------------------------ damping in the Y direction -------------------------
! call cpml_coeffs(ny, dy, dt, npoints_pml, d0_y, k_max, alpha_max, &
!             K_y, alpha_y, a_y, b_y, .FALSE.)

! call cpml_coeffs(ny, dy, dt, npoints_pml, d0_y, k_max, alpha_max, &
!             K_y_half, alpha_y_half, a_y_half, b_y_half, .TRUE.)

! ! ------------------------- damping in the Z direction-------------------------

! call cpml_coeffs(nz, dz, dt, npoints_pml, d0_z, k_max, alpha_max, &
!             K_z, alpha_z, a_z, b_z, .FALSE.)

! call cpml_coeffs(nz, dz, dt, npoints_pml, d0_z, k_max, alpha_max, &
!             K_z_half, alpha_z_half, a_z_half, b_z_half, .TRUE.)

  call loadcpml('kappax_cpml.dat', K_x)
  call loadcpml('alphax_cpml.dat', alpha_x)
  call loadcpml('acoefx_cpml.dat', a_x)
  call loadcpml('bcoefx_cpml.dat', b_x)
  
  call loadcpml('kappay_cpml.dat', K_y)
  call loadcpml('alphay_cpml.dat', alpha_y)
  call loadcpml('acoefy_cpml.dat', a_y)
  call loadcpml('bcoefy_cpml.dat', b_y)
  
  call loadcpml('kappaz_cpml.dat', K_z)
  call loadcpml('alphaz_cpml.dat', alpha_z)
  call loadcpml('acoefz_cpml.dat', a_z)
  call loadcpml('bcoefz_cpml.dat', b_z)
  
  call loadcpml('kappax_half_cpml.dat', K_x_half)
  call loadcpml('alphax_half_cpml.dat', alpha_x_half)
  call loadcpml('acoefx_half_cpml.dat', a_x_half)
  call loadcpml('bcoefx_half_cpml.dat', b_x_half)
  
  call loadcpml('kappay_half_cpml.dat', K_y_half)
  call loadcpml('alphay_half_cpml.dat', alpha_y_half)
  call loadcpml('acoefy_half_cpml.dat', a_y_half)
  call loadcpml('bcoefy_half_cpml.dat', b_y_half)
  
  call loadcpml('kappaz_half_cpml.dat', K_z_half)
  call loadcpml('alphaz_half_cpml.dat', alpha_z_half)
  call loadcpml('acoefz_half_cpml.dat', a_z_half)
  call loadcpml('bcoefz_half_cpml.dat', b_z_half)

  ! do i = 1,nz
  !   print *, K_z_half(i), alpha_z_half(i), a_z_half(i), b_z_half(i)
  ! enddo
  ! do i = 1,ny
  !   print *, K_y(i), alpha_y(i), a_y(i), b_y(i)
  ! enddo
    
    
! =============================== Forward Model ===============================
! initialize arrays
vx(:,:,:) = 0.d0
vy(:,:,:) = 0.d0
vz(:,:,:) = 0.d0

sigmaxx(:,:,:) = 0.d0
sigmayy(:,:,:) = 0.d0
sigmazz(:,:,:) = 0.d0
sigmaxy(:,:,:) = 0.d0
sigmaxz(:,:,:) = 0.d0
sigmayz(:,:,:) = 0.d0

! PML
memory_dvx_dx(:,:,:) = 0.d0
memory_dvx_dy(:,:,:) = 0.d0
memory_dvx_dz(:,:,:) = 0.d0

memory_dvy_dx(:,:,:) = 0.d0
memory_dvy_dy(:,:,:) = 0.d0
memory_dvy_dz(:,:,:) = 0.d0

memory_dvz_dx(:,:,:) = 0.d0
memory_dvz_dy(:,:,:) = 0.d0 
memory_dvz_dz(:,:,:) = 0.d0

memory_dsigmaxx_dx(:,:,:) = 0.d0
memory_dsigmayy_dy(:,:,:) = 0.d0
memory_dsigmazz_dz(:,:,:) = 0.d0

memory_dsigmaxy_dx(:,:,:) = 0.d0
memory_dsigmaxy_dy(:,:,:) = 0.d0
memory_dsigmaxz_dx(:,:,:) = 0.d0
memory_dsigmaxz_dz(:,:,:) = 0.d0
memory_dsigmayz_dy(:,:,:) = 0.d0
memory_dsigmayz_dz(:,:,:) = 0.d0

! Do it 
do it = 1,NSTEP
  !------------------------------------------------------------
  ! compute stress sigma and update memory variables for C-PML
  !------------------------------------------------------------
  ! Update in the x direction
  do k = 2,nz
    do j = 2,NY
      do i = 1,NX-1

        dvx_dx = (vx(i+1,j,k) - vx(i,j,k) ) / dx
        dvy_dx = (vy(i+1,j,k) - vy(i,j,k) ) / dx
        dvz_dx = (vz(i+1,j,k) - vz(i,j,k) ) / dx 
        dvy_dy = (vy(i,j,k) - vy(i,j-1,k) ) / dy
        dvx_dy = (vx(i,j,k) - vx(i,j-1,k) ) / dy
        dvz_dy = (vz(i,j,k) - vz(i,j-1,k) ) / dy
        dvz_dz = (vz(i,j,k) - vz(i,j,k-1) ) / dz
        dvx_dz = (vx(i,j,k) - vx(i,j,k-1) ) / dz
        dvy_dz = (vy(i,j,k) - vy(i,j,k-1) ) / dz
        
        memory_dvx_dx(i,j,k) = b_x_half(i) * memory_dvx_dx(i,j,k) + a_x_half(i) * dvx_dx
        memory_dvy_dx(i,j,k) = b_x_half(i) * memory_dvy_dx(i,j,k) + a_x_half(i) * dvy_dx
        memory_dvz_dx(i,j,k) = b_x_half(i) * memory_dvz_dx(i,j,k) + a_x_half(i) * dvz_dx
        memory_dvy_dy(i,j,k) = b_y(j) * memory_dvy_dy(i,j,k) + a_y(j) * dvy_dy
        memory_dvx_dy(i,j,k) = b_y(j) * memory_dvx_dy(i,j,k) + a_y(j) * dvx_dy
        memory_dvz_dy(i,j,k) = b_y(j) * memory_dvz_dy(i,j,k) + a_y(j) * dvz_dy
        memory_dvz_dz(i,j,k) = b_z(k) * memory_dvz_dz(i,j,k) + a_z(k) * dvz_dz
        memory_dvx_dz(i,j,k) = b_z(k) * memory_dvx_dz(i,j,k) + a_z(k) * dvx_dz
        memory_dvy_dz(i,j,k) = b_z(k) * memory_dvy_dz(i,j,k) + a_z(k) * dvy_dz

        dvx_dx = dvx_dx / K_x_half(i) + memory_dvx_dx(i,j,k)
        dvy_dx = dvy_dx / K_x_half(i) + memory_dvy_dx(i,j,k)
        dvz_dx = dvz_dx / K_x_half(i) + memory_dvz_dx(i,j,k)
        dvy_dy = dvy_dy / K_y(j) + memory_dvy_dy(i,j,k)
        dvx_dy = dvx_dy / K_y(j) + memory_dvx_dy(i,j,k)
        dvz_dy = dvz_dy / K_y(j) + memory_dvz_dy(i,j,k)
        dvz_dz = dvz_dz / K_z(k) + memory_dvz_dz(i,j,k)
        dvx_dz = dvx_dz / K_z(k) + memory_dvx_dz(i,j,k)
        dvy_dz = dvy_dz / K_z(k) + memory_dvy_dz(i,j,k)

        sigmaxx(i,j,k) = sigmaxx(i,j,k) + &
            ( ( ( c11(i+1,k) + c11(i,k) + c11(i,k-1) + c11(i+1,k-1))/4 ) * dvx_dx + &
              ( ( c12(i+1,k) + c12(i,k) + c12(i,k-1) + c12(i+1,k-1))/4 ) * dvy_dy + &
              ( ( c13(i+1,k) + c13(i,k) + c13(i,k-1) + c13(i+1,k-1))/4 ) * dvz_dz + &
              ( ( c14(i+1,k) + c14(i,k) + c14(i,k-1) + c14(i+1,k-1))/8 ) * (dvy_dz + dvz_dy) + &
              ( ( c15(i+1,k) + c15(i,k) + c15(i,k-1) + c15(i+1,k-1))/8 ) * (dvx_dz + dvz_dx) + &
              ( ( c16(i+1,k) + c16(i,k) + c16(i,k-1) + c16(i+1,k-1))/8 ) * (dvx_dy + dvz_dy) ) * dt

        ! Full 3D will need a gradient in the y-direction
        sigmayy(i,j,k) = sigmayy(i,j,k) + &
            ( ( ( c12(i+1,k) + c12(i,k) + c12(i,k-1) + c12(i+1,k-1))/4 ) * dvx_dx + &
              ( ( c22(i+1,k) + c22(i,k) + c22(i,k-1) + c22(i+1,k-1))/4 ) * dvy_dy + &
              ( ( c23(i+1,k) + c23(i,k) + c23(i,k-1) + c23(i+1,k-1))/4 ) * dvz_dz + &
              ( ( c24(i+1,k) + c24(i,k) + c24(i,k-1) + c24(i+1,k-1))/8 ) * (dvy_dz + dvz_dy) + &
              ( ( c25(i+1,k) + c25(i,k) + c25(i,k-1) + c25(i+1,k-1))/8 ) * (dvx_dz + dvz_dx) + &
              ( ( c26(i+1,k) + c26(i,k) + c26(i,k-1) + c26(i+1,k-1))/8 ) * (dvy_dx + dvx_dy) ) * dt

        sigmazz(i,j,k) = sigmazz(i,j,k) + &
            ( ( ( c13(i+1,k) + c13(i,k) + c13(i,k-1) + c13(i+1,k-1))/4 ) * dvx_dx + &
              ( ( c23(i+1,k) + c23(i,k) + c23(i,k-1) + c23(i+1,k-1))/4 ) * dvy_dy + &
              ( ( c33(i+1,k) + c33(i,k) + c33(i,k-1) + c33(i+1,k-1))/4 ) * dvz_dz + &
              ( ( c34(i+1,k) + c34(i,k) + c34(i,k-1) + c34(i+1,k-1))/8 ) * (dvy_dz + dvz_dy) + &
              ( ( c35(i+1,k) + c35(i,k) + c35(i,k-1) + c35(i+1,k-1))/8 ) * (dvx_dz + dvz_dx) + &
              ( ( c36(i+1,k) + c36(i,k) + c36(i,k-1) + c36(i+1,k-1))/8 ) * (dvy_dx + dvx_dy) ) * dt

      enddo
    enddo
  enddo

  ! Update sigmaxy, x-direction is full nodes
  do k = 2,nz
    do j = 1,ny-1
      do i = 2,nx

        dvx_dx = (vx(i,j,k) - vx(i-1,j,k)) / DX
        dvy_dx = (vy(i,j,k) - vy(i-1,j,k)) / DX
        dvz_dx = (vz(i,j,k) - vz(i-1,j,k)) / DX
        
        dvy_dy = (vy(i,j+1,k) - vy(i,j,k)) / DY
        dvx_dy = (vx(i,j+1,k) - vx(i,j,k)) / DY
        dvz_dy = (vz(i,j+1,k) - vz(i,j,k)) / DY
        
        dvz_dz = (vz(i,j,k) - vz(i,j,k-1)) / DZ
        dvx_dz = (vx(i,j,k) - vx(i,j,k-1)) / DZ
        dvy_dz = (vy(i,j,k) - vy(i,j,k-1)) / DZ

        memory_dvx_dx(i,j,k) = b_x(i) * memory_dvx_dx(i,j,k) + a_x(i) * dvx_dx
        memory_dvy_dx(i,j,k) = b_x(i) * memory_dvy_dx(i,j,k) + a_x(i) * dvy_dx
        memory_dvz_dx(i,j,k) = b_x(i) * memory_dvz_dx(i,j,k) + a_x(i) * dvz_dx
        memory_dvy_dy(i,j,k) = b_y_half(j) * memory_dvy_dy(i,j,k) + a_y_half(j) * dvy_dy
        memory_dvx_dy(i,j,k) = b_y_half(j) * memory_dvx_dy(i,j,k) + a_y_half(j) * dvx_dy
        memory_dvz_dy(i,j,k) = b_y_half(j) * memory_dvz_dy(i,j,k) + a_y_half(j) * dvz_dy
        memory_dvz_dz(i,j,k) = b_z_half(k) * memory_dvz_dz(i,j,k) + a_z_half(k) * dvz_dz
        memory_dvx_dz(i,j,k) = b_z_half(k) * memory_dvx_dz(i,j,k) + a_z_half(k) * dvx_dz
        memory_dvy_dz(i,j,k) = b_z_half(k) * memory_dvy_dz(i,j,k) + a_z_half(k) * dvy_dz

        dvx_dx = dvx_dx / K_x(i) + memory_dvx_dx(i,j,k)
        dvy_dx = dvy_dx / K_x(i) + memory_dvy_dx(i,j,k)
        dvy_dy = dvy_dy / K_y_half(j) + memory_dvy_dy(i,j,k)
        dvx_dy = dvx_dy / K_y_half(j) + memory_dvx_dy(i,j,k)
        dvz_dy = dvz_dy / K_y_half(j) + memory_dvz_dy(i,j,k)
        dvz_dz = dvz_dz / K_z_half(k) + memory_dvz_dz(i,j,k)
        dvy_dz = dvy_dz / K_z_half(k) + memory_dvy_dz(i,j,k)

        sigmaxy(i,j,k) = sigmaxy(i,j,k) + &
            ( ( ( c16(i,k) + c16(i-1,k) + c16(i,k-1) + c16(i-1,k-1))/4 ) * dvx_dx + &
              ( ( c26(i,k) + c26(i-1,k) + c26(i,k-1) + c26(i-1,k-1))/4 ) * dvy_dy + &
              ( ( c36(i,k) + c36(i-1,k) + c36(i,k-1) + c36(i-1,k-1))/4 ) * dvz_dz + &
              ( ( c46(i,k) + c46(i-1,k) + c46(i,k-1) + c46(i-1,k-1))/8 ) * (dvz_dy + dvy_dz) + &
              ( ( c56(i,k) + c56(i-1,k) + c56(i,k-1) + c56(i-1,k-1))/8 ) * (dvz_dx + dvx_dz) + &
              ( ( c66(i,k) + c66(i-1,k) + c66(i,k-1) + c66(i-1,k-1))/8 ) * (dvy_dx + dvx_dy) ) * dt

      enddo
    enddo
  enddo

  ! Update sigmaxz, z-direction is full nodes
  do k = 1,nz-1
    do j = 2,ny
      do i = 2,nx

        dvx_dx = (vx(i,j,k) - vx(i-1,j,k)) / DX
        dvy_dx = (vy(i,j,k) - vy(i-1,j,k)) / DX
        dvz_dx = (vz(i,j,k) - vz(i-1,j,k)) / DX
        dvy_dy = (vy(i,j,k) - vy(i,j-1,k)) / DY
        dvz_dy = (vz(i,j,k) - vz(i,j-1,k)) / DY
        dvx_dy = (vx(i,j,k) - vx(i,j-1,k)) / DY
        dvz_dz = (vz(i,j,k+1) - vz(i,j,k)) / DZ
        dvx_dz = (vx(i,j,k+1) - vx(i,j,k)) /DZ
        dvy_dz = (vy(i,j,k+1) - vy(i,j,k)) / DZ

        memory_dvx_dx(i,j,k) = b_x(i) * memory_dvx_dx(i,j,k) + a_x(i) * dvx_dx
        memory_dvy_dx(i,j,k) = b_x(i) * memory_dvy_dx(i,j,k) + a_x(i) * dvy_dx
        memory_dvz_dx(i,j,k) = b_x(i) * memory_dvz_dx(i,j,k) + a_x(i) * dvz_dx
        memory_dvy_dy(i,j,k) = b_y(j) * memory_dvy_dy(i,j,k) + a_y(j) * dvy_dy
        memory_dvx_dy(i,j,k) = b_y(j) * memory_dvx_dy(i,j,k) + a_y(j) * dvx_dy
        memory_dvz_dy(i,j,k) = b_y(j) * memory_dvz_dy(i,j,k) + a_y(j) * dvz_dy
        memory_dvz_dz(i,j,k) = b_z_half(k) * memory_dvz_dz(i,j,k) + a_z_half(k) * dvz_dz
        memory_dvx_dz(i,j,k) = b_z_half(k) * memory_dvx_dz(i,j,k) + a_z_half(k) * dvx_dz
        memory_dvy_dz(i,j,k) = b_z_half(k) * memory_dvy_dz(i,j,k) + a_z_half(k) * dvy_dz

        dvx_dx = dvx_dx / K_x(i) + memory_dvx_dx(i,j,k)
        dvy_dx = dvy_dx / K_x(i) + memory_dvy_dx(i,j,k)
        dvz_dx = dvz_dx / K_x(i) + memory_dvz_dx(i,j,k) 
        dvy_dy = dvy_dy / K_y(j) + memory_dvy_dy(i,j,k)
        dvx_dy = dvx_dy / K_y(j) + memory_dvx_dy(i,j,k)
        dvz_dy = dvz_dy / K_y(j) + memory_dvz_dy(i,j,k)
        dvz_dz = dvz_dz / K_z_half(k) + memory_dvz_dz(i,j,k)
        dvx_dz = dvx_dz / K_z_half(k) + memory_dvx_dz(i,j,k)
        dvy_dz = dvy_dz / K_z_half(k) + memory_dvy_dz(i,j,k)

        sigmaxz(i,j,k) = sigmaxz(i,j,k) + &
        ( ( ( c15(i,k+1) + c15(i,k) + c15(i-1,k) + c15(i-1,k+1))/4 ) * dvx_dx + &
          ( ( c25(i,k+1) + c25(i,k) + c25(i-1,k) + c25(i-1,k+1))/4 ) * dvy_dy + &
          ( ( c35(i,k+1) + c35(i,k) + c35(i-1,k) + c35(i-1,k+1))/4 ) * dvz_dz + &
          ( ( c45(i,k+1) + c45(i,k) + c45(i-1,k) + c45(i-1,k+1))/8 ) * ( dvx_dz + dvz_dx) + &
          ( ( c55(i,k+1) + c55(i,k) + c55(i-1,k) + c55(i-1,k+1))/8 ) * ( dvx_dz + dvz_dx) + &
          ( ( c56(i,k+1) + c56(i,k) + c56(i-1,k) + c56(i-1,k+1))/8 ) * ( dvx_dy + dvy_dx) ) * dt 

      enddo
    enddo

  !   ! update sigmayz, y-direction is full nodes
    do j = 1,ny-1
      do i = 1,nx-1
        
        dvx_dx = (vx(i+1,j,k) - vx(i,j,k)) / DX
        dvy_dx = (vy(i+1,j,k) - vy(i,j,k)) / DX
        dvz_dx = (vz(i+1,j,k) - vz(i,j,k)) / DX
        dvy_dy = (vy(i,j+1,k) - vy(i,j,k)) / DY
        dvx_dy = (vx(i,j+1,k) - vx(i,j,k)) / DY
        dvz_dy = (vz(i,j+1,k) - vz(i,j,k)) / DY
        dvz_dz = (vz(i,j,k+1) - vz(i,j,k)) / DZ
        dvx_dz = (vx(i,j,k+1) - vx(i,j,k)) / DZ 
        dvy_dz = (vy(i,j,k+1) - vy(i,j,k)) / DZ

        memory_dvx_dx(i,j,k) = b_x_half(i) * memory_dvx_dx(i,j,k) + a_x_half(i) * dvx_dx
        memory_dvy_dx(i,j,k) = b_x_half(i) * memory_dvy_dx(i,j,k) + a_x_half(i) * dvy_dx
        memory_dvz_dx(i,j,k) = b_x_half(i) * memory_dvz_dx(i,j,k) + a_x_half(i) * dvz_dx
        memory_dvy_dy(i,j,k) = b_y_half(j) * memory_dvy_dy(i,j,k) + a_y_half(j) * dvy_dy
        memory_dvx_dy(i,j,k) = b_y_half(j) * memory_dvx_dy(i,j,k) + a_y_half(j) * dvx_dy
        memory_dvz_dy(i,j,k) = b_y_half(j) * memory_dvz_dy(i,j,k) + a_y_half(j) * dvz_dy
        memory_dvz_dz(i,j,k) = b_z_half(k) * memory_dvz_dz(i,j,k) + a_z_half(k) * dvz_dz
        memory_dvx_dz(i,j,k) = b_z_half(k) * memory_dvx_dz(i,j,k) + a_z_half(k) * dvx_dz
        memory_dvy_dz(i,j,k) = b_z_half(k) * memory_dvy_dz(i,j,k) + a_z_half(k) * dvy_dz
        
        dvx_dx = dvx_dx / K_x_half(i) + memory_dvx_dx(i,j,k)
        dvy_dx = dvy_dx / K_x_half(i) + memory_dvy_dx(i,j,k)
        dvz_dx = dvz_dx / K_x_half(i) + memory_dvz_dx(i,j,k)
        dvy_dy = dvy_dy / K_y_half(j) + memory_dvy_dy(i,j,k)
        dvx_dy = dvx_dy / K_y_half(j) + memory_dvx_dy(i,j,k)
        dvz_dy = dvz_dy / K_y_half(j) + memory_dvz_dy(i,j,k)
        dvz_dz = dvz_dz / K_z_half(k) + memory_dvz_dz(i,j,k)
        dvx_dz = dvx_dz / K_z_half(k) + memory_dvx_dz(i,j,k)
        dvy_dz = dvy_dz / K_z_half(k) + memory_dvy_dz(i,j,k)
        
        sigmayz(i,j,k) = sigmayz(i,j,k)  + &
        ( ( ( c14(i+1,k) + c14(i,k) + c14(i,k+1) + c14(i+1,k+1))/4 ) * dvx_dx + &
          ( ( c24(i+1,k) + c24(i,k) + c24(i,k+1) + c24(i+1,k+1))/4 ) * dvy_dy + &
          ( ( c34(i+1,k) + c34(i,k) + c34(i,k+1) + c34(i+1,k+1))/4 ) * dvz_dz + &
          ( ( c44(i+1,k) + c44(i,k) + c44(i,k+1) + c44(i+1,k+1))/8 ) * ( dvy_dz + dvz_dy) + &
          ( ( c45(i+1,k) + c45(i,k) + c45(i,k+1) + c45(i+1,k+1))/8 ) * ( dvx_dz + dvz_dx) + &
          ( ( c46(i+1,k) + c46(i,k) + c46(i,k+1) + c46(i+1,k+1))/8 ) * ( dvy_dx + dvx_dy) ) * dt 

      enddo
    enddo
  enddo


!--------------------------------------------------------
! compute velocity and update memory variables for C-PML
!--------------------------------------------------------
  do k = 2,nz
    do j = 2,NY
      do i = 2,NX
        ! ds1/dx, ds6/dy, ds5,dz
        deltarho = (4 * rho(i,k) + rho(i-1,k) + rho(i,k-1) )/6
        
        dsigmaxx_dx = (sigmaxx(i,j,k) - sigmaxx(i-1,j,k) ) / dx
        dsigmaxy_dy = (sigmaxy(i,j,k) - sigmaxy(i,j-1,k) ) / dy
        dsigmaxz_dz = (sigmaxz(i,j,k) - sigmaxz(i,j,k-1) ) / dz

        memory_dsigmaxx_dx(i,j,k) = b_x(i) * &
                  memory_dsigmaxx_dx(i,j,k) + a_x(i) * dsigmaxx_dx
        memory_dsigmaxy_dy(i,j,k) = b_y(j) * &
                  memory_dsigmaxy_dy(i,j,k) + a_y(j) * dsigmaxy_dy
        memory_dsigmaxz_dz(i,j,k) = b_z(k) * &
                  memory_dsigmaxz_dz(i,j,k) + a_z(k) * dsigmaxz_dz

        dsigmaxx_dx = dsigmaxx_dx / K_x(i) + memory_dsigmaxx_dx(i,j,k)
        dsigmaxy_dy = dsigmaxy_dy / K_y(j) + memory_dsigmaxy_dy(i,j,k)
        dsigmaxz_dz = dsigmaxz_dz / K_z(k) + memory_dsigmaxz_dz(i,j,k) 

        vx(i,j,k) = vx(i,j,k) + &
            (dsigmaxx_dx + dsigmaxy_dy + dsigmaxz_dz) * &
            dt / deltarho !rho(i,k)

      enddo
    enddo

    do j = 1,ny-1
      do i = 1,nx-1
        ! ds6/dx, ds2/dy, ds4/dz
        deltarho = (4*rho(i,k) + rho(i+1,k) + rho(i,k-1) )/6

        dsigmaxy_dx = ( sigmaxy(i+1,j,k) - sigmaxy(i,j,k) ) / dx
        dsigmayy_dy = ( sigmayy(i,j+1,k) - sigmayy(i,j,k) ) / dy
        dsigmayz_dz = ( sigmayz(i,j,k) - sigmayz(i,j,k-1) ) / dz

        memory_dsigmaxy_dx(i,j,k) = b_x_half(i) * memory_dsigmaxy_dx(i,j,k) + a_x_half(i) * dsigmaxy_dx
        memory_dsigmayy_dy(i,j,k) = b_y_half(j) * memory_dsigmayy_dy(i,j,k) + a_y_half(j) * dsigmayy_dy
        memory_dsigmayz_dz(i,j,k) = b_z(k) * memory_dsigmayz_dz(i,j,k) + a_z(k) * dsigmayz_dz

        dsigmaxy_dx = dsigmaxy_dx / K_x_half(i) + memory_dsigmaxy_dx(i,j,k)
        dsigmayy_dy = dsigmayy_dy / K_y_half(j) + memory_dsigmayy_dy(i,j,k)
        dsigmayz_dz = dsigmayz_dz / K_z(k) + memory_dsigmayz_dz(i,j,k)

        vy(i,j,k) = vy(i,j,k) + &
            (dsigmaxy_dx + dsigmayy_dy + dsigmayz_dz) * &
            dt / deltarho !rho(i,k)

      enddo
    enddo
  enddo

  do k = 1,nz-1
    do j = 2,ny
      do i = 1,nx-1
        ! ds5/dx, ds4/dy, ds3/dz
        deltarho = ( rho(i+1,k) + rho(i,k+1) + 4*rho(i,k) )/6

        dsigmaxz_dx = ( sigmaxz(i+1,j,k) - sigmaxz(i,j,k) ) / dx
        dsigmayz_dy = ( sigmayz(i,j,k) - sigmayz(i,j-1,k) ) / dy
        dsigmazz_dz = ( sigmazz(i,j,k+1) - sigmazz(i,j,k) ) / dz
        
        memory_dsigmaxz_dx(i,j,k) = b_x_half(i) * memory_dsigmaxz_dx(i,j,k) + a_x_half(i) * dsigmaxz_dx
        memory_dsigmayz_dy(i,j,k) = b_y(j) * memory_dsigmayz_dy(i,j,k) + a_y(j) * dsigmayz_dy
        memory_dsigmazz_dz(i,j,k) = b_z_half(k) * memory_dsigmazz_dz(i,j,k) + a_z_half(k) * dsigmazz_dz

        dsigmaxz_dx = dsigmaxz_dx / K_x_half(i) + memory_dsigmaxz_dx(i,j,k)
        dsigmayz_dy = dsigmayz_dy / K_y(j) + memory_dsigmayz_dy(i,j,k)
        dsigmazz_dz = dsigmazz_dz / K_z_half(k) + memory_dsigmazz_dz(i,j,k)

        vz(i,j,k) = vz(i,j,k) + &
            (dsigmaxz_dx + dsigmayz_dy + dsigmazz_dz) * &
            dt / deltarho !rho(i,k)

      enddo
    enddo
  enddo

  vx(isource,jsource,ksource) = vx(isource,jsource,ksource) + srcx(it) * DT / rho(isource,ksource)
  vy(isource,jsource,ksource) = vy(isource,jsource,ksource) + srcy(it) * DT / rho(isource,ksource)
  vz(isource,jsource,ksource) = vz(isource,jsource,ksource) + srcz(it) * DT / rho(isource,ksource)

  ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
  vx(1,:,:) = 0.d0
  vy(1,:,:) = 0.d0
  vz(1,:,:) = 0.d0

  vx(:,1,:) = 0.d0
  vy(:,1,:) = 0.d0
  vz(:,1,:) = 0.d0
  
  vx(:,:,1) = 0.d0
  vy(:,:,1) = 0.d0
  vz(:,:,1) = 0.d0
  
  vx(NX,:,:) = 0.d0
  vy(NX,:,:) = 0.d0
  vz(NX,:,:) = 0.d0
  
  vx(:,NY,:) = 0.d0
  vy(:,NY,:) = 0.d0
  vz(:,NY,:) = 0.d0
  
  vx(:,:,NZ) = 0.d0
  vy(:,:,NZ) = 0.d0
  vz(:,:,NZ) = 0.d0

  ! check norm of velocity to make sure the solution isn't diverging
  velocnorm = maxval( sqrt(vx**2 + vy**2 + vz**2) )
  ! print *,'Time step # ',it,' out of ',NSTEP
  ! print *,'Time: ',(it-1)*DT,' seconds'
  print *,'Max vals for vx, vy, vz: ', maxval(vx), maxval(vy), maxval(vz)

  if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'

  ! Write the velocity values to an unformatted binary file
  call write_image(vx, nx, ny, nz, it, 'Vx')
  call write_image(vy, nx, ny, nz, it, 'Vy')
  call write_image(vz, nx, ny, nz, it, 'Vz')
  ! Now write the stress Values
  ! call write_image(sigmaxx, nx, ny, nz, it, 'S1')
  ! call write_image(sigmayy, nx, ny, nz, it, 'S2')
  ! call write_image(sigmazz, nx, ny, nz, it, 'S3')
  ! call write_image(sigmaxy, nx, ny, nz, it, 'S6')
  ! call write_image(sigmayz, nx, ny, nz, it, 'S4')
  ! call write_image(sigmaxz, nx, ny, nz, it, 'S5')

enddo   ! end of time loop

end subroutine seismic_cpml_25d


!==============================================================================
!==============================================================================
!==============================================================================
subroutine write_image(image_data, nx, ny, nz, it, channel)
! Write the 3D array out as single precision

implicit none

integer, parameter :: dp = kind(0.d0)
integer :: nx, ny, nz, it
real(kind=dp) :: image_data(nx, ny, nz)
character(len=2) :: channel
character(len=100) :: filename

WRITE (filename, "(a2, a1, i6.6, '.dat')" ) channel, '.', it

open(unit = 10, form = "unformatted", file = trim(filename) )
write(10) sngl(image_data)

close(unit = 10)

end subroutine write_image


end module seismicFDTD25d