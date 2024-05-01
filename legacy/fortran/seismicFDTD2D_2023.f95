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
! module seismicFDTD2d

!   implicit none
  
!   contains
    

    !============================================================================
    
    
    subroutine seismic_cpml_2d(nx, nz, dx, dz, npoints_pml, src, nstep)
    
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
      ! IMPORTANT : all our CPML codes work fine in single precision as well (which 
      !             is significantly faster). If you want you can thus force 
      !             automatic conversion to single precision at compile time or 
      !             change all the declarations and constants in the code from 
      !             real(kind=dp) to single.
      !
      ! INPUT
      !   im (INTEGER)
      !   nx, ny (INTEGER)
      !   c11, c12, c22, c66, rho (REAL)
      !   dx, dy (REAL)
      !   npoints_pml (INTEGER) - the thickness of the pml
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
        implicit none
    
        integer, parameter :: dp=kind(0.d0)
    
        ! total number of grid points in each direction of the grid
        integer :: nx
        integer :: nz
    
        ! thickness of the PML layer in grid points
        integer :: npoints_pml
        ! integer, dimension(nx,nz)
        real(kind=dp), dimension(nx,nz) :: c11, c13, c15, c33, c35, c55, rho
        real(kind=dp) :: deltarho
    
        ! total number of time steps
        integer :: nstep
    
        ! time step in seconds 
        real(kind=dp) :: DT
        real(kind=dp) :: dx, dz 
        ! parameters for the source
        real(kind=dp), parameter :: factor = 1.d7
    
        ! source
        integer,dimension(:) :: src
        integer :: isource, jsource
    
        ! value of PI
        real(kind=dp), parameter :: PI = 3.141592653589793238462643d0
    
        ! conversion from degrees to radians
        real(kind=dp), parameter :: DEGREES_TO_RADIANS = PI / 180.d0
    
        ! large value for maximum
        real(kind=dp), parameter :: HUGEVAL = 1.d+30
    
        ! velocity threshold above which we consider that the code became unstable
        real(kind=dp), parameter :: STABILITY_THRESHOLD = 1.d+25
    
        ! main arrays
        real(kind=dp), dimension(nx,nz) :: vx,vz,sigmaxx,sigmazz,sigmaxz
    
        ! power to compute d0 profile. Increasing this value allows for a larger dampening gradient in the PML
        real(kind=dp), parameter :: NPOWER = 2.d0
    
        ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
        real(kind=dp), parameter :: K_MAX_PML = 1.d1
        
        ! arrays for the memory variables
        ! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
        real(kind=dp), dimension(NX,NZ) :: &
            memory_dvx_dx, &
            memory_dvx_dz, &
            memory_dvz_dx, &
            memory_dvz_dz, &
            memory_dsigmaxx_dx, &
            memory_dsigmazz_dz, &
            memory_dsigmaxz_dx, &
            memory_dsigmaxz_dz
    
        real(kind=dp) :: &
            value_dvx_dx, &
            value_dvx_dz, &
            value_dvz_dx, &
            value_dvz_dz, &
            value_dsigmaxx_dx, &
            value_dsigmazz_dz, &
            value_dsigmaxz_dx, &
            value_dsigmaxz_dz
    
        ! 1D arrays for the damping profiles
        real(kind=dp), dimension(nx) :: K_x, alpha_x, a_x, b_x, &
                                    K_x_half, alpha_x_half, &
                                    a_x_half, b_x_half
        real(kind=dp), dimension(nz) ::K_z, alpha_z, a_z, b_z, &
                                    K_z_half, alpha_z_half, &
                                    a_z_half, b_z_half
    
        ! for the source
        real(kind=dp),dimension(nstep) :: srcx, srcz
    
        integer :: i,j,it
    
        real(kind=dp) :: velocnorm

        ! Name the f2py inputs 
        !f2py3 intent(in) :: nx, nx, dx, dx,
        !f2py3 intent(in) :: noints_pml, src, nstep
    
    
        ! -------------------- Load Stiffness Coefficients --------------------
    
        call material_rw('c11.dat', c11, .TRUE.)
        call material_rw('c13.dat', c13, .TRUE.)
        call material_rw('c15.dat', c15, .TRUE.)
        call material_rw('c33.dat', c33, .TRUE.)
        call material_rw('c35.dat', c35, .TRUE.)
        call material_rw('c55.dat', c55, .TRUE.)
        call material_rw('rho.dat', rho, .TRUE.)
    
        ! ------------------------ Assign some constants -----------------------
    
        isource = src(1)+npoints_pml
        jsource = src(2)+npoints_pml
    
        DT = minval( (/dx,dz/) )/ &
            (sqrt( 3.d0*( maxval( (/ c11/rho, c33/rho /) ) ) ) ) 

            ! ================================ LOAD SOURCE ================================
    
        call loadsource('seismicsourcex.dat', nstep, srcx)
        ! We are using the coordinate names x, Z but the math computes the source in 
        ! the x-z plane
        call loadsource('seismicsourcez.dat', nstep, srcz)
    
        ! -----------------------------------------------------------------------------
        !--- define profile of absorption in PML region

        ! Initialize PML 
        K_x(:) = 1.d0
        K_x_half(:) = 1.d0
        alpha_x(:) = 0.d0
        alpha_x_half(:) = 0.d0
        a_x(:) = 0.d0
        a_x_half(:) = 0.d0
        b_x(:) = 0.d0 
        b_x_half(:) = 0.d0

        K_z(:) = 1.d0
        K_z_half(:) = 1.d0 
        alpha_z(:) = 0.d0
        alpha_z_half(:) = 0.d0
        a_z(:) = 0.d0
        a_z_half(:) = 0.d0
        b_z(:) = 0.d0
        b_z_half(:) = 0.d0

      ! ------------------------------ Load the boundary ----------------------------        
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
    
        ! =============================================================================
    
    
        ! initialize arrays
        vx(:,:) = 0.d0
        vz(:,:) = 0.d0
        sigmaxx(:,:) = 0.d0
        sigmazz(:,:) = 0.d0
        sigmaxz(:,:) = 0.d0
    
        ! PML
        memory_dvx_dx(:,:) = 0.d0
        memory_dvx_dz(:,:) = 0.d0
        memory_dvz_dx(:,:) = 0.d0
        memory_dvz_dz(:,:) = 0.d0
        memory_dsigmaxx_dx(:,:) = 0.d0
        memory_dsigmazz_dz(:,:) = 0.d0
        memory_dsigmaxz_dx(:,:) = 0.d0
        memory_dsigmaxz_dz(:,:) = 0.d0
    
        !---
        !---  beginning of time loop
        !---
    
        do it = 1,NSTEP
          !------------------------------------------------------------
          ! compute stress sigma and update memory variables for C-PML
          !------------------------------------------------------------
          do j = 2,NZ
            do i = 1,NX-1
    
              value_dvx_dx = (vx(i+1,j) - vx(i,j)) / DX
              value_dvz_dz = (vz(i,j) - vz(i,j-1)) / DZ
              value_dvz_dx = (vz(i+1,j) - vz(i,j)) / DX
              value_dvx_dz = (vx(i,j) - vx(i,j-1)) / DZ

              memory_dvx_dx(i,j) = b_x_half(j) * memory_dvx_dx(i,j) + &
                                    a_x_half(i) * value_dvx_dx
              memory_dvz_dz(i,j) = b_z(j) * memory_dvz_dz(i,j) + &
                                    a_z(j) * value_dvz_dz
              memory_dvx_dz(i,j) = b_z_half(j) * memory_dvx_dz(i,j) + &
                                    a_z_half(j) * value_dvx_dz 
              memory_dvz_dx(i,j) = b_x(i) * memory_dvz_dx(i,j) + &
                                    a_x(i) * value_dvz_dx

              value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j)
              value_dvz_dz = value_dvz_dz / K_z(j) + memory_dvz_dz(i,j)
              value_dvz_dx = value_dvz_dx / K_x(i) + memory_dvz_dx(i,j)
              value_dvx_dz = value_dvx_dz / K_z_half(j) + memory_dvx_dz(i,j)
              
              sigmaxx(i,j) = sigmaxx(i,j) + &
                ( ( ( c11(i+1,j) + 2*c11(i,j) + c11(i,j-1) )/4) * value_dvx_dx + &
                  ( ( c13(i+1,j) + 2*c13(i,j) + c13(i,j-1) )/4) * value_dvz_dz + &
                  ( ( c15(i+1,j) + 2*c15(i,j) + c15(i,j-1) )/4) * &
                          (value_dvz_dx + value_dvx_dz) ) * DT
              sigmazz(i,j) = sigmazz(i,j) + &
                ( ( ( c13(i+1,j) + 2*c13(i,j) + c13(i,j-1) )/4) * value_dvx_dx + &
                  ( ( c33(i+1,j) + 2*c33(i,j) + c33(i,j-1) )/4) * value_dvz_dz + &
                  ( ( c35(i+1,j) + 2*c35(i,j) + c35(i,j-1) )/4) * &
                          (value_dvz_dx + value_dvx_dz) ) * DT
    
            enddo
          enddo
    
          do j = 1,NZ-1
            do i = 2,NX
    
              value_dvx_dx = (vx(i,j) - vx(i-1,j)) / DX
              value_dvz_dz = (vz(i,j+1) - vz(i,j)) / DZ
              value_dvz_dx = (vz(i,j) - vz(i-1,j)) / DX
              value_dvx_dz = (vx(i,j+1) - vx(i,j)) / DZ
              
              memory_dvx_dx(i,j) = b_x_half(i) * memory_dvx_dx(i,j) + &
                                    a_x_half(i) * value_dvx_dx
              memory_dvz_dz(i,j) = b_z(j) * memory_dvz_dz(i,j) + &
                                    a_z(j) * value_dvz_dz
              memory_dvx_dz(i,j) = b_z_half(j) * memory_dvx_dz(i,j) + &
                                    a_z_half(j) * value_dvx_dz 
              memory_dvz_dx(i,j) = b_x(i) * memory_dvz_dx(i,j) + &
                                    a_x(i) * value_dvz_dx
              
              value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j)
              value_dvz_dz = value_dvz_dz / K_z(j) + memory_dvz_dz(i,j)
              value_dvz_dx = value_dvz_dx / K_x(i) + memory_dvz_dx(i,j)
              value_dvx_dz = value_dvx_dz / K_z_half(j) + memory_dvx_dz(i,j)
    
              sigmaxz(i,j) = sigmaxz(i,j) + &
                ( ( (c15(i,j+1) + 2*c15(i,j) + c15(i-1,j))/4) * value_dvx_dx + &
                  ( (c35(i,j+1) + 2*c35(i,j) + c35(i-1,j))/4) * value_dvz_dz + &
                  ( (c55(i,j+1) + 2*c55(i,j) + c55(i-1,j))/4) * &
                          (value_dvz_dx + value_dvx_dz) ) * DT
    
            enddo
          enddo
    
        !--------------------------------------------------------
        ! compute velocity and update memory variables for C-PML
        !--------------------------------------------------------
          do j = 2,NZ
            do i = 2,NX
    
              deltarho = ( 2*rho(i,j) + rho(i-1,j) + rho(i,j-1) )/4
              value_dsigmaxx_dx = (sigmaxx(i,j) - sigmaxx(i-1,j)) / DX
              value_dsigmaxz_dz = (sigmaxz(i,j) - sigmaxz(i,j-1)) / DZ
    
              memory_dsigmaxx_dx(i,j) = b_x(i) * memory_dsigmaxx_dx(i,j) + &
                        a_x(i) * value_dsigmaxx_dx
              memory_dsigmaxz_dz(i,j) = b_z(j) * memory_dsigmaxz_dz(i,j) + &
                        a_z(j) * value_dsigmaxz_dz
    
              value_dsigmaxx_dx = value_dsigmaxx_dx / K_x(i) + &
                        memory_dsigmaxx_dx(i,j)
              value_dsigmaxz_dz = value_dsigmaxz_dz / K_z(j) + &
                        memory_dsigmaxz_dz(i,j)
    
              vx(i,j) = vx(i,j) + (value_dsigmaxx_dx + value_dsigmaxz_dz) * DT / rho(i,j)
    
            enddo
          enddo
    
          do j = 1,NZ-1
            do i = 1,NX-1
    
              deltarho = ( 2*rho(i,j) + rho(i+1,j) + rho(i,j+1) )/4
              value_dsigmaxz_dx = (sigmaxz(i+1,j) - sigmaxz(i,j)) / DX
              value_dsigmazz_dz = (sigmazz(i,j+1) - sigmazz(i,j)) / DZ
    
              memory_dsigmaxz_dx(i,j) = b_x_half(i) * memory_dsigmaxz_dx(i,j) + &
                        a_x_half(i) * value_dsigmaxz_dx
              memory_dsigmazz_dz(i,j) = b_z_half(j) * memory_dsigmazz_dz(i,j) + &
                        a_z_half(j) * value_dsigmazz_dz
    
              value_dsigmaxz_dx = value_dsigmaxz_dx / K_x_half(i) + memory_dsigmaxz_dx(i,j)
              value_dsigmazz_dz = value_dsigmazz_dz / K_z_half(j) + memory_dsigmazz_dz(i,j)
    
              vz(i,j) = vz(i,j) + (value_dsigmaxz_dx + value_dsigmazz_dz) * DT / deltarho
    
            enddo
          enddo
    
          ! Add the source term
          vx(isource,jsource) = vx(isource,jsource) + srcx(it) * DT / rho(isource,jsource)
          vz(isource,jsource) = vz(isource,jsource) + srcz(it) * DT / rho(isource,jsource)
    
          ! Dirichlet conditions (rigid boundaries) on the edges or at the 
          ! bottom of the PML layers
          vx(1,:) = 0.d0
          vx(NX,:) = 0.d0
    
          vx(:,1) = 0.d0
          vx(:,NZ) = 0.d0
    
          vz(1,:) = 0.d0
          vz(NX,:) = 0.d0
    
          vz(:,1) = 0.d0
          vz(:,NZ) = 0.d0
    
          ! print maximum of norm of velocity
          velocnorm = maxval(sqrt(vx**2 + vz**2))
          if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'
    
          call write_image2(vx, nx, nz, src, it, 'Vx')
          call write_image2(vz, nx, nz, src, it, 'Vz')
    
        enddo   ! end of time loop
    end subroutine seismic_cpml_2d
    
    
    !==============================================================================
    ! subroutine write_image2(image_data, nx, nz, src, it, channel)
    
    ! implicit none
    
    ! integer, parameter :: dp = kind(0.d0)
    ! integer :: nx, nz, it
    ! integer,dimension(2) :: src
    ! real(kind=dp) :: image_data(nx, nz)
    ! character(len=2) :: channel
    ! character(len=100) :: filename
    
    ! WRITE (filename, "(a2, a1, i6.6, a1, i0, a1, i0, a4)" ) &
    !                 channel,'.',it,'.', src(1),'.',src(2),'.dat'
    ! open(unit = 10, form = 'unformatted', file = trim(filename) )
    ! write(10) sngl(image_data)
    
    ! close(unit = 10)
    
    ! end subroutine write_image2
    
    ! end module seismicFDTD2d
    
