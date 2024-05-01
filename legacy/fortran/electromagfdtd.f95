module electromagfdtd

    use readwrite_routines

    implicit none
        
    contains
    
    ! =========================================================================
    subroutine electromag2(nx, nz, dx, dz, npoints_pml, src, nstep, OUTPUT_SINGLE)

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
        real(kind=dp), parameter :: mu0 = 4.0d0*pi*1.0d-7, eps0 = 8.85418782d-12

        ! typical relative permeability of common materials is close to unity but for
        ! the more specific case we can edit the following line to input permeability 
        ! as a 2D array 
        real(kind=dp), parameter :: mu = 1.d0

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
        ! character(len=6) :: src_type
        integer :: i,j, it

        real(kind=dp) :: velocnorm

        ! -------------------------------- PML parameters 
        ! 1D arrays for the damping profiles
        real(kind=dp), dimension(nx) :: K_x,alpha_x,a_x,b_x, &
                                        K_x_half, alpha_x_half,a_x_half,b_x_half
        real(kind=dp), dimension(nz) :: K_z,alpha_z,a_z,b_z, &
                                        K_z_half, alpha_z_half,a_z_half,b_z_half
                                        
        logical :: SINGLE
        logical, intent(in), optional :: OUTPUT_SINGLE 

        ! ------------------- Name the f2py inputs 
        !f2py3 intent(in) :: nx, nz, dx, dz,
        !f2py3 intent(in) :: noints_pml, src, nstep
        
        ! The default data output is single precision unless OUTPUT_SINGLE is 
        ! set to .FALSE.
        if (present(OUTPUT_SINGLE)) then 
            SINGLE = OUTPUT_SINGLE 
        else
            SINGLE = .TRUE.
        endif
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

            call write_image2(Ex, nx-1, nz, src, it, 'Ex', SINGLE)
            call write_image2(Ez, nx, nz-1, src, it, 'Ez', SINGLE)
        enddo
    end subroutine electromag2

    ! =========================================================================
    subroutine electromag2c(nx, nz, dx, dz, npoints_pml, src, nstep, OUTPUT_SINGLE)

        implicit none

        integer,parameter :: dp=kind(0.d0)

        ! total number of grid points in each direction of the grid
        integer :: nx, nz
        ! integer, dimension(:,:) :: rcx 

        ! thickness of the PML layer in grid points
        integer :: npoints_pml, nstep

        ! integer, dimension(nx,nz)
        complex(kind=dp), dimension(nx,nz) :: eps11, eps13, eps33, epsilonx, epsilonz
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
        real(kind=dp), parameter :: mu0 = 4.0d0*pi*1.0d-7, eps0 = 8.85418782d-12

        ! typical relative permeability of common materials is close to unity but for
        ! the more specific case we can edit the following line to input permeability 
        ! as a 2D array 
        real(kind=dp), parameter :: mu = 1.d0

        ! E-field threshold above which we consider that the code became unstable
        real(kind=dp), parameter :: STABILITY_THRESHOLD = 1.d+25

        ! main arrays
        complex(kind=dp), dimension(nx-1,nz) :: Ex
        complex(kind=dp), dimension(nx,nz-1) :: Ez
        complex(kind=dp), dimension(nx-1,nz-1) :: Hy

        ! we will compute the coefficients for the finite difference scheme 
        complex(kind=dp), dimension(nx, nz) :: caEx, cbEx
        complex(kind=dp), dimension(nx, nz) :: caEz, cbEz
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
        ! 1D arrays for the damping profiles
        real(kind=dp), dimension(nx) :: K_x,alpha_x,a_x,b_x, &
                                        K_x_half, alpha_x_half,a_x_half,b_x_half
        real(kind=dp), dimension(nz) :: K_z,alpha_z,a_z,b_z, &
                                        K_z_half, alpha_z_half,a_z_half,b_z_half
                                        
        ! Boolean flag to save as double precision or single precision 
        logical :: SINGLE
        logical, intent(in), optional :: OUTPUT_SINGLE 

        ! ------------------- Name the f2py inputs 
        !f2py3 intent(in) :: nx, nz, dx, dz,
        !f2py3 intent(in) :: noints_pml, src, nstep

        ! =============================================================================
        ! The default data output is single precision unless OUTPUT_SINGLE is 
        ! set to .FALSE.
        if (present(OUTPUT_SINGLE)) then 
            SINGLE = OUTPUT_SINGLE 
        else
            SINGLE = .TRUE.
        endif

        ! ----------------------- Load Permittivity Coefficients ----------------------

        call material_rwc('eps11.dat', eps11, .TRUE.)
        call material_rwc('eps13.dat', eps13, .TRUE.)
        call material_rwc('eps33.dat', eps33, .TRUE.) ! We will change y to z soon
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
            Ex(1,:) = complex(0.0d0, 0.0d0)
            Ex(nx-1,:) = complex(0.0d0, 0.0d0)
            Ex(:,1) = complex(0.0d0, 0.0d0)
            Ex(:,nz) = complex(0.0d0, 0.0d0)

            Ez(1,:) = complex(0.0d0, 0.0d0)
            Ez(nx,:) = complex(0.0d0, 0.0d0)
            Ez(:,1) = complex(0.0d0, 0.0d0)
            Ez(:,nz-1) = complex(0.0d0, 0.0d0)

            Hy(1,:) = complex(0.0d0, 0.0d0)
            Hy(nx-1,:) = complex(0.0d0, 0.0d0)
            Hy(:,1) = complex(0.0d0, 0.0d0)
            Hy(:,nz-1) = complex(0.0d0, 0.0d0)

            ! print maximum of norm of velocity
            velocnorm = maxval(abs(sqrt(Ex**2 + Ez**2)))
            if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'
            ! print *,'Max vals for Ex, Ey, Ez: ', maxval(REAL(Ex)), maxval(REAL(Ez))

            call write_image2c(Ex, nx-1, nz, src, it, 'Ex', SINGLE)
            call write_image2c(Ez, nx, nz-1, src, it, 'Ez', SINGLE)

        enddo   ! end of time loop


    end subroutine electromag2c

    ! =========================================================================
    subroutine electromag25(nx, ny, nz, dx, dy, dz, npoints_pml, src, nstep, OUTPUT_SINGLE)
        
        implicit none

        integer,parameter :: dp=kind(0.0d0)

        ! total number of grid points in each direction of the grid
        integer :: nx, ny, nz
        ! time step in seconds. decreasing the time step improves the pml attenuation
        ! but it should be inversely proportional to the center frequency of the 
        ! source frequency 
        integer :: npoints_pml, nstep

        ! source
        integer,dimension(:) :: src
        integer :: isource,jsource,ksource


        real(kind=dp), dimension(nx,nz) :: eps11, eps22, eps33, &
                                            eps12, eps13, eps23, &
                                            sig11, sig22, sig33, &
                                            sig12, sig13, sig23, &
                                            epsilonx, epsilony, epsilonz, &
                                            sigmax, sigmay, sigmaz

        real(kind=dp) :: DT, dx, dy, dz

        ! value of PI
        real(kind=dp), parameter :: PI = 3.141592653589793238462643d0

        ! speed of mother fuckin' light 
        real(kind=dp), parameter :: Clight = 2.9979458d+8

        ! permability and permittivity of free space 
        real(kind=dp), parameter :: mu0 = 4.0d0*pi*1.0d-7, eps0 = 8.85418782d-12

        ! typical relative permeability of common materials is close to unity but for
        ! the more specific case we can edit the following line to input permeability 
        ! as a 2D array 
        real(kind=dp), parameter :: mu = 1.0d0

        ! E-field threshold above which we consider that the code became unstable
        real(kind=dp), parameter :: STABILITY_THRESHOLD = 1.0d+25

        ! main arrays
        real(kind=dp), dimension(nx-1,ny,nz) :: Ex
        real(kind=dp), dimension(nx,ny-1,nz) :: Ey
        real(kind=dp), dimension(nx,ny,nz-1) :: Ez
        real(kind=dp), dimension(nx,ny-1,nz-1) :: Hx
        real(kind=dp), dimension(nx-1,ny,nz-1) :: Hy
        real(kind=dp), dimension(nx-1,ny-1,nz) :: Hz


        ! we will compute the coefficients for the finite difference scheme 
        real(kind=dp), dimension(nx, nz) :: caEx, cbEx, &
                                            caEy, cbEy, &
                                            caEz, cbEz

        !real(kind=dp), dimension(nx+1, ny+1) :: caEy, cbEy
        real(kind=dp) :: daHx, dbHx, daHy, dbHy, daHz, dbHz

        real(kind=dp) ::  dEx_dy, dEy_dx, &
                        dEy_dz, dEz_dy, &
                        dEz_dx, dEx_dz, &
                        dHx_dy, dHx_dz, &
                        dHy_dx, dHy_dz, &
                        dHz_dy, dHz_dx

        ! arrays for the memory variables
        ! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
        real(kind=dp),dimension(nx,ny,nz) :: memory_dEy_dx, memory_dEx_dy, &
                                            memory_dEz_dx, memory_dEx_dz, &
                                            memory_dEy_dz, memory_dEz_dy

        real(kind=dp),dimension(nx,ny,nz) :: memory_dHz_dx, memory_dHx_dz, &
                                            memory_dHy_dx, memory_dHx_dy, &
                                            memory_dHy_dz, memory_dHz_dy

        ! parameters for the source
        ! angle of source force clockwise with respect to vertical (Y) axis
        ! this will later be treated as an input
        real(kind=dp),dimension(nstep) :: srcx, srcy, srcz

        integer :: i,j,k,it

        real(kind=dp) :: velocnorm

        ! -------------------------------- PML parameters 

        ! 1D arrays for the damping profiles
        real(kind=dp),dimension(nx) ::  K_x,alpha_x,a_x,b_x, & 
                                        K_x_half, alpha_x_half,a_x_half,b_x_half
        real(kind=dp),dimension(ny) ::  K_y,alpha_y,a_y,b_y, &
                                        K_y_half, alpha_y_half,a_y_half,b_y_half
        real(kind=dp),dimension(nz) ::  K_z,alpha_z,a_z,b_z, &
                                        K_z_half, alpha_z_half,a_z_half,b_z_half

        ! integer :: npml_x,npml_y, npml_z
        logical :: SINGLE
        logical, intent(in), optional :: OUTPUT_SINGLE 

        ! Name the f2py inputs 
        !f2py3 intent(in) :: nx, ny, nz, dx, dy, dz,
        !f2py3 intent(in) :: noints_pml, src, nstep
                
        ! The default data output is single precision unless OUTPUT_SINGLE is 
        ! set to .FALSE.
        if (present(OUTPUT_SINGLE)) then 
            SINGLE = OUTPUT_SINGLE 
        else
            SINGLE = .TRUE.
        endif
        
        ! ------------------------ Load Permittivity Coefficients ------------------------
        ! Load Epsilon
        call material_rw('eps11.dat', eps11, .TRUE.)
        call material_rw('eps12.dat', eps12, .TRUE.)
        call material_rw('eps13.dat', eps13, .TRUE.)
        call material_rw('eps22.dat', eps22, .TRUE.)
        call material_rw('eps23.dat', eps23, .TRUE.)
        call material_rw('eps33.dat', eps33, .TRUE.)
        ! Load Sigma
        call material_rw('sig11.dat', sig11, .TRUE.)
        call material_rw('sig12.dat', sig12, .TRUE.)
        call material_rw('sig13.dat', sig13, .TRUE.)
        call material_rw('sig22.dat', sig22, .TRUE.)
        call material_rw('sig23.dat', sig23, .TRUE.)
        call material_rw('sig33.dat', sig33, .TRUE.)


        ! ------------------------ Assign some constants -----------------------
        ! Assign the source location indices
        isource = int(src(1)) + npoints_pml
        jsource = int(src(2)) + npoints_pml
        ksource = int(src(3)) + npoints_pml

        ! Define the 
        DT = minval( (/dx, dy, dz/) )/ ( 2.0d0 * Clight/ sqrt( minval( (/ eps11, eps22, eps33 /) ) ) )

        ! Compute the coefficients of the FD scheme. First scale the relative 
        ! permittivity and permeabilities to get the absolute values 
        epsilonx(:,:) = (eps11 + eps12 + eps13)*eps0 
        epsilony(:,:) = (eps12 + eps22 + eps23)*eps0
        epsilonz(:,:) = (eps13 + eps23 + eps33)*eps0
        sigmaX(:,:) = sig11 + sig12 + sig13
        sigmay(:,:) = sig12 + sig22 + sig23
        sigmaz(:,:) = sig13 + sig23 + sig33

        caEx(:,:) = ( 1.0d0 - sigmaX * dt / ( 2.0d0 * epsilonx ) ) / &
                    ( 1.0d0 + sigmaX * dt / (2.0d0 * epsilonx ) )
        cbEx(:,:) = (dt / epsilonx ) / ( 1.0d0 + sigmax * dt / ( 2.0d0 * epsilonx ) )

        caEy(:,:) = ( 1.0d0 - sigmay * dt / ( 2.0d0 * epsilony ) ) / &
                    ( 1.0d0 + sigmay * dt / (2.0d0 * epsilony ) )
        cbEy(:,:) = (dt / epsilony ) / ( 1.0d0 + sigmay * dt / ( 2.0d0 * epsilony ) )

        caEz(:,:) = ( 1.0d0 - sigmaz * dt / ( 2.0d0 * epsilonz ) ) / &
                    ( 1.0d0 + sigmaz * dt / (2.0d0 * epsilonz ) )
        cbEz(:,:) = (dt / epsilonz ) / ( 1.0d0 + sigmaz * dt / ( 2.0d0 * epsilonz ) )

        daHx = dt/(4.0d0*mu0*mu)
        dbHx = dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
        daHx = 1.0d0 ! (1-daHz)/(1+daHz) ! 

        daHy = dt/(4.0d0*mu0*mu)
        dbHy = dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
        daHy = 1.0d0 ! (1-daHz)/(1+daHz) ! 

        daHz = dt/(4.0d0*mu0*mu)
        dbHz = dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
        daHz = 1.0d0 ! (1-daHz)/(1+daHz) ! 


        ! ----------------------------------------------------------------------
        !---
        !--- program starts here
        !---

        ! ================================ LOAD SOURCE ================================

        call loadsource('electromagneticsourcex.dat', nstep, srcx)
        call loadsource('electromagneticsourcey.dat', nstep, srcy)
        call loadsource('electromagneticsourcez.dat', nstep, srcz)

        ! =============================================================================

        !--- define profile of absorption in PML region

        ! Initialize CPML damping variables
        K_x(:) = 1.0d0
        K_x_half(:) = 1.0d0
        alpha_x(:) = 0.0d0
        alpha_x_half(:) = 0.0d0
        a_x(:) = 0.0d0
        a_x_half(:) = 0.0d0
        b_x(:) = 0.0d0 
        b_x_half(:) = 0.0d0 

        K_y(:) = 1.0d0
        K_y_half(:) = 1.0d0
        alpha_y(:) = 0.0d0
        alpha_y_half(:) = 0.0d0
        a_y(:) = 0.0d0
        a_y_half(:) = 0.0d0
        b_y(:) = 0.d0
        K_z(:) = 1.0d0
        K_z_half(:) = 1.0d0
        alpha_z(:) = 0.0d0
        alpha_z_half(:) = 0.0d0
        a_z(:) = 0.0d0
        a_z_half(:) = 0.0d0

        ! ------------------------------ Load the boundary ----------------------------
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
        !   print *, K_z(i), alpha_z(i), a_z(i), b_z(i)
        ! enddo

        ! -----------------------------------------------------------------------------

        ! initialize arrays
        Ex(:,:,:) = 0.0d0
        Ey(:,:,:) = 0.0d0
        Ez(:,:,:) = 0.0d0

        Hx(:,:,:) = 0.0d0
        Hy(:,:,:) = 0.0d0
        Hz(:,:,:) = 0.0d0


        ! PML
        memory_dEx_dy(:,:,:) = 0.0d0
        memory_dEy_dx(:,:,:) = 0.0d0
        memory_dEx_dz(:,:,:) = 0.0d0
        memory_dEz_dx(:,:,:) = 0.0d0
        memory_dEz_dy(:,:,:) = 0.0d0
        memory_dEy_dz(:,:,:) = 0.0d0

        memory_dHz_dx(:,:,:) = 0.0d0
        memory_dHx_dz(:,:,:) = 0.0d0
        memory_dHz_dy(:,:,:) = 0.0d0
        memory_dHy_dz(:,:,:) = 0.0d0
        memory_dHx_dy(:,:,:) = 0.0d0
        memory_dHy_dx(:,:,:) = 0.0d0

        ! ---
        ! ---  beginning of time loop
        ! ---
        do it = 1,NSTEP
            !--------------------------------------------------------
            ! compute magnetic field and update memory variables for C-PML
            !--------------------------------------------------------
            ! Update Hx
            do k = 1,nz-1
                do i = 1,nx-1  
                    do j = 1,ny-1
                        ! Values needed for the magnetic field updates
                        dEz_dy = ( Ez(i,j,k) - Ez(i,j+1,k) )/dy
                        memory_dEz_dy(i,j,k) = b_y_half(j) * memory_dEz_dy(i,j,k) + a_y_half(j) * dEz_dy
                        dEz_dy = dEz_dy/ K_y_half(j) + memory_dEz_dy(i,j,k)

                        ! The rest of the equation needed for agnetic field updates
                        dEy_dz = ( Ey(i,j,k+1) - Ey(i,j,k) )/dz
                        memory_dEy_dz(i,j,k) = b_z_half(k) * memory_dEy_dz(i,j,k) + a_z_half(k) * dEy_dz
                        dEy_dz = dEy_dz/ K_z_half(k) + memory_dEy_dz(i,j,k)

                        ! Now update the Magnetic field
                        Hx(i,j,k) = daHx*Hx(i,j,k) + dbHx*( dEy_dz + dEz_dy )
                    enddo
                enddo  
            enddo

                ! Update Hy
            do k = 1,nz-1
                do i = 1,nx-1      
                    do j = 1,ny-1
                    
                        ! Values needed for the magnetic field updates
                        dEx_dz = ( Ex(i,j,k) - Ex(i,j,k+1) )/dz
                        memory_dEx_dz(i,j,k) = b_z(k) * memory_dEx_dz(i,j,k) + a_z(k) * dEx_dz
                        dEx_dz = dEx_dz/ K_z(k) + memory_dEx_dz(i,j,k)

                        ! The rest of the equation needed for agnetic field updates
                        dEz_dx = ( Ez(i+1,j,k) - Ez(i,j,k) )/dx
                        memory_dEz_dx(i,j,k) = b_x(i) * memory_dEz_dx(i,j,k) + a_x(i) * dEz_dx
                        dEz_dx = dEz_dx/ K_x(i) + memory_dEz_dx(i,j,k)

                        ! Now update the Magnetic field
                        Hy(i,j,k) = daHy*Hy(i,j,k) + dbHy*( dEz_dx + dEx_dz )

                    enddo
                enddo  
            enddo

                ! Update Hz
            do k = 2,nz-1
                do i = 1,nx-1      
                    do j = 1,ny-1
                        ! Values needed for the magnetic field updates
                        dEx_dy = ( Ex(i,j+1,k) - Ex(i,j,k) )/dy
                        memory_dEx_dy(i,j,k) = b_y(j) * memory_dEx_dy(i,j,k) + & 
                            a_y(j) * dEx_dy
                        dEx_dy = dEx_dy/ K_y(j) + memory_dEx_dy(i,j,k)

                        ! The rest of the equation needed for agnetic field updates
                        dEy_dx = ( Ey(i,j,k) - Ey(i+1,j,k) )/dx
                        memory_dEy_dx(i,j,k) = b_x(i) * memory_dEy_dx(i,j,k) + & 
                            a_x(i) * dEy_dx
                        dEy_dx = dEy_dx/ K_x(i) + memory_dEy_dx(i,j,k)

                        ! Now update the Magnetic field
                        Hz(i,j,k) = daHz*Hz(i,j,k) + dbHz*( dEy_dx + dEx_dy )
                    enddo
                enddo  
            enddo

            !--------------------------------------------------------
            ! compute electric field and update memory variables for C-PML
            !--------------------------------------------------------
                ! Compute the differences in the x-direction
            do k = 2,nz-1
                do i = 1,nx-1
                    do j = 2,ny-1  
                        ! Update the Ex field
                        dHz_dy = ( Hz(i,j,k) - Hz(i,j-1,k) )/dy
                        memory_dHz_dy(i,j,k) = b_y_half(j) * memory_dHz_dy(i,j,k) + & 
                            a_y_half(j) * dHz_dy
                        dHz_dy = dHz_dy/K_y_half(j) + memory_dHz_dy(i,j,k)

                        ! Changed from half to full node positions 
                        dHy_dz = ( Hy(i,j,k-1) - Hy(i,j,k) )/dz
                        memory_dHy_dz(i,j,k) = b_z(k) * memory_dHy_dz(i,j,k) + &
                            a_z(k) * dHy_dz
                        dHy_dz = dHy_dz/K_z(k) + memory_dHy_dz(i,j,k)
                        
                        Ex(i,j,k) = caEx(i,k)*Ex(i,j,k) + & 
                        cbEx(i,k)*(dHz_dy + dHy_dz) 
                    enddo
                enddo

                ! ! Compute the differences in the y-direction
                do i = 2,nx-1 
                    do j = 1,ny-1 
                        ! Update the Ey field
                        dHz_dx = ( Hz(i-1,j,k) - Hz(i,j,k) )/dx ! this is ny-1 length vector
                        memory_dHz_dx(i,j,k) = b_x_half(i) * memory_dHz_dx(i,j,k) + & 
                            a_x_half(i) * dHz_dx
                        dHz_dx = dHz_dx/K_x_half(i) + memory_dHz_dx(i,j,k)

                        dHx_dz = ( Hx(i,j,k) - Hx(i,j,k-1) )/dz ! this is ny-1 length vector
                        memory_dHx_dz(i,j,k) = b_z_half(k) * memory_dHx_dz(i,j,k) + &
                            a_z_half(k) * dHx_dz
                        dHx_dz = dHx_dz/K_z_half(k) + memory_dHx_dz(i,j,k)

                        Ey(i,j,k) = ( ( 4*caEy(i,k) + caEy(i-1,k) + caEy(i,k-1) )/6) * Ey(i,j,k) + & 
                        ( ( 4*cbEy(i,k) + cbEy(i-1,k) + cbEy(i,k-1) )/6 ) * & 
                        (dHz_dx + dHx_dz)
                    enddo
                enddo
            enddo 

                ! Compute the differences in the z-direction
            do k = 1,nz-1
                do i = 2,nx-1  
                    do j = 2,ny-1
                        ! Update the Ez field
                        dHx_dy = ( Hx(i,j-1,k) - Hx(i,j,k) )/dy
                        memory_dHx_dy(i,j,k) = b_y_half(j) * memory_dHx_dy(i,j,k) + &
                            a_y_half(j) * dHx_dy
                        dHx_dy = dHx_dy/K_y_half(j) + memory_dHx_dy(i,j,k)

                        dHy_dx = ( Hy(i,j,k) - Hy(i-1,j,k) )/dx
                        memory_dHy_dx(i,j,k) = b_x_half(i) * memory_dHy_dx(i,j,k) + &
                            a_x_half(i) * dHy_dx
                        dHy_dx = dHy_dx/K_x_half(i) + memory_dHy_dx(i,j,k)
                        
                        Ez(i,j,k) = ( ( 4*caEz(i,k) + caEz(i-1,k) + caEz(i,k+1) )/6 ) * Ez(i,j,k) + & 
                        ( ( 4*cbEz(i,k) + cbEz(i-1,k) + cbEz(i,k+1) )/6 ) * & 
                        (dHx_dy + dHy_dx)
                    enddo
                enddo
            enddo


            ! add the source (force vector located at a given grid point)
            Ex(isource,jsource,ksource) = Ex(isource,jsource,ksource) + srcx(it) * DT / eps11(isource,ksource)
            Ey(isource,jsource,ksource) = Ey(isource,jsource,ksource) + srcy(it) * DT / eps22(isource,ksource) 
            Ez(isource,jsource,ksource) = Ez(isource,jsource,ksource) + srcz(it) * DT / eps33(isource,ksource)
            
            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            Ex(1,:,:) = 0.0d0
            Ex(:,1,:) = 0.0d0
            Ex(:,:,1) = 0.0d0
            Ex(nx-1,:,:) = 0.0d0
            Ex(:,ny,:) = 0.0d0
            Ex(:,:,nz) = 0.0d0 

            Ey(1,:,:) = 0.0d0
            Ey(:,1,:) = 0.0d0
            Ey(:,:,1) = 0.0d0
            Ey(nx,:,:) = 0.0d0
            Ey(:,ny-1,:) = 0.0d0
            Ey(:,:,nz) = 0.0d0
            
            Ez(1,:,:) = 0.0d0
            Ez(:,1,:) = 0.0d0
            Ez(:,:,1) = 0.0d0
            Ez(nx,:,:) = 0.0d0
            Ez(:,ny,:) = 0.0d0
            Ez(:,:,nz-1) = 0.0d0
            
            Hx(1,:,:) = 0.0d0
            Hx(:,1,:) = 0.0d0
            Hx(:,:,1) = 0.0d0
            Hx(nx,:,:) = 0.0d0
            Hx(:,ny-1,:) = 0.0d0
            Hx(:,:,nz-1) = 0.0d0

            Hy(1,:,:) = 0.0d0
            Hy(:,1,:) = 0.0d0
            Hy(:,:,1) = 0.0d0
            Hy(nx-1,:,:) = 0.0d0
            Hy(:,ny,:) = 0.0d0
            Hy(:,:,nz-1) = 0.0d0
            
            Hz(1,:,:) = 0.0d0
            Hz(:,1,:) = 0.0d0
            Hz(:,:,1) = 0.0d0
            Hz(nx-1,:,:) = 0.0d0
            Hz(:,ny-1,:) = 0.0d0
            Hz(:,:,nz) = 0.0d0

            ! check norm of velocity to make sure the solution isn't diverging
            velocnorm = maxval(sqrt(Ex**2.0d0 + Ey**2.0d0 + Ez**2.0d0) )
            if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'
            print *,'Max vals for Ex, Ey, Ez: ', maxval(Ex), maxval(Ey), maxval(Ez)

            ! print *, maxval(Ex), maxval(Ey), maxval(Ez)
            call write_image3(Ex, nx-1, ny, nz, src, it, 'Ex', SINGLE)
            call write_image3(Ey, nx, ny-1, nz, src, it, 'Ey', SINGLE)
            call write_image3(Ez, nx, ny, nz-1, src, it, 'Ez', SINGLE)

        enddo   ! end of time loop

    end subroutine electromag25
    
    ! =========================================================================
    subroutine electromag25c(nx, ny, nz, dx, dy, dz, npoints_pml, src, nstep, OUTPUT_SINGLE)
        
        implicit none

        integer,parameter :: dp=kind(0.0d0)

        ! total number of grid points in each direction of the grid
        integer :: nx, ny, nz
        ! time step in seconds. decreasing the time step improves the pml attenuation
        ! but it should be inversely proportional to the center frequency of the 
        ! source frequency 
        integer :: npoints_pml, nstep

        ! source
        integer,dimension(:) :: src
        integer :: isource,jsource,ksource


        complex(kind=dp), dimension(nx,nz) :: eps11, eps22, eps33, &
                                            eps12, eps13, eps23, &
                                            epsilonx, epsilony, epsilonz
                                            
        real(kind=dp), dimension(nx, nz) :: sig11, sig22, sig33, &
                                            sig12, sig13, sig23, &
                                            sigmax, sigmay, sigmaz

        real(kind=dp) :: DT, dx, dy, dz

        ! value of PI
        real(kind=dp), parameter :: PI = 3.141592653589793238462643d0

        ! speed of mother fuckin' light 
        real(kind=dp), parameter :: Clight = 2.9979458d+8

        ! permability and permittivity of free space 
        real(kind=dp), parameter :: mu0 = 4.0d0*pi*1.0d-7, eps0 = 8.85418782d-12

        ! typical relative permeability of common materials is close to unity but for
        ! the more specific case we can edit the following line to input permeability 
        ! as a 2D array 
        real(kind=dp), parameter :: mu = 1.0d0

        ! E-field threshold above which we consider that the code became unstable
        real(kind=dp), parameter :: STABILITY_THRESHOLD = 1.0d+25

        ! main arrays
        complex(kind=dp), dimension(nx-1,ny,nz) :: Ex
        complex(kind=dp), dimension(nx,ny-1,nz) :: Ey
        complex(kind=dp), dimension(nx,ny,nz-1) :: Ez
        complex(kind=dp), dimension(nx,ny-1,nz-1) :: Hx
        complex(kind=dp), dimension(nx-1,ny,nz-1) :: Hy
        complex(kind=dp), dimension(nx-1,ny-1,nz) :: Hz


        ! we will compute the coefficients for the finite difference scheme 
        complex(kind=dp), dimension(nx, nz) :: caEx, cbEx, &
                                            caEy, cbEy, &
                                            caEz, cbEz

        !real(kind=dp), dimension(nx+1, ny+1) :: caEy, cbEy
        real(kind=dp) :: daHx, dbHx, daHy, dbHy, daHz, dbHz

        complex(kind=dp) ::  dEx_dy, dEy_dx, &
                        dEy_dz, dEz_dy, &
                        dEz_dx, dEx_dz, &
                        dHx_dy, dHx_dz, &
                        dHy_dx, dHy_dz, &
                        dHz_dy, dHz_dx

        ! arrays for the memory variables
        ! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
        ! complex(kind=dp),dimension(nx,ny,nz) :: memory_dEy_dx, memory_dEx_dy, &
        !                                     memory_dEz_dx, memory_dEx_dz, &
        !                                     memory_dEy_dz, memory_dEz_dy

        ! complex(kind=dp),dimension(nx,ny,nz) :: memory_dHz_dx, memory_dHx_dz, &
        !                                     memory_dHy_dx, memory_dHx_dy, &
        !                                     memory_dHy_dz, memory_dHz_dy

        ! parameters for the source
        ! angle of source force clockwise with respect to vertical (Y) axis
        ! this will later be treated as an input
        real(kind=dp),dimension(nstep) :: srcx, srcy, srcz

        integer :: i,j,k,it

        real(kind=dp) :: velocnorm

        ! -------------------------------- PML parameters 

        ! 1D arrays for the damping profiles
        real(kind=dp),dimension(nx) ::  K_x,alpha_x,a_x,b_x, & 
                                        K_x_half, alpha_x_half,a_x_half,b_x_half
        real(kind=dp),dimension(ny) ::  K_y,alpha_y,a_y,b_y, &
                                        K_y_half, alpha_y_half,a_y_half,b_y_half
        real(kind=dp),dimension(nz) ::  K_z,alpha_z,a_z,b_z, &
                                        K_z_half, alpha_z_half,a_z_half,b_z_half

        ! integer :: npml_x,npml_y, npml_z
        logical :: SINGLE
        logical, intent(in), optional :: OUTPUT_SINGLE 

        ! Name the f2py inputs 
        !f2py3 intent(in) :: nx, ny, nz, dx, dy, dz,
        !f2py3 intent(in) :: noints_pml, src, nstep
                
        ! The default data output is single precision unless OUTPUT_SINGLE is 
        ! set to .FALSE.
        if (present(OUTPUT_SINGLE)) then 
            SINGLE = OUTPUT_SINGLE 
        else
            SINGLE = .TRUE.
        endif
        
        ! ------------------------ Load Permittivity Coefficients ------------------------
        ! Load Epsilon
        call material_rwc('eps11.dat', eps11, .TRUE.)
        call material_rwc('eps12.dat', eps12, .TRUE.)
        call material_rwc('eps13.dat', eps13, .TRUE.)
        call material_rwc('eps22.dat', eps22, .TRUE.)
        call material_rwc('eps23.dat', eps23, .TRUE.)
        call material_rwc('eps33.dat', eps33, .TRUE.)
        ! Load Sigma
        call material_rw('sig11.dat', sig11, .TRUE.)
        call material_rw('sig12.dat', sig12, .TRUE.)
        call material_rw('sig13.dat', sig13, .TRUE.)
        call material_rw('sig22.dat', sig22, .TRUE.)
        call material_rw('sig23.dat', sig23, .TRUE.)
        call material_rw('sig33.dat', sig33, .TRUE.)


        ! ------------------------ Assign some constants -----------------------
        ! Assign the source location indices
        isource = int(src(1)) + npoints_pml
        jsource = int(src(2)) + npoints_pml
        ksource = int(src(3)) + npoints_pml

        ! Define the 
        DT = minval( (/dx, dy, dz/) )/ ( 2.0d0 * Clight/ sqrt( minval( (/ eps11, eps22, eps33 /) ) ) )

        ! Compute the coefficients of the FD scheme. First scale the relative 
        ! permittivity and permeabilities to get the absolute values 
        epsilonx(:,:) = (eps11 + eps12 + eps13)*eps0 
        epsilony(:,:) = (eps12 + eps22 + eps23)*eps0
        epsilonz(:,:) = (eps13 + eps23 + eps33)*eps0
        sigmaX(:,:) = sig11 + sig12 + sig13
        sigmay(:,:) = sig12 + sig22 + sig23
        sigmaz(:,:) = sig13 + sig23 + sig33

        caEx(:,:) = ( 1.0d0 - sigmaX * dt / ( 2.0d0 * epsilonx ) ) / &
                    ( 1.0d0 + sigmaX * dt / (2.0d0 * epsilonx ) )
        cbEx(:,:) = (dt / epsilonx ) / ( 1.0d0 + sigmax * dt / ( 2.0d0 * epsilonx ) )

        caEy(:,:) = ( 1.0d0 - sigmay * dt / ( 2.0d0 * epsilony ) ) / &
                    ( 1.0d0 + sigmay * dt / (2.0d0 * epsilony ) )
        cbEy(:,:) = (dt / epsilony ) / ( 1.0d0 + sigmay * dt / ( 2.0d0 * epsilony ) )

        caEz(:,:) = ( 1.0d0 - sigmaz * dt / ( 2.0d0 * epsilonz ) ) / &
                    ( 1.0d0 + sigmaz * dt / (2.0d0 * epsilonz ) )
        cbEz(:,:) = (dt / epsilonz ) / ( 1.0d0 + sigmaz * dt / ( 2.0d0 * epsilonz ) )

        daHx = dt/(4.0d0*mu0*mu)
        dbHx = dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
        daHx = 1.0d0 ! (1-daHz)/(1+daHz) ! 

        daHy = dt/(4.0d0*mu0*mu)
        dbHy = dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
        daHy = 1.0d0 ! (1-daHz)/(1+daHz) ! 

        daHz = dt/(4.0d0*mu0*mu)
        dbHz = dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
        daHz = 1.0d0 ! (1-daHz)/(1+daHz) ! 


        ! ----------------------------------------------------------------------
        !---
        !--- program starts here
        !---

        ! ================================ LOAD SOURCE ================================

        call loadsource('electromagneticsourcex.dat', nstep, srcx)
        call loadsource('electromagneticsourcey.dat', nstep, srcy)
        call loadsource('electromagneticsourcez.dat', nstep, srcz)

        ! =============================================================================

        !--- define profile of absorption in PML region

        ! Initialize CPML damping variables
        K_x(:) = 1.0d0
        K_x_half(:) = 1.0d0
        alpha_x(:) = 0.0d0
        alpha_x_half(:) = 0.0d0
        a_x(:) = 0.0d0
        a_x_half(:) = 0.0d0
        b_x(:) = 0.0d0 
        b_x_half(:) = 0.0d0 

        K_y(:) = 1.0d0
        K_y_half(:) = 1.0d0
        alpha_y(:) = 0.0d0
        alpha_y_half(:) = 0.0d0
        a_y(:) = 0.0d0
        a_y_half(:) = 0.0d0
        b_y(:) = 0.d0
        K_z(:) = 1.0d0
        K_z_half(:) = 1.0d0
        alpha_z(:) = 0.0d0
        alpha_z_half(:) = 0.0d0
        a_z(:) = 0.0d0
        a_z_half(:) = 0.0d0

        ! ------------------------------ Load the boundary ----------------------------
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
        !   print *, K_z(i), alpha_z(i), a_z(i), b_z(i)
        ! enddo

        ! -----------------------------------------------------------------------------

        ! initialize arrays
        Ex(:,:,:) = complex(0.0d0, 0.0d0)
        Ey(:,:,:) = complex(0.0d0, 0.0d0)
        Ez(:,:,:) = complex(0.0d0, 0.0d0)

        Hx(:,:,:) = complex(0.0d0, 0.0d0)
        Hy(:,:,:) = complex(0.0d0, 0.0d0)
        Hz(:,:,:) = complex(0.0d0, 0.0d0)


        ! PML
        dEx_dy(:,:,:) = complex(0.0d0, 0.0d0)
        dEy_dx(:,:,:) = complex(0.0d0, 0.0d0)
        dEx_dz(:,:,:) = complex(0.0d0, 0.0d0)
        dEz_dx(:,:,:) = complex(0.0d0, 0.0d0)
        dEz_dy(:,:,:) = complex(0.0d0, 0.0d0)
        dEy_dz(:,:,:) = complex(0.0d0, 0.0d0)

        dHz_dx(:,:,:) = complex(0.0d0, 0.0d0)
        dHx_dz(:,:,:) = complex(0.0d0, 0.0d0)
        dHz_dy(:,:,:) = complex(0.0d0, 0.0d0)
        dHy_dz(:,:,:) = complex(0.0d0, 0.0d0)
        dHx_dy(:,:,:) = complex(0.0d0, 0.0d0)
        dHy_dx(:,:,:) = complex(0.0d0, 0.0d0)

        ! ---
        ! ---  beginning of time loop
        ! ---
        do it = 1,NSTEP
            !--------------------------------------------------------
            ! compute magnetic field and update memory variables for C-PML
            !--------------------------------------------------------
            ! Update Hx
            do k = 1,nz-1
                do i = 1,nx-1  
                    do j = 1,ny-1
                        ! Values needed for the magnetic field updates
                        dEz_dy = ( Ez(i,j,k) - Ez(i,j+1,k) )/dy
                        dEz_dy(i,j,k) = b_y_half(j) * dEz_dy(i,j,k) + a_y_half(j) * dEz_dy
                        dEz_dy = dEz_dy/ K_y_half(j) + dEz_dy(i,j,k)

                        ! The rest of the equation needed for agnetic field updates
                        dEy_dz = ( Ey(i,j,k+1) - Ey(i,j,k) )/dz
                        dEy_dz(i,j,k) = b_z_half(k) * dEy_dz(i,j,k) + a_z_half(k) * dEy_dz
                        dEy_dz = dEy_dz/ K_z_half(k) + dEy_dz(i,j,k)

                        ! Now update the Magnetic field
                        Hx(i,j,k) = daHx*Hx(i,j,k) + dbHx*( dEy_dz + dEz_dy )
                    enddo
                enddo  
            enddo

                ! Update Hy
            do k = 1,nz-1
                do i = 1,nx-1      
                    do j = 1,ny-1
                    
                        ! Values needed for the magnetic field updates
                        dEx_dz = ( Ex(i,j,k) - Ex(i,j,k+1) )/dz
                        dEx_dz(i,j,k) = b_z(k) * dEx_dz(i,j,k) + a_z(k) * dEx_dz
                        dEx_dz = dEx_dz/ K_z(k) + dEx_dz(i,j,k)

                        ! The rest of the equation needed for agnetic field updates
                        dEz_dx = ( Ez(i+1,j,k) - Ez(i,j,k) )/dx
                        dEz_dx(i,j,k) = b_x(i) * dEz_dx(i,j,k) + a_x(i) * dEz_dx
                        dEz_dx = dEz_dx/ K_x(i) + dEz_dx(i,j,k)

                        ! Now update the Magnetic field
                        Hy(i,j,k) = daHy*Hy(i,j,k) + dbHy*( dEz_dx + dEx_dz )

                    enddo
                enddo  
            enddo

                ! Update Hz
            do k = 2,nz-1
                do i = 1,nx-1      
                    do j = 1,ny-1
                        ! Values needed for the magnetic field updates
                        dEx_dy = ( Ex(i,j+1,k) - Ex(i,j,k) )/dy
                        dEx_dy(i,j,k) = b_y(j) * dEx_dy(i,j,k) + & 
                            a_y(j) * dEx_dy
                        dEx_dy = dEx_dy/ K_y(j) + dEx_dy(i,j,k)

                        ! The rest of the equation needed for agnetic field updates
                        dEy_dx = ( Ey(i,j,k) - Ey(i+1,j,k) )/dx
                        dEy_dx(i,j,k) = b_x(i) * dEy_dx(i,j,k) + & 
                            a_x(i) * dEy_dx
                        dEy_dx = dEy_dx/ K_x(i) + dEy_dx(i,j,k)

                        ! Now update the Magnetic field
                        Hz(i,j,k) = daHz*Hz(i,j,k) + dbHz*( dEy_dx + dEx_dy )
                    enddo
                enddo  
            enddo

            !--------------------------------------------------------
            ! compute electric field and update memory variables for C-PML
            !--------------------------------------------------------
                ! Compute the differences in the x-direction
            do k = 2,nz-1
                do i = 1,nx-1
                    do j = 2,ny-1  
                        ! Update the Ex field
                        dHz_dy = ( Hz(i,j,k) - Hz(i,j-1,k) )/dy
                        dHz_dy(i,j,k) = b_y_half(j) * dHz_dy(i,j,k) + & 
                            a_y_half(j) * dHz_dy
                        dHz_dy = dHz_dy/K_y_half(j) + dHz_dy(i,j,k)

                        ! Changed from half to full node positions 
                        dHy_dz = ( Hy(i,j,k-1) - Hy(i,j,k) )/dz
                        dHy_dz(i,j,k) = b_z(k) * dHy_dz(i,j,k) + &
                            a_z(k) * dHy_dz
                        dHy_dz = dHy_dz/K_z(k) + dHy_dz(i,j,k)
                        
                        Ex(i,j,k) = caEx(i,k)*Ex(i,j,k) + & 
                        cbEx(i,k)*(dHz_dy + dHy_dz) 
                    enddo
                enddo

                ! ! Compute the differences in the y-direction
                do i = 2,nx-1 
                    do j = 1,ny-1 
                        ! Update the Ey field
                        dHz_dx = ( Hz(i-1,j,k) - Hz(i,j,k) )/dx ! this is ny-1 length vector
                        dHz_dx(i,j,k) = b_x_half(i) * dHz_dx(i,j,k) + & 
                            a_x_half(i) * dHz_dx
                        dHz_dx = dHz_dx/K_x_half(i) + dHz_dx(i,j,k)

                        dHx_dz = ( Hx(i,j,k) - Hx(i,j,k-1) )/dz ! this is ny-1 length vector
                        dHx_dz(i,j,k) = b_z_half(k) * dHx_dz(i,j,k) + &
                            a_z_half(k) * dHx_dz
                        dHx_dz = dHx_dz/K_z_half(k) + dHx_dz(i,j,k)

                        Ey(i,j,k) = ( ( 4*caEy(i,k) + caEy(i-1,k) + caEy(i,k-1) )/6) * Ey(i,j,k) + & 
                        ( ( 4*cbEy(i,k) + cbEy(i-1,k) + cbEy(i,k-1) )/6 ) * & 
                        (dHz_dx + dHx_dz)
                    enddo
                enddo
            enddo 

                ! Compute the differences in the z-direction
            do k = 1,nz-1
                do i = 2,nx-1  
                    do j = 2,ny-1
                        ! Update the Ez field
                        dHx_dy = ( Hx(i,j-1,k) - Hx(i,j,k) )/dy
                        dHx_dy(i,j,k) = b_y_half(j) * dHx_dy(i,j,k) + &
                            a_y_half(j) * dHx_dy
                        dHx_dy = dHx_dy/K_y_half(j) + dHx_dy(i,j,k)

                        dHy_dx = ( Hy(i,j,k) - Hy(i-1,j,k) )/dx
                        dHy_dx(i,j,k) = b_x_half(i) * dHy_dx(i,j,k) + &
                            a_x_half(i) * dHy_dx
                        dHy_dx = dHy_dx/K_x_half(i) + dHy_dx(i,j,k)
                        
                        Ez(i,j,k) = ( ( 4*caEz(i,k) + caEz(i-1,k) + caEz(i,k+1) )/6 ) * Ez(i,j,k) + & 
                        ( ( 4*cbEz(i,k) + cbEz(i-1,k) + cbEz(i,k+1) )/6 ) * & 
                        (dHx_dy + dHy_dx)
                    enddo
                enddo
            enddo


            ! add the source (force vector located at a given grid point)
            Ex(isource,jsource,ksource) = Ex(isource,jsource,ksource) + &
                                        srcx(it) * DT / eps11(isource,ksource)
            Ey(isource,jsource,ksource) = Ey(isource,jsource,ksource) + &
                                        srcy(it) * DT / eps22(isource,ksource) 
            Ez(isource,jsource,ksource) = Ez(isource,jsource,ksource) + &
                                        srcz(it) * DT / eps33(isource,ksource)
            
            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            Ex(1,:,:) = complex(0.0d0, 0.0d0)
            Ex(:,1,:) = complex(0.0d0, 0.0d0)
            Ex(:,:,1) = complex(0.0d0, 0.0d0)
            Ex(nx-1,:,:) = complex(0.0d0, 0.0d0)
            Ex(:,ny,:) = complex(0.0d0, 0.0d0)
            Ex(:,:,nz) = complex(0.0d0, 0.0d0)

            Ey(1,:,:) = complex(0.0d0, 0.0d0)
            Ey(:,1,:) = complex(0.0d0, 0.0d0)
            Ey(:,:,1) = complex(0.0d0, 0.0d0)
            Ey(nx,:,:) = complex(0.0d0, 0.0d0)
            Ey(:,ny-1,:) = complex(0.0d0, 0.0d0)
            Ey(:,:,nz) = complex(0.0d0, 0.0d0)
            
            Ez(1,:,:) = complex(0.0d0, 0.0d0)
            Ez(:,1,:) = complex(0.0d0, 0.0d0)
            Ez(:,:,1) = complex(0.0d0, 0.0d0)
            Ez(nx,:,:) = complex(0.0d0, 0.0d0)
            Ez(:,ny,:) = complex(0.0d0, 0.0d0)
            Ez(:,:,nz-1) = complex(0.0d0, 0.0d0)
            
            Hx(1,:,:) = complex(0.0d0, 0.0d0)
            Hx(:,1,:) = complex(0.0d0, 0.0d0)
            Hx(:,:,1) = complex(0.0d0, 0.0d0)
            Hx(nx,:,:) = complex(0.0d0, 0.0d0)
            Hx(:,ny-1,:) = complex(0.0d0, 0.0d0)
            Hx(:,:,nz-1) = complex(0.0d0, 0.0d0)

            Hy(1,:,:) = complex(0.0d0, 0.0d0)
            Hy(:,1,:) = complex(0.0d0, 0.0d0)
            Hy(:,:,1) = complex(0.0d0, 0.0d0)
            Hy(nx-1,:,:) = complex(0.0d0, 0.0d0)
            Hy(:,ny,:) = complex(0.0d0, 0.0d0)
            Hy(:,:,nz-1) = complex(0.0d0, 0.0d0)
            
            Hz(1,:,:) = complex(0.0d0, 0.0d0)
            Hz(:,1,:) = complex(0.0d0, 0.0d0)
            Hz(:,:,1) = complex(0.0d0, 0.0d0)
            Hz(nx-1,:,:) = complex(0.0d0, 0.0d0)
            Hz(:,ny-1,:) = complex(0.0d0, 0.0d0)
            Hz(:,:,nz) = complex(0.0d0, 0.0d0)

            ! check norm of velocity to make sure the solution isn't diverging
            velocnorm = maxval(sqrt(Ex**2.0d0 + Ey**2.0d0 + Ez**2.0d0) )
            if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'
            print *,'Max vals for Ex, Ey, Ez: ', maxval(Ex), maxval(Ey), maxval(Ez)

            ! print *, maxval(Ex), maxval(Ey), maxval(Ez)
            call write_image3c(Ex, nx-1, ny, nz, src, it, 'Ex', SINGLE)
            call write_image3c(Ey, nx, ny-1, nz, src, it, 'Ey', SINGLE)
            call write_image3c(Ez, nx, ny, nz-1, src, it, 'Ez', SINGLE)

        enddo   ! end of time loop

    end subroutine electromag25c

end module electromagfdtd