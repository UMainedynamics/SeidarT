module seismicfdtd
    
    use readwrite_routines
    
    implicit none
    
    contains

    ! =========================================================================    
    ! Computations are done in double precision and written to binary as single
    ! precision unless specified by the optional logical, OUTPUT_SINGLE.
    subroutine seismic2(nx, nz, dx, dz, npoints_pml, src, nstep, OUTPUT_SINGLE)

        ! 2D elastic finite-difference code in velocity and stress formulation
        ! with Convolutional-PML (C-PML) absorbing conditions for an 
        ! anisotropic medium
        !
        ! If using this program please give credit to the following: 
        !
        ! Dimitri Komatitsch, University of Pau, France, April 2007.
        ! Anisotropic implementation by Roland Martin and Dimitri Komatitsch, 
        ! University of Pau, France, April 2007.
        
        ! The second-order staggered-grid formulation of Madariaga (1976) and 
        ! Virieux (1986) is used:
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
        ! IMPORTANT : all our CPML codes work fine in single precision as well 
        ! (which is significantly faster). If you want you can thus force  
        ! automatic conversion to single precision at compile time or  change 
        ! all the declarations and constants in the code from  real(kind=dp) to
        ! single.
        !
        ! INPUT
        !   im (INTEGER)
        !   nx, ny (INTEGER)
        !   c11, c12, c22, c66, rho (REAL)
        !   dx, dy (REAL)
        !   npoints_pml (INTEGER) - the thickness of the pml
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
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

            ! source
            integer,dimension(:) :: src
            integer :: isource, jsource
        
            ! velocity threshold above which we consider that the code became unstable
            real(kind=dp), parameter :: STABILITY_THRESHOLD = 1.d+25
        
            ! main arrays
            real(kind=dp), dimension(nx,nz) :: vx,vz,sigmaxx,sigmazz,sigmaxz
                   
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
        
            real(kind=dp), dimension(nx,nz) :: gamma_x, gamma_z
            ! for the source
            real(kind=dp),dimension(nstep) :: srcx, srcz
        
            integer :: i,j,it
        
            real(kind=dp) :: velocnorm

            ! Boolean flag to save as double precision or single precision 
            logical :: SINGLE
            logical, intent(in), optional :: OUTPUT_SINGLE 

            ! Name the f2py inputs 
            !f2py3 intent(in) :: nx, nx, dx, dx,
            !f2py3 intent(in) :: noints_pml, src, nstep, OUTPUT_SINGLE

            ! The default data output is single precision unless OUTPUT_SINGLE is 
            ! set to .FALSE.
            if (present(OUTPUT_SINGLE)) then 
                SINGLE = OUTPUT_SINGLE 
            else
                SINGLE = .TRUE.
            endif
        
            ! -------------------- Load Stiffness Coefficients --------------------
        
            call material_rw('c11.dat', c11, .TRUE.)
            call material_rw('c13.dat', c13, .TRUE.)
            call material_rw('c15.dat', c15, .TRUE.)
            call material_rw('c33.dat', c33, .TRUE.)
            call material_rw('c35.dat', c35, .TRUE.)
            call material_rw('c55.dat', c55, .TRUE.)
            call material_rw('rho.dat', rho, .TRUE.)
            
            ! ------------------- Load Attenuation Coefficients --------------------
            call material_rw('gammas_x.dat', gamma_x, .TRUE.)
            call material_rw('gammas_z.dat', gamma_x, .TRUE.)
            
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
        
                vx(i,j) = vx(i,j)*(1 - gamma_x(i,j) ) + (value_dsigmaxx_dx + value_dsigmaxz_dz) * DT / rho(i,j)
        
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
        
                vz(i,j) = vz(i,j)*(1 - gamma_z(i,j) ) + (value_dsigmaxz_dx + value_dsigmazz_dz) * DT / deltarho
        
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
        
            call write_image2(vx, nx, nz, src, it, 'Vx', SINGLE)
            call write_image2(vz, nx, nz, src, it, 'Vz', SINGLE)
        
            enddo   ! end of time loop
    end subroutine seismic2
    
    ! =========================================================================
    subroutine seismic25(nx, ny, nz, dx, dy, dz, npoints_pml, src, nstep, OUTPUT_SINGLE)

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
        real(kind=dp) :: deltarho

        ! total number of time steps
        integer :: nstep

        ! time step in seconds. decreasing the time step improves the pml attenuation
        ! but it should be inversely proportional to the center frequency of the 
        ! source frequency 
        real(kind=dp) :: DT
        real(kind=dp) :: dx, dy, dz

        ! source
        integer,dimension(:) :: src
        integer :: isource, jsource, ksource

        ! velocity threshold above which we consider that the code became unstable
        real(kind=dp), parameter :: STABILITY_THRESHOLD = 1.d+25

        ! main arrays
        real(kind=dp), dimension(nx,ny,nz) :: vx,vy,vz, &
                                sigmaxx,sigmayy,sigmazz, &
                                sigmaxy,sigmaxz,sigmayz

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
        real(kind=dp), dimension(nx) :: K_x,alpha_x,a_x,b_x, &
                            K_x_half,alpha_x_half,a_x_half,b_x_half
        real(kind=dp), dimension(ny) :: K_y,alpha_y,a_y,b_y, &
                            K_y_half,alpha_y_half,a_y_half,b_y_half
        real(kind=dp), dimension(nz) :: K_z,alpha_z,a_z,b_z, &
                            K_z_half,alpha_z_half,a_z_half,b_z_half

        ! for the source
        real(kind=dp),dimension(nstep) :: srcx, srcy, srcz

        integer :: i,j,k,it

        real(kind=dp) :: velocnorm
        
        ! Boolean flag to save as double precision or single precision 
        logical :: SINGLE
        logical, intent(in), optional :: OUTPUT_SINGLE 
        
        ! Name the f2py inputs 
        !f2py3 intent(in) :: nx, ny, nz, dx, dy, dz,
        !f2py3 intent(in) :: noints_pml, src, nstep, OUTPUT_SINGLE

        ! The default data output is single precision unless OUTPUT_SINGLE is 
        ! set to .FALSE.
        if (present(OUTPUT_SINGLE)) then 
            SINGLE = OUTPUT_SINGLE 
        else
            SINGLE = .TRUE.
        endif
        
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

        ! ================================ LOAD SOURCE ================================

        call loadsource('seismicsourcex.dat', nstep, srcx)
        call loadsource('seismicsourcey.dat', nstep, srcy)
        call loadsource('seismicsourcez.dat', nstep, srcz)

        ! ==================================== PML ====================================
        ! Initialize PML 
        K_x(:) = 1.d0
        K_x_half(:) = 1.d0
        alpha_x(:) = 0.d0
        alpha_x_half(:) = 0.d0
        a_x(:) = 0.d0
        a_x_half(:) = 0.d0

        K_y(:) = 1.d0
        K_y_half(:) = 1.d0
        alpha_y(:) = 0.d0
        alpha_y_half(:) = 0.d0
        a_y(:) = 0.d0
        a_y_half(:) = 0.d0

        K_z(:) = 1.d0
        K_z_half(:) = 1.d0 
        alpha_z(:) = 0.d0
        alpha_z_half(:) = 0.d0
        a_z(:) = 0.d0
        a_z_half(:) = 0.d0

        ! ------------------------- Boundary Conditions -------------------------
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
            call write_image3(vx, nx, ny, nz, src, it, 'Vx', SINGLE)
            call write_image3(vy, nx, ny, nz, src, it, 'Vy', SINGLE)
            call write_image3(vz, nx, ny, nz, src, it, 'Vz', SINGLE)
            ! Now write the stress Values
            ! call write_image3(sigmaxx, nx, ny, nz, it, 'S1')
            ! call write_image3(sigmayy, nx, ny, nz, it, 'S2')
            ! call write_image3(sigmazz, nx, ny, nz, it, 'S3')
            ! call write_image3(sigmaxy, nx, ny, nz, it, 'S6')
            ! call write_image3(sigmayz, nx, ny, nz, it, 'S4')
            ! call write_image3(sigmaxz, nx, ny, nz, it, 'S5')

        enddo   ! end of time loop

    end subroutine seismic25

end module seismicfdtd
