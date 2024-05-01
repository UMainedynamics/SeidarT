module readwrite_routines
    implicit none
    ! integer, parameter :: dp = kind(0.d0)
    ! integer, parameter :: sp = kind(0e0)

    contains
    
    !==========================================================================
    subroutine loadsource(filename, N, srcfn)
        
        implicit none

        integer,parameter :: dp = kind(0.d0)
        character(len=*) :: filename
        integer :: N
        real(kind=dp),dimension(N) :: srcfn
        
        open(unit = 13, form="unformatted", file = trim(filename))
        read(13) srcfn
        
        close(unit = 13)

    end subroutine loadsource

    !==========================================================================
    subroutine loadcpml(filename, image_data)

        implicit none

        integer,parameter :: dp = kind(0.d0)
        character(len=*) :: filename
        real(kind=dp),dimension(:) :: image_data

        open(unit = 13, form="unformatted", file = trim(filename), access='stream')
        read(13) image_data
        close(unit = 13)
    end subroutine loadcpml
    
    !==========================================================================
    subroutine permittivity_write(im, mlist, npoints_pml, nx, nz) 
        ! STIFFNESS_ARRAYS takes a matrix containing the material integer 
        ! identifiers and creates the same size array for each independent 
        ! coefficient of the stiffness matrix along with a density matrix. 
        ! Since we ae using PML boundaries, we will extend the the boundary 
        ! values through the PML region.
        ! 
        ! INPUT 
        !   im (INTEGER)  
        !   mlist (REAL)
        !   eps11(i,j), sig11(i,j), eps22(i,j), sig22, (REAL) -
        !   npoints_pml (INTEGER) - the 
        !   
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        implicit none 
        
        integer :: nx,nz
        integer,parameter :: dp = kind(0.d0)
        integer,dimension(nx,nz) :: im
        integer :: i, j, npoints_pml
        real(kind=dp), dimension(:,:) :: mlist
        real(kind=dp), dimension(2*npoints_pml+nx,2*npoints_pml+nz) :: &
                eps11, eps22, eps33, &
                eps12, eps13, eps23, &
                sig11, sig22, sig33, &
                sig12, sig13, sig23
        
        !f2py3 intent(in):: im, mlist, npoints_pml, nx, nz
        
        ! Allocate space for permittivity and conductivity values
        eps11(:,:) = 0.d0
        eps12(:,:) = 0.d0
        eps13(:,:) = 0.d0
        eps22(:,:) = 0.d0
        eps23(:,:) = 0.d0
        eps33(:,:) = 0.d0
        sig11(:,:) = 0.d0
        sig12(:,:) = 0.d0
        sig13(:,:) = 0.d0
        sig22(:,:) = 0.d0
        sig23(:,:) = 0.d0
        sig33(:,:) = 0.d0
        
        do i=npoints_pml+1,nx + npoints_pml
            do j=npoints_pml+1,nz + npoints_pml
                eps11(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 2)
                eps12(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),3)
                eps13(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),4)
                eps22(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 5)
                eps23(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),6)
                eps33(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 7)
                
                sig11(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 8) 
                sig12(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),9)
                sig13(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),10)
                sig22(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 11)
                sig23(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),12)
                sig33(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 13)
            end do
        end do
        
        ! Extend the boundary values of the stiffnesses into the PML region
        do i = 1,npoints_pml+1
            ! top and bottom
            eps11( i, : ) = eps11(npoints_pml+1,:)
            eps22( i, : ) = eps22(npoints_pml+1,:)
            eps33( i, : ) = eps33(npoints_pml+1,:)
            eps12( i, : ) = eps12(npoints_pml+1,:)
            eps13( i, : ) = eps13(npoints_pml+1,:)
            eps23( i, : ) = eps23(npoints_pml+1,:)

            eps11( nx+npoints_pml-1+i, : ) = eps11(nx+npoints_pml-1,:)
            eps22( nx+npoints_pml-1+i, : ) = eps22(nx+npoints_pml-1,:)
            eps33( nx+npoints_pml-1+i, : ) = eps33(nx+npoints_pml-1,:)
            eps12( nx+npoints_pml-1+i, : ) = eps12(nx+npoints_pml-1,:)
            eps13( nx+npoints_pml-1+i, : ) = eps13(nx+npoints_pml-1,:)
            eps23( nx+npoints_pml-1+i, : ) = eps23(nx+npoints_pml-1,:)
            
            sig11( i, : ) = sig11(npoints_pml+1,:)
            sig22( i, : ) = sig22(npoints_pml+1,:)
            sig33( i, : ) = sig33(npoints_pml+1,:)
            sig12( i, : ) = sig12(npoints_pml+1,:)
            sig13( i, : ) = sig13(npoints_pml+1,:)
            sig23( i, : ) = sig23(npoints_pml+1,:)

            sig11( nx+npoints_pml-1+i, : ) = sig11(nx+npoints_pml-1,:)
            sig22( nx+npoints_pml-1+i, : ) = sig22(nx+npoints_pml-1,:)
            sig33( nx+npoints_pml-1+i, : ) = sig33(nx+npoints_pml-1,:)
            sig13( nx+npoints_pml-1+i, : ) = sig12(nx+npoints_pml-1,:)
            sig13( nx+npoints_pml-1+i, : ) = sig13(nx+npoints_pml-1,:)
            sig23( nx+npoints_pml-1+i, : ) = sig23(nx+npoints_pml-1,:)
            
            !!!!!  ! left and right
            eps11( :, i ) = eps11(:, npoints_pml+1)
            eps22( :, i ) = eps22(:, npoints_pml+1)
            eps33( :, i ) = eps33(:, npoints_pml+1)
            eps12( :, i ) = eps12(:, npoints_pml+1)
            eps13( :, i ) = eps13(:, npoints_pml+1)
            eps23( :, i ) = eps23(:, npoints_pml+1)

            eps11( :, nz+npoints_pml-1+i ) = eps11(:,nz+npoints_pml-1)    
            eps22( :, nz+npoints_pml-1+i ) = eps22(:,nz+npoints_pml-1)
            eps33( :, nz+npoints_pml-1+i ) = eps33(:,nz+npoints_pml-1)
            eps12( :, nz+npoints_pml-1+i ) = eps12(:,nz+npoints_pml-1)    
            eps13( :, nz+npoints_pml-1+i ) = eps13(:,nz+npoints_pml-1)
            eps23( :, nz+npoints_pml-1+i ) = eps23(:,nz+npoints_pml-1)
            
            sig11( :, i ) = sig11(:, npoints_pml+1)
            sig22( :, i ) = sig22(:, npoints_pml+1)
            sig33( :, i ) = sig33(:, npoints_pml+1)
            sig12( :, i ) = sig11(:, npoints_pml+1)
            sig13( :, i ) = sig13(:, npoints_pml+1)
            sig23( :, i ) = sig33(:, npoints_pml+1)
            
            sig11( :, nz+npoints_pml-1+i ) = sig11(:,nz+npoints_pml-1)    
            sig22( :, nz+npoints_pml-1+i ) = sig22(:,nz+npoints_pml-1)
            sig33( :, nz+npoints_pml-1+i ) = sig33(:,nz+npoints_pml-1)
            sig12( :, nz+npoints_pml-1+i ) = sig12(:,nz+npoints_pml-1)    
            sig13( :, nz+npoints_pml-1+i ) = sig13(:,nz+npoints_pml-1)
            sig23( :, nz+npoints_pml-1+i ) = sig23(:,nz+npoints_pml-1)
        end do 

        ! Write each of the matrices to file
        call material_rw('eps11.dat', eps11, .FALSE.)
        call material_rw('eps12.dat', eps12, .FALSE.)
        call material_rw('eps13.dat', eps13, .FALSE.)
        call material_rw('eps22.dat', eps22, .FALSE.)
        call material_rw('eps23.dat', eps23, .FALSE.)
        call material_rw('eps33.dat', eps33, .FALSE.)
        call material_rw('sig11.dat', sig11, .FALSE.)
        call material_rw('sig12.dat', sig12, .FALSE.)
        call material_rw('sig13.dat', sig13, .FALSE.)
        call material_rw('sig22.dat', sig22, .FALSE.)
        call material_rw('sig23.dat', sig23, .FALSE.)
        call material_rw('sig33.dat', sig33, .FALSE.)

    end subroutine permittivity_write
    
    !==========================================================================
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
    
    !==========================================================================
    subroutine attenuation_write(im, alist, npoints_pml, nx, nz, cpmlvalue, seismic) 
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
        real(kind=dp),dimension(:,:) :: alist
        real(kind=dp) :: cpmlvalue
        logical :: seismic
        real(kind=dp),dimension(2*npoints_pml+nx,2*npoints_pml+nz) :: gamma_x, gamma_y, gamma_z
                                                                     
        

        !f2py3 intent(in) :: im, alist, npoints_pml, nx, nz, gradient
        
        ! call material_rw('rho.dat', rho, .TRUE. )
        gamma_x(:,:) = cpmlvalue 
        gamma_y(:,:) = cpmlvalue 
        gamma_z(:,:) = cpmlvalue 
        
        !Assign between the PML regions
        do i = npoints_pml+1, nx+npoints_pml
            do j = npoints_pml+1, nz+npoints_pml
                gamma_x(i,j) = alist( im(i-npoints_pml,j-npoints_pml), 2)!*dt/rho(i,j)
                gamma_y(i,j) = alist( im(i-npoints_pml,j-npoints_pml), 3)!*dt/rho(i,j)
                gamma_z(i,j) = alist( im(i-npoints_pml,j-npoints_pml), 4)!*dt/rho(i,j)
            enddo
        enddo
        
        ! Write each of the matrices to file
        if ( seismic ) then
            call material_rw('gammas_x.dat', gamma_x, .FALSE.)
            call material_rw('gammas_y.dat', gamma_y, .FALSE.)
            call material_rw('gammas_z.dat', gamma_z, .FALSE.)
        else
            call material_rw('gammae_x.dat', gamma_x, .FALSE.)
            call material_rw('gammae_y.dat', gamma_y, .FALSE.)
            call material_rw('gammae_z.dat', gamma_z, .FALSE.)
        end if    

    end subroutine attenuation_write
    
    !==========================================================================
    subroutine material_rw(filename, image_data, readfile)

        implicit none
        
        integer,parameter :: dp = kind(0.d0)
        character(len=*) :: filename
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

    ! =========================================================================
    subroutine write_image2(image_data, nx, nz, src, it, channel, SINGLE)
    
        implicit none

        integer, parameter :: dp = kind(0.d0)
        integer :: nx, nz, it
        integer,dimension(2) :: src
        real(kind=dp) :: image_data(nx, nz)
        character(len=2) :: channel
        character(len=100) :: filename
        logical :: SINGLE

        ! WRITE (filename, "(a2, i6.6, '.dat')" ) channel, it
        WRITE (filename, "(a2, a1, i6.6, a1, i0, a1, i0, a1, a4)" ) &
                    channel,'.', it,'.', src(1),'.', src(2), '.','.dat'
        open(unit = 10, form = 'unformatted', file = trim(filename) )

        if (SINGLE) then
            write(10) sngl(image_data)
        else
            write(10) image_data 
        end if 


        close(unit = 10)

    end subroutine write_image2
    
    ! =========================================================================
    subroutine write_image3(image_data, nx, ny, nz, src, it, channel, SINGLE)
    
        implicit none
    
        integer, parameter :: dp = kind(0.d0)
        integer :: nx, ny, nz, it
        integer,dimension(3) :: src
        real(kind=dp) :: image_data(nx, ny, nz)
        character(len=2) :: channel
        character(len=80) :: filename
        logical :: SINGLE
        
        WRITE (filename, "(a2, a1, i6.6, a1, i0, a1, i0, a1, i0, a4)" ) &
                        channel,'.',it,'.', src(1),'.',src(2),'.',src(3),'.dat'
        
        open(unit = 10, form = 'unformatted', file = trim(filename) )
        
        if (SINGLE) then
            write(10) sngl(image_data)
        else
            write(10) image_data 
        end if 
        
        close(unit = 10)

    end subroutine write_image3

    ! ---------------------------------------------------------------------
    ! subroutine write_image3c(image_data, nx, ny, nz, src, it, channel, SINGLE)
    
    !     implicit none
    
    !     integer, parameter :: dp = kind(0.d0)
    !     integer, parameter :: sp = kind(1e0)
    !     integer :: nx, ny, nz, it
    !     integer,dimension(3) :: src
    !     complex(kind=dp) :: image_data(nx, ny, nz)
    !     real(kind=sp) :: real_part(nx,ny,nz), imag_part(nx,ny,nz)
    !     character(len=2) :: channel
    !     character(len=80) :: filename
    !     logical :: SINGLE
        
    !     WRITE (filename, "(a2, a1, i6.6, a1, i0, a1, i0, a1, i0, a4)" ) &
    !                     channel,'.',it,'.', src(1),'.',src(2),'.',src(3),'.dat'
        
    !     open(unit = 10, form = 'unformatted', file = trim(filename) )
        
    !     if (SINGLE) then
    !         real_part = real(image_data, kind=sp)
    !         imag_part = aimag(image_data)
    !         write(10) real_part, imag_part
    !     else
    !         write(10) image_data 
    !     end if 
        
    !     close(unit = 10)

    ! end subroutine write_image3c
    
    !==========================================================================
    subroutine permittivity_write_c(im, mlist, npoints_pml, nx, nz) 
        ! STIFFNESS_ARRAYS takes a matrix containing the material integer 
        ! identifiers and creates the same size array for each independent 
        ! coefficient of the stiffness matrix along with a density matrix. 
        ! Since we ae using PML boundaries, we will extend the the boundary 
        ! values through the PML region.
        ! 
        ! INPUT 
        !   im (INTEGER)  
        !   mlist (REAL)
        !   eps11(i,j), sig11(i,j), eps22(i,j), sig22, (REAL) -
        !   npoints_pml (INTEGER) - the 
        !   
        ! sigma is related to the complex permittivity value, but it is not a 
        ! complex valued number. The variable mlist is a complex valued array 
        ! but for sigma values, the complex component is 0. 
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        implicit none 
        
        integer :: nx,nz
        integer,parameter :: dp = kind(0.d0)
        integer,dimension(nx,nz) :: im
        integer :: i, j, npoints_pml
        complex(kind=dp), dimension(:,:) :: mlist
        complex(kind=dp), dimension(2*npoints_pml+nx,2*npoints_pml+nz) :: &
                eps11, eps22, eps33, &
                eps12, eps13, eps23
        real(kind=dp), dimension(2*npoints_pml+nx,2*npoints_pml+nz) :: &
                sig11, sig22, sig33, &
                sig12, sig13, sig23
        
        !f2py3 intent(in):: im, mlist, npoints_pml, nx, nz
        
        ! Allocate space for permittivity and conductivity values
        eps11(:,:) = 0.d0
        eps12(:,:) = 0.d0
        eps13(:,:) = 0.d0
        eps22(:,:) = 0.d0
        eps23(:,:) = 0.d0
        eps33(:,:) = 0.d0
        sig11(:,:) = 0.d0
        sig12(:,:) = 0.d0
        sig13(:,:) = 0.d0
        sig22(:,:) = 0.d0
        sig23(:,:) = 0.d0
        sig33(:,:) = 0.d0
        
        do i=npoints_pml+1,nx + npoints_pml
            do j=npoints_pml+1,nz + npoints_pml
                eps11(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 2)
                eps12(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),3)
                eps13(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),4)
                eps22(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 5)
                eps23(i,j) = mlist( im(i-npoints_pml, j-npoints_pml),6)
                eps33(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 7)
                
                sig11(i,j) = abs(mlist( im(i-npoints_pml,j-npoints_pml), 8) )
                sig12(i,j) = abs(mlist( im(i-npoints_pml, j-npoints_pml),9) )
                sig13(i,j) = abs(mlist( im(i-npoints_pml, j-npoints_pml),10) )
                sig22(i,j) = abs(mlist( im(i-npoints_pml,j-npoints_pml), 11) )
                sig23(i,j) = abs(mlist( im(i-npoints_pml, j-npoints_pml),12) )
                sig33(i,j) = abs(mlist( im(i-npoints_pml,j-npoints_pml), 13) )
            end do
        end do
        
        ! Extend the boundary values of the stiffnesses into the PML region
        do i = 1,npoints_pml+1
            ! top and bottom
            eps11( i, : ) = eps11(npoints_pml+1,:)
            eps22( i, : ) = eps22(npoints_pml+1,:)
            eps33( i, : ) = eps33(npoints_pml+1,:)
            eps12( i, : ) = eps12(npoints_pml+1,:)
            eps13( i, : ) = eps13(npoints_pml+1,:)
            eps23( i, : ) = eps23(npoints_pml+1,:)

            eps11( nx+npoints_pml-1+i, : ) = eps11(nx+npoints_pml-1,:)
            eps22( nx+npoints_pml-1+i, : ) = eps22(nx+npoints_pml-1,:)
            eps33( nx+npoints_pml-1+i, : ) = eps33(nx+npoints_pml-1,:)
            eps12( nx+npoints_pml-1+i, : ) = eps12(nx+npoints_pml-1,:)
            eps13( nx+npoints_pml-1+i, : ) = eps13(nx+npoints_pml-1,:)
            eps23( nx+npoints_pml-1+i, : ) = eps23(nx+npoints_pml-1,:)
            
            sig11( i, : ) = sig11(npoints_pml+1,:)
            sig22( i, : ) = sig22(npoints_pml+1,:)
            sig33( i, : ) = sig33(npoints_pml+1,:)
            sig12( i, : ) = sig12(npoints_pml+1,:)
            sig13( i, : ) = sig13(npoints_pml+1,:)
            sig23( i, : ) = sig23(npoints_pml+1,:)

            sig11( nx+npoints_pml-1+i, : ) = sig11(nx+npoints_pml-1,:)
            sig22( nx+npoints_pml-1+i, : ) = sig22(nx+npoints_pml-1,:)
            sig33( nx+npoints_pml-1+i, : ) = sig33(nx+npoints_pml-1,:)
            sig13( nx+npoints_pml-1+i, : ) = sig12(nx+npoints_pml-1,:)
            sig13( nx+npoints_pml-1+i, : ) = sig13(nx+npoints_pml-1,:)
            sig23( nx+npoints_pml-1+i, : ) = sig23(nx+npoints_pml-1,:)
            
            !!!!!  ! left and right
            eps11( :, i ) = eps11(:, npoints_pml+1)
            eps22( :, i ) = eps22(:, npoints_pml+1)
            eps33( :, i ) = eps33(:, npoints_pml+1)
            eps12( :, i ) = eps12(:, npoints_pml+1)
            eps13( :, i ) = eps13(:, npoints_pml+1)
            eps23( :, i ) = eps23(:, npoints_pml+1)

            eps11( :, nz+npoints_pml-1+i ) = eps11(:,nz+npoints_pml-1)    
            eps22( :, nz+npoints_pml-1+i ) = eps22(:,nz+npoints_pml-1)
            eps33( :, nz+npoints_pml-1+i ) = eps33(:,nz+npoints_pml-1)
            eps12( :, nz+npoints_pml-1+i ) = eps12(:,nz+npoints_pml-1)    
            eps13( :, nz+npoints_pml-1+i ) = eps13(:,nz+npoints_pml-1)
            eps23( :, nz+npoints_pml-1+i ) = eps23(:,nz+npoints_pml-1)
            
            sig11( :, i ) = sig11(:, npoints_pml+1)
            sig22( :, i ) = sig22(:, npoints_pml+1)
            sig33( :, i ) = sig33(:, npoints_pml+1)
            sig12( :, i ) = sig11(:, npoints_pml+1)
            sig13( :, i ) = sig13(:, npoints_pml+1)
            sig23( :, i ) = sig33(:, npoints_pml+1)
            
            sig11( :, nz+npoints_pml-1+i ) = sig11(:,nz+npoints_pml-1)    
            sig22( :, nz+npoints_pml-1+i ) = sig22(:,nz+npoints_pml-1)
            sig33( :, nz+npoints_pml-1+i ) = sig33(:,nz+npoints_pml-1)
            sig12( :, nz+npoints_pml-1+i ) = sig12(:,nz+npoints_pml-1)    
            sig13( :, nz+npoints_pml-1+i ) = sig13(:,nz+npoints_pml-1)
            sig23( :, nz+npoints_pml-1+i ) = sig23(:,nz+npoints_pml-1)
        end do 

        ! Write each of the matrices to file
        call material_rwc('eps11.dat', eps11, .FALSE.)
        call material_rwc('eps12.dat', eps12, .FALSE.)
        call material_rwc('eps13.dat', eps13, .FALSE.)
        call material_rwc('eps22.dat', eps22, .FALSE.)
        call material_rwc('eps23.dat', eps23, .FALSE.)
        call material_rwc('eps33.dat', eps33, .FALSE.)
        call material_rw('sig11.dat', sig11, .FALSE.)
        call material_rw('sig12.dat', sig12, .FALSE.)
        call material_rw('sig13.dat', sig13, .FALSE.)
        call material_rw('sig22.dat', sig22, .FALSE.)
        call material_rw('sig23.dat', sig23, .FALSE.)
        call material_rw('sig33.dat', sig33, .FALSE.)

    end subroutine permittivity_write_c
    
    !==========================================================================
    ! subroutine loadcpml_c(filename, image_data)

    !     implicit none

    !     integer,parameter :: dp = kind(0.d0)
    !     character(len=*) :: filename
    !     complex(kind=dp),dimension(:) :: image_data

    !     open(unit = 13, form="unformatted", file = trim(filename), access='stream')
    !     read(13) image_data
    !     close(unit = 13)
    ! end subroutine loadcpml_c
    
    !==========================================================================
    ! subroutine stiffness_write_c(im, mlist, npoints_pml, nx, nz, gradient) 
    !     ! STIFFNESS_ARRAYS takes a matrix containing the material integer identifiers 
    !     ! and creates the same size array for each independent coefficient of the 
    !     ! stiffness matrix along with a density matrix. Since we ae using PML
    !     ! boundaries, we will extend the the boundary values through the PML region.
    !     ! 
    !     ! INPUT 
    !     !   im (INTEGER)  
    !     !   mlist (COMPLEX)
    !     !   c11(i,j), c12(i,j), c22(i,j), c66 (COMPLEX)
    !     !   rho(i,j) (REAL)
    !     !   npoints_pml (INTEGER) - the 
    !     !   
    !     ! rho is not complex valued, but it is contained in the variable 'mlist' as a
    !     ! complex number; r' + j r'' where r'' = 0
    !     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !     implicit none 

    !     integer,parameter :: dp = kind(0.d0)
    !     integer :: nx, nz
    !     integer,dimension(nx,nz) :: im
    !     integer :: i, j, npoints_pml
    !     complex(kind=dp),dimension(:,:) :: mlist
    !     real(kind=dp),dimension(nx,nz) :: gradient
    !     complex(kind=dp),dimension(2*npoints_pml+nx,2*npoints_pml+nz) :: c11,c12,c13,&
    !                                                                 c14,c15,c16, &
    !                                                                 c22,c23,c24,c25,c26,&
    !                                                                 c33,c34,c35,c36, &
    !                                                                 c44,c45,c46, &
    !                                                                 c55,c56,c66
    !     real(kind=dp),dimension(2*npoints_pml+nx, 2*npoints_pml*nz) :: rho

    !     !f2py3 intent(in) :: im, mlist, npoints_pml, nx, nz, gradient

    !     c11(:,:) = 0.d0 
    !     c12(:,:) = 0.d0 
    !     c13(:,:) = 0.d0
    !     c14(:,:) = 0.d0 
    !     c15(:,:) = 0.d0 
    !     c16(:,:) = 0.d0 
    !     c22(:,:) = 0.d0 
    !     c23(:,:) = 0.d0 
    !     c24(:,:) = 0.d0 
    !     c25(:,:) = 0.d0 
    !     c26(:,:) = 0.d0 
    !     c33(:,:) = 0.d0 
    !     c34(:,:) = 0.d0 
    !     c35(:,:) = 0.d0 
    !     c36(:,:) = 0.d0 
    !     c44(:,:) = 0.d0 
    !     c45(:,:) = 0.d0 
    !     c46(:,:) = 0.d0 
    !     c55(:,:) = 0.d0 
    !     c56(:,:) = 0.d0 
    !     c66(:,:) = 0.d0 
    !     rho(:,:) = 0.d0 

    !     !Assign between the PML regions
    !     do i = npoints_pml+1, nx+npoints_pml
    !     do j = npoints_pml+1, nz+npoints_pml
    !         c11(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 2)
    !         c12(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 3)
    !         c13(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 4)
    !         c14(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 5)
    !         c15(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 6)
    !         c16(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 7)
    !         c22(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 8)
    !         c23(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 9)
    !         c24(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 10)
    !         c25(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 11)
    !         c26(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 12)
    !         c33(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 13)
    !         c34(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 14)
    !         c35(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 15)
    !         c36(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 16)
    !         c44(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 17)
    !         c45(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 18)
    !         c46(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 19)
    !         c55(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 20)
    !         c56(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 21)
    !         c66(i,j) = mlist( im(i-npoints_pml,j-npoints_pml), 22)
    !         rho(i,j) = abs(mlist( im(i-npoints_pml,j-npoints_pml), 23) )
    !     enddo
    !     enddo

    !     rho(npoints_pml+1:nx+npoints_pml,npoints_pml+1:nz+npoints_pml) = &
    !         rho(npoints_pml+1:nx+npoints_pml,npoints_pml+1:nz+npoints_pml)*gradient

    !     ! Extend the boundary values of the stiffnesses into the PML region
    !     do i = 1,npoints_pml+1
    !     ! top 
    !     c11( i, :) = c11(npoints_pml+1,:)
    !     c12( i, :) = c12(npoints_pml+1,:)
    !     c13( i, :) = c13(npoints_pml+1,:)
    !     c14( i, :) = c14(npoints_pml+1,:)
    !     c15( i, :) = c15(npoints_pml+1,:)
    !     c16( i, :) = c16(npoints_pml+1,:)
    !     c22( i, :) = c22(npoints_pml+1,:)
    !     c23( i, :) = c23(npoints_pml+1,:)
    !     c24( i, :) = c24(npoints_pml+1,:)
    !     c25( i, :) = c25(npoints_pml+1,:)
    !     c26( i, :) = c26(npoints_pml+1,:)
    !     c33( i, :) = c33(npoints_pml+1,:)
    !     c34( i, :) = c34(npoints_pml+1,:)
    !     c35( i, :) = c35(npoints_pml+1,:)
    !     c36( i, :) = c36(npoints_pml+1,:)
    !     c44( i, :) = c44(npoints_pml+1,:)
    !     c45( i, :) = c45(npoints_pml+1,:)
    !     c46( i, :) = c46(npoints_pml+1,:)
    !     c55( i, :) = c55(npoints_pml+1,:)
    !     c56( i, :) = c56(npoints_pml+1,:)
    !     c66( i, :) = c66(npoints_pml+1,:)
    !     rho( i, :) = rho(npoints_pml+1,:)

    !     ! bottom
    !     c11( nx+npoints_pml-1+i, :) = c11(nx+npoints_pml-1,:)
    !     c12( nx+npoints_pml-1+i, :) = c12(nx+npoints_pml-1,:)
    !     c13( nx+npoints_pml-1+i, :) = c13(nx+npoints_pml-1,:)
    !     c14( nx+npoints_pml-1+i, :) = c14(nx+npoints_pml-1,:)
    !     c15( nx+npoints_pml-1+i, :) = c15(nx+npoints_pml-1,:)
    !     c16( nx+npoints_pml-1+i, :) = c16(nx+npoints_pml-1,:)
    !     c22( nx+npoints_pml-1+i, :) = c22(nx+npoints_pml-1,:)
    !     c23( nx+npoints_pml-1+i, :) = c23(nx+npoints_pml-1,:)
    !     c24( nx+npoints_pml-1+i, :) = c24(nx+npoints_pml-1,:)
    !     c25( nx+npoints_pml-1+i, :) = c25(nx+npoints_pml-1,:)
    !     c26( nx+npoints_pml-1+i, :) = c26(nx+npoints_pml-1,:)
    !     c33( nx+npoints_pml-1+i, :) = c33(nx+npoints_pml-1,:)
    !     c34( nx+npoints_pml-1+i, :) = c34(nx+npoints_pml-1,:)
    !     c35( nx+npoints_pml-1+i, :) = c35(nx+npoints_pml-1,:)
    !     c36( nx+npoints_pml-1+i, :) = c36(nx+npoints_pml-1,:)
    !     c44( nx+npoints_pml-1+i, :) = c44(nx+npoints_pml-1,:)
    !     c45( nx+npoints_pml-1+i, :) = c45(nx+npoints_pml-1,:)
    !     c46( nx+npoints_pml-1+i, :) = c46(nx+npoints_pml-1,:)
    !     c55( nx+npoints_pml-1+i, :) = c55(nx+npoints_pml-1,:)
    !     c56( nx+npoints_pml-1+i, :) = c56(nx+npoints_pml-1,:)
    !     c66( nx+npoints_pml-1+i, :) = c66(nx+npoints_pml-1,:)
    !     rho( nx+npoints_pml-1+i, :) = rho(nx+npoints_pml-1,:)

    !     ! left 
    !     c11( :, i) = c11(:, npoints_pml+1)
    !     c12( :, i) = c12(:, npoints_pml+1)
    !     c13( :, i) = c13(:, npoints_pml+1)
    !     c14( :, i) = c14(:, npoints_pml+1)
    !     c15( :, i) = c15(:, npoints_pml+1)
    !     c16( :, i) = c16(:, npoints_pml+1)
    !     c22( :, i) = c22(:, npoints_pml+1)
    !     c23( :, i) = c23(:, npoints_pml+1)
    !     c24( :, i) = c24(:, npoints_pml+1)
    !     c25( :, i) = c25(:, npoints_pml+1)
    !     c26( :, i) = c26(:, npoints_pml+1)
    !     c33( :, i) = c33(:, npoints_pml+1)
    !     c34( :, i) = c34(:, npoints_pml+1)
    !     c35( :, i) = c35(:, npoints_pml+1)
    !     c36( :, i) = c36(:, npoints_pml+1)
    !     c44( :, i) = c44(:, npoints_pml+1)
    !     c45( :, i) = c45(:, npoints_pml+1)
    !     c46( :, i) = c46(:, npoints_pml+1)
    !     c55( :, i) = c55(:, npoints_pml+1)
    !     c56( :, i) = c56(:, npoints_pml+1)
    !     c66( :, i) = c66(:, npoints_pml+1)
    !     rho( :, i) = rho(:, npoints_pml+1)

    !     ! right
    !     c11( :, nz+npoints_pml-1+i) = c11(:,nz+npoints_pml-1)
    !     c12( :, nz+npoints_pml-1+i) = c12(:,nz+npoints_pml-1)
    !     c13( :, nz+npoints_pml-1+i) = c13(:,nz+npoints_pml-1)      
    !     c14( :, nz+npoints_pml-1+i) = c14(:,nz+npoints_pml-1)      
    !     c15( :, nz+npoints_pml-1+i) = c15(:,nz+npoints_pml-1)      
    !     c16( :, nz+npoints_pml-1+i) = c16(:,nz+npoints_pml-1)      
    !     c22( :, nz+npoints_pml-1+i) = c22(:,nz+npoints_pml-1)
    !     c23( :, nz+npoints_pml-1+i) = c23(:,nz+npoints_pml-1)
    !     c24( :, nz+npoints_pml-1+i) = c24(:,nz+npoints_pml-1)
    !     c25( :, nz+npoints_pml-1+i) = c25(:,nz+npoints_pml-1)
    !     c26( :, nz+npoints_pml-1+i) = c26(:,nz+npoints_pml-1)
    !     c33( :, nz+npoints_pml-1+i) = c33(:,nz+npoints_pml-1)
    !     c34( :, nz+npoints_pml-1+i) = c34(:,nz+npoints_pml-1)
    !     c35( :, nz+npoints_pml-1+i) = c35(:,nz+npoints_pml-1)
    !     c36( :, nz+npoints_pml-1+i) = c36(:,nz+npoints_pml-1)
    !     c44( :, nz+npoints_pml-1+i) = c44(:,nz+npoints_pml-1)
    !     c45( :, nz+npoints_pml-1+i) = c45(:,nz+npoints_pml-1)
    !     c46( :, nz+npoints_pml-1+i) = c46(:,nz+npoints_pml-1)
    !     c55( :, nz+npoints_pml-1+i) = c55(:,nz+npoints_pml-1)      
    !     c56( :, nz+npoints_pml-1+i) = c56(:,nz+npoints_pml-1)      
    !     c66( :, nz+npoints_pml-1+i) = c66(:,nz+npoints_pml-1)
    !     rho( :, nz+npoints_pml-1+i) = rho(:,nz+npoints_pml-1)

    !     end do 

    !     ! Write each of the matrices to file
    !     call material_rwc('c11.dat', c11, .FALSE.)
    !     call material_rwc('c12.dat', c12, .FALSE.)
    !     call material_rwc('c13.dat', c13, .FALSE.)
    !     call material_rwc('c14.dat', c14, .FALSE.)
    !     call material_rwc('c15.dat', c15, .FALSE.)
    !     call material_rwc('c16.dat', c16, .FALSE.)
    !     call material_rwc('c22.dat', c22, .FALSE.)
    !     call material_rwc('c23.dat', c23, .FALSE.)
    !     call material_rwc('c24.dat', c24, .FALSE.)
    !     call material_rwc('c25.dat', c25, .FALSE.)
    !     call material_rwc('c26.dat', c26, .FALSE.)
    !     call material_rwc('c33.dat', c33, .FALSE.)
    !     call material_rwc('c34.dat', c34, .FALSE.)
    !     call material_rwc('c35.dat', c35, .FALSE.)
    !     call material_rwc('c36.dat', c36, .FALSE.)
    !     call material_rwc('c44.dat', c44, .FALSE.)
    !     call material_rwc('c45.dat', c45, .FALSE.)
    !     call material_rwc('c46.dat', c46, .FALSE.)
    !     call material_rwc('c55.dat', c55, .FALSE.)
    !     call material_rwc('c56.dat', c56, .FALSE.)
    !     call material_rwc('c66.dat', c66, .FALSE.)
    !     call material_rw('rho.dat', rho, .FALSE. )

    ! end subroutine stiffness_write_c
    
    ! ---------------------------------------------------------------------
    subroutine write_image2c(image_data, nx, nz, src, it, channel, SINGLE)
    
        implicit none

        integer, parameter :: dp = kind(0.d0)
        integer, parameter :: sp = kind(1e0)
        integer :: nx, nz, it
        integer,dimension(2) :: src
        complex(kind=dp) :: image_data(nx, nz)
        real(kind=sp) :: real_part(nx, nz), imag_part(nx, nz)
        character(len=2) :: channel
        character(len=100) :: filename
        logical :: SINGLE

        ! WRITE (filename, "(a2, i6.6, '.dat')" ) channel, it
        WRITE (filename, "(a2, a1, i6.6, a1, i0, a1, i0, a1, a4)" ) &
                    channel,'.', it,'.', src(1),'.', src(2), '.','.dat'
        
        open(unit = 10, form = 'unformatted', file = trim(filename) )
        if (SINGLE) then
            real_part = real(image_data, kind = sp)
            imag_part = aimag(image_data)
            write(10) real_part, imag_part
        else
            write(10) image_data 
        end if 
        

        close(unit = 10)

    end subroutine write_image2c
    
    !==========================================================================
    subroutine material_rwc(filename, image_data, readfile)

        implicit none
        
        integer,parameter :: dp = kind(0.d0)
        character(len=*) :: filename
        complex(kind=dp),dimension(:,:) :: image_data
        logical :: readfile
        
        open(unit = 13, form="unformatted", file = trim(filename))
        
        if ( readfile ) then
            read(13) image_data
        else
            write(13) image_data
        endif
        
        close(unit = 13)

    end subroutine material_rwc
    
end module readwrite_routines
