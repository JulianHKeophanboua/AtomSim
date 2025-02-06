
program Console2 
    use iso_fortran_env
    implicit none
    
    ! Variables
    real(kind=REAL64) :: t_0 = 0 ! Initial Time [s]
    real(kind=REAL64) :: t = 100.0D-12 ! Total Time [s]; in order of picoseconds
    real(kind=REAL64) :: dt = 0.1D-15 ! Timestep [s]; in order of femtoseconds
    
    real(kind=REAL64) :: m = (40.0D-3) / (6.02D23) ! Argon Mass [kg]
    
! Went inside subroutine also (rather than making it an argument)
    real(kind=REAL64) :: LJ_E = 0.0103 / (6.242D18) ! [J]
    real(kind=REAL64) :: LJ_S = 3.405D-10 ! Sigma [m]; in order of 10^-10 (Angstroms)
  !----------------------------------------------------------------------------------  
    
    real(kind=REAL64) :: r ! Lattice Parameter [m]; in order of 10^-10 (Angstroms)
    real(kind=REAL64), dimension(50) :: x_0 ! Initial X Position [m]
    real(kind=REAL64), dimension(50) :: y_0 ! Initial Y Position [m]
    real(kind=REAL64), dimension(50) :: vx_0 ! Initial X Velocity [m/s]
    real(kind=REAL64), dimension(50) :: vy_0 ! Initial Y Velocity [m/s]
    integer :: i
    integer :: j = 0
    integer :: k = 0
    
    real(kind=REAL64), dimension(50) :: Fx_old
    real(kind=REAL64), dimension(50) :: Fy_old
    real(kind=REAL64), dimension(50) :: Fx_new
    real(kind=REAL64), dimension(50) :: Fy_new
    
    real(kind=REAL64) :: KE = 0
    real(kind=REAL64), dimension(50) :: PE_A
    real(kind=REAL64) :: PE_tot = 0
    
    r = (2.0**(1.0/6.0)) * LJ_S
    
    do i = 1, SIZE(y_0)
        y_0(i) = r*j
        x_0(i) = r*k
        j = j + 1
        
        if (j == 10) then
            k = k + 1
            j = 0
        end if
    end do
    
    
    !do i = 1, SIZE(x_0)
    !    print *, i, x_0(i), y_0(i)
    !end do



    call find_forces(x_0, y_0, Fx_old, Fy_old, PE_A)

    !print *, PE_A
    
    
    open(1, file = "AtomSim_HW2_Energy.csv", status = "new")
    open(2, file = "AtomSim_HW2_Position.csv", status = "new")    
    write(2,'(*(g0,","))') "X0 [A]", "Y0 [A]"    
    
    do i = 1, SIZE(x_0)
        write(2,'(*(g0,","))') (x_0(i)*(1.0D10)), (y_0(i)*(1.0D10))
    end do
    
    write(1,'(*(g0,","))') "Time [ps]", "KE [eV]", "PE [eV]", "Total Energy [eV]"
   
    ! Body of Console1
    do while (t_0 <= t)
        
        ! Find Kinetic Energy (of System)
        do i = 1, SIZE(vx_0)
            KE = KE + SQRT(((0.5 * m * (vx_0(i)**2.0))**2.0) + ((0.5 * m * (vy_0(i)**2.0))**2.0))
        end do
        
        ! Find Potential Energy (of System)
        do i = 1, SIZE(PE_A)
            PE_tot = PE_tot + PE_A(i)
        end do
        
        write(1,'(*(g0,","))') (t_0*(1.0D12)), (KE*(6.242D18)), (PE_tot*(6.242D18)), ((KE*(6.242D18)) + (PE_tot*(6.242D18)))
        
        ! Reset Variables (Each one only depends on the current time and not the last iteration)
        KE = 0
        PE_tot = 0
        PE_A = 0
        Fx_new = 0
        Fy_new = 0
        
        ! Find New Positions
        do i = 1, SIZE(x_0)
            x_0(i) = x_0(i) + (dt*vx_0(i)) + (((dt**2.0)*Fx_old(i)) / (2.0*m))
            y_0(i) = y_0(i) + (dt*vy_0(i)) + (((dt**2.0)*Fy_old(i)) / (2.0*m))
        end do
        
        !Calculate New Forces + PE
        call find_forces(x_0, y_0, Fx_new, Fy_new, PE_A)
        
        ! Find New Velocities
        do i = 1, SIZE(vx_0)
            vx_0(i) = vx_0(i) + ((dt*(Fx_old(i) + Fx_new(i))) / (2.0*m))
            vy_0(i) = vy_0(i) + ((dt*(Fy_old(i) + Fy_new(i))) / (2.0*m))
        end do
        
        ! Copy New to Old
        do i = 1, SIZE(Fx_old)
            Fx_old(i) = Fx_new(i)
            Fy_old(i) = Fy_new(i)
        end do
        
        t_0 = t_0 + dt
    end do

        write(2,'(*(g0,","))') "Xf [A]", "Yf [A]"    
    do i = 1, SIZE(x_0)
        write(2,'(*(g0,","))') (x_0(i)*(1.0D10)), (y_0(i)*(1.0D10))
    end do
    
end program Console2

    
    
subroutine find_forces(x, y, Fx_new, Fy_new, PE_A)
    use iso_fortran_env
    implicit none
    real(kind=REAL64), dimension(50), intent(in) :: x
    real(kind=REAL64), dimension(50), intent(in) :: y
    real(kind=REAL64), dimension(50), intent(out) :: Fx_new
    real(kind=REAL64), dimension(50), intent(out) :: Fy_new
    real(kind=REAL64), dimension(50), intent(out) :: PE_A
    integer :: i
    integer :: j
    real(kind=REAL64) :: r_ij ! [m]
    real(kind=REAL64) :: Fx ! Placeholders for singular calculation
    real(kind=REAL64) :: Fy ! Placeholders for singular calculation
    real(kind=REAL64) :: Fx_sum = 0 ! Placeholders for summation
    real(kind=REAL64) :: Fy_sum = 0 ! Placeholders for summation
    real(kind=REAL64) :: LJ_E = 0.0103 / (6.242D18) ! [J]
    real(kind=REAL64) :: LJ_S = 3.405D-10 ! Sigma [m]; in order of 10^-10 (Angstroms)
    real(kind=REAL64) :: r ! Lattice Parameter [m]
    r = (2.0**(1.0/6.0)) * LJ_S
    

    do i = 1, SIZE(x)
        do j = 1, SIZE(x)
            if (j /= i) then
                
                !print *, "I:", i, "J:", j
                !print *, "XI:", x(i), "XJ:", x(j), "YI:", y(i), "YJ", y(j)
                r_ij = SQRT(((x(i) - x(j))**2.0) + ((y(i) - y(j))**2.0))
                !print *, "R1", r_ij, i, j

                ! Implementation of Potential Range            
                !if (r_ij <=  SQRT((r**2.0)+(r**2.0))) then
                !    !print *, r_ij, SQRT((r**2.0)+(r**2.0)), i, j
                !    !print *, "PE1", PE_A(i), i, j
                !    PE_A(i) = PE_A(i) + (0.5*(4.0 * LJ_E * (((LJ_S / r_ij)**12.0) - ((LJ_S / r_ij)**6.0))))
                !    !print *, "PE2", PE_A(i), i, j
                !    
                !! Forces between atoms w/i potential range
                !    Fx = 24.0 * ((LJ_E * (LJ_S**6.0) * (x(i) - x(j))) / (r_ij**8.0)) * (1.0 - (2.0 * ((LJ_S / r_ij)**6.0)))
                !    Fy = 24.0 * ((LJ_E * (LJ_S**6.0) * (y(i) - y(j))) / (r_ij**8.0)) * (1.0 - (2.0 * ((LJ_S / r_ij)**6.0)))
                !
                !    Fx_sum = Fx_sum + Fx
                !    Fy_sum = Fy_sum + Fy
                !end if
                
                ! W/o Potential Range
                PE_A(i) = PE_A(i) + (0.5*(4.0 * LJ_E * (((LJ_S / r_ij)**12.0) - ((LJ_S / r_ij)**6.0))))
                
                ! Forces from all other atoms acting on element i
                Fx = 24.0 * ((LJ_E * (LJ_S**6.0) * (x(i) - x(j))) / (r_ij**8.0)) * (1.0 - (2.0 * ((LJ_S / r_ij)**6.0)))
                Fy = 24.0 * ((LJ_E * (LJ_S**6.0) * (y(i) - y(j))) / (r_ij**8.0)) * (1.0 - (2.0 * ((LJ_S / r_ij)**6.0)))
                ! print *, "Fx", Fx, "Fy", Fy, "J", j, "I", i
                Fx_sum = Fx_sum + Fx
                Fy_sum = Fy_sum + Fy
                !print *, Fx, Fx_sum, j
                
                ! Other way of calculating PE (ignores previous bonds)
                !if (j>i) then
                    !print *, "R2", r_ij, i, j
                    !print *, "PE1", PE_A(i), i, j
                    !PE_A(i) = PE_A(i) + (4.0 * LJ_E * (((LJ_S / r_ij)**12.0) - ((LJ_S / r_ij)**6.0)))
                    !print *, "PE2", PE_A(i), i, j
                !end if
             !print *, "PE", PE_A(i), i, j
            
            
            end if
            
        end do
        !print *, PE_A(i)
        Fx_new(i) = -Fx_sum
        Fy_new(i) = -Fy_sum
        
        Fx_sum = 0
        Fy_sum = 0
    end do
    

end subroutine find_forces
    
