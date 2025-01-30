!  Console1.f90 
!
!  FUNCTIONS:
!  Console1 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Console1
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program Console1

    implicit none

    ! Variables - Parameter keyword = values cannot be changed --> that is why a new copy of the variable is made
    
    real, parameter :: k = 200.0 !Spring constant [N/m]
    real, parameter :: x_0 = 0 !Equilibrium Position [m]
    real, parameter :: m_1  = 1.0 !Mass 1 [kg]
    real, parameter :: m_2 = 1.0 !Mass 2 [kg]
    real, parameter :: v_10 = 0 !Initial Velocity 1 [m/s]
    real, parameter :: v_20 = 0 !Initial Velocity 2 [m/s]
    real, parameter :: x_10 = -0.05 !Initial Position 1 [m]
    real, parameter :: x_20 = 0.05 !Initial Position 2 [m]
    real, parameter :: t = 1 !Time [s]
    real, parameter :: dt = 0.01 !Timestep [s]
    real :: t_0 = 0 !Initial Time [s]
    real :: f_0 = k*(x_20 - x_10 - x_0) !Initial Force
    real :: f_old
    real :: f_new
    real :: x_old1
    real :: x_old2
    real :: v_old1
    real :: v_old2
    f_old = f_0
    x_old1 = x_10
    x_old2 = x_20
    v_old1 = v_10
    v_old2 = v_20
    
    ! Body of Console1
    open(1, file = 'AtomSim_HW1.csv', status = 'new')
        
        write (1,*) "Time", ",", "Position 1", ",", "Position 2", ",", "Velocity 1", ",", "Velocity 2"
        
        do while (t_0 <= t)
            !write (1, *) t_0,x_old1,x_old2,v_old1,v_old2
            write(1,'(*(g0,","))') t_0,x_old1,x_old2,v_old1,v_old2
            print *, "Time", t_0
            print *, "Positions (1|2)", x_old1, x_old2
            print *, "Velocities (1|2)", v_old1, v_old2
            !print *, "Before", x_old1, x_old2, v_old1, v_old2
            !print *, "Old x", x_old1, x_old2
            ! Position @ current time
            x_old1 = x_old1 + (dt*v_old1) + ((dt**2)*f_old/(2*m_1))
            x_old2 = x_old2 + (dt*v_old2) - ((dt**2)*f_old/(2*m_2))
            !print *, "New x", x_old1, x_old2
            ! New Force
            f_new = k * (x_old2 - x_old1 - x_0)
            !print *, "New f", f_new
            ! New Velocities
            !print *, "Old v", v_old1, v_old2
            v_old1 = v_old1 + ((dt*(f_old+f_new))/(2*m_1))
            v_old2 = v_old2 - ((dt*(f_old+f_new))/(2*m_2))
            !print *, "New v", v_old1, v_old2
            ! Copy New to Old
            !print *, "Old f", f_old
            f_old = f_new
            !print *, "Copy f", f_old
            t_0 = t_0 + dt
            print *, "After", x_old1, x_old2, v_old1, v_old2
        end do

    close(1)
    
end program Console1

