!==================================================================================================================================! 
!===================================================================================================================== EM_constants
!==================================================================================================================================! 
    module EM_constants    
    
    implicit none      
!
! A few constants:  If you change these to something bad, you are truly a demon child and should be banished from the planet :)
!
    real(8),    parameter :: PI  = 3.141592653589793d0
    real(8),    parameter :: PI2 = 2d0*pi        
    real(8),    parameter :: MU0 = 1.256637061435917d-6 
    real(8),    parameter :: EPSILON = 0d0 ! 8.854187817620391d-12  ! if this is non-zero, there can be problems when
                                                                    ! kx^2 == omega^2*mu*epsilon
    
    complex(8), parameter :: IC  = (0d0,1d0)     
    real(8),    parameter :: rad2deg = 57.295779513082323d0
    real(8),    parameter :: deg2rad =  0.017453292519943d0
    integer,    parameter :: eps(6) = (/ 1,2,3,1,2,3 /)  
            
    end module EM_constants