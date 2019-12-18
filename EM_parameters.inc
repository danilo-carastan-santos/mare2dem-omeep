!
! CSEM1D_parameters.inc
!
! Stores magic numbers for CSEM and MT data types
!

!
! Data Type Index Dictionary
!
! These are the codes used in the data file parameter table.
! 
!
! CSEM Data:
! 
! This format is compatible with 2D and 3D CSEM modeling.
! 
!   CSEM Electric Data, Real/Imag:
    integer, parameter :: indRealEx = 1
    integer, parameter :: indImagEx = 2
    integer, parameter :: indRealEy = 3
    integer, parameter :: indImagEy = 4
    integer, parameter :: indRealEz = 5
    integer, parameter :: indImagEz = 6

    
!   CSEM Magnetic Data, Real/Imag:
    integer, parameter :: indRealBx = 11
    integer, parameter :: indImagBx = 12
    integer, parameter :: indRealBy = 13
    integer, parameter :: indImagBy = 14
    integer, parameter :: indRealBz = 15
    integer, parameter :: indImagBz = 16
    
!   CSEM Electric Data, Amp/Phs:
    integer, parameter :: indAmpEx = 21
    integer, parameter :: indPhsEx = 22
    integer, parameter :: indAmpEy = 23
    integer, parameter :: indPhsEy = 24
    integer, parameter :: indAmpEz = 25
    integer, parameter :: indPhsEz = 26
    
!   CSEM Electric Data, log10 Amp:
    integer, parameter :: indLog10AmpEx = 27
    integer, parameter :: indLog10AmpEy = 28
    integer, parameter :: indLog10AmpEz = 29
    
!   CSEM Magnetic Data, Amp/Phs:
    integer, parameter :: indAmpBx = 31
    integer, parameter :: indPhsBx = 32
    integer, parameter :: indAmpBy = 33
    integer, parameter :: indPhsBy = 34
    integer, parameter :: indAmpBz = 35
    integer, parameter :: indPhsBz = 36
    
!   CSEM Magnetic Data, log10 Amp:
    integer, parameter :: indLog10AmpBx = 37
    integer, parameter :: indLog10AmpBy = 38
    integer, parameter :: indLog10AmpBz = 39


!   CSEM Horizontal Polarization Ellipses:
    integer, parameter :: iPEmax = 41
    integer, parameter :: iPEmin = 42
    
    integer, parameter :: iPBmax = 43
    integer, parameter :: iPBmin = 44
    

!   
! Magnetotelluric Data:  
!
! This format is compatible with 2D and 3D MT modeling codes.
!
!   Apparent Resistivity and Phase:
    integer, parameter :: indRhoZXX  = 101   ! reserved for 3D
    integer, parameter :: indPhsZXX  = 102   ! reserved for 3D
    integer, parameter :: indRhoZXY  = 103   ! 
    integer, parameter :: indPhsZXY  = 104   !  
    integer, parameter :: indRhoZYX  = 105   !
    integer, parameter :: indPhsZYX  = 106   ! 
    integer, parameter :: indRhoZYY  = 107   ! reserved for 3D
    integer, parameter :: indPhsZYY  = 108   ! reserved for 3D
   
    ! log10(apparent resistivity)
    integer, parameter :: indlog10RhoZXX  = 121 ! reserved for 3D  
    integer, parameter :: indlog10RhoZXY  = 123   
    integer, parameter :: indlog10RhoZYX  = 125   
    integer, parameter :: indlog10RhoZYY  = 127 ! reserved for 3D  

        
!   Impedances: 
    integer, parameter :: indRealZXX  = 111   ! reserved for 3D
    integer, parameter :: indImagZXX  = 112   ! reserved for 3D
    integer, parameter :: indRealZXY  = 113   ! 
    integer, parameter :: indImagZXY  = 114   ! 
    integer, parameter :: indRealZYX  = 115   !
    integer, parameter :: indImagZYX  = 116   ! 
    integer, parameter :: indRealZYY  = 117   ! reserved for 3D
    integer, parameter :: indImagZYY  = 118   ! reserved for 3D 
    
! Magnetic Tippers:
    integer, parameter :: indRealMZY  = 133  
    integer, parameter :: indImagMZY  = 134  

! To do:  other transfer functions ( marine Jz tipper, interstation E-E, H-H and hybrid E-H), etc

