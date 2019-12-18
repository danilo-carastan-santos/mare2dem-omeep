module binaryTree
    
    implicit none 
    
    public  :: insert_by_value, insert_by_index, get_data, print_tree, delete_tree, replace_value, get_length
    public  ::                  insert_by_index_z, get_data_z, print_tree_z, delete_tree_z, replace_value_z, get_length_z  
    private :: walk_data, walk_tree, walk_data_z, walk_tree_z
    
    type, public :: node 
        integer :: index
        real(8) :: value 
        type(node), pointer :: left  => null(), &
                               right => null()
    end type node 
 
     type, public :: node_z 
        integer :: index
        complex(8) :: value 
        type(node_z), pointer :: left  => null(), &
                                 right => null()
    end type node_z 
       
    integer :: ict 
    
    contains
    
!----------------------------------------------------------------------
! Insert by value:
!----------------------------------------------------------------------  
    recursive subroutine insert_by_value(t, index, value)
        type(node), pointer :: t  
        integer, intent(in) :: index
        real(8), intent(in) :: value
        
        ! If (sub)tree is empty, put number at root 
        if (.not. associated(t)) then
            allocate (t)
            t%index = index
            t%value = value
            
        ! Otherwise, insert into correct subtree 
        else if (value < t%value) then
            call insert_by_value(t%left, index, value) 
        else
            call insert_by_value(t%right, index, value)
        end if
        
    end subroutine insert_by_value
    
!----------------------------------------------------------------------
! Insert by index:
!---------------------------------------------------------------------- 
    recursive subroutine insert_by_index(t, index, value)
        type(node), pointer :: t  
        integer, intent(in) :: index
        real(8), intent(in) :: value
        
        ! If (sub)tree is empty, put number at root 
        if (.not. associated(t)) then
            allocate (t)
            t%index = index
            t%value = value
            
        ! Otherwise, insert into correct subtree 
        else if (index < t%index) then
            call insert_by_index(t%left, index, value) 
        elseif (index == t%index) then 
            t%value = t%value + value   
        else
            call insert_by_index(t%right, index, value)
        end if
        
    end subroutine insert_by_index
    
    recursive subroutine insert_by_index_z(t, index, value)
        type(node_z), pointer :: t  
        integer, intent(in) :: index
        complex(8), intent(in) :: value
        
        ! If (sub)tree is empty, put number at root 
        if (.not. associated(t)) then
            allocate (t)
            t%index = index
            t%value = value
            
        ! Otherwise, insert into correct subtree 
        else if (index < t%index) then
            call insert_by_index_z(t%left, index, value) 
        elseif (index == t%index) then 
            t%value = t%value + value      
        else
            call insert_by_index_z(t%right, index, value)
        end if
        
    end subroutine insert_by_index_z
    
!----------------------------------------------------------------------
! Get data
!----------------------------------------------------------------------  
    subroutine get_data(t,indices,values) ! gets all indices and values in order

        type(node), pointer, intent(in) :: t ! A tree
        integer, dimension(:), intent(out) :: indices
        real(8), dimension(:), intent(out) :: values
        
        ict = 0
        call walk_data(t,indices,values) 

    end subroutine get_data
    
    
    recursive subroutine walk_data(t,indices,values) 
    
        type(node), pointer, intent(in) :: t ! A tree
        integer, dimension(:), intent(inout) :: indices
        real(8), dimension(:), intent(inout) :: values
        
        if (associated(t)) then
            
            call walk_data(t%left,indices,values)
            ict = ict + 1
            indices(ict) = t%index
            values(ict)  = t%value
            call walk_data(t%right,indices,values)
         
        end if    
    
    end subroutine walk_data

    subroutine get_data_z(t,indices,values) ! gets all indices and values in order

        type(node_z), pointer, intent(in)  :: t ! A tree
        integer, dimension(:), intent(out) :: indices
        complex(8), dimension(:), intent(out) :: values
        
        ict = 0
        call walk_data_z(t,indices,values) 

    end subroutine get_data_z
    
    
    recursive subroutine walk_data_z(t,indices,values) 
    
        type(node_z), pointer, intent(in) :: t ! A tree
        integer, dimension(:), intent(inout) :: indices
        complex(8), dimension(:), intent(inout) :: values
        
        if (associated(t)) then
            
            call walk_data_z(t%left,indices,values)
            ict = ict + 1
            indices(ict) = t%index
            values(ict)  = t%value
            call walk_data_z(t%right,indices,values)
         
        end if    
    
    end subroutine walk_data_z
    
        
!----------------------------------------------------------------------
! Get length:
!----------------------------------------------------------------------      
     function get_length(t)  result(length)
     
        type(node), pointer :: t ! A tree
        integer :: length
        
        ict = 0
        call walk_tree(t)
        length = ict
        

    end function get_length
    
    recursive subroutine walk_tree(t)    
        type(node), pointer :: t ! A tree
        
        if (associated(t)) then
            call walk_tree(t%left) 
            call walk_tree(t%right)
            ict = ict + 1
        end if
    
    end subroutine walk_tree
 
      function get_length_z(t)  result(length)
     
        type(node_z), pointer :: t ! A tree
        integer :: length
        
        ict = 0
        call walk_tree_z(t)
        length = ict
        

    end function get_length_z
    
    recursive subroutine walk_tree_z(t)    
        type(node_z), pointer :: t ! A tree
        
        if (associated(t)) then
            call walk_tree_z(t%left) 
            call walk_tree_z(t%right)
            ict = ict + 1
        end if
    
    end subroutine walk_tree_z
       
!----------------------------------------------------------------------
! Print:
!----------------------------------------------------------------------  
    recursive subroutine print_tree(t) ! Print tree in  order

        type(node), pointer :: t ! A tree

        if (associated(t)) then
            call print_tree(t%left) 
            write(*,*) t%index, t%value
            call print_tree(t%right)
        end if

    end subroutine print_tree
 
    recursive subroutine print_tree_z(t) ! Print tree in  order

        type(node_z), pointer :: t ! A tree

        if (associated(t)) then
            call print_tree_z(t%left) 
            write(*,*) t%index, t%value
            call print_tree_z(t%right)
        end if

    end subroutine print_tree_z

!----------------------------------------------------------------------
! Replace value at a specific index
!----------------------------------------------------------------------
     recursive subroutine replace_value(t, index, value)
      
        type(node), pointer :: t  
        integer, intent(in) :: index
        real(8), intent(in) :: value
        
        if (.not. associated(t)) then
            write(*,*) 'error replacing value at index: ',index
            write(*,*) 'I did not find that index...'
            return 
        elseif (index < t%index) then
            call replace_value(t%left, index, value) 
        elseif (index == t%index)  then
            t%value = value   
        else
            call replace_value(t%right, index, value)
        end if
        
    end subroutine replace_value
    
     recursive subroutine replace_value_z(t, index, value)
      
        type(node_z), pointer :: t  
        integer, intent(in) :: index
        complex(8), intent(in) :: value
        
        if (.not. associated(t)) then
            write(*,*) 'error replacing value at index: ',index
            write(*,*) 'I did not find that index...'
            return 
        elseif (index < t%index) then
            call replace_value_z(t%left, index, value) 
        elseif (index == t%index)  then
            t%value = value   
        else
            call replace_value_z(t%right, index, value)
        end if
        
    end subroutine replace_value_z
           
!----------------------------------------------------------------------
! Delete tree:
!----------------------------------------------------------------------
     recursive subroutine delete_tree(t) ! delete tree in  order

        type(node), pointer :: t ! A tree

        if (associated(t)) then
        
            ! Delete left and right subtrees
            call delete_tree(t%left) 
            call delete_tree(t%right)
        
            deallocate(t)
        
        end if

    end subroutine delete_tree   
    
    recursive subroutine delete_tree_z(t) ! delete tree in  order

        type(node_z), pointer :: t ! A tree

        if (associated(t)) then
        
            ! Delete left and right subtrees
            call delete_tree_z(t%left) 
            call delete_tree_z(t%right)
        
            deallocate(t)
        
        end if

    end subroutine delete_tree_z   
      
end module binaryTree