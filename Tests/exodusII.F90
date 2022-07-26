module exodusIIF90_interface
implicit none
#include <exodusII.inc>

	interface
	   subroutine exppv(idexo, time_step, var_type, var_index, &
						obj_id, start_index, num_entities,     &
						var_vals, ierr)	
		INTEGER idexo           ! (R)
		INTEGER time_step       ! (R)
		INTEGER var_type        ! (R)
		INTEGER var_index       ! (R)
		INTEGER obj_id          ! (R)
		INTEGER start_index     ! (R)
		INTEGER num_entities    ! (R)
		REAL    var_vals(*)     ! (W)
		INTEGER ierr            ! (W)
		end subroutine exppv
	end interface
end module exodusIIF90_interface

module exodusIIF90
	use exodusIIF90_interface, &
	    ex_put_partial_var => exppv
end module exodusIIF90
