module randomnum
	implicit none
	
	contains
		subroutine randomnum_int_srange(lowerbound, upperbound, output)						!generates a uniform randomnumber /INT/ in a given range
			implicit none
			
			integer, intent(in) :: lowerbound, upperbound								!two input variables /int/		
			integer :: output													!the output variable /int/
			real(kind = 8) :: rn												!interim variable for random number generation
			
			call random_seed()
			call random_number(rn)
			
			rn = rn * (upperbound - lowerbound) + lowerbound							!standardization in a given range :: real number :: expect a truncation error
			output = nint(rn)													!integer format conversion
			
		end subroutine randomnum_int_srange
		
end module randomnum
