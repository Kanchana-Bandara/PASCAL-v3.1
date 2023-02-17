!#################################################################################################################################################################################################################
!									A SELF-ADAPTIVE INDIVIDUAL-BASE MODEL FOR COPEPOD BEHAVIOR AND LIFE HISTORY - BANDARA K. 2019							 			     #
!											MODEL BASE ADOPTED FROM BANDARA ET AL. 2018 & BANDARA ET AL. 2019										 			     #
!#################################################################################################################################################################################################################
!														MODULE FILE	:: SEEDING												 			     	     	     #
!#################################################################################################################################################################################################################

module seeding																!outputs pseudorandom numbers into a array of desired dimensions (currently supports 1D)
	implicit none
	
	contains
		subroutine seeding_uniform(i1, i2, o1) 										!takes the end vector and the length of random numbers as inputs from a uniform prob. distribution
			implicit none
			
			integer, intent(in) :: i1, i2										!input variables : i1 = NStart, i2 = NMax :: see main program
			integer :: i													!iteration tracker
			real (kind = 8), dimension(i2) :: o1									!output vector to NMax
			real (kind = 8) :: RN												!variable to save random numbers iteratively
			
			do i = 1, i1
				call random_number(RN)
				o1(i) = RN
			end do
			
			if(i1 < i2) then													!fills the rest of the vector (i2-i1) with a uniform value of -1 for easy indexing
				o1(i1:i2) = -1.00
			end if
			
		end subroutine seeding_uniform
		
		subroutine seeding_uniform_irange(i1, i2, i3, i4, o1)							!generates pseudorandom numbers in a given range (currently supports 1D)
			implicit none						
			
			integer, intent(in) :: i1, i2, i3, i4									!input variables :: i1 =  NStart; i2 = lower bound; i3 = upper bound, i4 = NMax
			integer :: i													!iteration tracker
			real (kind = 8), dimension(i4) :: o1_raw, o1								!output vector to NMax
			real (kind = 8) :: RN												!variable to save random numbers iteratively
			
			do i = 1, i4													!this iterates over NMax (i4) to avoid range scaling issues
				call random_number(RN)
				o1_raw(i) = RN
			end do
			
			o1 = o1_raw * (i3 - i2) + i2											!scaling to a given range :: i2 = min, i3 = max
			
			if(i1 < i4) then													!fills the rest of the vector (i2-i1) with a uniform value of -1 for easy indexing
				o1(i1:i4) = -1.00
			end if
			
		end subroutine seeding_uniform_irange
		
		subroutine seeding_seedvalint(i1, i2, i3, i4, o1)								!assigns integer values for seeding (non-random; currently supports 1D)
			implicit none
			
			integer, intent(in) :: i1, i2, i3, i4									!i1 = NStart, i2 = NMax, i3 = seedval, i4 = default val, o1 = output vector
			integer, dimension(i2) :: o1											!output vector NMax 
			
			o1(1:i1) = i3													!seeds value to the seeded elements with a length of NStart
			o1(i1 + 1: i2) = i4												!assigns default value to the rest of the array from (NStart + 1 : NMAx)
			
		end subroutine seeding_seedvalint
		
		subroutine seeding_seedvalreal(i1, i2, i3, i4, o1)								!assigns real values for seeding (non-random; currently supports 1D)
			implicit none
			
			integer, intent(in) :: i1, i2										!i1 = NStart, i2 = NMax
			real, intent(in) :: i3, i4											!i3 = seedval, i4 = default val, o1 = output vector
			real (kind = 8), dimension(i2) :: o1									!output vector NMax 
			
			o1(1:i1) = i3													!seeds value to the seeded elements with a length of NStart
			o1(i1 + 1: i2) = i4												!assigns default value to the rest of the array from (NStart + 1 : NMAx)
			
		end subroutine seeding_seedvalreal
		
end module seeding
