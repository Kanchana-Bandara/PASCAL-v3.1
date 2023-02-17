!#################################################################################################################################################################################################################
!									A SELF-ADAPTIVE INDIVIDUAL-BASE MODEL FOR COPEPOD BEHAVIOR AND LIFE HISTORY - BANDARA K. 2019							 			     #
!											MODEL BASE ADOPTED FROM BANDARA ET AL. 2018 & BANDARA ET AL. 2019										 			     #
!#################################################################################################################################################################################################################
!													MODULE FILE	:: IMPORTING												 			     	     	     #
!#################################################################################################################################################################################################################

module importfile
	implicit none
	
	contains
		subroutine importfile_2darr(i1, o1)															!imports a 2D array :: change dimensions accordingly
			implicit none
			
			integer :: i, j
			real (kind = 8), dimension(8760, 500) :: o1													!2D array for storing incoming data
			character (len = 46) :: i1																!gives a warning if the lengths of allocation and string lengths differ
			
			open(unit = 2000, file = i1)																!arbitrary unit number :: recycled :: +1for 1darr
			read(2000 ,*) o1
			close(2000)
			
			!testprint
			print *, "5x10 testprint of the 2D datafile: "
			
			do j = 1, 10
				print *, (o1(j, i), i = 1, 5)
			end do
			print *, " "
									
		end subroutine importfile_2darr
		
		subroutine importfile_1darr(i1, o1)
			implicit none
			
			real (kind = 8), dimension(8760) :: o1														!1D array for storing incoming data 
			character (len = 46) :: i1
			
			open(unit = 2001, file = i1)																!arbitrary unit number :: recycled :: 11for 2darr
			read(2001 ,*) o1
			close(2001)
			
			!testprint
			print *, "1x5 testprint of the 1D datafile: "
			print *, o1(1:5)
			print *, " "
						
		end subroutine importfile_1darr
		
		subroutine importfile_1darr_int(i1, o1)
			implicit none
			
			integer, dimension(8760) :: o1
			character(len = 46) :: i1
			
			open(unit = 2001, file = i1)
			read(2001 ,*) o1
			close(2001)
			
			!testprint
			print *, "1x5 testprint of the 1D datafile: "
			print *, o1(1:5)
			print *, " "
			
		end subroutine importfile_1darr_int
		
end module importfile
