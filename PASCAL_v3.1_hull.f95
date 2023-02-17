!################################################################################################################################################################################################################
!									A SELF-ADAPTIVE INDIVIDUAL-BASE MODEL FOR COPEPOD BEHAVIOR AND LIFE HISTORY - BANDARA K. 2019							 			    #
!											MODEL BASE ADOPTED FROM BANDARA ET AL. 2018 & BANDARA ET AL. 2019										 			    #
!################################################################################################################################################################################################################
!																BASE MODEL FILE													 			    #
!################################################################################################################################################################################################################

program IBSM1D_v1
	!1. modules																				!modules used in the model
	use importfile
	use simulation
	use randomnum
	use csv_file
	
	implicit none
	
	!2. initialize variables and constants
	!2.1 essentials
	integer, parameter :: NMax = 1000000															!system-limited carrying capacity of the population
	integer :: ind = 0																		!individual ID (for iteration)
	logical :: LoopingCondition = .true.															!logical variable to proceed iteration
	logical :: Seeding = .true.																	!logical variable to enable seeding
	logical :: Spawning = .false.																	!logical variable to enable seeding of spawned eggs
	integer :: CurrentSeedSize = 10																!no.of initial seeds launched in each iteration
	integer :: CurrentSpawnSize = 0																!no.of spawned seeds launched in each iteration
	integer :: i = 1 																			!universal iterator-1
	integer :: j = 1																			!universal iterator-2
	integer :: RNI = -1																		!universal random number destination /int/
	real(kind = 8) :: RN = -1.00																	!universal random number destination /real-8b/
	character(len = 8) :: filename																!filename for the output file
	character(len = 36) :: directory = "/home/kanchana/Modelling/M1/Outputs/"									!output directory name
	character(len = 4) :: extension = ".txt"															!output file extension
	character(len = 48) :: output																!output name concatenation
	integer :: progress = 0																		!iteration progress
	integer :: popsize = 0																		!current population size
	integer, dimension(12) :: months = (/1, 768, 1440, 2184, 2904, 3648, 4368, 5112, 5856, 6576, 7320, 8040/)				!30-d breaks of the year
	
	!2.2 evolvable parameters
	real(kind = 8), dimension(NMax) :: GAP = -1.00														!growth allocation parameter
	real(kind = 8), dimension(NMax) :: SDP = -1.00														!seasonal descent parameter
	real(kind = 8), dimension(NMax) :: SAP = -1.00														!seasonal ascent parameter
	real(kind = 8), dimension(NMax) :: BSP = -1.00														!body size parameter

	real(kind = 8) :: CurrentGAP = -1.00															!GAP for currently processed individual 
	real(kind = 8) :: CurrentSDP = -1.00															!SDP for currently processed individual	
	real(kind = 8) :: CurrentSAP = -1.00															!SAP for currently processed individual
	real(kind = 8) :: CurrentBSP = -1.00															!BSP for currently processed individual
	
	!2.3 state variables
	integer, dimension(NMax) :: DevStage = -1														!the developmental stage (egg = 1; adult = 13) :: storage array
	integer, dimension(NMax) :: Fecundity = -1														!no.of eggs produced per lifetime :: storage array
	integer, dimension(NMax) :: LDS = 0															!life and death status (0 = dead; 1 = alive) :: storage array
	real(kind = 8), dimension(NMax) :: CWeight = -1.00					 								!structural mass (carbon, ug) :: storage array
	real(kind = 8), dimension(NMax) :: SWeight = -1.00													!mass of the energy reserve (carbon, ug) :: storage array

	integer :: CurrentDevStage = -1																!developmental stage of a given individual at a given time
	integer :: CurrentFecundity = -1																!fecundity of a given individual at a given time
	real(kind = 8) :: CurrentCWeight = -1.00															!structural mass of a given individual at a given time
	real(kind = 8) :: CurrentSWeight = -1.00															!energy reserve mass of a given individual at a given time
	
	!2.4 performance tracers
	integer, dimension(NMax) :: TSD = -1															!time of seasonal descent : storage array 
	integer, dimension(NMax) :: TSA = -1															!time of seasonal ascent : storage array
	integer, dimension(NMax) :: OSP = -1															!onset of egg production : storage array
	real(kind = 8), dimension(NMax) :: SSM = -1.00														!size at osp : storage array
	real(kind = 8), dimension(NMax) :: CWSD = -1.00													!structural mass at seasonal descent
	real(kind = 8), dimension(NMax) :: SWSD = -1.00													!reserve mass at seasonal descent
	real(kind = 8), dimension(NMax) :: Fitness = -1.00													!the no. of viable offspring produced : storage array
	
	!2.5 time-related variables
	integer :: CurrentPing = 0																	!current time, in pings
	integer :: CurrentYear	= 1																	!current time, in years
	integer, parameter :: EndYear = 50																!termination year
	integer, dimension(NMax) :: TB = -1															!birth time :: storage array
	
	!2.6 position-related variables
	integer, dimension(NMax) :: Depth = -1															!vertical position : storage array
	integer, parameter :: MaxDepth = 500															!maximum simulated depth
	integer, parameter :: MinDepth = 1																!minimum simulated depth
	integer :: CurrentDepth	 = -1																	!vertical position of a given individual at a given time
	integer :: MaxDistance = -1																	!the potential maximum distance can be searched by an individual
	integer :: MFloor = -1																		!deepest depth bin of the search
	integer :: MCeiling = -1																	!shallowest depth bin of the search
	integer :: SDistance = -1																	!the search range :: MFloor - MCeiling
	integer :: VShift	 = -1																		!vertical distance shift between time pings
	real(kind = 8), allocatable, dimension(:) :: PVisit												!probability list of visiting a depth bin in search range	
	integer, dimension(8760, 500) :: VPos = 0														!annual vertical distribution of the population
	integer :: CurrentVPos = 0																	!vpos of a given individual at a given time
	integer, dimension(1) :: ML = -1																!placeholder for maxloc() function output
	
	!2.7 variables related to physical submodel
	real(kind = 8), dimension(8760, 500) :: Temp = -1.00													!matrix for temperature (oC)
	real(kind = 8), dimension(8760, 500) :: Fcon = -1.00													!matrix for food concentration (mgChl-a. m-3)
	real(kind = 8), dimension(8760) :: Sird = -1.00													!vector for surface irradiance (umol photones m-2 sec-1)
	integer, dimension(8760) :: MLD = -1															!vector for mixed layer depth (m)
	
	real(kind = 8) :: CurrentTemp = -1.00															!temperature at a given time in a given depth bin
	real(kind = 8) :: CurrentFcon = -1.00															!food concentration at a given time in a given depth
	real(kind = 8) :: CurrentSird = -1.00															!surface irradiance at a given time in a given depth bin
	
	real(kind = 8), allocatable, dimension(:) :: CurrentTempRange											!allocatable array for temperature at the a given search range
	real(kind = 8), allocatable, dimension(:) :: CurrentFconRange											!allocatable array for food concentration at the a given search range
	
	!2.8 variables related to development and growth submodel
	integer, dimension(NMax) :: StageDuration = -1														!stage durations (h) :: storage array
	integer :: CurrentStageDuration = -1															!stage duration of a given individual at a given time
	real(kind = 8), dimension(NMax) :: DTTrack = -1.00													!traced development time :: storage array
	real(kind = 8) :: CurrentDevTime = -1.00															!estimated dev.time for a given individual at a given time
	real(kind = 8) :: MeanDevTime = -1.00															!geometric mean of lifetime dev.time of a given individual at a given time
	real(kind = 8) :: DTProd = -1.00																!product of CurrentDevTime and DTTrack of a given individual for geom.mean.est
	real(kind = 8) :: CMM = -1.00																!stage specific critical molting mass
	real(kind = 8) :: CWeightCeiling = -1.00															!maximum possible weight gained at a given body mass trajectory
	real(kind = 8), dimension(12) :: WMin = (/0.23, 0.23, 0.23, 0.49, 0.70, 0.95, 1.34, 2.94, 6.35, 13.96, 36.58, 141.84/)		!lower trajectory of critical molting masses
	real(kind = 8), dimension(12) :: WMax = (/0.23, 0.23, 0.23, 0.63, 1.01, 1.49, 2.26, 5.51, 12.82, 30.20, 84.36, 349.34/)	!upper trajectory of critical molting masses
	
	real(kind = 8), parameter :: AsimCoef = 0.60														!assimilation coefficient
	real(kind = 8), dimension(NMax) :: MaxCWeight = -1.00												!maximum structural mass reached by a given individual : storage array
	real(kind = 8) :: CurrentMaxCWeight = -1.00														!MaxCWeight of a given individual at a given time
	real(kind = 8) :: CurrentGGRate = -1.00															!estimated gross growth rate (assimilation rate; ug C ind-1 h-1)
	real(kind = 8) :: CurrentNGRate = -1.00															!estimated net growth rate (ug C ind-1 h-1)
	real(kind = 8), allocatable, dimension(:) :: PGP													!allocatable array for percieved growth potential at the a given search range
	real(kind = 8), parameter :: MaxPG = 1710.72														!maximum percieved growth potential :: Fcon x Temp (K) = 285.15 x 6
	real(kind = 8) :: MaxRG = -1.00																!theoreitical maximum growth rate for a given body mass 
	
	!2.9 variables related to reproductive submodel
	real(kind = 8), dimension(12, NMax) :: Pairing	 = -1.00												!2D array to store Ma, Pa & Da values of evolvable parameters
	integer, dimension(NMax) :: NEggs	= 0															!no.of eggs produced at each time ping : storage array
	integer :: CurrentNEggs = -1																	!no. of eggs produced by a given female at a given time
	integer, parameter :: FecundityCeiling = 1000														!max.fecundity allowed for female
	character(len = 1), dimension(NMax) :: Gender = "U"													!gender of each individual : storage array
	integer :: pairID	 = 1																		!tracking variable for pair selection, crossover and mutation
	real(kind = 8), dimension(4) :: RNCX = -1.00														!random numbers for crossover operation (for 4 x evolvable parameters)
	real(kind = 8), dimension(4) :: CXProb = -1.00														!random numbers for crossover decision (for 4 x evolvable parameters)
	real(kind = 8), dimension(4) :: MProb = -1.00														!random numbers for mutation decision (for 4 x evolvable parameters)
	real(kind = 8), dimension(4) :: RNM = -1.00														!random numbers for mutation operation(for 4 x evolvable parameters)
	real(kind = 8) :: CXTemp = -1.00																!temporary placeholder for crossover value
	real(kind = 8) :: MTemp = -1.00																!temporary placeholder for mutation value
	real(kind = 8), parameter :: GenderTreshold = 0.50													!gender selection treshold (50:50 male:female)
	real(kind = 8), parameter :: CXTreshold = 0.70														!crossover treshold (70% per evol.par)
	real(kind = 8), parameter :: MTreshold = 0.10														!mutation treshold (10% per evol.par)
	real(kind = 8), dimension(NMax) :: RWeight = -1.00													!mass allocated to reproductive development : storage array
	real(kind = 8) :: CurrentRWeight = -1.00															!reproductive allocation for a given individual at a given time
	character(len = 1), dimension(NMax) :: ISState = "U"												!insemination state :: U = undefined, I = inseminated, N = non-inseminated
	character(len = 1) :: CurrentISState = "U"														!insemination state of a given female/male? at a given time
	integer, dimension(NMax) :: IS_F = 0															!insemination array for females
	integer, dimension(NMax) :: IS_M = 0															!insemination array for males
	integer :: NISF = 0																		!number of inseminated females
	integer :: NISM = 0																		!number of inseminated males
	integer, dimension(NMax) :: SelectedMate	 = -1														!the id of randomly selected mate : storage array
	integer :: FID = -1																		!individual id of the female in mate selection iteration
	integer :: MID = -1																		!individual id of the selected mate
	integer :: NEPF = 0																		!number of egg producing females
	
	!2.10 variables related to survival submodel
	real(kind = 8) :: CurrentBMRate = -1.00															!estimated basal metabolic rate
	real(kind = 8) :: CurrentAMRate = -1.00															!estimated active metabolic rate
	real(kind = 8) :: CWLoss = -1.00																!catabolized structural mass compared to max.structural mass
	integer, dimension(NMax) :: Age = -1															!age of the animal : storage array
	integer :: CurrentAge																		!age of an animal at a given time
	real(kind = 8), dimension(8760, 500) :: VPR = -1.00													!matrix for temperature and visual predation risk values (dim.less)
	real(kind = 8) :: CurrentStarvationRisk = 0.00														!starvation risk of a given individual at a given time
	real(kind = 8) :: CurrentMortalityRisk = 0.00														!sum of all sources of mortality at a given time for a given individual
	real(kind = 8) :: CurrentVPRisk = 0.00															!visual predation risk of a given individual at a given time
	real(kind = 8) :: CurrentNVPRisk = 0.00															!other sources of mortality of a given individual at a given time
	real(kind = 8) :: CurrentDDMRisk = 0.00															!density-dependent mortality risk :: based on artificial pop. ceiling
	real(kind = 8) :: VPRScalar = 0.1																!scalar for visual predation risk
	real(kind = 8), allocatable, dimension(:) :: CurrentVPRiskRange										!allocatable array for visual predation risk in a given search range
	real(kind = 8), allocatable, dimension(:) :: PVPR													!allocatable array for pervieved visual predation risk
	
	character(len = 1), dimension(NMax) :: DiapState = "U"												!diapause state : storage array
	character(len = 1) :: CurrentDiapState = "U"														!diap.state of a given ind.at t. :: "U"-undf, "A"-active, "R"-ready,"D"-diap.
	integer, dimension(NMax) :: DiapDepth = -1														!overwintering depth : sorage array
	integer :: CurrentDiapDepth = -1																!overwintering depth of a given individual at a given time during diapause
	
	integer, dimension(NMax) :: LCL = -1															!life cycle length : storage array
	integer :: CurrentLCL = -1																	!prescribed LCL of a given individual
	real(kind = 8) :: LCLProb = -1.00																!probability of reproducing within the same calendar year
	
	real(kind = 8) :: SurvivalProb																!probability of survival for a given individual from time t to t+1
	
	!2.11 variables related to I/O functionality
	integer :: n_sd = 0																		!count of no seasonal decents per year
	integer :: n_sa = 0																		!count of no of seasonal ascents per year
	integer :: n_sm = 0																		!count of sexually mature adults
	integer :: sum_tsd = 0																		!aggregated timing of seasonal descent
	integer :: sum_tsa = 0																		!aggregated timing of seasonal ascent
	real(kind = 8) :: sum_cwsd = 0.00																!mean structural mass at seasonal descent
	real(kind = 8) :: sum_swsd = 0.00																!mean reserve mass at seasonal descent
	real(kind = 8) :: sum_cwsm = 0.00																!aggregated size of sexual maturity
	real(kind = 8), dimension(EndYear, 5) :: IO_sdsa = 0.00												!output array of seasonal descent and ascent properties
	integer, dimension(8760) :: log_popsize = 0														!hourly annual population size log
	integer, dimension(8760) :: log_neggs = 0														!hourly annual log of egg production
	integer :: NM = 0																			!no.of males
	integer :: NF = 0																			!no.of females
	integer, dimension(8760) :: log_nmales = 0														!hourly annual log of male:female sex-ratio 
	integer, dimension(8760) :: log_nfemales = 0														!hourly annual log of male:female sex-ratio
	integer, dimension(8760, 13) :: log_devstage = 0													!hourly annual log of developmental stage composition
	integer :: devstagecount = 0																	!updated value of devstage count
	
	
	!3. generate model environment
	print "(a42)", "[READING DISC] the system is loading files"
	
	call importfile_2darr("/home/kanchana/Modelling/M1/DataFiles/Temp.txt", Temp)									!imports temperature 2D array :: 46 char variable name
	call importfile_2darr("/home/kanchana/Modelling/M1/DataFiles/Fcon.txt", Fcon)									!imports food.con 2D array :: 46 char variable name
	call importfile_1darr("/home/kanchana/Modelling/M1/DataFiles/Sird.txt", Sird)									!imports surf.irad 1D array :: 46 char variable name
	call importfile_1darr_int("/home/kanchana/Modelling/M1/DataFiles/MLdp.txt", MLD)								!imports mixed layer depth /int/ to 1D array :: 46 char variable name
	call importfile_2darr("/home/kanchana/Modelling/M1/DataFiles/VIRS.txt", VPR)									!imports vertical pred.risk 2D array :: 46 char variable name
	
	print "(a7)", "[READY]"
	print "(a2, i3, a3, 1X, a18, i8)", "[ ", progress, "% ]", "population size : ", popsize							!prints progress and population size
	
	!4. simulation
	!4.1 time-level iteration
	do while(LoopingCondition)																	!main loop iterates over time until LoopCondition = .true.
		!4.1.1 timecoding and ping mapping
		CurrentPing = CurrentPing + 1
		
		if(CurrentPing == 8761) then																!at the end of a calendar year;
			CurrentYear = CurrentYear + 1															!year is updated
			CurrentPing = 1																	!ping counter is reset
			VPos = 0																		!set the vertical distibution tracker to zero
		end if																					
		
		if(CurrentYear > 1)then																	!seeding occurs only in the 1st year
			
			if(popsize <= 100)then
				Seeding = .true.
				CurrentSeedSize = 10
			else
				Seeding = .false.
				CurrentSeedSize = 0
			end if

		else
			Seeding = .true.
			CurrentSeedSize = 10																!this has to be reset here because seeding loop rewrites this value
		end if
		
		CurrentSpawnSize = sum(NEggs)																!sets the number of eggs to sprawn
		
		if(CurrentSpawnSize > 0)then																!sets whether or not to initiate spawning
			Spawning = .true.
		else
			Spawning = .false.
		end if
		
		!4.1.2 seeding												
		if(Seeding)then																		!seeding occurs if Seeding = .true.
			goto 100																		!goto seeding codeblock :: ref.address #100
		else
			goto 101																		!goto no seeding codeblock :: ref.address #101
		end if
		
		100 do i = 1, size(LDS)																!searches for space in the population to seed
			
			if(LDS(i) == 0)then																!if there is a space in the population; does the following;
				CurrentSeedSize = CurrentSeedSize - 1												!one seed less to be seeded
				TB(i) = CurrentPing															!time of birth is set to current ping
				call random_seed()
				call random_number(RN)															!random number for GAP
				GAP(i) = RN																	!growth alloc. param. selected from uniformly random pdf
				LCLProb = 1.00/ (1.00 + exp((0.30 - RN) / 0.10))										!life cycle length prob. is a function of GAP
				call random_seed()
				call random_number(RN)															!life cycle length prob. checker
				
				if(RN > LCLProb)then
					LCL(i) = 0																!multyple generations per year (0-year)
				else
					LCL(i) = 1																!one generation per year (1-year)
				end if
				
				call random_seed()
				call random_number(RN)
				SDP(i) = RN																	!seas.des.param. selected from uniformly random pdf
				call random_seed()
				call random_number(RN)
				SAP(i) = RN																	!seas.asc.param. selected from uniformly random pdf
				call random_seed()
				call random_number(RN)
				BSP(i) = RN																	!body size param. selected from uniformly random pdf
				call random_seed()
				call random_number(RN)															!gender prob. checker
					
				if(RN <= GenderTreshold)then
					Gender(i) = "F"															!generates a female at a given (0.5) probability
				else
					Gender(i) = "M"															!generates a male at a given (0.5) probability
				end if
				
				Depth(i) = 1																!sets a basic current depth
				NEggs(i) = 0																!sets the number of eggs produced
				RWeight(i) = 0.00																!sets the mass allocated to reproductive development
				ISState(i) = "N"																!sets the insemination status
				StageDuration(i) = 0															!sets the stage duration
				DTTrack(i) = 1.00																!sets the tracking variable for development time est.
				CWeight(i) = 0.23																!sets the initial weight for 0.23 ugC for C.fin
				SWeight(i) = 0.00																!no stores to begin with (embroyonic stores disregarded)
				MaxCWeight(i) = 0.23															!sets the maximum structural mass
				DevStage(i) = 1																!starts at embryonic stage (= 1)
				Fecundity(i) = 0																!has no fecundity until adulthood
				LDS(i) = 1																	!animal is alive to begin with
				DiapState(i) = "A"															!sets current diapause state :: all active at the begining
				call randomnum_int_srange(300, 450, RNI)												!randomly selects diapause depth (*** LB: 300, UB: 450 ***)
				DiapDepth(i) = RNI
				TSD(i) = 0; CWSD(i) = 0.23; SWSD(i) = 0.00; TSA(i) = 0; SSM(i) = 0.00; OSP(i) = 0; Fitness(i) = 0.00		!performance tracers
					
				if(CurrentSeedSize <= 0)then														!seeding ends when all seeds have been seeded
					exit
				end if
				
			else
				cycle
			end if
			
		end do
		
		101 if(Spawning)then
			goto 102
		else
			goto 103
		end if
		
		!4.1.3 spawning
		102 pairID = 1																		!initiate sprawing iteration at pair 1
		
		do i = 1, size(LDS)
		
			if(LDS(i) == 0)then																!if there is a space in the population;
				CurrentSpawnSize = CurrentSpawnSize - 1												!one individual less to be spawned
				TB(i) = CurrentPing															!simple assignment of birth time
				call random_seed()
				call random_number(CXProb)														!generates crossover probabilities :: all param (= 4)
				call random_seed()
				call random_number(RNCX)														!generates crossover beta values :: (Radcliff, 1991):: all param (= 4)
				call random_seed()
				call random_number(MProb)														!generates mutation probabilities :: all param (= 4)
				call random_seed()
				call random_number(RNM)														!generates mutation-based random replacements:: all param (= 4)
						
				do j = 1, 4																	!generates new values for #4 evo.params following CX and mutation 
					if(CXProb(j) <= CXTreshold)then												!70% chance for crossover
						CXTemp = RNCX(j) * Pairing(j, pairID) + (1 - RNCX(j)) * Pairing(j + 4, pairID)			!blend crossover value
					else
						CXTemp = Pairing(j, pairID)												!no crossover :: direct Ma-value inheritance
					end if
							
					if(MProb(j) <= MTreshold)then
						MTemp = RNM(j)														!mutate with random replacement [0-1]
					else
						MTemp = CXTemp														!no mutation
					end if
							
					Pairing(j + 8, pairID) = MTemp												!new parameter value for the spawned individual
				end do
						
				GAP(i) = Pairing(9, pairID)														!growth alloc. param. from cx-m operations
				LCLProb = 1.00/ (1.00 + exp((0.30 - Pairing(9, pairID)) / 0.10))								!life cycle length prob. based on the value for GAP
				call random_seed()
				call random_number(RN)															!life cycle prob.cutoff
				
				if(RN > Pairing(9, pairID))then
					LCL(i) = 0																!multyple generations per year (0-year)
				else
					LCL(i) = 1																!one generation per year (1-year)
				end if
				
				SDP(i) = Pairing(10, pairID)														!seas.des.param. from cx-m operations
				SAP(i) = Pairing(11, pairID)														!seas.as.param. from cx-m operations
				BSP(i) = Pairing(12, pairID)														!b.size.param. from cx-m operations
						
				call random_seed()
				call random_number(RN)															!generate gender probability
						
				if(RN <= GenderTreshold)then
					Gender(i) = "F"															!generates a female at a given (0.5) probability
				else
					Gender(i) = "M"															!generates a male at a given (0.5) probability
				end if
				
				Depth(i) = 1																!sets a basic depth		
				StageDuration(i) = 0															!sets the stage duration
				DTTrack(i) = 1.00																!sets the tracking variable for development time est.
				CWeight(i) = 0.23																!sets the initial weight for 0.23 ugC for C.fin
				SWeight(i) = 0.00																!no stores to begin with (embroyonic stores disregarded)
				MaxCWeight(i) = 0.23															!sets the maximum structural mass
				DevStage(i) = 1																!starts at embryonic stage (= 1)
				Fecundity(i) = 0																!has no fecundity until adulthood
				LDS(i) = 1																	!animal is alive to begin with
				DiapState(i) = "A"															!sets current diapause state :: all active at the begining
				call randomnum_int_srange(300, 450, RNI)												!randomly selects diapause depth
				DiapDepth(i) = RNI
				TSD(i) = 0; CWSD(i) = 0.23; SWSD(i) = 0.00; TSA(i) = 0; SSM(i) = 0.00; OSP(i) = 0; Fitness(i) = 0.00		!performance tracers
				RWeight(i) = 0.00																!sets the mass allocated to reproductive development
				ISState(i) = "N"																!sets the insemination status
				
				NEggs(pairID) = NEggs(pairID) - 1													!sprawned eggs are removed from the system
				
				if(NEggs(pairID) == 0)then														!when eggs per given female is zero, move to next pair
					pairID = pairID + 1
				end if
						
				if(CurrentSpawnSize <= 0)then
					exit
				end if
					
			else
				cycle
			end if
			
		end do
			
		103 CurrentSpawnSize = 0																!force-resets the value for next iteration
		NEggs = 0																			!force-resets the value for next iteration
		Pairing = -1.00																		!force-resets the value for next iteration
		NISF = 0																			!force-resets the value for next iteration
		NISM = 0																			!force-resets the value for next iteration
		NEPF = 0																			!force-resets the value for next iteration
		IS_F = 0																			!force-resets the value for next iteration
		IS_M = 0																			!force-resets the value for next iteration
		NM = 0																			!force-resets the value for next iteration
		NF = 0																			!force-resets the value for next iteration
		
		popsize = sum(LDS)																	!estimated population size at the begining of time t
		
		!4.1.4 individual-level iteration
		do ind = 1, NMax																		!sub loop that iterates over NMax individuals
			
			if(LDS(ind) == 0)then																! <<< checks if a live individual is available in each slot
				cycle																		!if NA (i.e. LDS(i) == 0) the loop skips to the next iteration
			else
				CurrentAge = Age(ind)															! <<< selects current age of the animal
				CurrentAge = CurrentAge + 1														!ping age
				Age(ind) = CurrentAge															! >>> updates age in storage array
				CurrentDevStage = DevStage(ind)													! <<< selects current developmental stage
				CurrentDiapState = DiapState(ind)													! <<< selects the active/ready/diapause state
				CurrentCWeight = CWeight(ind)														! <<< selects structural mass
				
				devstagecount = log_devstage(CurrentPing, CurrentDevStage) + 1								!counts the developmental stages in the population 
				log_devstage(CurrentPing, CurrentDevStage) = devstagecount									!updates developmental stage count to the tracking variable
		
				!4.1.4.1 depth selection
				if(CurrentDiapState == "A" .or. CurrentDiapState == "T")then								!depth selection for active individuals pre (A) and post (T) diapause
				
					if(CurrentDevStage <= 3)then													!depth selection for egg, NI and NII
						call randomnum_int_srange(1, MLD(CurrentPing), CurrentDepth)						!upto NIII, the depth is dependent on the turbulent mixing (*** 50 m max ***)
						Depth(ind) = CurrentDepth												! >>> update vertical position in storage list
					else 																	!depth selection for NIII onwards
						MaxDistance = nint(5.22866 * CurrentCWeight ** 0.48616)							!maximum distance that can be searched by an individual
						CurrentDepth = Depth(ind)												! <<< presently occupied depth
					
						if(MaxDistance < MinDepth)then											!lower limit correction in case of rounding errors generating negative vals.
							MaxDistance = 1
						end if
					
						if(CurrentDepth + MaxDistance > MaxDepth)then									!upper limit correction to avoid overshooting maximum depth
							MFloor = MaxDepth	
						else	
							MFloor = CurrentDepth + MaxDistance
						end if
						
					
						if(CurrentDepth - MaxDistance < MinDepth)then									!lower limit correction to avoid undershooting minimum depth
							MCeiling = MinDepth
						else
							MCeiling = CurrentDepth - MaxDistance
						end if
					
						SDistance = MFloor - MCeiling	+ 1											!the search distance; i.e. the dimension of the search parameter arrays
					
						allocate(CurrentTempRange(SDistance))										!allocatable arrays for depth computation :: length = sdistance ::temperature
						allocate(CurrentFconRange(SDistance))										! :: food concentration
						allocate(CurrentVPRiskRange(SDistance))										! :: est.visual predation risk
						allocate(PGP(SDistance))												! :: percieved growth potential
						allocate(PVPR(SDistance))												! :: percieved visual predation risk
						allocate(PVisit(SDistance))												! :: probability of visiting
					
						CurrentTempRange = Temp(CurrentPing, MCeiling:MFloor)								! <<< extract data :: Temperature range
						CurrentTempRange = CurrentTempRange + 273.15									!convert to Kelvin to avoid zero-issues
						CurrentFconRange = Fcon(CurrentPing, MCeiling:MFloor)								! <<< extract data :: Food concentration range
						CurrentVPRiskRange = VPR(CurrentPing, MCeiling:MFloor)							! <<< extract data :: Vis.pred.risk range
						PGP = (CurrentTempRange * CurrentFconRange) / MaxPG								!percived growth potential
						PVPR = CurrentVPRiskRange * (1.00 / (1.00 + exp((125.00 - CurrentCWeight) / 25.00))) * VPRScalar
						call random_seed()
						call random_number(RN)													!probability checker for PVisit
					
						do i = 1, size(PVPR)													!evaluate visiting probability of all cells
					
							if(RN >= PVPR(i))then												!visit only if the PVPR is lesser than a uniform rand.no.
								PVisit(i) = PGP(i) 											!2nd eval. using perceived g.p :: because PGP(i) x 1.00 = PGP(i)
							else
								PVisit(i) = 0.00												!zero probability of visiting if uniform ran.no. < PVPR
							end if
						
						end do
						
						ML = maxloc(PVisit)													!saves the maxloc() output (i.e. 1-element array) to placeholder
						CurrentDepth = ML(1) + MCeiling											!visit the cell with max.g.potential
						VShift = abs(CurrentDepth - Depth(ind))										! <<< estimates the distance travelled
						Depth(ind) = CurrentDepth												! >>> update vertical position in storage array
					
						deallocate(CurrentTempRange)												!deallocate depth computation variables
						deallocate(CurrentFconRange)
						deallocate(CurrentVPRiskRange)
						deallocate(PGP)
						deallocate(PVPR)
						deallocate(PVisit)
					
					end if

				elseif(CurrentDiapState == "D")then												!depth selection for diapausing individuals
					CurrentDiapDepth = DiapDepth(ind)												! <<< selects randomly pre-selected diapause depth
					call randomnum_int_srange(CurrentDiapDepth - 25, CurrentDiapDepth + 25, CurrentDepth)			!random movement within the 'diapausing habitat' (*** +/- 25 m ***)
					Depth(ind) = CurrentDepth													! >>> update vertical position in storage array
	
				else																		!depth selection for diapause-ready individuals ("R")
					MaxDistance = nint(5.22866 * CurrentCWeight ** 0.48616)								!maximum distance that can be searched by an individual
					CurrentDepth = Depth(ind)													! <<< presently occupied depth
					CurrentDiapDepth = DiapDepth(ind)												! <<< randomly pre-selected diapause depth
					
					if(MaxDistance < MinDepth)then												!lower limit correction in case of rounding errors generating < 0 vals.
						MaxDistance = 1	
					end if
					
					if(CurrentDepth + MaxDistance > MaxDepth)then										!upper limit correction to avoid overshooting maximum depth
						MFloor = MaxDepth	
					else
						MFloor = CurrentDepth + MaxDistance	
					end if
					
					if(CurrentDepth - MaxDistance < MinDepth)then										!lower limit correction to avoid undershooting mimimum depth
						MCeiling = MinDepth	
					else
						MCeiling = CurrentDepth - MaxDistance
					end if
					
					if(CurrentDepth == CurrentDiapDepth)then											!depth of the diapause-ready individual equates the diapause depth 
						CurrentDepth = CurrentDiapDepth											!animal is now at prescribed diapause depth
						VShift = abs(CurrentDepth - Depth(ind))										! <<< estimates the distance travelled
						Depth(ind) = CurrentDepth												! >>> update vertical position in storage array
						CurrentDiapState = "D"													!switch from diapause-ready mode ("R") to diapause mode ("D")
						DiapState(ind) = CurrentDiapState											! >>> update diapause state in the storage array
					elseif(CurrentDepth > CurrentDiapDepth)then										!depth of diapause-ready individual exceed the prescribed diapause depth
						
						if(MCeiling > CurrentDiapDepth)then											!ceiling of migration lies below diapause depth
							CurrentDepth = MCeiling
							VShift = abs(CurrentDepth - Depth(ind))									! <<< estimates the distance travelled
							Depth(ind) = CurrentDepth											! >>> update vertical position in storage array
						else																!ceiling of migration depth lies above or equal diapause depth
							CurrentDepth = CurrentDiapDepth										!animal is now at prescribed diapause depth
							VShift = abs(CurrentDepth - Depth(ind))									! <<< estimates the distance travelled
							Depth(ind) = CurrentDepth											! >>> update vertical position in storage array
							CurrentDiapState = "D"												!switch from diapause-ready mode ("R") to diapause mode ("D")
							DiapState(ind) = CurrentDiapState										! >>> update diapause state in the storage array	
						end if
					
					else																	!depth of diapause-ready individual is shallower than the diapause depth
						
						if(MFloor < CurrentDiapDepth)then											!floor of migration lies above diapause depth 
							CurrentDepth = MFloor
							VShift = abs(CurrentDepth - Depth(ind))									! <<< estimates the distance travelled
							Depth(ind) = CurrentDepth											!update tracking variable Depth
						else																! >>> update vertical position in storage array
							CurrentDepth = CurrentDiapDepth
							VShift = abs(CurrentDepth - Depth(ind))									! <<< estimates the distance travelled
							Depth(ind) = CurrentDepth											! >>> update vertical position in storage array
							CurrentDiapState = "D"												!switch from diapause-ready mode ("R") to diapause mode ("D")
							DiapState(ind) = CurrentDiapState										! >>> update diapause state in the storage array
						end if
					
					end if
				
				end if
				
				CurrentVPos = VPos(CurrentPing, CurrentDepth) + 1										!extract and increment the population size of current depth at current time
				VPos(CurrentPing, CurrentDepth) = CurrentVPos											! >>> update VPos in storage array
				
				!4.1.5 simulate growth & development, survival and reproduction			
				if(CurrentDevStage == 1)then														!simulates growth & development and survival of embryonic stages
					CurrentStageDuration = StageDuration(ind) + 1										! <<< selects and updates stage duration
					CurrentTemp = Temp(CurrentPing, CurrentDepth)										! <<< subsets and select depth- and time-specific temperature
					
					call simulation_E(CurrentTemp, CurrentCWeight, CurrentDevTime, CurrentBMRate)					!estimates the current development time and loss of mass (degrowth)
					
					DTProd = DTTrack(ind) * CurrentDevTime											! <<< calculate the product for geometric mean estimation
					MeanDevTime = DTProd ** (1.00 / CurrentStageDuration)									!geometric mean of the development time :: development history
					CurrentCWeight = CurrentCWeight - CurrentBMRate										!estimate current structural mass after loss	
					
					if(CurrentCWeight <= 0.00)then
						CurrentCWeight = 0.00													!structural masses cannot drop below this lower treshold
						CWeight(ind) = CurrentCWeight												! >>> update structural mass in the storage array
						CurrentStarvationRisk = 1.00												!mortality risk peaks due to starvation
					else
						CWeight(ind) = CurrentCWeight												! >>> update structural mass in the storage array
						CurrentStarvationRisk = 0.00												!mortality is zero due to non-starvation
					end if
					
					if(CurrentStageDuration >= MeanDevTime)then
						CurrentDevStage = CurrentDevStage + 1										!molts to N1 (stage = 2)
						DevStage(ind) = CurrentDevStage											! >>> update developmental stage in storage array
						StageDuration(ind) = 0													! >>> reset stage duration in storage array
						DTTrack(ind) = 1.00													! >>> reset dttrack in the storage array
					else
						StageDuration(ind) = CurrentStageDuration										! >>> updates stage duration in storage array
						DTTrack(ind) = DTProd													! >>> update dttrack in the storage array
					end if
					
				elseif(CurrentDevStage == 2 .or. CurrentDevStage == 3)then									!simulates growth & development and survival of non-feeding nauplii
					CurrentStageDuration = StageDuration(ind) + 1										! <<< selects and updates stage duration
					CurrentTemp = Temp(CurrentPing, CurrentDepth)										! <<< subsets and select depth- and time-specific temperature
					
					call simulation_NFS(CurrentDevStage, CurrentTemp, CurrentCWeight, CurrentDevTime, CurrentBMRate)	!estimates the current development time and loss of mass (degrowth)
					
					DTProd = DTTrack(ind) * CurrentDevTime											! <<< calculate the product for geometric mean estimation
					MeanDevTime = DTProd ** (1.00 / CurrentStageDuration)									!geometric mean of the development time :: development history
					CurrentCWeight = CurrentCWeight - CurrentBMRate										!estimate current structural mass after loss
					
					if(CurrentCWeight <= 0.00)then												!structural masses cannot drop below this lower treshold
						CurrentCWeight = 0.00													
						CWeight(ind) = CurrentCWeight												! >>> update structural mass in the storage array
						CurrentStarvationRisk = 1.00												!mortality risk peaks due to starvation
					else
						CWeight(ind) = CurrentCWeight												! >>> update structural mass in the storage array
						CurrentStarvationRisk = 0.00												!no mortality risk due to lack of starvation
					end if
					
					if(CurrentStageDuration >= MeanDevTime)then
						CurrentDevStage = CurrentDevStage + 1										!molts to N1 (stage = 2) or NII (stage = 3)
						DevStage(ind) = CurrentDevStage											! >>> update developmental stage in storage array
						StageDuration(ind) = 0													! >>> reset stage duration in storage array
						DTTrack = 1.00														! >>> reset dttrack in the storage array
					else
						StageDuration(ind) = CurrentStageDuration										! >>> updates stage duration in storage array	
						DTTrack(ind) = DTProd													! >>> updates dttrack in storage array
					end if
					
				elseif(CurrentDevStage > 3 .and. CurrentDevStage < 11)then									!simulates growth & development and survival of feeding stages (NIII-CIII)	
					CurrentTemp = Temp(CurrentPing, CurrentDepth)										! <<< subsets and selects depth- and time-specific temperature
					CurrentFcon = Fcon(CurrentPing, CurrentDepth)										! <<< subsets and selects depth- and time-specific food concentration
					CurrentMaxCWeight = MaxCWeight(ind)												! <<< selects the maximum structural mass achieved
					CurrentBSP = BSP(ind)														! <<< selects the value for body-size parameter for current individual
								
					call simulation_FS(CurrentTemp, CurrentFcon, CurrentCWeight, CurrentGGRate, CurrentBMRate)		!processes growth of feeding and non-storing stages
					
					CurrentAMRate = (real(VShift) / real(MaxDistance)) * 1.5 * CurrentBMRate					!estimate active metabolic rate, i.e. swimming cost
					CurrentNGRate = (AsimCoef * CurrentGGRate) - (CurrentBMRate + CurrentAMRate)					!net growth potential
					
					if(CurrentNGRate >= 0.00)then													!if there is surplus growth; do the following
						
						CurrentCWeight = CurrentCWeight + CurrentNGRate									!assign net growth potential to current str.mass
						CWeight(ind) = CurrentCWeight												! >>> update the structural mass in storage array
						
						if(CurrentCWeight > CurrentMaxCWeight)then
							CurrentMaxCWeight = CurrentCWeight										!sets a new maximum structural mass
							MaxCWeight(ind) = CurrentMaxCWeight										! >>> update the new maximum structural mass in storage array
						end if
						
						CurrentStarvationRisk = 0.00												!no mortality risk due to lack of starvation
					else
						CurrentCWeight = CurrentCWeight + CurrentNGRate									!CurrentNGRate is a negative value
						CWeight(ind) = CurrentCWeight												!>>> update the structural mass in storage array
						CWLoss = (CurrentMaxCWeight - CurrentCWeight)  / CurrentMaxCWeight					!estimate propotional loss of structural mass
						
						if(CWLoss <= 0.10)then
							CurrentStarvationRisk = 0.00											!tolerant to modest loss of structural mass
						elseif(CWLoss > 0.10 .and. CWLoss <= 0.50)then
							CurrentStarvationRisk = 0.01 / (1.00 + exp((0.30 - CWLoss) / 0.05))				!estimate mortality risk via starvation model
						else
							CurrentStarvationRisk = 1.00											!certain death
						end if
						
					end if
					
					CMM = WMin(CurrentDevStage) + (WMax(CurrentDevStage) - WMin(CurrentDevStage)) * CurrentBSP		!estimate stage-specific critical molting mass
					
					if(CurrentCWeight >= CMM)then
						CurrentDevStage = CurrentDevStage + 1										!molts into the next developmental stage
						DevStage(ind) = CurrentDevStage											! >>> update developmental stage in storage array
					end if
					
				elseif(CurrentDevStage == 11 .or. CurrentDevStage == 12)then								!simulate growth & development and survival of potential overwintering stages
					
					if(CurrentDiapState == "A")then												!if the animal is in pre-diapause (1-year) or no-diapause (0-year) state
						CurrentLCL = LCL(ind)													! <<< selects the life cycle length
						CurrentTemp = Temp(CurrentPing, CurrentDepth)									! <<< subsets and selects depth and time specific temperature
						CurrentFcon = Fcon(CurrentPing, CurrentDepth)									! <<< subsets and selects depth and time specific food concentration
						CurrentMaxCWeight = MaxCWeight(ind)											! <<< selects maximum structural mass attained
						CurrentSWeight = SWeight(ind)												! <<< selects the energy reserve mass
						CurrentGAP = GAP(ind)													! <<< selects the growth allocation parameter
						CurrentSDP = SDP(ind)													! <<< selects the seasonal descent parameter
						CurrentBSP = BSP(ind)													! <<< selects the body size parameter
						
						call simulation_PDS(CurrentTemp, CurrentFcon, CurrentCWeight, CurrentSWeight, CurrentGGRate, CurrentBMRate)
						
						CurrentAMRate = (real(VShift) / real(MaxDistance)) * 1.5 * CurrentBMRate				!estimate active metabolic rate, i.e. swimming cost
						CurrentNGRate = (AsimCoef * CurrentGGRate) - (CurrentBMRate + CurrentAMRate)				!net growth rate
						
						if(CurrentNGRate > 0.00)then												!when there is surplus net growth
							
							if(CurrentLCL == 0)then												!no limitation to structural growth for 0-yr group :: direct development
							
								if(CurrentSWeight / CurrentCWeight >= 1.00)then							!when reserves are full;
									CurrentCWeight = CurrentCWeight + CurrentNGRate						!all assimilated carbon to build-up str. mass
								else														!when the reserves are not full;
									CurrentSWeight = CurrentSWeight + (CurrentNGRate * CurrentGAP)			!allocate individual-specific amount of surplus aquisition to reserve
									CurrentCWeight = CurrentCWeight + (CurrentNGRate * (1.00 - CurrentGAP))		!rest (if any :: where GAP = 1.00) to strutural mass 
								end if				
							
							else															!limits the structural growth at the onset of female size :: no direct devel. 
								CWeightCeiling = WMin(12) + (WMax(12) - WMin(12)) * CurrentBSP				!estimated female size (j = 13; in array = 13 - 1)
								
								if(CurrentCWeight >= CWeightCeiling)then								!when animal tries to grow beyond the estimated size of female
									CurrentSWeight = CurrentSWeight + CurrentNGRate						!entire surplus aquisition to store build-up
								else														!when animal is below the estimated size of female
									CurrentSWeight = CurrentSWeight + (CurrentNGRate * CurrentGAP)			!allocate ind.-specific amount of surplus assim. to reserve
									CurrentCWeight = CurrentCWeight + (CurrentNGRate * (1.00 - CurrentGAP))		!rest (if any :: where GAP = 1.00) to str. mass
								end if
								
							end if
							
							if(CurrentCWeight > CurrentMaxCWeight)then								!only when a new height of structural mass is reached;
								CurrentMaxCWeight = CurrentCWeight									!that maximum structural mass is incremented
								MaxCWeight(ind) = CurrentMaxCWeight									! >>> update maximum structural mass in storage array
							end if	

							CurrentStarvationRisk = 0.00											!no starvation risk
							
							CWeight(ind) = CurrentCWeight											! >>> update structural mass in storage array
							SWeight(ind) = CurrentSWeight											! >>> update reserve mass in storage array
							
						elseif(CurrentNGRate == 0.00)then											!when there is zero net growth	
							CurrentStarvationRisk = 0.00											!no starvation risk
						else																!when there is negative net growth;
							
							if(CurrentSWeight >= abs(CurrentNGRate))then								!when reserves can balance the negative growth;	
								CurrentSWeight = CurrentSWeight + CurrentNGRate							!reserves mobilized :: CurrentNGRate is a negative value
								CurrentStarvationRisk = 0.00										!no mortality risk due to no starvation
								SWeight(ind) = CurrentSWeight										! >>> update reserve mass in storage array	
							else															!when reserves cannot balance the potential degrowth;
								CurrentCWeight = CurrentCWeight + CurrentSWeight + CurrentNGRate	 			!mobilize reserves if any :: CurrentNGRate is negative
								CurrentSWeight = 0.00											!no reserves after mobilization of any available reserves
								CWLoss = (CurrentMaxCWeight - CurrentCWeight)  / CurrentMaxCWeight			!estimated loss of carbon as a fraction to max. str. mass
								
								if(CWLoss <= 0.10)then
									CurrentStarvationRisk = 0.00									!tolerant to modest loss of structural mass
								elseif(CWLoss > 0.10 .and. CWLoss <= 0.50)then
									CurrentStarvationRisk = 0.01 / (1.00 + exp((0.30 - CWLoss) / 0.05))		!estimate mortality risk via starvation model
								else
									CurrentStarvationRisk = 1.00									!certain death
								end if

								CWeight(ind) = CurrentCWeight										! >>> update structural mass in storage array
								SWeight(ind) = CurrentSWeight										! >>> update reserve mass in storage array
							end if
		
						end if	
						
						CMM = WMin(CurrentDevStage) + (WMax(CurrentDevStage) - WMin(CurrentDevStage)) * CurrentBSP	!estimated stage-specific critical molting mass
							
						if(CurrentLCL == 1)then													!when the animal has a 1-year life cycle;
							
							if(CurrentDevStage == 11)then											!when the animal is in CIV stage;
	
								if(CurrentCWeight >= CMM)then										!when current structural mass exceeds the critical molting mass;
									CurrentDevStage = CurrentDevStage + 1							!molts to CV
									DevStage(ind) = CurrentDevStage								! >>> update developmental stage in storage array
								elseif(CurrentSWeight / CurrentCWeight >= CurrentSDP)then					!when the animal meets the diapause condition
									CurrentDiapState = "R"										!animal is now ready to enter diapause
									DiapState(ind) = CurrentDiapState								! >>> update diapause state in storage array
									SWSD(ind) = CurrentSWeight									! >>> update reserve mass at seasonal descent in storage array
									CWSD(ind) = CurrentCWeight									! >>> update structural mass at seasonal descent in storage array
									TSD(ind) = CurrentPing										! >>> update time of seasonal descent in storage array
									n_sd = n_sd + 1											!count the seasonal descents
									sum_tsd = sum_tsd + CurrentPing								!aggregate seasonal descent timings
									sum_cwsd = sum_cwsd + CurrentCWeight							!aggregate structural mass at seasonal descent
									sum_swsd = sum_swsd + CurrentSWeight							!aggregate reserve mass at seasonal descent
								end if
							
							else															!when the animal is in CV stage;
	
								if(CurrentSWeight / CurrentCWeight >= CurrentSDP)then						!when the animal meets the diapause condition
									CurrentDiapState = "R"										!animal is now ready to enter diapause
									DiapState(ind) = CurrentDiapState								! >>> update diapause state in storage array
									SWSD(ind) = CurrentSWeight									! >>> update reserve mass at seasonal descent in storage array
									CWSD(ind) = CurrentCWeight									! >>> update structural mass at seasonal descent in storage array
									TSD(ind) = CurrentPing										! >>> update time of seasonal descent in storage array
									n_sd = n_sd + 1											!count the seasonal descents
									sum_tsd = sum_tsd + CurrentPing								!aggregate seasonal descent timings
									sum_cwsd = sum_cwsd + CurrentCWeight							!aggregate structural mass at seasonal descent
									sum_swsd = sum_swsd + CurrentSWeight							!aggregate reserve mass at seasonal descent
								end if
									
							end if
						else																!when the animal has a 0-year life cycle
						
							if(CurrentCWeight >= CMM)then
								CurrentDevStage = CurrentDevStage + 1								!molts to CV
								DevStage(ind) = CurrentDevStage									! >>> update developmental stage in storage array
								
								if(CurrentDevStage == 13)then
									n_sm = n_sm + 1											!count the number reaching sexual maturity
									sum_cwsm = sum_cwsm + CurrentCWeight							!aggregate the structural mass of sexually mature adults
								end if
								
							end if
							
						end if
						
					elseif(CurrentDiapState == "R")then											!if the animal is in a diapause-ready state;
						CurrentTemp = Temp(CurrentPing, CurrentDepth)									! <<< subset and select time and depth specific temperature
						CurrentMaxCWeight = MaxCWeight(ind)											! <<< select maximum structutal mass attained
						CurrentSWeight = SWeight(ind)												! <<< selects energy reserve mass
						
						call simulation_DRS(CurrentTemp, CurrentCWeight, CurrentSWeight, CurrentBMRate)			!estimates degrowth of diapause-ready stages
						
						CurrentAMRate = (real(VShift) / real(MaxDistance)) * 1.5 * CurrentBMRate				!estimate active metabolic rate, i.e. swimming cost
						
						if(CurrentSWeight >= CurrentAMRate + CurrentBMRate)then							!when energy reserves are sufficient to balance metabolic losses
							CurrentSWeight = CurrentSWeight - (CurrentAMRate + CurrentBMRate)					!reserves are mobilized to balance metabolic demands
							CurrentStarvationRisk = 0.00											!no mortality risk due to no starvation
							SWeight(ind) = CurrentSWeight											! >>> update reserve mass in storage array
						else																!when energy reserves not sufficient to balance metabolic demands
							CurrentCWeight = (CurrentCWeight + CurrentSWeight) - (CurrentAMRate + CurrentBMRate)	!structural mass is mobilized to meet metabolic demands
							CurrentSWeight = 0.00												!whatever the remaining stores are mobilized to balance metabolic demands
							CWLoss = (CurrentMaxCWeight - CurrentCWeight)  / CurrentMaxCWeight				!estimated loss of carbon as a fraction to max. str. mass
							
							if(CWLoss <= 0.10)then
								CurrentStarvationRisk = 0.00										!tolerant to modest loss of structural mass
							elseif(CWLoss > 0.10 .and. CWLoss <= 0.50)then
								CurrentStarvationRisk = 0.01 / (1.00 + exp((0.30 - CWLoss) / 0.05))			!estimate mortality risk via starvation model
							else
								CurrentStarvationRisk = 1.00										!certain death
							end if
							
							CWeight(ind) = CurrentCWeight											! >>> update structural mass in storage array
							SWeight(ind) = CurrentSWeight											! >>> update reserve mass in storage array
						end if
						
					elseif(CurrentDiapState == "D")then											!if the animal is in a diapause state;
						CurrentTemp = Temp(CurrentPing, CurrentDepth)									! <<< subset and select time and depth specific temperature
						CurrentMaxCWeight = MaxCWeight(ind)											! <<< select maximum structutal mass attained
						CurrentSWeight = SWeight(ind)												! <<< selects energy reserve mass
						CurrentSAP = SAP(ind)													! <<< selects seasonal ascent parameter
						
						call simulation_DS(CurrentTemp, CurrentCWeight, CurrentSWeight, CurrentBMRate)			!estimates degrowth of stages in diapause ("D")
						
						CurrentSWeight = CurrentSWeight - CurrentBMRate									!metabolic demands are balanced by energy reserves
						
						if(CurrentSWeight < 0.00)then												!in an extremely rare case if one time ping causes negative reserves;
							CurrentSWeight = 0.00												!reserves are reset to zero
						end if
						
						SWeight(ind) = CurrentSWeight												! >>> update energy reserve mass in storage array
						CurrentStarvationRisk = 0.00												!no mortality risk in diapause due to no starvation
						
						if(CurrentSWeight <= SWSD(ind) * CurrentSAP)then								!when the animal meets the diapause termination condition;
							CurrentDiapState = "T"												!the diapause is now terminated :: animal is now active ("T"-state)
							DiapState(ind) = CurrentDiapState										! >>> update diapause state in storage array
							TSA(ind) = CurrentPing												! >>> update time of seasonal ascent in storage array
							n_sa = n_sa + 1													!count the seasonal descents
							sum_tsa = sum_tsa + CurrentPing										!aggregate seasonal descent timings			
						end if
						
					else																	!if the animal is in the diapause terminated state;
						CurrentTemp = Temp(CurrentPing, CurrentDepth)									! <<< subsets and selects depth and time specific temperature
						CurrentFcon = Fcon(CurrentPing, CurrentDepth)									! <<< subsets and selects depth and time specific food concentration
						CurrentMaxCWeight = MaxCWeight(ind)											! <<< selects maximum structural mass attained
						CurrentSWeight = SWeight(ind)												! <<< selects the energy reserve mass
						CurrentBSP = BSP(ind)													! <<< selects the body size parameter
						
						call simulation_PDS(CurrentTemp, CurrentFcon, CurrentCWeight, CurrentSWeight, CurrentGGRate, CurrentBMRate)
						
						CurrentAMRate = (real(VShift) / real(MaxDistance)) * 1.5 * CurrentBMRate				!estimate active metabolic rate, i.e. swimming cost
						CurrentNGRate = (AsimCoef * CurrentGGRate) - (CurrentBMRate + CurrentAMRate)				!net growth rate
						
						if(CurrentNGRate >= 0.00)then												!when there is a net surplus growth or zero growth;
							CurrentCWeight = CurrentCWeight + CurrentNGRate								!all surplus aquisition is channeled to structural growtgh
							CurrentStarvationRisk = 0.00											!no mortality risk due to no starvation
							
							if(CurrentCWeight > CurrentMaxCWeight)then								!only when a new height of structural mass is reached;
								CurrentMaxCWeight = CurrentCWeight									!that maximum structural mass is incremented
								MaxCWeight(ind) = CurrentMaxCWeight									! >>> update maximum structural mass in storage array
							end if
								
							CWeight(ind) = CurrentCWeight											! >>> update structural mass in storage array
						else
							
							if(CurrentSWeight >= abs(CurrentNGRate))then								!when reserves can balance the metabolic demands;
								CurrentSWeight = CurrentSWeight + CurrentNGRate							!reserves are mobilized to meet the metabolic costs
								CurrentStarvationRisk = 0.00										!no mortality risk due to no starvation
								
								SWeight(ind) = CurrentSWeight										! >>> update reserve mass in storage array
							else
								CurrentCWeight = CurrentCWeight + CurrentSWeight + CurrentNGRate				!structural mass mobilized to meet metabolic demands :: CurrentNGRate is -ve
								CurrentSWeight = 0.00											!reserves are fully spent
								CWLoss = (CurrentMaxCWeight - CurrentCWeight)  / CurrentMaxCWeight			!estimated loss of carbon as a fraction to max. str. mass
								
								if(CWLoss <= 0.10)then
									CurrentStarvationRisk = 0.00									!tolerant to modest loss of structural mass
								elseif(CWLoss > 0.10 .and. CWLoss <= 0.50)then
									CurrentStarvationRisk = 0.01 / (1.00 + exp((0.30 - CWLoss) / 0.05))		!estimate mortality risk via starvation model
								else
									CurrentStarvationRisk = 1.00									!certain death
								end if
								
								CWeight(ind) = CurrentCWeight										! >>> update structural mass in storage array
								SWeight(ind) = CurrentSWeight										! >>> update reserve mass in storage array
							end if
							
						end if
						
						CMM = WMin(CurrentDevStage) + (WMax(CurrentDevStage) - WMin(CurrentDevStage)) * CurrentBSP	!estimated stage-specific critical molting mass
						
						if(CurrentCWeight >= CMM)then
							CurrentDevStage = CurrentDevStage + 1									!molts to CV or Adult
							DevStage(ind) = CurrentDevStage										! >>> update developmental stage in storage array
						end if
						
					end if			
						
				else																		!simulates growth, development and reproduction of adult stages		
					CurrentTemp = Temp(CurrentPing, CurrentDepth)										! <<< subsets and selects depth and time specific temperature
					CurrentFcon = Fcon(CurrentPing, CurrentDepth)										! <<< subsets and selects depth and time specific food concentration
					CurrentMaxCWeight = MaxCWeight(ind)												! <<< selects maximum structural mass attained
					CurrentRWeight = RWeight(ind)													! <<< selects the allocated mass for reproductive development
					CurrentSWeight = SWeight(ind)													! <<< selects the mass of energy reserves
					CurrentISState = ISState(ind)													! <<< selects the insemination state
					CurrentFecundity = Fecundity(ind)												! <<< selects the no. of eggs produced
					
					if(Gender(ind) == "M")then
						NM  = NM + 1														!add to a male tracker
					elseif(Gender(ind) == "F")then
						NF = NF + 1															!add to a female tracker
					end if
						
					call simulation_PDS(CurrentTemp, CurrentFcon, CurrentCWeight, CurrentSWeight, CurrentGGRate, CurrentBMRate)
						
					CurrentAMRate = (real(VShift) / real(MaxDistance)) * 1.5 * CurrentBMRate					!estimate active metabolic rate, i.e. swimming cost
					CurrentNGRate = (AsimCoef * CurrentGGRate) - (CurrentBMRate + CurrentAMRate)					!net growth rate
						
					if(CurrentNGRate >= 0.00)then													!if there is surplus growth;
					
						if(Gender(ind) == "F")then												!if the individual is a female;
							
							if(CurrentISState == "N")then											!if the female is not inseminated;
								NISF = NISF + 1												!add to the count of number of non-inseminated females
								IS_F(NISF) = ind												! >>> update individual ID in non-inseminated male list
							else															!if the female is inseminated;
								CurrentRWeight = CurrentRWeight + CurrentNGRate							!all surplus aquisition channeled to egg production
								
								if(CurrentRWeight >= 0.23)then									!when the surplus reproductive allocation is greater than unit egg mass;
									NEPF = NEPF + 1											!counts the number of egg-producing females
									CurrentNEggs = nint((CurrentRWeight) / 0.23)						!eggs are produced
									CurrentFecundity = CurrentFecundity + CurrentNEggs					!no. of eggs produced by the current time ping
									CurrentRWeight = CurrentRWeight - (0.23 * CurrentNEggs)				!remainder of reproductive allocation
									Fecundity(ind) = CurrentFecundity								! >>> update egg production estimate in storage array
									RWeight(ind) = CurrentRWeight									! >>> update reproductive allocation in storage array
									
									MID = SelectedMate(ind)										!selected mate id
									Pairing(1, NEPF) = GAP(ind)									!set GAP 'aleele' from mother
									Pairing(2, NEPF) = SDP(ind)									!set SDP 'aleele' from mother
									Pairing(3, NEPF) = SAP(ind)									!set SAP 'aleele' from mother
									Pairing(4, NEPF) = BSP(ind)									!set BSP 'aleele' from mother
									Pairing(5, NEPF) = GAP(MID)									!set GAP 'aleele' from father
									Pairing(6, NEPF) = SDP(MID)									!set SDP 'aleele' from father
									Pairing(7, NEPF) = SAP(MID)									!set SAP 'aleele' from father
									Pairing(8, NEPF) = BSP(MID)									!set BSP 'aleele' from father
									NEggs(NEPF) = CurrentNEggs									!set no. of eggs produced (align with Pairing)
								else
									RWeight(ind) = CurrentRWeight									! >>> update reproductive allocation in storage array
								end if
								
							end if
							
						else																!if the individual is a male
							NISM = NISM + 1													!add to the count of number of inseminate-ready males
							IS_M(NISM) = ind													! >>> update individual ID in inseminate-ready male list
						end if
						
						CurrentStarvationRisk = 0.00												!no mortality risk due to lack of starvation
						
					else																	!if there is no surplus net growth;
					
						if(CurrentSWeight >= abs(CurrentNGRate))then									!if reserves can balance the metabolic demands
							CurrentSWeight = CurrentSWeight + CurrentNGRate								!energy demands balanced by reserves :: CurrentNGRate is a negative val.
							CurrentStarvationRisk = 0.00											!no mortality risk due to no starvation
							SWeight(ind) = CurrentSWeight											! >>> update reserve mass in storage array
						else																!if reserves cannot balance metabolic demands
							CurrentCWeight = CurrentCWeight + CurrentSWeight + CurrentNGRate					!structural mass mobilized to meet metabolic demands
							CurrentSWeight = 0.00												!mobilize whatever the reserves available to counter starvation
							CWLoss = (CurrentMaxCWeight - CurrentCWeight)  / CurrentMaxCWeight				!estimated loss of carbon as a fraction to max. str. mass
							
							if(CWLoss <= 0.10)then
								CurrentStarvationRisk = 0.00										!tolerant to modest loss of structural mass
							elseif(CWLoss > 0.10 .and. CWLoss <= 0.50)then
								CurrentStarvationRisk = 0.01 / (1.00 + exp((0.30 - CWLoss) / 0.05))			!estimate mortality risk via starvation model
							else
								CurrentStarvationRisk = 1.00										!certain death
							end if
							
							CWeight(ind) = CurrentCWeight											! >>> update structural mass in storage array
							SWeight(ind) = CurrentSWeight											! >>> update reserve mass in storage array
						end if
						
						if(Gender(ind) == "M")then												!if the individual is a male;
							NISM = NISM + 1													!add to the count of number of inseminate-ready males
							IS_M(NISM) = ind													! >>> update individual ID in inseminate-ready male list
						else																!if the individual is a female;
							
							if(CurrentISState == "N")then											!if the female is not inseminated;
								NISF = NISF + 1												!add to the count of number of non-inseminated females
								IS_F(NISF) = ind												! >>> update individual ID in non-inseminated female list
							end if
							
						end if
							
					end if
											
				end if
				
				!4.1.6 mortality risk and population management
				CurrentNVPRisk = 0.00001														!fixed non-visual predation risk (10% of max.visual pred.risk)
				CurrentVPRisk = VPR(CurrentPing, CurrentDepth) * (1.00 / (1.00 + exp((125.00 - CurrentCWeight) / 25.00))) * VPRScalar	!size-, depth-, time- and scale-dependent visual predation risk
				CurrentDDMRisk = 0.01 / (1.00 + exp((0.85 * NMax - popsize) / (0.025 * NMax)))					!estimated density dependent mortality :: no size-dependency			
				CurrentMortalityRisk = CurrentNVPRisk + CurrentVPRisk + CurrentStarvationRisk + CurrentDDMRisk			!sum of all sources of mortality			
				
				call random_seed()
				call random_number(SurvivalProb)													!survival probability from time t to t+1
				
				if(SurvivalProb < CurrentMortalityRisk .or. CurrentFecundity > FecundityCeiling .or. CurrentAge > 17520)then
					LDS(ind) = 0															!LDS set to 0(death) :: all other placeholders reset at seeding/spawning
				end if
								
			end if
			
		end do
		
		!4.2 mate selection and spawning elements
		if(NISF > 0 .and. NISM > 0)then															!if there are either 1 or more ins.-ready males and non-inseminated females
			
			do i = 1, NISF																	!do the following for each inseminate-ready female
				call randomnum_int_srange(1, NISM, RNI)												!randomly select a mate :: *** no encounter rates ***
				FID = IS_F(i)																!individual ID of the female
				MID = IS_M(RNI)																!individual ID of the selected male
				SelectedMate(FID) = MID															! >>> update selected mate in storage array		
				ISState(FID) = "I"															!change the insemination state of the mated female
			end do
			
		end if
		
		CurrentSpawnSize = sum(NEggs)																!the number of eggs to spawn at the onset of next time ping
		
		!4.3 I/O communications
		if(any(months == CurrentPing))then		
			progress = (real(CurrentYear)/real(EndYear)) * 100											!completed iterations as a %
			!popsize = sum(VPos(CurrentPing, :))													!estimates the current population size :: method-1
			!popsize = sum(LDS)																!same as above :: method-2
			print "(a2, i3, a3, 1X, a18, i8)", "[ ", progress, "% ]", "population size : ", popsize					!prints progress and population size
		end if
		
		!4.4 data logging
		log_neggs(CurrentPing) = CurrentSpawnSize														!log the current egg production
		!log_popsize(CurrentPing) = popsize															!log the current population size
		log_nmales(CurrentPing) = NM																!log the current no. of males
		log_nfemales(CurrentPing) = NF															!log the current no. of females
		
		if(CurrentPing == 8760)then																!at the end of each year;
			!4.4.1 I/0 type - 1 :: 2D-hourly vertical position log 
			write(filename, "(A5, I3)") "vpos_", CurrentYear											!initialize filename
			output = directory//filename//extension													!concatenate
			open(unit = 1000, file = output)														!create empty file for output
			call csv_write_integer_2d(1000, VPos)													!write csv file
			close(1000)																		!close unit	
			
			!4.4.2 I/O type -2 :: 1D-hourly annual logs
			!write(filename, "(A5, I3)") "posz_", CurrentYear
			!output = directory//filename//extension
			!open(unit = 1002, file = output)
			!call csv_write(1002, log_popsize, .false.)
			!close(1002)
			
			write(filename, "(A5, I3)") "eprd_", CurrentYear
			output = directory//filename//extension
			open(unit = 1003, file = output)
			call csv_write(1003, log_neggs, .false.)
			close(1003)
			
			write(filename, "(A5, I3)") "nmle_", CurrentYear
			output = directory//filename//extension
			open(unit = 1004, file = output)
			call csv_write(1004, log_nmales, .false.)
			close(1004)
			
			write(filename, "(A5, I3)") "nfem_", CurrentYear
			output = directory//filename//extension
			open(unit = 1005, file = output)
			call csv_write(1005, log_nfemales, .false.)
			close(1005)
			
			write(filename, "(A5, I3)") "ndev_", CurrentYear
			output = directory//filename//extension
			open(unit = 1006, file = output)
			call csv_write_integer_2d(1006, log_devstage)
			
			!4.4.3 I/O type - 2 :: once in a lifetime event logs
			IO_sdsa(CurrentYear, 1) = real(sum_tsd) / real(n_sd)											!calculate and update mean annual timing of seasonal descent
			IO_sdsa(CurrentYear, 2) = sum_cwsd / real(n_sd)												!calculate and update mean annual structural mass of seasonal descent
			IO_sdsa(CurrentYear, 3) = sum_swsd / real(n_sd)												!calculate and update mean annual reserve mass of seasonal descent
			IO_sdsa(CurrentYear, 4) = real(sum_tsa) / real(n_sa)											!calculate and update mean annual timing of seasonal ascent
			IO_sdsa(CurrentYear, 5) = sum_cwsm / real(n_sm)												!calculate and update mean annual size of sexual maturity
			
			n_sd = 0; n_sa = 0; n_sm = 0; sum_tsd = 0; sum_tsa = 0; sum_cwsd = 0.00; sum_swsd = 0.00; sum_cwsm = 0.00		!reset values for next year's logging
			log_nmales = 0; log_nfemales = 0; log_popsize = 0; NM = 0; NF = 0; log_devstage = 0						!reset values for next iteration
		end if																				

		!4.5. termination
		if(CurrentYear >= EndYear)then															!breaking condition for main loop
			print "(a17)", "[WRITING TO DISC]"
			open(unit = 1001, file = "/home/kanchana/Modelling/M1/Outputs/AnnualLog.txt")
			call csv_write_dble_2d(1001, IO_sdsa)
			close(1001)		
			print "(a10)", "[COMPLETE]"
			
			exit
		end if
		 
	end do

end program IBSM1D_v1
