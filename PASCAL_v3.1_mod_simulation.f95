module simulation
	implicit none
	
	contains
		subroutine simulation_E(temperature, bodymass, devtime, estmrate)												!estimates dev.time and loss of mass :: embroyo
			implicit none
			
			real(kind = 8), intent(in) :: temperature, bodymass														!inputs = temperature & bodymass
			real(kind = 8) :: devtime, estmrate																!outputs = devtime & estmrate
			real(kind = 8) :: mrate_m, mrate_t																	!interim :: temp and mass specific metabolic rates
			integer, parameter :: dp1 = 595																	!development parameter-1 :: C.fin (Campbell et al. 2001)
			real(kind = 8), parameter :: dp2 = 9.11																!development parameter-2 :: C.fin (Campbell et al. 2001)
			real(kind = 8), parameter :: M_masscoef = 0.0008487														!mass coef. of metabolism (all stages)
			real(kind = 8), parameter :: M_massexpo = 0.7503														!mass exponent of metabolism
			real(kind = 8), parameter :: M_tempcoef = 1.2956														!temperature coefficient of metabolism
			real(kind = 8), parameter :: M_tempexpo = 0.1170														!temperature exponent of metabolism
			
			devtime = 24 * (dp1 * (dp2 + temperature) ** (-2.05))														!estimates development time in hours
			mrate_m = M_masscoef * bodymass ** M_massexpo															!mass-specific metabolic rate at -2C
			mrate_t = mrate_m * M_tempcoef * exp(M_tempexpo * temperature) * 0.5											!temperature-adjusted metabolic rate
			estmrate = mrate_t																			!assign and export value
			
		end subroutine simulation_E
		
		subroutine simulation_NFS(devstage, temperature, bodymass, devtime, estmrate)										!estimates dev.time and loss of mass :: non-feeding NI & NII
			implicit none
			
			real(kind = 8), intent(in) :: temperature, bodymass														!inputs = temperature & bodymass
			integer, intent(in) :: devstage																	!input = dev.stage
			real(kind = 8) :: devtime, estmrate																!outputs = devtime & estmrate
			real(kind = 8) :: mrate_m, mrate_t																	!interim :: temp and mass specific metabolic rates
			integer, parameter :: dp1_a = 388																	!development parameter-1 :: C.fin (Campbell et al. 2001)
			integer, parameter :: dp1_b = 581																	!development parameter-1 :: C.fin (Campbell et al. 2001)
				
			real(kind = 8), parameter :: dp2 = 9.11																!development parameter-2 :: C.fin (Campbell et al. 2001)
			real(kind = 8), parameter :: M_masscoef = 0.0008487														!mass coef. of metabolism (all stages)
			real(kind = 8), parameter :: M_massexpo = 0.7503														!mass exponent of metabolism
			real(kind = 8), parameter :: M_tempcoef = 1.2956														!temperature coefficient of metabolism
			real(kind = 8), parameter :: M_tempexpo = 0.1170														!temperature exponent of metabolism
			
			if(devstage == 2)then
				devtime = 24 * (dp1_a * (dp2 + temperature) ** (-2.05))												!estimates development time in hours :: NI	
			else
				devtime = 24 * (dp1_b * (dp2 + temperature) ** (-2.05))												!estimates development time in hours :: NII
			end if
			
			mrate_m = M_masscoef * bodymass ** M_massexpo															!mass-specific metabolic rate at -2C
			mrate_t = mrate_m * M_tempcoef * exp(M_tempexpo * temperature) * 0.5											!temperature-adjusted metabolic rate
			estmrate = mrate_t																			!assign and export value
			
		end subroutine simulation_NFS
		
		subroutine simulation_FS(temperature, foodcon, bodymass, estgrate, estmrate)											!estimates growth of feeding, non-storing NIII-NVI & CI-CIII
			implicit none
			
			real(kind = 8), intent(in) :: temperature, foodcon, bodymass												!inputs = temp., foodcon. and structural mass 
			real(kind = 8) :: estgrate, estmrate																!outputs = growth rate and metabolic rate 
			real(kind = 8), parameter :: FC = 30.00																!food conversion factor
			real(kind = 8), parameter :: G_masscoef = 0.009283														!mass coef. of assimilation (all stages)
			real(kind = 8), parameter :: G_massexpo = 0.7524														!mass expo. of assimilation
			real(kind = 8), parameter :: G_tempcoef = 1.2382														!temperature coef. of assimilation
			real(kind = 8), parameter :: G_tempexpo = 0.0966														!temperature expo. of assimilation
			real(kind = 8), parameter :: M_masscoef = 0.0008487														!mass coef. of metabolism (all stages)
			real(kind = 8), parameter :: M_massexpo = 0.7503														!mass exponent of metabolism
			real(kind = 8), parameter :: M_tempcoef = 1.2956														!temperature coefficient of metabolism
			real(kind = 8), parameter :: M_tempexpo = 0.1170														!temperature exponent of metabolism
			real(kind = 8) :: mrate_m, mrate_t																	!interim :: temp and mass specific metabolic rates
			real(kind = 8) :: GRate_mass, GRate_temp, GRate, gp1														!interim :: mass, temp and food-adjusted assimilation rates
			
			GRate_mass = G_masscoef * bodymass ** G_massexpo														!mass-adjusted assimilation rate
			GRate_temp = GRate_mass * G_tempcoef * exp(G_tempexpo * temperature)											!temperature-scaled assimilation rate
			gp1 = 0.3 * bodymass ** (-0.138)																	!scaling parameter :: parameter 'd' in Bandara et al. 2019
			GRate = (GRate_temp * gp1 * foodcon * FC) / (1 + gp1 * foodcon * FC)											!food-adjusted assimilation rate
			estgrate = GRate																				!assign and export value-1
			
			mrate_m = M_masscoef * bodymass ** M_massexpo															!mass-specific metabolic rate at -2C
			mrate_t = mrate_m * M_tempcoef * exp(M_tempexpo * temperature)												!temperature-adjusted metabolic rate
			estmrate = mrate_t																			!assign and export value-2
			
		end subroutine simulation_FS
		
		subroutine simulation_PDS(temperature, foodcon, bodymass, reservemass, estgrate, estmrate)								!estimates growth of feeding, storing CIV, CV
			implicit none
			
			real(kind = 8), intent(in) :: temperature, foodcon, bodymass, reservemass										!inputs = temp., foodcon. reserve and structural mass 
			real(kind = 8) :: estgrate, estmrate																!outputs = growth rate and metabolic rate 
			real(kind = 8), parameter :: FC = 30.00																!food conversion factor
			real(kind = 8), parameter :: G_masscoef = 0.009283														!mass coef. of assimilation (all stages)
			real(kind = 8), parameter :: G_massexpo = 0.7524														!mass expo. of assimilation
			real(kind = 8), parameter :: G_tempcoef = 1.2382														!temperature coef. of assimilation
			real(kind = 8), parameter :: G_tempexpo = 0.0966														!temperature expo. of assimilation
			real(kind = 8), parameter :: M_masscoef = 0.0008487														!mass coef. of metabolism (all stages)
			real(kind = 8), parameter :: M_massexpo = 0.7503														!mass exponent of metabolism
			real(kind = 8), parameter :: M_tempcoef = 1.2956														!temperature coefficient of metabolism
			real(kind = 8), parameter :: M_tempexpo = 0.1170														!temperature exponent of metabolism
			real(kind = 8) :: mrate_m, mrate_t, totalmass															!interim :: temp and mass specific metabolic rates
			real(kind = 8) :: GRate_mass, GRate_temp, GRate, gp1														!interim :: mass, temp and food-adjusted assimilation rates
			
			GRate_mass = G_masscoef * bodymass ** G_massexpo														!mass-adjusted assimilation rate
			GRate_temp = GRate_mass * G_tempcoef * exp(G_tempexpo * temperature)											!temperature-scaled assimilation rate
			gp1 = 0.3 * bodymass ** (-0.138)																	!scaling parameter :: parameter 'd' in Bandara et al. 2019
			GRate = (GRate_temp * gp1 * foodcon * FC) / (1 + gp1 * foodcon * FC)											!food-adjusted assimilation rate
			estgrate = GRate																				!assign and export value-1

			totalmass = bodymass + reservemass																	!estimated total mass
			mrate_m = M_masscoef * totalmass ** M_massexpo															!mass-specific metabolic rate at -2C
			mrate_t = mrate_m * M_tempcoef * exp(M_tempexpo * temperature)												!temperature-adjusted metabolic rate
			estmrate = mrate_t																			!assign and export value-2			
			
		end subroutine simulation_PDS
		
		subroutine simulation_DRS(temperature, bodymass, reservemass, estmrate)
			implicit none
			
			real(kind = 8), intent(in) :: temperature, bodymass, reservemass												!inputs = temp., reserve and structural mass
			real(kind = 8) :: estmrate																		!output = est.metabolic rate
			real(kind = 8) :: mrate_m, mrate_t, totalmass															!interim :: temp and mass specific metabolic rates
			real(kind = 8), parameter :: M_masscoef = 0.0008487														!mass coef. of metabolism (all stages)
			real(kind = 8), parameter :: M_massexpo = 0.7503														!mass exponent of metabolism
			real(kind = 8), parameter :: M_tempcoef = 1.2956														!temperature coefficient of metabolism
			real(kind = 8), parameter :: M_tempexpo = 0.1170														!temperature exponent of metabolism
			
			totalmass = bodymass + reservemass																	!estimated total mass
			mrate_m = M_masscoef * totalmass ** M_massexpo															!mass-specific metabolic rate at -2C
			mrate_t = mrate_m * M_tempcoef * exp(M_tempexpo * temperature)												!temperature-adjusted metabolic rate
			estmrate = mrate_t
			
		end subroutine simulation_DRS
		
		subroutine simulation_DS(temperature, bodymass, reservemass, estmrate)
			implicit none
			
			real(kind = 8), intent(in) :: temperature, bodymass, reservemass												!inputs = temp., reserve and structural mass
			real(kind = 8) :: estmrate																		!output = est.metabolic rate
			real(kind = 8), parameter :: M_masscoef = 0.00068177														!mass coef. of metabolism (all stages)
			real(kind = 8), parameter :: M_massexpo = 0.7503														!mass exponent of metabolism
			real(kind = 8), parameter :: M_tempcoef = 1.1216														!temperature coefficient of metabolism
			real(kind = 8), parameter :: M_tempexpo = 0.0875														!temperature exponent of metabolism
			real(kind = 8) :: mrate_m, mrate_t, totalmass															!interim :: temp and mass specific metabolic rates
			
			totalmass = bodymass + reservemass																	!estimated total mass
			mrate_m = M_masscoef * totalmass ** M_massexpo															!mass-specific metabolic rate at -2C
			mrate_t = mrate_m * M_tempcoef * exp(M_tempexpo * temperature)												!temperature-adjusted metabolic rate
			estmrate = mrate_t * 0.25																		!metabolic rate reduces by 75% during diapause
			
		end subroutine simulation_DS	
		
		subroutine simulation_TMGR(bodymass, estgrate)
			implicit none
			
			real(kind = 8), intent(in) :: bodymass																!inputs = structural mass 
			real(kind = 8) :: estgrate																		!outputs = growth rate
			real(kind = 8), parameter :: FC = 30.00																!food conversion factor
			real(kind = 8), parameter :: temperature = 12.00														!maximum possible temperature
			real(kind = 8), parameter :: foodcon = 6.00															!maximum possible food concentration
			real(kind = 8), parameter :: G_masscoef = 0.009283														!mass coef. of assimilation (all stages)
			real(kind = 8), parameter :: G_massexpo = 0.7524														!mass expo. of assimilation
			real(kind = 8), parameter :: G_tempcoef = 1.2382														!temperature coef. of assimilation
			real(kind = 8), parameter :: G_tempexpo = 0.0966														!temperature expo. of assimilation
			real(kind = 8) :: GRate_mass, GRate_temp, GRate, gp1														!interim :: mass, temp and food-adjusted assimilation rates
			
			GRate_mass = G_masscoef * bodymass ** G_massexpo														!mass-adjusted assimilation rate
			GRate_temp = GRate_mass * G_tempcoef * exp(G_tempexpo * temperature)											!temperature-scaled assimilation rate
			gp1 = 0.3 * bodymass ** (-0.138)																	!scaling parameter :: parameter 'd' in Bandara et al. 2019
			GRate = (GRate_temp * gp1 * foodcon * FC) / (1 + gp1 * foodcon * FC)											!food-adjusted assimilation rate
			estgrate = GRate																				!assign and export value-1

		end subroutine simulation_TMGR
		
end module simulation
