## Calibration model for COVID-19 outbreaks in Laos  from 11 April 2021 to 12 May 2022
## Treat the whole population as a non age structure population of Laos size 7,389,060 (2021)
## Because of novel coronavirus, a initial proportion susceptible =100% 
## consider four ASC(ascertainment rates): ASC1[30,60] ASC2[61,249]ASC3[250,350]ASC4
## Target: using the confirmed case, deaths and recoveries data to infer the COVID-19 activity and its severity 


###############################
CContact=c("Google mobility","GovResponseIndex","GeometricMean","ArithmeticMean")
## first of all; select the combined contact model: 1 (Google mobility); 2(Oxford Gov Response Index) 3(Geometricl mean) 4(Arithmetic mean)
i_CC =1
EndDate="12/05/2022"				#date up to which the data will be used for model calibration


################################ 
## Read the whole nation data of COVID-19 outbreaks in Laos
Read_in_DATA<-function(Policy) {
 Data_20Jan2021_15Octl2022 <- read.csv(file = "Laos2020-01-22to2022-10-15.csv", head = TRUE)  

 WorkPeriod = which(Data_20Jan2021_15Octl2022$Date_reported=="04/04/2021"):which(Data_20Jan2021_15Octl2022$Date_reported==EndDate);  
###################  =========== Oxford Government pCOVID-19 policy data ============= ############################
 if(Policy=="GovernmentResponseIndex")  LA_PolicyIndex=Data_20Jan2021_15Octl2022$GovernmentResponseIndex[WorkPeriod]
 if(Policy=="StringencyIndex")          LA_PolicyIndex=Data_20Jan2021_15Octl2022$StringencyIndex[WorkPeriod]
 if(length(which(is.na(LA_PolicyIndex)))>0) LA_PolicyIndex[which(is.na(LA_PolicyIndex))] =LA_PolicyIndex[min(which(is.na(LA_PolicyIndex)))-1];	

############ ========== Google mobility data ============= ####################################
 retail_and_recreation= Data_20Jan2021_15Octl2022$retail.and.recreation[WorkPeriod]
 grocery_and_pharmacy = Data_20Jan2021_15Octl2022$grocery.and.pharmacy[WorkPeriod] 
 parks                = Data_20Jan2021_15Octl2022$parks[WorkPeriod]
 transit_stations     = Data_20Jan2021_15Octl2022$transit.stations[WorkPeriod]
 workplaces           = Data_20Jan2021_15Octl2022$workplaces[WorkPeriod]
 residential          = Data_20Jan2021_15Octl2022$residential[WorkPeriod]
 overallaverage       = (retail_and_recreation+grocery_and_pharmacy+parks+transit_stations+workplaces+residential)/6;  # composite mobility
 if(length(which(is.na(overallaverage)))>0) overallaverage[which(is.na(overallaverage))]=overallaverage[min(which(is.na(overallaverage)))-1];	

####   avergae over one week
  retail_and_recreationW= retail_and_recreation
  grocery_and_pharmacyW = grocery_and_pharmacy
  parksW                = parks
  transit_stationsW     = transit_stations   
  workplacesW           = workplaces
  residentialW          = residential

for(i in 8:length(WorkPeriod)) {							#forward 7-day averages
  retail_and_recreationW[i]= sum(retail_and_recreation[(i-7):i])/7
  grocery_and_pharmacyW[i] = sum(grocery_and_pharmacy[(i-7):i])/7
  parksW[i]                = sum(parks[(i-7):i])/7
  transit_stationsW[i]     = sum(transit_stations[(i-7):i])/7   
  workplacesW[i]           = sum(workplaces[(i-7):i])/7
  residentialW[i]          = sum(residential[(i-7):i])/7
}
 LA_overallaverageW       = (retail_and_recreationW+grocery_and_pharmacyW+parksW+transit_stationsW+workplacesW+residentialW)/6
 if(length(which(is.na(LA_overallaverageW)))>0) LA_overallaverageW[which(is.na(LA_overallaverageW))]=LA_overallaverageW[min(which(is.na(LA_overallaverageW)))-1];	

 WorkPeriod2 = which(Data_20Jan2021_15Octl2022$Date_reported=="11/04/2021"):which(Data_20Jan2021_15Octl2022$Date_reported==EndDate);  	

##Outbreak data
 LA_DateDay        = Data_20Jan2021_15Octl2022$Date_reported[WorkPeriod2]
 LA_Cases_local    = Data_20Jan2021_15Octl2022$NewCases_local[WorkPeriod2] 
 LA_Cases_imported = Data_20Jan2021_15Octl2022$NewCases_imported[WorkPeriod2]
 LA_Death          = Data_20Jan2021_15Octl2022$New_deaths[WorkPeriod2]	
 LA_Recovery       = Data_20Jan2021_15Octl2022$Recovery[WorkPeriod2]   

 WorkPeriod3 = which(Data_20Jan2021_15Octl2022$Date_reported=="29/03/2021"):which(Data_20Jan2021_15Octl2022$Date_reported==EndDate);  	#2 weeks delay in effect of vaccination
 LA_sing_vaccine= Data_20Jan2021_15Octl2022$Atleast1dose[WorkPeriod3]            
 LA_Full_vaccine= Data_20Jan2021_15Octl2022$fullvaccine[WorkPeriod3] 

 DATA=list(LA_DateDay,LA_Cases_local,LA_Cases_imported,LA_Death,LA_Recovery,LA_overallaverageW,LA_PolicyIndex,LA_sing_vaccine,LA_Full_vaccine)
 return(DATA)
}
##############################################################################################################




########################################################################################################
##   SEEIIR epidemic dynamics to generate TID (daily no of infections) given model parameter THETA    ##
########################################################################################################
SEEIIRprocess<-function(OldP) {        
 IOnsetday<- rep(0,Dayrun); 					  		# daily infections
 SuscepM  <- rep(1,Dayrun); 					  		# daily susceptibility
 RepNM    <- rep(1,Dayrun); 					  		# Rt 

 Bt_omicrom    = OldP[13]*deltaT;						# transmission coeff of Omicron variant
 T_brk_omicron = OldP[16];							# middlepoint in transition from delta to omicron
 ve1_o= OldP[17];    								# vaccine efficacy for at leat 1 dose after the changepoint t_ve 
 ve2_o= OldP[18];									# vaccine efficacy for full vaccinated after the changepoint t_ve  

 C_delta  =  OldP[19];								# transition speed parameter from Alpha to Delta
 C_omicron=  OldP[20];								# transition speed parameter from Delta to Omicron

 psi_r <- OldP[1];                        			  	# initial growth rate/day from 11 April 2021 before the control measures 
 beta  <- psi_r*((psi_r*dL/2+1)^2)/(1-1/(psi_r*dI/2+1)^2)
 Bt_alpha    <- beta*deltaT;          	   			   	# transmission rate per unit of time (deltaT) for infection with alpha varinat
 I0    <- OldP[4]; 			  		   	   		# initial seeds on 11 April 2021
 T_brk_delta = floor(OldP[5])       					# midddle point of transition from alpha to delta
 Bt_delta    <- OldP[6]*deltaT;   						# transmission rate after T_brk_delta

 sigma_delta    <- OldP[9]; 		sigma_omicron = OldP[10]      # susceptibility changes along with variants of concern
        	
 btt<-rep(Bt_alpha,Dayrun);			

#t_start_delta =86;   								#number of days from 11 April to 6 July 2021 (Delta first reported in Laos), 
#t_start_omicron =296;   							#number of days from 11 April to 01 Feb 2022 (Omicron first reported in Laos), 

 t_mid_delta   =t_start_delta+T_brk_delta;   				#T_brk_delta is assumed as the half-width of transition from Alpha to Delta
 t_mid_omicron =t_start_omicron+T_brk_omicron;   			#T_brk_omicron for transition from Delta to Omicron

 btt[t_start_delta:Dayrun]  <-Bt_alpha +(Bt_delta- Bt_alpha)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)*C_delta));		# sigmoidal function  delta emerges /overtakes alpha
 btt[t_start_omicron:Dayrun]<-Bt_delta+(Bt_omicrom-Bt_delta)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)*C_omicron));	# sigmoidal function omicron emerges /overtakes delta

 VE1=rep(ve1_a,Dayrun);				VE2=rep(ve2_a,Dayrun);	
 VE1[t_start_delta:Dayrun] = ve1_a+(ve1_d-ve1_a)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)*C_delta));
 VE2[t_start_delta:Dayrun] = ve2_a+(ve2_d-ve2_a)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)*C_delta));

 VE1[t_start_omicron:Dayrun] = ve1_d+(ve1_o-ve1_d)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)*C_omicron));
 VE2[t_start_omicron:Dayrun] = ve2_d+(ve2_o-ve2_d)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)*C_omicron));
	

# relative susceptibility change with replacement of Alpha, Delta and Omicron
 sigma =rep(sigma_alpha,Dayrun);	
 sigma[t_start_delta:Dayrun]   = sigma_alpha+(sigma_delta  -sigma_alpha)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)*C_delta));
 sigma[t_start_omicron:Dayrun] = sigma_delta+(sigma_omicron-sigma_delta)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)*C_omicron));


# waning rate change with the replacement of Alpha, Delta and Omicron	
 OG =rep(og0_a,Dayrun);	
 OG[t_start_delta:Dayrun]   = og0_a+(og0_d  -og0_a)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)*C_delta));
 OG[t_start_omicron:Dayrun] = og0_d+(og0_o-og0_d)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)*C_omicron));

# reduced transmissibility change with the replacement of Alpha, Delta and Omicron	after 1 dose vaccine and fully vaccinated
 epsilon2_o =OldP[11];     	epsilon1_o =epsilon2_o;		
 epsilon1_a = epsilon2_a;  	epsilon1_d = epsilon2_d;	
 Epsilon1 =rep(epsilon1_a,Dayrun);	
 Epsilon1[t_start_delta:Dayrun]   = epsilon1_a+(epsilon1_d  -epsilon1_a)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)*C_delta));
 Epsilon1[t_start_omicron:Dayrun] = epsilon1_d+(epsilon1_o-epsilon1_d)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)*C_omicron));

 Epsilon2 =rep(epsilon2_a,Dayrun);	
 Epsilon2[t_start_delta:Dayrun]   = epsilon2_a+(epsilon2_d  -epsilon2_a)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)*C_delta));
 Epsilon2[t_start_omicron:Dayrun] = epsilon2_d+(epsilon2_o-epsilon2_d)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)*C_omicron));


p<-rep(0,17);								#17 compartments

## initial conditions
p[4] <- I0/(1+gm/(gm+psi_r));						# I1
p[5] <- I0-p[4];								# I2
p[3] <- p[4]*(gm+psi_r)/sm;						# E2
p[2] <- p[3]*(sm+psi_r)/sm;						# E1
p[6] <- dV1[1];								# those having at least 1 dose of vaccine
p[11]<- dV2[1];								# those having full vaccine 
p[1] <- Npop*Suscep-p[2]-p[3]-p[4]-p[5]-p[6]-p[11];  		# assuming the proportion Suscep (=100%)
######################

## Runge-Kutta 4th Order to solve the model equations ##
 K1<- rep(0,length(p))
 K2<- rep(0,length(p))
 K3<- rep(0,length(p))
 K4<- rep(0,length(p))

 TIday <-0;        tday<-1;    daycount <- 0;
 IOnsetday[tday]= I0;
 SuscepM[tday] = (p[1]+(1-VE1[tday])*p[6]+(1-VE2[tday])*p[11]+sigma[tday]*p[17])/Npop;	# susceptibility 
 RepNM[tday]   = SuscepM[tday]*(btt[tday]/deltaT)*dI*ccontact[i_CC,tday]*((p[4]+p[5])+Epsilon1[tday]*(p[9]+p[10])+Epsilon2[tday]*(p[14]+p[15]))/(p[4]+p[5]+p[9]+p[10]+p[14]+p[15])

for(i in 1:(daytimes*(Dayrun-1))) {
  ModelP =c(btt[1+tday]*ccontact[i_CC,1+tday],VE1[1+tday],VE2[1+tday],sigma[1+tday],Epsilon1[1+tday],Epsilon2[1+tday],OG[1+tday]);					

  K1<- RKFourstep(ModelP,p); 
  K2<- RKFourstep(ModelP,p+K1/2); 
  K3<- RKFourstep(ModelP,p+K2/2); 
  K4<- RKFourstep(ModelP,p+K3)

  pp <-p + (K1+2*K2+2*K3+K4)/6;	 

###############     New infections    ##      #############
  TIday <-TIday + (pp[3]+ pp[8]+ pp[13])*sm;  							#new infections
   daycount<-daycount+1

  if(daycount==daytimes) {  										# accumulating locall cases to a day
    tday<-tday+1;												# count the number of days
    IOnsetday[tday]<-TIday;      TIday<-0;  							# renew
    if(tday<=Nday) pp[4] <- pp[4]+Cases_imported[tday-1];					# imported cases enter into the trasnmission dynamics  
# --- move the Vaccinated to the their groups occurring daily ---
    pp[1] <- max(0,pp[1]    - dV1[tday]);
    pp[6] <- max(0,pp[6]    + dV1[tday]- dV2[tday]);
    pp[11]<- min(Npop,pp[11]+ dV2[tday]);
    SuscepM[tday] = (pp[1]+(1-VE1[tday])*pp[6]+(1-VE2[tday])*pp[11]+sigma[tday]*pp[17])/Npop;	# susceptibility 
    RepNM[tday]   = SuscepM[tday]*(btt[tday]/deltaT)*dI*ccontact[i_CC,tday]*((pp[4]+pp[5])+Epsilon1[tday]*(pp[9]+pp[10])+Epsilon2[tday]*(pp[14]+pp[15]))/(pp[4]+pp[5]+pp[9]+pp[10]+pp[14]+pp[15])
    	daycount <-0;     	  	   	
  }  		

 p<-pp;    						#Renew for the next step    
}						

 return(list(IOnsetday,SuscepM,RepNM));  		# return time series of infection cases, susceptibility
}
#######################################################################################################################


########## calculate K1,K2,K3,K4 for Runge-Kutta 4th order numerical for solving the differential equatuions ##########
RKFourstep<-function(MP,p){  ##V1 V2 to W
   pp<-rep(0,17)
   BETA=MP[1]; VE1=MP[2];  VE2=MP[3];  SIGMA = MP[4]; EPSILON1= MP[5]; EPSILON2= MP[6];  OG0= MP[7];  OG1 = 2*OG0;  OG2= OG0;    # Limm1=Limm0/2

 LAMBDA = BETA*((p[4]+p[5])+EPSILON1*(p[9]+p[10])+EPSILON2*(p[14]+p[15]))/Npop
 pp[1] <-  -p[1]*LAMBDA;									# S 
 pp[6] <-  -p[6]*LAMBDA*(1-VE1)-OG1*p[6];						      # vaccinated with at least one dose (V1)
 pp[11]<- -p[11]*LAMBDA*(1-VE2)-OG2*p[11];						# fully vaccinated (V2)

 pp[2]<- p[1]*LAMBDA -sm*p[2];								# E1
 pp[3]<- sm*p[2] -sm*p[3];									# E2
 pp[4]<- sm*p[3] -gm*p[4];									# I1
 pp[5]<- gm*p[4] -gm*p[5];									# I2

 pp[7] <- p[6]*LAMBDA*(1-VE1) -sm*p[7];							# E1^V1
 pp[8] <- sm*p[7] -sm*p[8];									# E2^V1
 pp[9] <- sm*p[8] -gm*p[9];									# I1^V1
 pp[10]<- gm*p[9] -gm*p[10];									# I2^V1

 pp[12]<- p[11]*LAMBDA*(1-VE2)+p[17]*LAMBDA*SIGMA -sm*p[12];			# E1^V2
 pp[13]<- sm*p[12] -sm*p[13];									# E2^V2
 pp[14]<- sm*p[13] -gm*p[14];									# I1^V2
 pp[15]<- gm*p[14] -gm*p[15];									# I2^V2

 pp[16]<- gm*(p[5]+p[10]+p[15])-OG0*p[16];						# Recovered and fully immune
 pp[17]<- OG0*p[16]+OG1*p[6]+OG2*p[11]-p[17]*LAMBDA*SIGMA;    			# Waned and partly immune V1,V2 and R return to W once losing immunity

  return(pp)
}       
#########################################################################################################################



############################ Observational model  #####################################
 Reporting<-function(TID,Parameters) {
#  infections:    TIDaily[day]   these include both symptomatic and asymptomatic cases  

#  XiOnsetC ---- distribution of the delay from symptom onset to confirmation    
#  XiCDeath ---  distribution of the delay from confirmation to Death
#  XiCRecovery---distribution of the delay from confirmation to Recovery

 ASC_1 <-Parameters[2];						# Ascertainment rate among infections before changepoint t12
 ASC_2 <-Parameters[24];					# Ascertainment rate among infections btwn changepoint t12 and t23
 ASC_3 <-Parameters[25];					# Ascertainment rate among infections btwn changepoint t23 and t34
 ASC_4 <-Parameters[15];					# Ascertainment rate among infections after changepoint t34
 t12   <-Parameters[23];					# changepoint in confirmation rate among inmfections
 t23   <-Parameters[12];					# changepoint in confirmation rate among inmfections
 t34   <-Parameters[14];					# changepoint in confirmation rate among inmfections

ASCs <- rep(ASC_1,Dayrun)					# continuous changing Ascertainment rate 
ASCs[floor(t12):floor(t23-1)] = ASC_2
ASCs[floor(t23):floor(t34-1)] = ASC_3
ASCs[floor(t34):Dayrun]       = ASC_4

 CFR_alpha    <-Parameters[3];				# case fatality rate given confirmation 
 C_delta      = Parameters[19];
 C_omicron    = Parameters[20];
 CFR_delta    = Parameters[21];
 CFR_omicron  = Parameters[22];
 T_brk_delta  = floor(Parameters[5]);      
 T_brk_omicron= Parameters[16];		

 t_mid_delta  = t_start_delta+T_brk_delta;   		#T_brk_delta is assumed as the half-width of transition from Alpha to Delta
 t_mid_omicron= t_start_omicron+T_brk_omicron;   	#T_brk_omicron is assumed as the half-width of transition from Delta to Omicron

 CFRs<-rep(CFR_alpha,Dayrun)
 CFRs[t_start_delta:Dayrun]  <-CFR_alpha +(CFR_delta- CFR_alpha)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)*C_delta));		# sigmoidal function  delta emerges /overtakes alpha
 CFRs[t_start_omicron:Dayrun]<-CFR_delta+(CFR_omicron-CFR_delta)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)*C_omicron));	# sigmoidal function omicron emerges /overtakes delta

##-----------  Delay from Illness to confirmation/Hospital admission to death/Recovery --------##
##1. Estimation of delay from Symptom onset to Confirmation/Hospitalization 
 meanOnsetConf<-mOC;					
 varOnsetConf <-vOC;
 rateOnsetConf <-meanOnsetConf/varOnsetConf;      
 shapeOnsetConf<-rateOnsetConf*meanOnsetConf;      				

XiOnsetC<-rep(0,45);          				 	## probability of delay from symptom onset to hospital admission in each of 45 days
for(i in 1:45) { XiOnsetC[i]<-pgamma(i,shapeOnsetConf,rateOnsetConf)-pgamma(i-1,shapeOnsetConf,rateOnsetConf); }

##2 Estimation of delay from Hospital to Death 
 meanCDeath <-mCD; 		
 varCDeath  <-vCD;
 rateCDeath <-meanCDeath/varCDeath;
 shapeCDeath<-rateCDeath*meanCDeath; 			
 
XiCDeath<-rep(0,45);          					## probability of delay from Hospital to Deaths in each of 45 days
for(i in 1:45) { XiCDeath[i]<-pgamma(i,shapeCDeath,rateCDeath)-pgamma(i-1,shapeCDeath,rateCDeath); }

##3 Estimation of delay from Hospital to RecoveryDeath 
 meanCRecovery <-mCR; 		
 varCRecovery  <-vCR;
 rateCRecovery <-meanCRecovery/varCRecovery;
 shapeCRecovery<-rateCRecovery*meanCRecovery; 			
 
XiCRecovery<-rep(0,45);          					## probability of delay from Hospital to Recovery in each of 45 days
for(i in 1:45) { XiCRecovery[i]<-pgamma(i,shapeCRecovery,rateCRecovery)-pgamma(i-1,shapeCRecovery,rateCRecovery); }
######----------------------------------------------------------------------------------------########


### three independent observations #
 CONFM     <- rep(0,Dayrun);    	 			# daily Confirmed cases 
 DeathM    <- rep(0,Dayrun);    	 			# daily Deaths 
 RecoveryM <- rep(0,Dayrun);    	 			# daily Recoveries

 for(j in 1:Dayrun) {     					# obtain the number of hospitalizations  from infections
   CONFM[j]<-0;         			
    for(jH in 0:(length(XiOnsetC)-1)) {   if(jH<j) CONFM[j]<-CONFM[j]+ ASCs[j-jH]*TID[j-jH]*XiOnsetC[1+jH];   }   
   }

 CONFMtemp<-CONFM;    						# save for the recovery and death processes

 for(j in 1:Dayrun) {     					# obtain the number of Deaths from HOS
   DeathM[j]<-0;   RecoveryM[j]<-0;       			
   for(jD in 0:(length(XiCDeath)-1)) {   
	if(jD<j) { DeathM[j]  <-DeathM[j]   + CFRs[j-jD]*CONFMtemp[j-jD]*XiCDeath[1+jD];   
		  RecoveryM[j]  <-RecoveryM[j]+ (1-CFRs[j-jD])*CONFMtemp[j-jD]*XiCRecovery[1+jD];   
		  CONFMtemp[j-jD]<-max(0,CONFMtemp[j-jD]-CFRs[j-jD]*CONFMtemp[j-jD]*XiCDeath[1+jD]-(1-CFRs[j-jD])*CONFMtemp[j-jD]*XiCRecovery[1+jD]);    ## remove those died or recoveried from confirmations
	}
   }
 }

 return (list(CONFM,DeathM,RecoveryM));
}


############################################################################################################
### calculate negative log NB likelihood for outbreak data:Confirmastion,Death and cumulative recoveries ###
LLikelihood<-function(HOSM,DeathM,RecoveryM,eta1,eta2) {  		
  Recovery_last=max(which(Recovery>0))   
  LLike =dpois(x =sum(round(Recovery[1:Recovery_last])), lambda =sum(RecoveryM[1:Recovery_last]),log=T);		# poisson density for cumulative cases to Recovery_last
 
  aa1<-log(eta1);      			ab1 <- log(1.-1./eta1);	
  aa2<-log(eta2);      			ab2 <- log(1.-1./eta2);	
	
  rr1<- HOSM[Day_report1]/(eta1-1.0);
  rr2<- DeathM[Day_report2]/(eta2-1.0);		
  
  LLike<- LLike -(LogGmMHOSplus1 + LogGmMDeathplus1) +ab1*sum(Cases_local[Day_report1]) + ab2*sum(Death[Day_report2]);
  LLike<- LLike + sum(lgamma(Cases_local[Day_report1]+rr1)-lgamma(rr1)-rr1*aa1)+ sum(lgamma(Death[Day_report2]+rr2)-lgamma(rr2)-rr2*aa2);

  return (-LLike);								# negative log likelihood
}



############# Sample new parameters with Normal walk around the old parameters for HM algoritham of MCMC ################
SampleParameter<-function(OldPara,stDev,Ich) { 	
 Para = OldPara; 			Para[Ich] <- -10;

 if(Ich!=17&Ich!=18) { while((Para[Ich]<para_a[Ich])|(Para[Ich]>para_b[Ich])) { Para[Ich] = rnorm(1)*stDev[Ich] + OldPara[Ich];  }  }

 if(Ich==17) { while((Para[Ich]<para_a[Ich])|(Para[Ich]>para_b[Ich])|(Para[Ich]>Para[Ich+1])) { Para[Ich] = rnorm(1)*stDev[Ich] + OldPara[Ich];  }  }    #Ve1
 if(Ich==18) { while((Para[Ich]<para_a[Ich])|(Para[Ich]>para_b[Ich])|(Para[Ich]<Para[Ich-1])) { Para[Ich] = rnorm(1)*stDev[Ich] + OldPara[Ich];  }  } 	  #Ve2 to guarantee Ve2>=Ve1

 return(Para);
}



############## Sample new parameters from priori distribution for Bayesian framework ###################
SampleParameterSeeds<-function() {
 NewPara<-rep(-1,ndim);
 for(i in 1:ndim) NewPara[i] <- runif(1,min=PI_shape[i],max=PI_rate[i]); ## 
 return (NewPara);
}


#########------------------------------------------------------------------------------#############
###    Draw the time series of samplings and its distribution                                   ####
Draw_Sample_Density<-function(Name_Parameter,MCMC_Sample,Prior_L_U_MLE) {
 if(Name_Parameter=="negative log likelihood") {
  plot(1:length(MCMC_Sample),MCMC_Sample,col="pink",xlab="#Sample",ylab=Name_Parameter,lty=5,type="b",cex=.4,pch=20,
  main=paste(round(median(MCMC_Sample),2),"[",round(min(MCMC_Sample),2),",",round(max(MCMC_Sample),2),"]")); 

  hist(MCMC_Sample, xlab=Name_Parameter,ylab="Number of simulations",main=paste("AIC:",round(Prior_L_U_MLE[1],2),"DIC_G:",round(Prior_L_U_MLE[2],2),"DIC_S:",round(Prior_L_U_MLE[3],2),sep=""))

 } else { 
##----------------------- MCMC samples ---------------------------##
 Qant_MCMC<-quantile(MCMC_Sample,probs=c(0.025,0.5,0.975));
 plot(1:length(MCMC_Sample),MCMC_Sample,col="pink",xlab="#Sample",ylab=Name_Parameter,lty=5,type="b",cex=.4,pch=20,
 main=paste("MLE=",round(Prior_L_U_MLE[3],3))); 

##-----------------estimate density and produce plot   -----------##
  hist(MCMC_Sample, xlab=Name_Parameter,ylab="Number of simulations",main=paste(round(Qant_MCMC[2],3)," (",round(Qant_MCMC[1],3),",",round(Qant_MCMC[3],3),")",sep=""))
  abline(v=Prior_L_U_MLE[3],lty=2,col="green")
  abline(v=Prior_L_U_MLE[1],lty=2,col="red");    abline(v=Prior_L_U_MLE[2],lty=2,col="red")
 }
}
#-----------------------------------------------------------------------------------------------------#


##----------- Transparent colors    ## Mark Gardener 2015  ## www.dataanalytics.org.uk
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color

## Get RGB values for named color
rgb.val <- col2rgb(color)

## Make new color using input color as base and alpha set by transparency
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255,
             alpha = (100 - percent) * 255 / 100,
             names = name)

## Save the color
invisible(t.col)
}
##------------ END
 mycol2 <- t_col("green", perc = 90, name = "lt.green")



####=======================================================================================#########
## MCMC sampling to calibrate the model to observed data of Cases_local, Deaths and Recoveries #####
MCMC_Sample<-function(Cases_imported,Cases_local,Death,Recovery) {

### parameters for running MCMC programme ########################################################
Loop	<-  10000;    	
Lburnin<- 800000;     
IOP  <-  300;        					# steps for output  "thinning of sampling"
Iadp <-  50;      					# steps for adaptive change of stdeviation of ndim parameters
##############################################################################################

mixingthresholdUB<-0.40; 					# The upper boundary  for the accepted ratio
mixingthresholdLB<-0.12; 					# The lower boundary  for the accepted ratio
Increasingfactor <- 1.2;  					# increase the stDev if accepted ratio is lower than thresholdUB
Decreasingfactor<- 0.8;  					# decrease the stDev if accepted ratio is higher than thresholdLB


########################################################################################################
################################### Variable for storing the epidemic curves #######################################
INFD  <- matrix(0,nrow=Dayrun,ncol=Loop);     		# Storage for mean nos of infections over Dayrun days 
HOSD  <- matrix(0,nrow=Dayrun,ncol=Loop);    	 	# Storage for actual nos of confirmation/HOSpitalizationsover Dayrun days
DeathD<- matrix(0,nrow=Dayrun,ncol=Loop);     		# Storage for actual  nos of Deaths over Dayrun days 
RecoveryD<- matrix(0,nrow=Dayrun,ncol=Loop);    	# Storage for actual  nos of Recoverys over Dayrun days 
SuscepD  <- matrix(0,nrow=Dayrun,ncol=Loop);    	# Storage for susceptibility over Dayrun days   
RepND <- matrix(0,nrow=Dayrun,ncol=Loop);    		# Storage for Reproduction number Rt over Dayrun days   
 
TIDaily      <- rep(0,Dayrun);     				# Daily new number of infections
HOSDaily     <- rep(0,Dayrun);     				# Daily new number of confirmed cases/hospitalizations due to local transmisison
DeathDaily   <- rep(0,Dayrun);     				# Daily new number of Deaths
RecoveryDaily<- rep(0,Dayrun);     				# Daily new number of Recoverys


i_TS		<-0;							# record select iteration
logAlph	<-0;
OldP 		<-rep(0,ndim);			NewP <-rep(0,ndim);    		# sample of parameters 
OldQ 		<-rep(0,ndim);			NewQ <-rep(0,ndim);     	# the inverse of the proposal function
PriorOld	<-rep(0,ndim);			PriorNew<-rep(0,ndim);  	# prior density 
NameSample<-c("psi_r","asc1","CFRa","I0","T_brk_d","beta_d","eta_Case","eta_Death","Sigma_delta","Sigma_omicron","epsilon2_o","t23","Beta_o","t34","asc4","T_brk_o","ve1_o","ve2_o","CFRd","CFRo","t12","asc2","acs3","LLikelihood");

SPsi_r <-rep(0,Loop);	SR   <-rep(0,Loop);		SR2  <-rep(0,Loop);	SBeta_a<-rep(0,Loop);     	SBeta_d<-rep(0,Loop);   	SSigma_o<-rep(0,Loop);    
SASC1 <-rep(0,Loop);  	St12<-rep(0,Loop); 		SASC2 	<-rep(0,Loop);  	St23<-rep(0,Loop);  	SASC3 	<-rep(0,Loop);  	
SCFRa<-rep(0,Loop);	SCFRd<-rep(0,Loop);		SCFRo<-rep(0,Loop);	SI0<-rep(0,Loop);			ST_brk_delta<-rep(0,Loop); 	SSigma_d <-rep(0,Loop); 	 	 	
Seps2_o <-rep(0,Loop);	SBeta_o<-rep(0,Loop);		St34  <-rep(0,Loop);	SASC4 <-rep(0,Loop);		ST_brk_omicron<-rep(0,Loop);	
Sve1_o <-rep(0,Loop);	Sve2_o <-rep(0,Loop);		SC_delta <-rep(0,Loop);	SC_omicron <-rep(0,Loop);

SETA1  <-rep(0,Loop);	SETA2<-rep(0,Loop);		SLL    <-rep(0,Loop)

OldP<-SampleParameterSeeds(); 					

   Outcome  <-SEEIIRprocess(OldP);  			      ## transmission dynamics model to generate newly infections
   TIDaily = Outcome[[1]];

ObsOutput<-Reporting(TIDaily,OldP);  				# Observational model: disease reporting to generate confirmed cases (local), Deaths and REcoverys
   HOSDaily     <-ObsOutput[[1]];		     			#from infection to hospital admissions (Conf/HOS series)  ##
   DeathDaily   <-ObsOutput[[2]];  					#from hospital to Deaths (Death series) ##
   RecoveryDaily<-ObsOutput[[3]];  					#from hospital to Deaths (Death series) ##
 	
### To cover a wide dispersion in observed data Negative Binomial distribution is chosen
  prevLogLike<- LLikelihood(HOSDaily,DeathDaily,RecoveryDaily+1e-10,OldP[7],OldP[8]); 	#calculate neg log NB likelihood of outbreak data 

for(i in 1:length(stDev)) {   					# for the use in calculation of acceptance rate
 OldQ[i] = pnorm((para_b[i]-OldP[i])/stDev[i])-pnorm((para_a[i]-OldP[i])/stDev[i]);      
 PriorOld[i]<-dunif(OldP[i],para_a[i],para_b[i]);   		
}


##  ----------------------------------- Markov Chains -------------------------------------------
print("psi_r  ASC1 CFR_a I0 T_brk_delta Beta_d  eta1  eat2 Sigma_d Sigam_o eps2_o T_brk23 Beta_o T_brk34 ASC4 T_brk_omicron ve1_o ve2_o C_d C_o CFR_d CFR_o T_brk12 ASC2 ASC3 LL");
j<-0;    									# total accepted steps
while(i_TS<Loop) {
numberAccepted<-rep(0,ndim);					    	# empty the accepted number box used to monitor the accepted ratio
for(jj in 1:Iadp) {    							#adaptive
   if(i_TS==Loop) break;
   for(ii in 1:IOP) {   						#output   "thinning"

      for(i1 in SampleSet) {  #1:ndim) {
         NewP<-SampleParameter(OldP,stDev,i1);   		
	 
		if((i1<7)|(i1>8)) {					# for dispersion used in model comparison, no need to generate the SIR time series again
               Outcome  <-SEEIIRprocess(NewP);  		# transmission dynamics model to generate newly infections
               TIDaily = Outcome[[1]]
		   ObsOutput <-Reporting(TIDaily,NewP);		# Observational model:disease reporting to generate Confirmatioin/Hospitalizations,Deaths and Recoveries
		   HOSDaily     <-ObsOutput[[1]];		     	#from infection to hospital admissions (Conf/HOS series)  ##
  		   DeathDaily   <-ObsOutput[[2]];  			#from hospital to Deaths (Death series) ##
  		   RecoveryDaily<-ObsOutput[[3]];  			#from hospital to Deaths (Death series) ##
		}

		postLogLike<- LLikelihood(HOSDaily,DeathDaily,RecoveryDaily+1e-10,NewP[7],NewP[8]); 		
		NewQ[i1] = pnorm((para_b[i1]-NewP[i1])/stDev[i1])-pnorm((para_a[i1]-NewP[i1])/stDev[i1]);   	
		PriorNew[i1]<- dunif(NewP[i1],para_a[i1],para_b[i1]);							  	      
		logAlph = -(postLogLike - prevLogLike);
															   			
            if(runif(1)<=(PriorNew[i1]/PriorOld[i1])*exp(logAlph)*(OldQ[i1]/NewQ[i1])) { 				
              prevLogLike = postLogLike;
              OldP[i1]    =NewP[i1];
              OldQ[i1]    =NewQ[i1];
		  PriorOld[i1]=PriorNew[i1];
              numberAccepted[i1]<- numberAccepted[i1]+1;
		  j <- j+1; 						
            }
       }   									# ndim parameters renew in turm
   }       									#   IOP
if(j>Lburnin) {    
	i_TS <- i_TS+1;

# generate the 'actual numbers of infections, confirmeation/hospitalizations, Deaths, Recovery from the means + ETA by sampling from NB distribution
 INFD[,i_TS]     <-rnbinom(Dayrun,mu=TIDaily[1:Dayrun],size=TIDaily[1:Dayrun]/(OldP[7]-1));
 HOSD[,i_TS]     <-rnbinom(Dayrun,mu=HOSDaily[1:Dayrun],size=HOSDaily[1:Dayrun]/(OldP[7]-1));
 DeathD[,i_TS]   <-rnbinom(Dayrun,mu=DeathDaily[1:Dayrun],size=DeathDaily[1:Dayrun]/(OldP[8]-1));
 RecoveryD[,i_TS]<-rpois(Dayrun,RecoveryDaily[1:Dayrun]);

 SuscepD[,i_TS]= Outcome[[2]];  		# susceptibility
 RepND[,i_TS]  = Outcome[[3]];  		# Reproduction number
	
 SPsi_r[i_TS]<-OldP[1];	
 SASC1[i_TS]     <-OldP[2]; 		SCFRa[i_TS]<-OldP[3];  		SI0[i_TS]  <-OldP[4];   	ST_brk_delta[i_TS]<-OldP[5];  	SBeta_d[i_TS]<-OldP[6];     
 SETA1[i_TS]   <-OldP[7]; 		SETA2[i_TS] <-OldP[8];		SSigma_d[i_TS] <-OldP[9]; 	SSigma_o[i_TS]   <-OldP[10];	
 Seps2_o[i_TS]   <-OldP[11];		St23[i_TS]<-OldP[12];		SBeta_o[i_TS]<-OldP[13];	St34[i_TS]<-OldP[14];			SASC4[i_TS]<-OldP[15]; 
 ST_brk_omicron[i_TS]<-OldP[16];	Sve1_o[i_TS]<-OldP[17];		Sve2_o[i_TS]<-OldP[18]; 	SC_delta[i_TS]<-OldP[19];		SC_omicron[i_TS]<-OldP[20]; 
 SCFRd[i_TS]<-OldP[21];  		SCFRo[i_TS]<-OldP[22];  	
 St12[i_TS]<-OldP[23];  		SASC2[i_TS]<-OldP[24]; 		SASC3[i_TS]<-OldP[25]; 

 SLL[i_TS]  <-postLogLike;
print(paste("psi_r:",round(SPsi_r[i_TS],2)," ASC1:",round(100*SASC1[i_TS],1)," CFRa:",round(100*SCFRa[i_TS],1)," CFRd:",round(100*SCFRd[i_TS],1)," CFRo:",round(100*SCFRo[i_TS],1)," I0:",round(SI0[i_TS])," T_b1:",round(ST_brk_delta[i_TS])," bt_d:",round(SBeta_d[i_TS],3)," Sigma_d:",round(100*SSigma_d[i_TS],1)," Sigma_o:",round(100*SSigma_o[i_TS],1)," eps2_o:",round(100*Seps2_o[i_TS])," T_23:",round(St23[i_TS])," bt_o:",round(SBeta_o[i_TS],3)," T_34:",round(St34[i_TS])," ASC4:",round(100*SASC4[i_TS])," T_b2:",round(ST_brk_omicron[i_TS])," ve1_o:",round(100*Sve1_o[i_TS])," ve2_o:",round(100*Sve2_o[i_TS])," C_d:",round(SC_delta[i_TS],1)," C_o:",round(SC_omicron[i_TS],1)," T_12:",round(St12[i_TS])," ASC2:",round(100*SASC2[i_TS],1)," ASC3:",round(100*SASC3[i_TS],1)," Like=",round(SLL[i_TS],1)," j=",j," i_TS=",i_TS," #Acc=",numberAccepted,sep=""));
}   ### j

} 		# Iadp  

  for(i0 in 1:ndim) {
     ratio = numberAccepted[i0]/(IOP*Iadp);
     if (ratio>mixingthresholdUB) { stDev[i0] = stDev[i0]*Increasingfactor; }
     else if (ratio<mixingthresholdLB) { stDev[i0]=stDev[i0]*Decreasingfactor; }
     print(paste("i0=",i0," #Accepted=",numberAccepted[i0]," ratio=",ratio," SD=",stDev[i0]));
  }

} ## while

 Output <-list(SPsi_r,SASC1,SCFRa,SI0,ST_brk_delta,SBeta_d,SETA1,SETA2,SSigma_d,SSigma_o,Seps2_o,St23,SBeta_o,St34,SASC4,ST_brk_omicron,Sve1_o,Sve2_o,SC_delta,SC_omicron,SCFRd,SCFRo,St12,SASC2,SASC3,SLL,SuscepD,RepND,INFD,HOSD,DeathD,RecoveryD)
 return(Output)
}



####################################################################################################
WeekAverage<-function(ODailyData) {
  L_data = length(ODailyData)
  WeeklyAv= rep(0,L_data);
  WeeklyAv[1:3] =ODailyData[1:3]; 

  for(i in 4:(L_data-3)) WeeklyAv[i] = mean(ODailyData[(i-3):(i+3)]); 	#central average
   WeeklyAv[(L_data-2):(L_data)] =ODailyData[(L_data-2):(L_data)]; 
 return(WeeklyAv)
}


############################################     Mains   ##############################################
##### Data of overall_numbers for Confirmation/Hospitals,deaths and others obtained       	  #####
## DATA ###########----------------------------------------------------------##########################
 Npop <- 7389060;  							# 2020 population size of Laos  

# load the data of Loas## the first two case reported on 24 March 2020
# choose GOV Policy index
Policy="GovernmentResponseIndex";   	#"StringencyIndex" 
LA_data=Read_in_DATA(Policy) 

# set up the length of running epidemic  from 11 APril 2021 to 12 May 2022 				
 DateDay        = LA_data[[1]]; 			
 Cases_local_O  = LA_data[[2]]; 	Cases_local=WeekAverage(Cases_local_O)
 Cases_imported = LA_data[[3]]; 
 Death_O        = LA_data[[4]];   	Death =WeekAverage(Death_O)
 Recovery_O     = LA_data[[5]]; 	Recovery_O[which(is.na(Recovery_O))]=0;	Recovery=WeekAverage(Recovery_O)

 Nday = length(DateDay);							# set up the length of running
###---  ---  	---   	&&&   	****** 	&&&  		--- 	--- 	----    ##
 Npred <- 0; #				      		# number of days projected into future  # set to 0 as we only calibrate data to 12 May 2022
 Dayrun<-Nday+Npred;   						#model epidemic starts from 11/04/2021

 rel_contact = (LA_data[[6]]+100)/100;    		# Google mobility# transform into change relative to baseline activities
 GovPolicy   = (100-LA_data[[7]])/100;			# Oxford Government policy data 

# combined contact models 4 models #Google mobility #Oxford Gov Response Index #geometricl mean #arithmetic mean
#  ccontact=rbind(rel_contact[1:Nday],GovPolicy[1:Nday],sqrt(rel_contact[1:Nday]*GovPolicy[1:Nday]),0.5*(rel_contact[1:Nday]+GovPolicy[1:Nday]))
# extend to 4 week in future  implying that the current control measures are maintained for 4 weeks in future 

rel_contact1 =c(rel_contact,rep(rel_contact[length(rel_contact)],28))
GovPolicy1   =c(GovPolicy,rep(GovPolicy[length(GovPolicy)],28))
  ccontact=rbind(rel_contact1[1:Dayrun],GovPolicy1[1:Dayrun],sqrt(rel_contact1[1:Dayrun]*GovPolicy1[1:Dayrun]),0.5*(rel_contact1[1:Dayrun]+GovPolicy1[1:Dayrun]))


 StartDate<-"days from 11/04/2021";   			# reporting date of index cases

 Day_report1<- 5:length(Cases_local);       		# which(HOS!=0)
 DateDay =as.Date(DateDay, format="%d/%m/%Y") 
 Day_report2<- 5:length(Death);    				#which(Death!=0);    #


###--- generate sum of log(gamma(Death[j]+1)) ... up to day Nday for use in calculation of Negative Binomial likelihood --- ###
 LogGmMHOSplus1  <-sum(lgamma(Cases_local[Day_report1]+1)); 
 LogGmMDeathplus1<-sum(lgamma(Death[Day_report2]+1));

#-------------- daily number of the vaccinated -----------------------
 dV1<- rep(0,Dayrun);  											#these having one dose of vaccine
 dV2<- rep(0,Dayrun);											#those fully vaccinated

 accV1 = LA_data[[8]] 
 dV1[1:length(accV1)]= c(accV1[1],diff(accV1));							#transform into daily number of people vaccinated
 if(Dayrun>length(accV1)) dV1[(1+length(accV1)):Dayrun]=mean(dV1[(length(accV1)-7):length(accV1)]);   	#daily no of vaccinations = daily mean of the past week

 accV2 = LA_data[[9]]   
 dV2[1:length(accV2)]= c(accV2[1],diff(accV2));							#transform into daily number of people vaccinated
 if(Dayrun>length(accV1)) dV2[(1+length(accV2)):Dayrun]=mean(dV2[(length(accV2)-7):length(accV2)]);   	#daily no of vaccinations = daily mean of the past week
 dV1[dV1<=0]=0;
 dV2[dV2<=0]=0;

Suscep<- 1;    				# 2019-nCov is a novel virus and 100% susceptibility is assumed

### parameters for running the transmission model	
daytimes  <-1; #12;   					# parts of day times
deltaT    <- 1.0/daytimes;      			# time step chosen for calculation 


## parameters of transmission dynamics:
#  delay data taken from the estimates in China (Verity et al 2020)
 mOC=4.75;   		 vOC=4.7;    			#mean and variance of delay from symptom onset to confirmation
 mCD=19.5;	   		 vCD=22.7;  			#mean and variance of delay from confirmation to death
 mCR=15.4;  		 vCR=15.5;	  			#mean and variance of delay from confirmation to recovery

##incubation period, infectious period of symptomatic case  from Andre's summary
dL <-5.2;						## incubation (latent) period: 5 days 	based on Qun Li et al 2020; Adam Kucharski et al 2020
sm <-(2/dL)/daytimes ;  			## 1/sm    progression rate from E1-->E2 and E2--> I1
dI <- 3.5;				      	## infectious period  a guess  		based on Jon Read et al 2020; Adam Kucharski et al 2020
gm <-(2/dI)/daytimes ;  			## 1/gm    progression rate from I1 --> I2 I2-->R

##### initial immunity (at week 1 after vaccination) and the duration of immunity for three variants of concern of SARS-CoV-2 virus
ve1_a  = 0.72;					# VE for at leat 1 dose before the changepoint   (Alpha) from Hall et al 2021 
ve2_a  = 0.86;					# VE for full vaccinated before the changepoint 

ve1_d = (0.58+0.43)/2;    			# VE for at leat 1 dose after the changepoint (Delta) from Pouwels et al 2021 Nat Med (Pfizer-BioNTech+AstraZeneca)/2
ve2_d = (0.67+0.82)/2;				# VE for full vaccinated after the changepoint 


### Alpha
L_imm_a =11.6*30.5; #0.7; #1.3;		# overal mean of 6.71+16.39 for fully immune   changed date: 19 Jan 2023
og0_a = (1/L_imm_a)/daytimes;			# loss rate of immunity of the recovered 
og1_a = (1/L_imm_a)/daytimes;			# loss rate of immunity of the vaccinated 1 dose
og2_a = (1/L_imm_a)/daytimes;			# loss rate of immunity of the fully vaccinated  

### Delta
L_imm_d =9.42*30.5;					# overall mea of 11.43+7.78+5.35+4.82+13.39+10.64+12.58 see the file: Decay of vaccine effectiveness  
og0_d = (1/L_imm_d)/daytimes;			# loss rate of immunity of the recovered 
og1_d = (1/L_imm_d)/daytimes;			# loss rate of immunity of the vaccinated 1 dose
og2_d = (1/L_imm_d)/daytimes;			# loss rate of immunity of the fully vaccinated  

### Omicron
L_imm_o =3.22*30.5;				# overal mean of 4.02+2.86+2.78 the file: Decay of vaccine effectiveness
og0_o = (1/L_imm_o)/daytimes;			# loss rate of immunity of the recovered 
og1_o = (1/L_imm_o)/daytimes;			# loss rate of immunity of the vaccinated 1 dose
og2_o = (1/L_imm_o)/daytimes;			# loss rate of immunity of the fully vaccinated  

sigma_alpha =0.16;				# relative susceptible of the wanned from vaccination or natural infection  1-0.84 from Hall et al; 1-0.805  from Hansen et al 2021

### reduction in transmission for the infections that were vaccoinated before Eyre et al 2021 (Oxford U)
epsilon2_a= 0.40;		#average over BNT162b2(0.32) and ChAdOx1(0.48)  that is, compared with the unvaccoinated, infectivity reduces from 1 to epsilon. 
epsilon2_d= 0.63;		#average over BNT162b2(0.50) and ChAdOx1(0.76)   # no data for epsilon_o, to be estimated


 t_start_delta   =86;   								#86: number of days from 11 April to 6 July 2021 (Delta first reported in Laos), 
 t_start_omicron =296;   								#296: number of days from 11 April to 01 Feb 2022 (Omicron first reported in Laos), 

ndim<-25

## ---------- hyperpareameter for prior distribution of model parameters -----------------
# we assume uninformative priors uniform [PI_shape,PI_rate]
cfr0 = max(0.0001,sum(Death)/sum(Cases_local)); T_H=Nday/2
#        Psi_r, ASC1, CFRa,  I0, T_brk_d,beta_d,eta1,eta2,sigma_d,sigma_o,eps2_o,t23,beta_o,t34,ASC4,T_brk_o ve1_o,ve2_o,C_d, C_o,  CFRd,   CFRo,  t12, ASc2, ASC3
PI_shape<-c(.04,.10,cfr0/2,1500.,55.0, 0.70, 35.,  1.01, 0.20,  0.20,   0.63,   90.0,   1.50,   280.,   0.10, 15,   0.17, 0.20, 0.07,0.10, cfr0/2, cfr0/3, 40.0,  0.02, 0.05);  	
PI_rate <-c(.10,.30,cfr0*2,2500, 75.0, 1.20, 70.,  1.10, 0.35,  0.35,   1.00,   155.,   3.00,   340.,   0.25, 35,   0.36, 0.63, 1.00,1.00, cfr0*2, cfr0*1, 60.0,  0.10, 0.10);  

para_a<- c(.02,.01,cfr0/4,100., 20.0, 0.50, 5.0,  1.005,0.13,  0.13,   0.25,   71.,    0.70,   250.,   0.01,  5.,  0.05, 0.10, 0.05,0.05, cfr0/4, cfr0/5, 30.0,  0.01, 0.01);  #lower boundary for the proposal of normal walk
para_b<- c(0.4,.95,cfr0*4,5000, 100., 4.00, 200,  2.00, 0.70,  0.70,   1.00,   249,    10.0,   360.,   0.95,  60,  0.50, 0.75, 2.00,2.00, cfr0*6, cfr0*4, 70.0,  0.95, 0.95);  #upper boundary for the proposal of normal walk
stDev  <-c(.02,.05,cfr0,  250., 7.00, 0.15, 12.,  0.01, 0.10,  0.10,   0.12,   20.,    0.60,   60.0,   0.04,  7.,  0.10, 0.12, 0.30,0.30, cfr0,   cfr0,   10.0,  0.02, 0.02);  #SD for sampling new parameters 

SampleSet=c(1:11,12,13:ndim);		

CalibrationResult=MCMC_Sample(Cases_imported,Cases_local,Death,Recovery);
NameSample<-c("psi_r","ASC1","CFRa","I0","T_brk_delta","beta_d","eta_Case","eta_Death","Sigma_d","Sigma_o","eps2_o","t23","beta_o","t34","ASC4","T_brk_o","ve1_o","ve2_o","C_d","C_o","CFRd","CFRo","t12","ASC2","ASC3","LLikelihood");
 SamplePsi_r   =  CalibrationResult[[1]];			#sample of the model parameters
 SampleASC1      =  CalibrationResult[[2]];
 SampleCFRa=  CalibrationResult[[3]];
 SampleI0      =  CalibrationResult[[4]];
 SampleT_brk_delta  =  CalibrationResult[[5]];		Tbt1 =floor(SampleT_brk_delta)
 SampleBeta_d  =  CalibrationResult[[6]];
 SampleETA1    =  CalibrationResult[[7]];
 SampleETA2    =  CalibrationResult[[8]];
 SampleSigma_d =  CalibrationResult[[9]];		
 SampleSigma_o =  CalibrationResult[[10]];
 Sampleeps2_o  =  CalibrationResult[[11]];
 Samplet23  =  CalibrationResult[[12]];
 SampleBeta_o  =  CalibrationResult[[13]];
 Samplet34 =  CalibrationResult[[14]];
 SampleASC4     =  CalibrationResult[[15]];
 SampleT_brk_omicron =  CalibrationResult[[16]];
 Sampleve1_o   =  CalibrationResult[[17]];
 Sampleve2_o   =  CalibrationResult[[18]];
 SampleC_delta =  CalibrationResult[[19]];
 SampleC_omicron=  CalibrationResult[[20]];

 SampleCFRd    =  CalibrationResult[[21]];
 SampleCFRo    =  CalibrationResult[[22]];

 Samplet12  =  CalibrationResult[[23]];
 SampleASC2     =  CalibrationResult[[24]];
 SampleASC3     =  CalibrationResult[[25]];
 SampleLL      =  CalibrationResult[[26]];

 SuscepDSample = CalibrationResult[[27]];
 SampleR  = CalibrationResult[[28]];
 INFDSample     = CalibrationResult[[29]];		#samples of the epidmic
 HOSDSample     = CalibrationResult[[30]];
 DeathDSample   = CalibrationResult[[31]];
 RecoveryDSample= CalibrationResult[[32]];

 SampleBeta_a <- SamplePsi_r*((SamplePsi_r*dL/2+1)^2)/(1-1/(SamplePsi_r*dI/2+1)^2);     	#transmission coeff before T_brk_delta

i_1=1; 
LS =length(SampleLL);   					# select the sample range
S_USED2= i_1:LS
## selection the MLE(Maximum Likelihood Estimations) of model parameters   postLogLike is negative log likelihood
 MaxLL=min(SampleLL[S_USED2])
 i_MLE= which(SampleLL==MaxLL)

 MaxPsi_r   <- SamplePsi_r[i_MLE];  
 MaxASC1 <- SampleASC1[i_MLE];     				Maxt12 <- Samplet12[i_MLE];  	MaxASC2 <- SampleASC2[i_MLE];     MaxASC3 <- SampleASC3[i_MLE];    
 MaxCFRa<- SampleCFRa[i_MLE];					MaxCFRd<- SampleCFRd[i_MLE];			MaxCFRo<- SampleCFRo[i_MLE];	
 MaxI0      <- SampleI0[i_MLE];	
 MaxT_brk_delta<- SampleT_brk_delta[i_MLE];		MaxBeta_a   <- SampleBeta_a[i_MLE];		MaxBeta_d   <- SampleBeta_d[i_MLE];	
 MaxSigma_d <- SampleSigma_d[i_MLE];			MaxSigma_o<- SampleSigma_o[i_MLE];			
 MaxETA1    <- SampleETA1[i_MLE];				MaxETA2    <- SampleETA2[i_MLE];	
 MaxR       <- SampleR[,i_MLE]/SuscepDSample[,i_MLE];		
 Maxeps2_o   = Sampleeps2_o[i_MLE];				Maxt23 = Samplet23[i_MLE];	 
 MaxBeta_o = SampleBeta_o[i_MLE];				Maxt34 = Samplet34[i_MLE];	MaxASC4 = SampleASC4[i_MLE];
 MaxT_brk_omicron = SampleT_brk_omicron[i_MLE];		Maxve1_o= Sampleve1_o[i_MLE];			Maxve2_o= Sampleve2_o[i_MLE];
 MaxC_delta= SampleC_delta[i_MLE];				MaxC_omicron = SampleC_omicron[i_MLE];	 		


###### Calculate the Deviance information Criterion (DIC) to measure the model fitting  #############
 Devbar<- mean(2*SampleLL);
 DevV  <- var(2*SampleLL);

   p_DG  <- DevV/2;			#Gelman et al (2004)'s definition of effective no of model pareameters
   DIC_G <- p_DG +Devbar;		#Gelman's definition of DIC
###-------------------------------------------------------------------------------------##############
#====== calculate the Deviance at the mean of sample
mPsi_r=mean(SamplePsi_r[S_USED2]);   			mASC1=mean(SampleASC1[S_USED2]);				mCFRa=mean(SampleCFRa[S_USED2])
mI0=mean(SampleI0[S_USED2]);					mT_brk_delta=mean(SampleT_brk_delta[S_USED2]);		mBeta_d=mean(SampleBeta_d[S_USED2])
mETA1=mean(SampleETA1[S_USED2]);				mETA2=mean(SampleETA2[S_USED2]);				mSigma_d=mean(SampleSigma_d[S_USED2])
mSigma_o=mean(SampleSigma_o[S_USED2]);			meps2_o=mean(Sampleeps2_o[S_USED2]);			mt23=mean(Samplet23[S_USED2]); 	mASC3=mean(SampleASC3[S_USED2])
mBeta_o=mean(SampleBeta_o[S_USED2]);			mt34=mean(Samplet34[S_USED2]);			mASC4=mean(SampleASC4[S_USED2])
mT_brk_omicron=mean(SampleT_brk_omicron[S_USED2]);	mve1_o=mean(Sampleve1_o[S_USED2]);				mve2_o=mean(Sampleve2_o[S_USED2])
mC_delta=mean(SampleC_delta[S_USED2]);			mC_omicron=mean(SampleC_omicron[S_USED2]);		mCFRd=mean(SampleCFRd[S_USED2])
mCFRo=mean(SampleCFRo[S_USED2]);				mt12=mean(Samplet12[S_USED2]);			mASC2=mean(SampleASC2[S_USED2])

NewP	<-c(mPsi_r,mASC1,mCFRa,mI0,mT_brk_delta,mBeta_d,mETA1,mETA2,mSigma_d,mSigma_o,meps2_o,mt23,mBeta_o,mt34,mASC4,mT_brk_omicron,mve1_o,mve2_o,mC_delta,mC_omicron,mCFRd,mCFRo,mt12,mASC2,mASC3)

Outcome  <-SEEIIRprocess(NewP);  			      	# transmission dynamics model to generate newly infections
TIDm = Outcome[[1]]
ObsOutput <-Reporting(TIDm,NewP);					# Observational model:disease reporting to generate Hospitalizations,and Deaths

HOSDm     <-ObsOutput[[1]];  		     				#from symptomatic infections to hospital (HOS time series)  ## 
DeathDm   <-ObsOutput[[2]];    						#from hospital to Deaths (Death series)
RecoveryDm<-ObsOutput[[3]];  						#from hospital to Recoverys (Recovery series)


  Dev_at_mean <-2*LLikelihood(HOSDm,DeathDm,RecoveryDm,mETA1,mETA2); 									# Deviance = 2* neg log likelihood
#=== Spiegelhalter et al (2002)'s definition of the effective no of model parameters ====#
  p_DS  <- Devbar -Dev_at_mean
  DIC_S <- p_DS +Devbar

#====== AIC  ====================
  AIC = 2*length(SampleSet)+2*min(SampleLL)



#######################################################  outputs     #########################################################
#### --------------------------------------------------  outpus of epidemics  ------------------------------------------- ####
###   reproduce the maximum likely epidemic   ######
TIDMAX   <- rep(0,Dayrun);     	 					# Total infections each day for MLE epdiemic
HOSDMAX  <- rep(0,Dayrun);     	 					# Hospitalizations each day
DeathDMAX<- rep(0,Dayrun);     	 					# Total Deaths each day for MLE epdiemic
RecoveryDMAX<- rep(0,Dayrun);     	 					# Total Recoverys each day for MLE epdiemic
NewP	<-c(MaxPsi_r,MaxASC1,MaxCFRa,MaxI0,MaxT_brk_delta,MaxBeta_d,MaxETA1,MaxETA2,MaxSigma_d,MaxSigma_o,Maxeps2_o,Maxt23,MaxBeta_o,Maxt34,MaxASC4,MaxT_brk_omicron,Maxve1_o,Maxve2_o,MaxC_delta,MaxC_omicron,MaxCFRd,MaxCFRo,Maxt12,MaxASC2,MaxASC3)

Outcome  <-SEEIIRprocess(NewP);  			      		# transmission dynamics model to generate newly infections
TIDMAX = Outcome[[1]]
ObsOutput <-Reporting(TIDMAX,NewP);						# Observational model:disease reporting to generate Hospitalizations,and Deaths

HOSDMAX     <-ObsOutput[[1]];  		     				#from symptomatic infections to hospital (HOS time series)  ## 
DeathDMAX   <-ObsOutput[[2]];    						#from hospital to Deaths (Death series)
RecoveryDMAX<-ObsOutput[[3]];  						#from hospital to Recoverys (Recovery series)

 SuscepDMax =Outcome[[2]];  							# susceptibility
 ReprDMax =Outcome[[3]];  							# reproduction number
###########################################################################################################
INF.quant   <-apply(INFDSample[,S_USED2],1,quantile,probs=c(0.025,0.5,0.975)); 			# generate the 95% CI and median for infections
SusD.quant  <-apply(SuscepDSample[,S_USED2],1,quantile,probs=c(0.025,0.5,0.975));  		# generate the 95% CI and median for Susceptibility
ReprD.quant  <-apply(SampleR[,S_USED2],1,quantile,probs=c(0.025,0.5,0.975));  			# generate the 95% CI and median for effective reproduction number
HOSD.quant  <-apply(HOSDSample[,S_USED2],1,quantile,probs=c(0.025,0.5,0.975));  			# generate the 95% CI and median for Hospitalizations

DeathD.quant   <-apply(DeathDSample[,S_USED2],1,quantile,probs=c(0.025,0.5,0.975));   		# generate the 95% CI and median for Deaths
RecoveryD.quant<-apply(RecoveryDSample[,S_USED2],1,quantile,probs=c(0.025,0.5,0.975));   	# generate the 95% CI and median for Recoverys


####===========================################# estimates of model parameters #############----==================================================###
OutputFileName<-paste("MCMCParameters-",Nday,today,"_","CContact_",CContact[i_CC],".pdf",sep="")
pdf(OutputFileName);
par(mfrow=c(3,2))
 Draw_Sample_Density("Psi_r",SamplePsi_r[S_USED2],c(para_a[1],para_b[1],MaxPsi_r)) 								# inital growth rate psi_r ----------------------------#
 Draw_Sample_Density("R0(Alpha)",SampleBeta_a[S_USED2]*dI,c(0,0,MaxBeta_a*dI)) 								# basic reproductive number R -------------------------#
 Draw_Sample_Density("Initial seeds(I0)",SampleI0[S_USED2],c(para_a[4],para_b[4],MaxI0))							# initial seeds I0 ------------------------------------#
 Draw_Sample_Density("Beta_Alpha",SampleBeta_a[S_USED2],c(0,0,MaxBeta_a)) 									# Beta_a before T_brk_delta ---------------------------#
 Draw_Sample_Density("Beta_Delta",SampleBeta_d[S_USED2],c(para_a[6],para_b[6],MaxBeta_d)) 						# Beta_d between two turning points -------------------#
 Draw_Sample_Density("beta_omicron",SampleBeta_o[S_USED2],c(para_a[13],para_b[13],MaxBeta_o)) 						# Beta_omicron after T_brk_delta   --------------------#
 Draw_Sample_Density("midpoint from alpha to delta",SampleT_brk_delta[S_USED2],c(para_a[5],para_b[5],MaxT_brk_delta)) 		# turning point T_brk_delta in transmission -----------#
 Draw_Sample_Density("midpoint from delta to omicron",SampleT_brk_omicron[S_USED2],c(para_a[16],para_b[16],MaxT_brk_omicron)) # changepoint T_brk_omicron in transmission -----------#
 Draw_Sample_Density("Sigma_delta(suscep)",SampleSigma_d[S_USED2],c(para_a[9],para_b[9],MaxSigma_d)) 					# susceptibility to delta variant of SARS-CoV-2 virus--#
 Draw_Sample_Density("Sigma_omicron(suscep)",SampleSigma_o[S_USED2],c(para_a[10],para_b[10],MaxSigma_o)) 				# susceptibility to omicron        --------------------#
 Draw_Sample_Density("ASC1(%)",100*SampleASC1[S_USED2],100*c(para_a[2],para_b[2],MaxASC1)) 						# ascertainment rate 1 among Infections  --------------#
 Draw_Sample_Density("TurningPoint1",Samplet12[S_USED2],c(para_a[23],para_b[23],Maxt12)) 							# 1st changepoint in ASC    ---------------------------#
 Draw_Sample_Density("ASC2(%)",100*SampleASC2[S_USED2],100*c(para_a[24],para_b[24],MaxASC2)) 						# ascertainment rate 2 among Infections  --------------#
 Draw_Sample_Density("TurningPoint2",Samplet23[S_USED2],c(para_a[12],para_b[12],Maxt23)) 							# 2nd changepoint in ASC -#
 Draw_Sample_Density("ASC3(%)",100*SampleASC3[S_USED2],100*c(para_a[25],para_b[25],MaxASC3)) 						# ascertainment rate 3  among Infections  -------------#
 Draw_Sample_Density("TurningPoint3",Samplet34[S_USED2],c(para_a[14],para_b[14],Maxt34)) 							# 3rd changepoint in ASC    ---------------------------#
 Draw_Sample_Density("ASC4(%)",100*SampleASC4[S_USED2],100*c(para_a[15],para_b[15],MaxASC4)) 						# ascertainment rate 4    -----------------------------#
 Draw_Sample_Density("ETA_Case",SampleETA1[S_USED2],c(para_a[7],para_b[7],MaxETA1)) 							# dispersion of confirmation/Hospitalizations   -------#
 Draw_Sample_Density("ETA_Death",SampleETA2[S_USED2],c(para_a[8],para_b[8],MaxETA2)) 							# dispersion of deaths among Hospitalizations   -------#
 Draw_Sample_Density("CFR_Alpha(%)",100*SampleCFRa[S_USED2],100*c(para_a[3],para_b[3],MaxCFRa)) 					# proportion of Death among confirmed cases(CFRa)(Alpha)-#
 Draw_Sample_Density("CFR_delta(%)",100*SampleCFRd[S_USED2],100*c(para_a[21],para_b[21],MaxCFRd)) 					# CFR_delta    ----------------------------------------#
 Draw_Sample_Density("CFR_omicron(%)",100*SampleCFRo[S_USED2],100*c(para_a[22],para_b[22],MaxCFRo)) 					# CFR_omicron   ---------------------------------------#
 Draw_Sample_Density("IFR(%)",100*SampleASC1[S_USED2]*SampleCFRa[S_USED2],100*c(para_a[2]*para_a[3],para_b[2]*para_b[3],MaxASC1*MaxCFRa)) # IFR proportion of Deaths among infections#
 Draw_Sample_Density("rel_infectivity(%)_full doses(omicron)",100*Sampleeps2_o[S_USED2],100*c(para_a[11],para_b[11],Maxeps2_o)) # rel infectivity of 2 doses vaccinated individuals (epsilon2_Omicron)--------#
 Draw_Sample_Density("ve(%)_single dose(omicron)",100*Sampleve1_o[S_USED2],100*c(para_a[17],para_b[17],Maxve1_o)) 		# 1 dose vaccine efficacy omicron  (ve1_omicron) ------#
 Draw_Sample_Density("ve(%)_full doses(omicron)",100*Sampleve2_o[S_USED2],100*c(para_a[18],para_b[18],Maxve2_o)) 			# vaccine efficacy omicron (ve2_omicron)    -----------#
 Draw_Sample_Density("C_delta",SampleC_delta[S_USED2],c(para_a[19],para_b[19],MaxC_delta)) 						# C_delta, width of transition of delta  --------------#
 Draw_Sample_Density("C_omicron",SampleC_omicron[S_USED2],c(para_a[20],para_b[20],MaxC_omicron)) 					# C_omicron, width of transition of omicron -----------#
 Draw_Sample_Density("negative log likelihood",SampleLL[S_USED2],c(AIC,DIC_G,DIC_S))   							# negative log likelihood among samples and DIC  ------#

dev.off();



#####  ---------------------  distribution of epidmeiological statuses ------------ ############
pdf("EpiFitSTATUS.pdf");
 par(mfcol=c(4,1),mar=numeric(4),oma=c(4,4,.5,5.5),mgp=c(2,0.6,0))
##### -----------------------1. daily effective reproduction number -----------------------------------------------------------
plot(as.Date(DateDay),ReprD.quant[2,],col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Rt",lty=1,type="b",ylim=c(0,1.2*max(ReprD.quant[3,]))); #,main =paste("")); 
 abline(h=1,col="black",lty=5);
 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(ReprD.quant[1,],rev(ReprD.quant[3,])),col="grey",border=NA); 

 abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);		
 abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2);	
 lines(as.Date(DateDay),ReprD.quant[2,],lty="solid",col="black",lwd=2.5)
 legend("topleft","A) Effective reproduction numer",bty="n",cex=1.35);
 axis(2,at=seq(0,2.0,by=.5),labels=seq(0,2.0,by=.5),cex.axis=1.4,col="blue",col.lab="blue",las=2) #axis(4L,cex.axis=1.15);  
 mtext("Rt", side=2, line=2.5, cex=0.9,las=0, col="black")
 box(); 

##### -----------------------2. susceptiblity -----------------------------------------------------------
plot(as.Date(DateDay),SusD.quant[2,],col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Daily susceptibility",lty=1,type="b",ylim=c(0,1.2*max(SusD.quant[2,]))); #,main =paste("")); 

 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(SusD.quant[1,],rev(SusD.quant[3,])),col="grey",border=NA); 
 abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);		
 abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2);	
 lines(as.Date(DateDay),SusD.quant[2,],lty="solid",col="black",lwd=2.5)
 legend("topleft","B) Susceptibility",bty="n",cex=1.35);
 axis(2,at=seq(0,1.0,by=.2),labels=seq(0,1.0,by=.2),cex.axis=1.4,col="blue",col.lab="blue",las=2); #axis(2L,cex.axis=1.15); 
 mtext("Susceptibility", side=2, line=2.5, cex=0.9,las=0, col="black")
 legend("bottomleft", legend=c("mid-point from Alpha to Delta", "mid-point from Delta to Omicron"),col=c("red", "purple"), lty = 1:2,bty="n", cex=1.1)
 box(); 


##### --------------------3 daily vaccine effectiveness (VE1-- 1 dose, VE2 -- full vaccinated) -------------------------max(MaxVE1,MaxVE2)-------------------
  t_mid_delta   =t_start_delta+median(SampleT_brk_delta);   				#T_brk_delta is assumed as the half-width of transition from Alpha to Delta
  t_mid_omicron =t_start_omicron+median(SampleT_brk_omicron);   				#T_brk_omicron is assumed as the half-width of transition from Delta to Omicron

 MaxVE1=rep(ve1_a,Dayrun);				MaxVE2=rep(ve2_a,Dayrun);	
 MaxVE1[t_start_delta:Dayrun] = ve1_a+(ve1_d-ve1_a)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)/8));
 MaxVE2[t_start_delta:Dayrun] = ve2_a+(ve2_d-ve2_a)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)/8));

 MaxVE1[t_start_omicron:Dayrun] = ve1_d+(median(Sampleve1_o)-ve1_d)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)/8));
 MaxVE2[t_start_omicron:Dayrun] = ve2_d+(median(Sampleve2_o)-ve2_d)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)/8));

plot(as.Date(DateDay),MaxVE2,col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Daily Vaccine Effectiveness",lty=1,type="b",ylim=c(0,1.2),lwd=.5); #,main =paste("full(black) 1dose(blue)")); 

 abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);		
 abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2);	
 lines(as.Date(DateDay), MaxVE1, col = "black",  type = "l");		#lty = 3,
 legend("topleft","C) Vaccine effectiveness and vaccinations",bty="n",cex=1.35);
 axis(2,at=seq(0,1.0,by=.2),labels=seq(0,1.0,by=.2),cex.axis=1.4,col="blue",col.lab="blue",las=2)  #axis(2L,cex.axis=1.15); 
 mtext("Vaccine effectiveness", side=2, line=2.5, cex=0.9,las=0, col="black")
box(); 

##### --------------------4 cumulative vaccine effectiveness (VE1, VE2) --------------------------------------------------------
accV_max = max(accV1,accV2)
rel_V    =1/Npop 
 axis(4,at=seq(0,1.0,by=.2),labels=seq(0,100,by=20),cex.axis=1.4,col="blue",las=2)
 mtext("Vaccinations(%)", side=4, line=3.0, cex.lab=0.9,las=0, col="blue")
 lines(as.Date(DateDay),rel_V*accV2[1:Dayrun],col="blue",lty=5,type="b")
 lines(as.Date(DateDay), rel_V*accV1[1:Dayrun], col = "blue", type = "l"); 


##### ---------------------5 weekly average google mobility data --------------------------------------------------------
  rel_contact =(LA_data[[6]]+100)/100;  	## Google mobility
  GovPolicy  =(100-LA_data[[7]])/100;	## Gov policy "GovernmentResponseIndex";  

plot(as.Date(DateDay),rel_contact[1:length(DateDay)],col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Relative Mobility",lty=5,type="b",ylim=c(0,max(1.2,rel_contact))); 
 round(cor(rel_contact[1:length(DateDay)],GovPolicy[1:length(DateDay)]),3)

abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);			
abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2);	
abline(h=1,col="black",lty=6);	
legend("topleft","D) Relative mobility and Government response index",bty="n",cex=1.35); 		
axis.Date(1, at=seq(min(as.Date(DateDay)), max(as.Date(DateDay)), by="1 mon"), format="%d/%m/%Y",cex.axis=1.2); 
 axis(2,at=seq(0,1.0,by=.2),labels=seq(0,1.0,by=.2),cex.axis=1.4,col="blue",col.lab="blue",las=2); ##axis(2L,cex.axis=1.15);  
 mtext("Relative mobility", side=2, line=2.5, cex=0.9,las=0, col="black")
box();

##### ---------------------6 weekly average OXford government policy data --
 axis(4,at=seq(0,1.0,by=.2),labels=seq(0,1.0,by=.2),cex.axis=1.4,col="blue",col.lab="blue",las=2)
 mtext("Gov response index", side=4, line=3.0, cex.lab=0.9,las=0, col="blue")
lines(as.Date(DateDay),GovPolicy[1:length(DateDay)],col="blue",lty=5,type="b")
 mtext("Dates",side=1,outer=TRUE,line=2.2)
dev.off();



###-----------------------------------------------------------------------------------------------------####
#ascertainment rate and breaking points
ASC_rate1=quantile(SampleASC1[S_USED2],probs=c(0.025,0.5,0.975));  ASC_rate2= quantile(SampleASC2[S_USED2],probs=c(0.025,0.5,0.975));
ASC_rate3=quantile(SampleASC3[S_USED2],probs=c(0.025,0.5,0.975)); ASC_rate4= quantile(SampleASC4[S_USED2],probs=c(0.025,0.5,0.975));
T_bk1 = ceiling(Samplet12[i_MLE])
T_bk2 = ceiling(Samplet23[i_MLE])
T_bk3 = ceiling(Samplet34[i_MLE])
ASC_RATE= matrix(0,ncol=Nday,nrow=3)

ASC_RATE[,        1:T_bk1]=ASC_rate1; 
ASC_RATE[,(1+T_bk1):T_bk2]=ASC_rate2; 
ASC_RATE[,(1+T_bk2):T_bk3]=ASC_rate3; 
ASC_RATE[,(1+T_bk3):Nday] =ASC_rate4; 


pdf("EpiCurveFit.pdf");
 par(mfcol=c(4,1),mar=numeric(4),oma=c(4,4,.5,5.5),mgp=c(2,0.6,0));   #par(mfrow=c(2,2))
##### ----------------------- infections ---------daily incidence--------------------------------------------------
plot(as.Date(DateDay),INF.quant[2,],col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Daily number of infections",lty=1,type="b",ylim=c(0,1.15*max(INF.quant[3,]))); 
 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(INF.quant[1,],rev(INF.quant[3,])),col="grey",border=NA); 

abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);		
abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2)
   lines(as.Date(DateDay),INF.quant[2,],lty="solid",col="black",lwd=2.5)

legend("topleft","A) Infections",bty="n",cex=1.35);
max_infection =max(INF.quant[,1:Nday])
 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(max_infection*ASC_RATE[1,]/.6,rev(max_infection*ASC_RATE[3,]/.6)),col=mycol2,border=NA); 
   lines(as.Date(DateDay),max_infection*ASC_RATE[2,]/.6,lty="solid",col="green",lwd=2.5)
 axis(4,at=seq(0,max_infection,by=max_infection/6),labels=seq(0,60,by=10),cex.axis=1.4,col="blue",col.lab="blue",las=2); ##axis(2L,cex.axis=1.15);  
 mtext("Ascertainment rate(%)", side=4, line=2.5, cex=0.9,las=0, col="green")

axis(2L,cex.axis=1.2);  box(); 
 mtext("Incidence", side=2, line=2.5, cex=0.9,las=0, col="black")


##### -----------------confirmed cases ---------------------------------------------
plot(as.Date(DateDay),HOSD.quant[2,],col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Daily no. of confirmed cases",ylim=c(0,max(Cases_local,Cases_imported,HOSD.quant[3,])),cex=1.0,lty=1,type="b");
 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(HOSD.quant[1,],rev(HOSD.quant[3,])),col="grey",border=NA); 

abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);		
abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2);	

points(as.Date(DateDay),Cases_local,cex=.8,pch=17,col="blue")
points(as.Date(DateDay[which(Cases_imported>0)]),Cases_imported[which(Cases_imported>0)],cex=.8,pch=17,col="red")
   lines(as.Date(DateDay),HOSD.quant[2,],lty="solid",col="black",lwd=2.5)
legend("topleft","B) Confirmed infections",bty="n",cex=1.35);
axis(2L,cex.axis=1.2); 
 mtext("confirmed cases", side=2, line=2.5, cex=0.9,las=0, col="black")
legend("left", legend=c("mid-point from Alpha to Delta", "mid-point from Delta to Omicron"),col=c("red", "purple"), lty = 1:2,bty="n", cex=1.1)

##### ----------------------- Death -----------------------------------------------------------
plot(as.Date(DateDay),DeathD.quant[2,],col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Daily no. of deaths",lty=1,type="b",ylim=c(0,max(Death,DeathD.quant[3,]))); 
 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(DeathD.quant[1,],rev(DeathD.quant[3,])),col="grey",border=NA); 

abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);		
abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2);	

points(as.Date(DateDay[which(Death>0)]),Death[which(Death>0)],cex=.8,pch=17,col="blue") 
   lines(as.Date(DateDay),DeathD.quant[2,],lty="solid",col="black",lwd=2.5)
legend("topleft","C) Deaths",bty="n",cex=1.35);
 mtext("Deaths", side=2, line=2.5, cex=0.9,las=0, col="black")
axis(2L,cex.axis=1.2); 

##### ----------------------- Recovery -----------------------------------------------------------
plot(as.Date(DateDay),RecoveryD.quant[2,],col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Daily no. of recoveries",lty=1,type="b",ylim=c(0,max(Recovery,RecoveryD.quant[3,]))); 
 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(RecoveryD.quant[1,],rev(RecoveryD.quant[3,])),col="grey",border=NA); 

abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);		
abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2);	

points(as.Date(DateDay[1:max(which(Recovery>0))]),Recovery[1:max(which(Recovery>0))],,cex=.8,pch=17,col="blue")
   lines(as.Date(DateDay),RecoveryD.quant[2,],lty="solid",col="black",lwd=2.5)
legend("topleft","D) Recoveries",bty="n",cex=1.35);
axis.Date(1, at=seq(min(as.Date(DateDay)), max(as.Date(DateDay)), by="1 mon"), format="%d/%m/%Y",cex.axis=1.2);
axis(2L,cex.axis=1.2); 
 mtext("Dates",side=1,cex.lab=1.25,outer=TRUE,line=2.2);  mtext("The recovered", side=2, line=2.5, cex=0.9,las=0, col="black")
dev.off();


save(SamplePsi_r,SampleI0,SampleASC1,SampleASC2,Samplet12,SampleASC3,Samplet23,SampleCFRa,SampleCFRd,SampleCFRo,SampleT_brk_delta,SampleBeta_d,SampleSigma_d,SampleSigma_o,Sampleeps2_o,SampleBeta_o,Samplet34,SampleASC4,SampleT_brk_omicron,Sampleve1_o,Sampleve2_o,SampleC_delta,SampleC_omicron,SampleETA1,SampleETA2,SampleLL,file ="MPEstimates.RData")
save(SampleR,SuscepDSample,INFDSample,HOSDSample,DeathDSample,RecoveryDSample, file ="EPCurves.RData")
