# Counterfactual analysis based on the estimates of model parameters obtained from LaosCOVID-19_Variants_Vacc_Mobility.R

#here we simulate the epicurve based on the sampled values of model parameters and assess the impact of control measures and vaccination
i_CMV =0   
# 0-- baseline (Vaccination & NPIs);  
# 1-- No NPIs (vaccination only);     
# 2-- No Vaccination (NPI only);             
# 3-- No NPIs and No Vaccination          
# 4-- Ve=0.5*VE  5--Ve=1   
# 6-- permenant immunity 
# 7-- Ve=1 & permenant immunity
# 8-- No NPIs;Ve=1 
# 9-- No NPIs;permenant immunity 
# 10--No NPIs; Ve=1 & permenant immunity 
# 11--No NPIs; Ve=0.5Ve


# Select the combined contact model: 1 (Google mobility); 2(Oxford Gov Response Index) 3(Geometricl mean) 4(Arithmetic mean)

if(i_CMV==1|i_CMV==3|i_CMV>7) { i_CC =5		#counterfacual situation: mixing rate remains unchanged	
} else {			i_CC =1;			#baseline situation: mixing rate varied according to Google mobility data
}
EndDate="12/05/2022"
###############################
CContact=c("Google mobility","GovResponseIndex","GeometricMean","ArithmeticMean","Counterfactual")

#######     Read the whole nation data of Laos    #################
Read_in_DATA<-function(Policy) {
 Data_20Jan2021_15Octl2022 <- read.csv(file = "Laos2020-01-22to2022-10-15.csv", head = TRUE)  

 WorkPeriod = which(Data_20Jan2021_15Octl2022$Date_reported=="04/04/2021"):which(Data_20Jan2021_15Octl2022$Date_reported==EndDate);  	#forward average of mobility
###################  =========== Oxford Government pCOVID-19 policy data ============= ############################
 if(Policy=="GovernmentResponseIndex")  LA_PolicyIndex=Data_20Jan2021_15Octl2022$GovernmentResponseIndex[WorkPeriod]
 if(Policy=="StringencyIndex")          LA_PolicyIndex=Data_20Jan2021_15Octl2022$StringencyIndex[WorkPeriod]
 if(length(which(is.na(LA_PolicyIndex)))>0) LA_PolicyIndex[which(is.na(LA_PolicyIndex))] =LA_PolicyIndex[min(which(is.na(LA_PolicyIndex)))-1];	#the last 14 days are assumed to be the same as that of last day 15

############ ========== Google mobility data ============= ####################################
 retail_and_recreation= Data_20Jan2021_15Octl2022$retail.and.recreation[WorkPeriod]
 grocery_and_pharmacy = Data_20Jan2021_15Octl2022$grocery.and.pharmacy[WorkPeriod] 
 parks                = Data_20Jan2021_15Octl2022$parks[WorkPeriod]
 transit_stations     = Data_20Jan2021_15Octl2022$transit.stations[WorkPeriod]
 workplaces           = Data_20Jan2021_15Octl2022$workplaces[WorkPeriod]
 residential          = Data_20Jan2021_15Octl2022$residential[WorkPeriod]
 overallaverage       = (retail_and_recreation+grocery_and_pharmacy+parks+transit_stations+workplaces+residential)/6;  # composite mobility
 if(length(which(is.na(overallaverage)))>0) overallaverage[which(is.na(overallaverage))]=overallaverage[min(which(is.na(overallaverage)))-1];	# the last 3 days are assumed to be the same as that of last day 4

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
 if(length(which(is.na(LA_overallaverageW)))>0) LA_overallaverageW[which(is.na(LA_overallaverageW))]=LA_overallaverageW[min(which(is.na(LA_overallaverageW)))-1];	#the last 3 days are assumed to be the same as that of last day 4

 WorkPeriod2 = which(Data_20Jan2021_15Octl2022$Date_reported=="11/04/2021"):which(Data_20Jan2021_15Octl2022$Date_reported==EndDate);  	#in view of delayed action, Policy started from 04/04/2021 start from 11 Apr 2021  day 101 from 1 Jan 2021
##Outbreak data
 LA_DateDay        = Data_20Jan2021_15Octl2022$Date_reported[WorkPeriod2]
 LA_Cases_local    = Data_20Jan2021_15Octl2022$NewCases_local[WorkPeriod2] 
 LA_Cases_imported = Data_20Jan2021_15Octl2022$NewCases_imported[WorkPeriod2]
 LA_Death          = Data_20Jan2021_15Octl2022$New_deaths[WorkPeriod2]	
 LA_Recovery       = Data_20Jan2021_15Octl2022$Recovery[WorkPeriod2]   

 WorkPeriod3 = which(Data_20Jan2021_15Octl2022$Date_reported=="29/03/2021"):which(Data_20Jan2021_15Octl2022$Date_reported==EndDate);  	#14 days since vacciona started from 21/03/2021 start from 11 Apr 2021  day 101 from 1 Jan 2021
if(i_CMV==2|i_CMV==3) {		# choices for no vaccinations
 LA_sing_vaccine=  rep(0,length(WorkPeriod3)); #         
 LA_Full_vaccine=  rep(0,length(WorkPeriod3)); # 
} else {
 LA_sing_vaccine= Data_20Jan2021_15Octl2022$Atleast1dose[WorkPeriod3];         
 LA_Full_vaccine= Data_20Jan2021_15Octl2022$fullvaccine[WorkPeriod3];  
}

 return(list(LA_DateDay,LA_Cases_local,LA_Cases_imported,LA_Death,LA_Recovery,LA_overallaverageW,LA_PolicyIndex,LA_sing_vaccine,LA_Full_vaccine))
}
##############################################################################################################




########################################################################################################
##   SEEIIR epidemic dynamics to generate TID (daily no of infections) given model parameter THETA    ##
########################################################################################################
SEEIIRprocess<-function(OldP) {        
 IOnsetday<- rep(0,Dayrun); 					  		#daily infections
 SuscepM  <- rep(1,Dayrun); 					  		#daily susceptibility
 RepNM    <- rep(1,Dayrun); 					  		#daily Rt 

 Bt_omicrom    = OldP[13]*deltaT;						# transmission coeff of Omicron variant
 T_brk_omicron = OldP[16];							# middlepoint in transition from delta to omicron
 ve1_o= OldP[17];    								# VE for at leat 1 dose after the changepoint  
 ve2_o= OldP[18];									# VE for full vaccinated after the changepoint  

 C_delta  =  OldP[19];  							# transition speed parameter from Alpha to Delta
 C_omicron  =  OldP[20];							# transition speed parameter from Delta to Omicron

 psi_r <- OldP[1];                        			  	# initial growth rate per day from 11 April 2021 before the control measures 
 beta  <- psi_r*((psi_r*dL/2+1)^2)/(1-1/(psi_r*dI/2+1)^2)
 Bt_alpha    <- beta*deltaT;          	   			   	# transmission rate per unit of time DeltaT for alpha variant
 I0    <- OldP[4];   		   	   					# initial seeds 
 T_brk_delta = floor(OldP[5])       					# middle point of transition from alpha to delta
 Bt_delta    <- OldP[6]*deltaT;   						# transmission rate for delta

 sigma_delta    <- OldP[9]; 		sigma_omicron = OldP[10]      # susceptibility change along with variants of concern that circulate
        	
 btt<-rep(Bt_alpha,Dayrun);			

 t_mid_delta   =t_start_delta+T_brk_delta;   				#T_brk_delta is assumed as the half-width of transition from Alpha to Delta
 t_mid_omicron =t_start_omicron+T_brk_omicron;   			#T_brk_omicron is assumed as the half-width of transition from Delta to Omicron

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

### test halving the VE  also changing the immunity period 
 if(i_CMV==4|i_CMV==11) {   VE1 =0.5*VE1;    			VE2 =0.5*VE2; 						}
 if(i_CMV==5|i_CMV==8)  {   VE1 =rep(1,length(VE1));    	VE2 =rep(1,length(VE2)); 				}
 if(i_CMV==6|i_CMV==9)  {   OG  =rep(0,length(OG)) 										}	
 if(i_CMV==7|i_CMV==10) {   VE1 =rep(1,length(VE1));    	VE2 =rep(1,length(VE2)); OG =rep(0,length(OG)) 	}	
 
# reduced transmissibility change with the replacement of Alpha, Delta and Omicron	after 1 dose vaccine and fully vaccinated
 epsilon2_o =OldP[11];     		epsilon1_o =epsilon2_o;		
 epsilon1_a = epsilon2_a;  		epsilon1_d = epsilon2_d;	
 Epsilon1 =rep(epsilon1_a,Dayrun);	
 Epsilon1[t_start_delta:Dayrun]   = epsilon1_a+(epsilon1_d  -epsilon1_a)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)*C_delta));
 Epsilon1[t_start_omicron:Dayrun] = epsilon1_d+(epsilon1_o-epsilon1_d)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)*C_omicron));

 Epsilon2 =rep(epsilon2_a,Dayrun);	
 Epsilon2[t_start_delta:Dayrun]   = epsilon2_a+(epsilon2_d  -epsilon2_a)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)*C_delta));
 Epsilon2[t_start_omicron:Dayrun] = epsilon2_d+(epsilon2_o-epsilon2_d)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)*C_omicron));


p<-rep(0,17)

##seeding -- initial conditions: 
p[4] <- I0/(1+gm/(gm+psi_r));						# I1
p[5] <- I0-p[4];								# I2
p[3] <- p[4]*(gm+psi_r)/sm;						# E2
p[2] <- p[3]*(sm+psi_r)/sm;						# E1
p[6] <- dV1[1];								# those having at least 1 dose of vaccine
p[11]<- dV2[1];								# those having full vaccine 
p[1] <- Npop*Suscep-p[2]-p[3]-p[4]-p[5]-p[6]-p[11];  		# assuming the proportion Suscep(=100%)
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
  TIday <-TIday + (pp[3]+ pp[8]+ pp[13])*sm;  
   daycount<-daycount+1

  if(daycount==daytimes) {  										# accumulating locall cases to a day
    tday<-tday+1;												# count the number of days
    IOnsetday[tday]<-TIday;      TIday<-0;  							
    if(tday<=Nday) pp[4] <- pp[4]+Cases_imported[tday-1];					# imported cases enter into the trasnmission dynamics  
# -------move the Vaccinated to the their groups --------------------
    pp[1] <- max(0,pp[1]    - dV1[tday]);
    pp[6] <- max(0,pp[6]    + dV1[tday]- dV2[tday]);
    pp[11]<- min(Npop,pp[11]+ dV2[tday]);
    SuscepM[tday] = (pp[1]+(1-VE1[tday])*pp[6]+(1-VE2[tday])*pp[11]+sigma[tday]*pp[17])/Npop;	# susceptibility 
    RepNM[tday]   = SuscepM[tday]*(btt[tday]/deltaT)*dI*ccontact[i_CC,tday]*((pp[4]+pp[5])+Epsilon1[tday]*(pp[9]+pp[10])+Epsilon2[tday]*(pp[14]+pp[15]))/(pp[4]+pp[5]+pp[9]+pp[10]+pp[14]+pp[15])
    	daycount <-0;     	  	   	
  }  		####

p<-pp;    						 
}							

 return(list(IOnsetday,SuscepM,RepNM));  					# return time series of infection cases, susceptibility
}
#######################################################################################################################


############ calculate K1,K2,K3,K4 for Runge-Kutta 4th order numerical for soving the differential equatuions ##########
RKFourstep<-function(MP,p){  
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
 pp[17]<- OG0*p[16]+OG1*p[6]+OG2*p[11]-p[17]*LAMBDA*SIGMA;    			# Waned and partly immune   V2 and R return to W once losing immunity

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



## Transparent colors    ## Mark Gardener 2015  ## www.dataanalytics.org.uk
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
invisible(t.col)		## Save the color
}  ## END
 mycol2 <- t_col("green", perc = 90, name = "lt.green")


#-----------------------------------------------------------------------------------------------------#
## Draw the distribution of sampled paramneter                                     ####
Plot_samples<-function(SamplePARA,MLE_PARA,LABEL) {
  QantPARA <- quantile(SamplePARA,probs=c(0.025,0.5,0.975));
  hist(SamplePARA, xlab=LABEL,ylab="Number of simulations",main=paste(round(QantPARA[2],3)," (",round(QantPARA[1],3),",",round(QantPARA[3],3),")",sep=""))
  abline(v=MLE_PARA,lty=2,col="green")
  #legend("topright",c(paste("median=",round(QantPARA[2],3)),paste("95%CI=[",round(QantPARA[1],3),",",round(QantPARA[3],3) ,"]")),cex=0.8);
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
### Data of overall_numbers for Confirmation/Hospitals,deaths, vaccnation and others obtained  	   ####
## DATA ###########----------------------------------------------------------##########################
 Npop <- 7389060;  							## 2020 population size of Laos  

# load the data of Loas## the first two case reported on 24 March 2020
# choose your GOV Policy index
Policy="GovernmentResponseIndex";   #"StringencyIndex" 
LA_data=Read_in_DATA(Policy) 
				
 DateDay        = LA_data[[1]]; 						# set up the length of running epidemic  from 11 APril 2021 to 12 May 2022 
 Cases_local_O  = LA_data[[2]]; 	Cases_local=WeekAverage(Cases_local_O)
 Cases_imported = LA_data[[3]]; 
 Death_O        = LA_data[[4]];   	Death=WeekAverage(Death_O)
 Recovery_O     = LA_data[[5]]; 	Recovery_O[which(is.na(Recovery_O))]=0;	Recovery=WeekAverage(Recovery_O)

 Nday = length(DateDay);							# set up the length of running
###---  ---  	---   	&&&   	****** 	&&&  		--- 	--- 	----    ##
 Npred <- 0; #				      				# number of days projected into future
 Dayrun<-Nday+Npred;   								#model epidemic starts from 11/04/2021

  rel_contact = (LA_data[[6]]+100)/100;    		## Google mobility# transform into change relative to baseline activities
  GovPolicy   = (100-LA_data[[7]])/100;			## Oxford Government policy data 

## combined contact models 4 models #Google mobility #Oxford Gov Response Index #geometricl mean #arithmetic mean
#  ccontact=rbind(rel_contact[1:Nday],GovPolicy[1:Nday],sqrt(rel_contact[1:Nday]*GovPolicy[1:Nday]),0.5*(rel_contact[1:Nday]+GovPolicy[1:Nday]))
# extend to 4 week in future  implying that the current control measures are maintained for 4 weeks in future 
rel_contact1 =c(rel_contact,rep(rel_contact[length(rel_contact)],28))
GovPolicy1   =c(GovPolicy,rep(GovPolicy[length(GovPolicy)],28))
 ccontact=rbind(rel_contact1[1:Dayrun],GovPolicy1[1:Dayrun],sqrt(rel_contact1[1:Dayrun]*GovPolicy1[1:Dayrun]),0.5*(rel_contact1[1:Dayrun]+GovPolicy1[1:Dayrun]),rep(1,Dayrun))

 StartDate<-"days from 11/04/2021";   		# reporting date of index cases

 Day_report1<- 5:length(Cases_local);       	# which(HOS!=0)
 DateDay =as.Date(DateDay, format="%d/%m/%Y") 
 Day_report2<- 5:length(Death);    			#which(Death!=0);    #

###--- generate sum of log(gamma(Death[j]+1)) ... up to day Nday for use in calculation of Negative Binomial likelihood --- ###
 LogGmMHOSplus1  <-sum(lgamma(Cases_local[Day_report1]+1)); 
 LogGmMDeathplus1<-sum(lgamma(Death[Day_report2]+1));

#-------------- daily number of the vaccinated -----------------------
 dV1<- rep(0,Dayrun);  											#these having one dose of vaccine
 dV2<- rep(0,Dayrun);											#those fully vaccinated

 accV1 = LA_data[[8]] 
 dV1[1:length(accV1)]= c(accV1[1],diff(accV1));						#transform into daily number of people vaccinated
 if(Dayrun>length(accV1)) dV1[(1+length(accV1)):Dayrun]=mean(dV1[(length(accV1)-7):length(accV1)]);   	#daily no of vaccinations = daily mean of the past week

 accV2 = LA_data[[9]]   
 dV2[1:length(accV2)]= c(accV2[1],diff(accV2));						#transform into daily number of people vaccinated
 if(Dayrun>length(accV1)) dV2[(1+length(accV2)):Dayrun]=mean(dV2[(length(accV2)-7):length(accV2)]);   	#daily no of vaccinations = daily mean of the past week
 dV1[dV1<=0]=0;
 dV2[dV2<=0]=0;

 Suscep<- 1;    						# 2019-nCov is a novel virus and 100% susceptibility is assumed

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
ve1_a  = 0.72 #-0.15;				# vaccine efficacy for at leat 1 dose before the changepoint t_ve   (Alpha) from Hall et al 2021 
ve2_a  = 0.86 #-0.15;				# vaccine efficacy for full vaccinated before the changepoint t_ve  over infections

ve1_d = (0.58+0.43)/2 #-0.15;    		# vaccine efficacy for at leat 1 dose after the changepoint t_ve (Delta) from Pouwels et al 2021 Nat Med (Pfizer-BioNTech+AstraZeneca)/2
ve2_d = (0.67+0.82)/2 #-0.15;			# vaccine efficacy for full vaccinated after the changepoint t_ve  over infections

### Alpha
L_imm_a =11.6*30; 				# overal mean of 6.71+16.39 for fully immune   changed date: 19 Jan 2023
og0_a = (1/L_imm_a)/daytimes;			# loss rate of immunity of the recovered 
og1_a = (1/L_imm_a)/daytimes;			# loss rate of immunity of the vaccinated 1 dose
og2_a = (1/L_imm_a)/daytimes;			# loss rate of immunity of the fully vaccinated  

### Delta
L_imm_d =9.42*30;					# overall mea of 11.43+7.78+5.35+4.82+13.39+10.64+12.58 see the file: Decay of vaccine effectiveness  
og0_d = (1/L_imm_d)/daytimes;			# loss rate of immunity of the recovered 
og1_d = (1/L_imm_d)/daytimes;			# loss rate of immunity of the vaccinated 1 dose
og2_d = (1/L_imm_d)/daytimes;			# loss rate of immunity of the fully vaccinated  

### Omicron
L_imm_o =3.22*30;					# overal mean of 4.02+2.86+2.78 the file: Decay of vaccine effectiveness
og0_o = (1/L_imm_o)/daytimes;			# loss rate of immunity of the recovered 
og1_o = (1/L_imm_o)/daytimes;			# loss rate of immunity of the vaccinated 1 dose
og2_o = (1/L_imm_o)/daytimes;			# loss rate of immunity of the fully vaccinated  

sigma_alpha =0.16;				# relative susceptible of the wanned from vaccination or natural infection  1-0.84 from Hall et al; 1-0.805  from Hansen et al 2021

### reduction in transmission for the infections that were vaccoinated before Eyre et al 2021 (Oxford U)
epsilon2_a= 0.40;					#average over BNT162b2(0.32) and ChAdOx1(0.48)  that is, compared with the unvaccoinated, infectivity reduces from 1 to epsilon. 
epsilon2_d= 0.63;					#average over BNT162b2(0.50) and ChAdOx1(0.76)  # no data for epsilon_o, to be estimated

t_start_delta   =86;   				#86: number of days from 11 April to 6 July 2021 (Delta first reported in Laos), 
t_start_omicron =296;   			#296: number of days from 11 April to 01 Feb 2022 (Omicron first reported in Laos), 

			

#######################################################  outputs     #########################################################
#------------------------input of esimates of model parametrs--------------------------------# 
 load("MPEstimates.RData")
 LLoop=length(SampleLL);	i_MLE = which.min(SampleLL)		

################################### Variable for storing the epidemic curves #######################################
INFD  <- matrix(0,nrow=Dayrun,ncol=LLoop);     	 # Storage for mean nos of infections  
HOSD  <- matrix(0,nrow=Dayrun,ncol=LLoop);     	 # Storage for nos of Confirmation/HOSpitalizations
DeathD<- matrix(0,nrow=Dayrun,ncol=LLoop);     	 # Storage for nos of Deaths
RecoveryD<- matrix(0,nrow=Dayrun,ncol=LLoop);    # Storage for nos of Recoverys
SuscepD  <- matrix(0,nrow=Dayrun,ncol=LLoop);    # Storage for susceptibility 
RepND <- matrix(0,nrow=Dayrun,ncol=LLoop);    	 # Storage for Reproduction number Rt
 
TIDaily      <- rep(0,Dayrun);     			 # Daily new number of infections
HOSDaily     <- rep(0,Dayrun);     			 # Daily new number of Confirmation/Hospitalizations due to local transmisison
DeathDaily   <- rep(0,Dayrun);     			 # Daily new number of Deaths
RecoveryDaily<- rep(0,Dayrun);     			 # Daily new number of Recoverys

S_USED= 1:LLoop

 SampleBeta_a <- SamplePsi_r*((SamplePsi_r*dL/2+1)^2)/(1-1/(SamplePsi_r*dI/2+1)^2);     	#transmission coeff before T_brk_delta
for(iL in S_USED) {
#### --------------------------------------------------  outpus of epidemics  ------------------------------------------- ####
NewP	<-c(SamplePsi_r[iL],SampleASC1[iL],SampleCFRa[iL],SampleI0[iL],SampleT_brk_delta[iL],SampleBeta_d[iL],SampleETA1[iL],SampleETA2[iL],SampleSigma_d[iL],SampleSigma_o[iL],Sampleeps2_o[iL],Samplet23[iL],SampleBeta_o[iL],Samplet34[iL],SampleASC4[iL],SampleT_brk_omicron[iL],Sampleve1_o[iL],Sampleve2_o[iL],SampleC_delta[iL],SampleC_omicron[iL],SampleCFRd[iL],SampleCFRo[iL],Samplet12[iL],SampleASC2[iL],SampleASC3[iL])
  
Outcome  <-SEEIIRprocess(NewP);  			      ## transmission dynamics model to generate newly infections
 TIDaily = Outcome[[1]]
 ObsOutput <-Reporting(TIDaily,NewP);					# Observational model:disease reporting to generate Hospitalizations,and Deaths

 HOSDaily      <-ObsOutput[[1]];  		     				#from symptomatic infections to hospital (HOS time series)  ## 
 DeathDaily    <-ObsOutput[[2]];    					#from hospital to Deaths (Death series)
 RecoveryDaily <-ObsOutput[[3]];  						#from hospital to Recoverys (Recovery series)

### generate the 'actual numbers of infections, hospitalizations, Deaths, Recovery from the means + ETA by sampling from NB distribution
 INFD[,iL]  <- rnbinom(Dayrun,mu=TIDaily[1:Dayrun],size=TIDaily[1:Dayrun]/(NewP[7]-1));
 HOSD[,iL]  <- rnbinom(Dayrun,mu=HOSDaily[1:Dayrun],size=HOSDaily[1:Dayrun]/(NewP[7]-1));

 DeathD[,iL]<- rnbinom(Dayrun,mu=DeathDaily[1:Dayrun],size=DeathDaily[1:Dayrun]/(NewP[8]-1));
 RecoveryD[,iL]<-rpois(Dayrun,RecoveryDaily[1:Dayrun]);

 SuscepD[,iL]= Outcome[[2]];  	# susceptibility
 RepND[,iL]  = Outcome[[3]];  	# Reproduction number

}  #iL

T_a= round(t_start_delta+SampleT_brk_delta)
T_d= round(t_start_omicron+SampleT_brk_omicron)

 ### Calculate IFR for three variants of SARS-CoV-2 virus
  Alpha_ifr=Delta_ifr=Omicron_ifr=rep(0,dim(INFD)[2])
  for(i_ifr in 1:dim(INFD)[2]) {
     Alpha_ifr[i_ifr]   =sum(DeathD[1:T_a[i_ifr],i_ifr])/sum(INFD[1:T_a[i_ifr],i_ifr]); 
     Delta_ifr[i_ifr]   =sum(DeathD[(1+T_a[i_ifr]):T_d[i_ifr],i_ifr])/sum(INFD[(1+T_a[i_ifr]):T_d[i_ifr],i_ifr]); 
     Omicron_ifr[i_ifr] =sum(DeathD[(1+T_d[i_ifr]):Nday,i_ifr])/sum(INFD[(1+T_d[i_ifr]):Nday,i_ifr]); 
  }
   QantAlpha_ifr <- quantile(Alpha_ifr,probs=c(0.025,0.5,0.975));
   QantDelta_ifr <- quantile(Delta_ifr,probs=c(0.025,0.5,0.975));
   QantOmicron_ifr <- quantile(Omicron_ifr,probs=c(0.025,0.5,0.975));
   ifr_a=QantAlpha_ifr[2]; ifr_d=QantDelta_ifr[2]; ifr_o=QantOmicron_ifr[2]; 
         DeathTotal<-round(quantile(colSums(DeathD),probs=c(0.025,0.5,0.975)));  

if(i_CMV==3) { 	Inf_baseline =colSums(INFD);  Death_baseline =colSums(DeathD);  						# counterfactual situatiion: No NPIs and No vaccine
} else {       	Inf_avoid 	=100*(colSums(INFD)-Inf_baseline)/Inf_baseline; 
			Death_avoid =100*(colSums(DeathD)-Death_baseline)/Death_baseline;
         		#QuantInf_avoid  <-round(quantile(Inf_avoid,probs=c(0.025,0.5,0.975)),2);  
         		#QuantDeath_avoid<-round(quantile(Death_avoid,probs=c(0.025,0.5,0.975)),2);  
}			# proportion (%) of infections avoided compared to counterfactual situation

###########################################################################################################
 INF.quant   <-apply(INFD,1,quantile,probs=c(0.025,0.5,0.975)); 			# generate the 95% CI and median for infections
 SusD.quant  <-apply(SuscepD,1,quantile,probs=c(0.025,0.5,0.975));  		# generate the 95% CI and median for Susceptibility
 ReprD.quant  <-apply(RepND,1,quantile,probs=c(0.025,0.5,0.975));  		# generate the 95% CI and median for effective reproduction number
 HOSD.quant  <-apply(HOSD,1,quantile,probs=c(0.025,0.5,0.975));  			# generate the 95% CI and median for Hospitalizations

 DeathD.quant   <-apply(DeathD,1,quantile,probs=c(0.025,0.5,0.975));   		# generate the 95% CI and median for Deaths
 RecoveryD.quant<-apply(RecoveryD,1,quantile,probs=c(0.025,0.5,0.975));   	# generate the 95% CI and median for Recoverys


#####  ---------------------  distribution of Infections and model fitting ------------ ############
QantT_brk_delta   <-quantile(t_start_delta+SampleT_brk_delta[S_USED],probs=c(0.025,0.5,0.975));
QantT_brk_omicron <-quantile(t_start_omicron+SampleT_brk_omicron[S_USED],probs=c(0.025,0.5,0.975));



pdf(paste("Fig5_",i_CMV,".pdf"));
 par(mfcol=c(4,1),mar=numeric(4),oma=c(4,4,.5,5.5),mgp=c(2,0.6,0))
##### -----------------------1. daily effective reproduction number -----------------------------------------------------------
plot(as.Date(DateDay),ReprD.quant[2,],col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Rt",lty=1,type="b",ylim=c(0,2.3)); #1.4*max(ReprD.quant[3,]))); 
 abline(h=1,col="black",lty=5);
 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(ReprD.quant[1,],rev(ReprD.quant[3,])),col="grey",border=NA); 

 abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);		
 abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2);	
 lines(as.Date(DateDay),ReprD.quant[2,],lty="solid",col="black",lwd=2.5)
 legend("topleft","a",bty="n",col="dark",cex=1.6); 
 axis(2,at=seq(0,2.0,by=.5),labels=seq(0,2.0,by=.5),cex.axis=1.4,col="blue",col.lab="blue",las=2) #axis(4L,cex.axis=1.15);  
 mtext("Rt", side=2, line=2.5, cex=0.9,las=0, col="black")
 box(); 

##### -----------------------. transmission coefficients -----------------------------------------------------------
  t_mid_delta   =t_start_delta+median(SampleT_brk_delta);   				#T_brk_delta is assumed as the half-width of transition from Alpha to Delta
  t_mid_omicron =t_start_omicron+median(SampleT_brk_omicron);   				#T_brk_omicron is assumed as the half-width of transition from Delta to Omicron

  QantS_Beta_a <- quantile(SampleBeta_a,probs=c(0.025,0.5,0.975));    
  QantS_Beta_d <- quantile(SampleBeta_d,probs=c(0.025,0.5,0.975));
  QantS_Beta_o <- quantile(SampleBeta_o,probs=c(0.025,0.5,0.975));
  Beta=matrix(0,nrow=3,ncol=Dayrun);				
for(i_ci in 1:3) {
 Beta[i_ci,1:t_start_delta] = QantS_Beta_a[i_ci]
 Beta[i_ci,t_start_delta:Dayrun] = QantS_Beta_a[i_ci]+(QantS_Beta_d[i_ci]-QantS_Beta_a[i_ci])/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)*median(SampleC_delta)));
 Beta[i_ci,t_start_omicron:Dayrun] = QantS_Beta_d[i_ci]+(QantS_Beta_o[i_ci]-QantS_Beta_d[i_ci])/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)*median(SampleC_omicron)));
}

 Max_Beta = max(Beta)
 axis(4,at=seq(0,2.0,by=.4),labels=seq(0,3.5,by=.4*(3.5/2)),cex.axis=1.4,col="green",col.lab="green",las=2)
 mtext("Transmission coeff", side=4, line=3.0, cex.lab=0.9,las=0, col="green") 
 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c((2/3.5)*Beta[1,],rev((2/3.5)*Beta[3,])),col=mycol2,border=NA); 
lines(as.Date(DateDay),(2/3.5)*Beta[2,],col="green",lty=5,type="b")


##### --------------------3 daily vaccine effectiveness (VE1-- 1 dose, VE2 -- full vaccinated) -------------------------max(MaxVE1,MaxVE2)-------------------
  QantS_ve1_o <- quantile(Sampleve1_o,probs=c(0.025,0.5,0.975));    QantS_ve2_o <- quantile(Sampleve2_o,probs=c(0.025,0.5,0.975));
 MaxVE1=matrix(ve1_a,nrow=3,ncol=Dayrun);				MaxVE2=matrix(ve2_a,nrow=3,ncol=Dayrun);	
for(i_ci in 1:3) {
 MaxVE1[i_ci,t_start_delta:Dayrun] = ve1_a+(ve1_d-ve1_a)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)*median(SampleC_delta)));
 MaxVE2[i_ci,t_start_delta:Dayrun] = ve2_a+(ve2_d-ve2_a)/(1+exp(-(t_start_delta:Dayrun-t_mid_delta)*median(SampleC_omicron)));

 MaxVE1[i_ci,t_start_omicron:Dayrun] = ve1_d+(QantS_ve1_o[i_ci]-ve1_d)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)*median(SampleC_delta)));
 MaxVE2[i_ci,t_start_omicron:Dayrun] = ve2_d+(QantS_ve2_o[i_ci]-ve2_d)/(1+exp(-(t_start_omicron:Dayrun-t_mid_omicron)*median(SampleC_omicron)));
}


plot(as.Date(DateDay),MaxVE2[2,],col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Daily Vaccine Effectiveness",lty=1,type="b",ylim=c(0,1.2),lwd=.5);  
 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(MaxVE2[1,],rev(MaxVE2[3,])),col="grey",border=NA); 
 abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);		
 abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2);	
 lines(as.Date(DateDay),MaxVE2[2,],lty="solid",col="black",lwd=2.5)

 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(MaxVE1[1,],rev(MaxVE1[3,])),col="grey",border=NA); 
 lines(as.Date(DateDay), MaxVE1[2,], col = "black",  type = "l");		

 legend("topleft","b",bty="n",col="dark",cex=1.6);	
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
#### Google mobility
  rel_contact =ccontact[i_CC,]   #(LA_data[[6]]+100)/100;  
## Gov policy "GovernmentResponseIndex";  
  GovPolicy  = (100-LA_data[[7]])/100

plot(as.Date(DateDay),rel_contact[1:length(DateDay)],col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Relative Mobility",lty=5,type="b",ylim=c(0,max(1.2,rel_contact))); 

abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);			
abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2);	
abline(h=1,col="black",lty=6);	
legend("topleft","c",bty="n",col="dark",cex=1.6); #legend("topleft","C) Relative mobility and susceptibility",bty="n",cex=1.35); 		

 axis(2,at=seq(0,1.0,by=.2),labels=seq(0,1.0,by=.2),cex.axis=1.4,col="blue",col.lab="blue",las=2); 
 mtext("Relative mobility", side=2, line=2.5, cex=0.9,las=0, col="black")
box();

##### -----------------------2. susceptiblity -----------------------------------------------------------
 axis(4,at=seq(0,1.0,by=.2),labels=seq(0,1.0,by=.2),cex.axis=1.4,col="green",col.lab="green",las=2)
 mtext("Susceptibility", side=4, line=3.0, cex.lab=0.9,las=0, col="green") 
 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(SusD.quant[1,],rev(SusD.quant[3,])),col=mycol2,border=NA); 
lines(as.Date(DateDay),SusD.quant[2,],col="green",lty=5,type="b")
legend("bottomleft", legend=c("mid-point from Alpha to Delta", "mid-point from Delta to Omicron"),col=c("red", "purple"), lty = 1:2,bty="n", cex=1.0)


##### ---------------------6 weekly average OXford government policy data --
# axis(4,at=seq(0,1.0,by=.2),labels=seq(0,1.0,by=.2),cex.axis=1.4,col="blue",col.lab="blue",las=2)
# mtext("Gov response index", side=4, line=3.0, cex.lab=0.9,las=0, col="blue")
#lines(as.Date(DateDay),GovPolicy[1:length(DateDay)],col="blue",lty=5,type="b")

##### ----------------------- infections ---------daily incidence--------------------------------------------------
plot(as.Date(DateDay),INF.quant[2,],col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Daily number of infections",lty=1,type="b",ylim=c(0,130000)); #1.15*max(INF.quant[3,]))); 
#axis.Date(1, at=seq(min(as.Date(DateDay)), max(as.Date(DateDay)), by="1 mon"), format="%d-%m-%Y");
 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(INF.quant[1,],rev(INF.quant[3,])),col="grey",border=NA); 

abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);		
abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2)
   lines(as.Date(DateDay),INF.quant[2,],lty="solid",col="black",lwd=2.5)

CIN1=round(100*cumsum(INF.quant[1,])[Dayrun]/Npop,1)
CIN2=round(100*cumsum(INF.quant[2,])[Dayrun]/Npop,1)
CIN3=round(100*cumsum(INF.quant[3,])[Dayrun]/Npop,1)
legend("topleft",paste("d"),bty="n",col="dark",cex=1.6); 		#legend("topleft",paste("D) Incidence"),bty="n",cex=1.35); 
legend("top",c(paste("Prevelance: ",CIN2,"%(",CIN1,",",CIN3,"%)",sep=""),paste("Deaths: ",DeathTotal[2],"(",DeathTotal[1],",",DeathTotal[3],")",sep="")),bty="n",cex=1.2); 

axis.Date(1, at=seq(min(as.Date(DateDay)), max(as.Date(DateDay)), by="1 mon"), format="%d/%m/%Y",cex.axis=1.2); 
axis(2L,cex.axis=1.2);  box();  mtext("Daily infections", side=2, line=2.5, cex=0.9,las=0, col="black")

 mtext("Dates",side=1,outer=TRUE,line=2.2)
dev.off();



#ascertainment rate and breaking points
ASC_rate1=quantile(SampleASC1[S_USED],probs=c(0.025,0.5,0.975));  ASC_rate2= quantile(SampleASC2[S_USED],probs=c(0.025,0.5,0.975));
ASC_rate3=quantile(SampleASC3[S_USED],probs=c(0.025,0.5,0.975));  ASC_rate4= quantile(SampleASC4[S_USED],probs=c(0.025,0.5,0.975));
T_bk1 = ceiling(Samplet12[i_MLE])
T_bk2 = ceiling(Samplet23[i_MLE])
T_bk3 = ceiling(Samplet34[i_MLE])
ASC_RATE= matrix(0,ncol=Nday,nrow=3)

ASC_RATE[,        1:T_bk1]=ASC_rate1; 
ASC_RATE[,(1+T_bk1):T_bk2]=ASC_rate2; 
ASC_RATE[,(1+T_bk2):T_bk3]=ASC_rate3; 
ASC_RATE[,(1+T_bk3):Nday]=ASC_rate4; 


pdf(paste("EpiFitFig3_",i_CMV,".pdf"));
 if(i_CMV!=3) {
   par(mfrow=c(3,2))
   Plot_samples(Inf_avoid[S_USED],Inf_avoid[i_MLE],"Infections avoided (%)")		#
   Plot_samples(Death_avoid[S_USED],Death_avoid[i_MLE],"Deaths avoided (%)")	#
 }

 par(mfcol=c(4,1),mar=numeric(4),oma=c(4,4,.5,5.5),mgp=c(2,0.6,0));   #par(mfrow=c(2,2))
##### ----------------------- infections ---------daily incidence--------------------------------------------------
plot(as.Date(DateDay),INF.quant[2,],col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Daily number of infections",lty=1,type="b",ylim=c(0,1.15*max(INF.quant[3,])));  
 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(INF.quant[1,],rev(INF.quant[3,])),col="grey",border=NA); 

abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);		
abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2) 
   lines(as.Date(DateDay),INF.quant[2,],lty="solid",col="black",lwd=2.5)

ipeak=which.max(INF.quant[2,1:median(T_d)])

legend("topleft","a",bty="n",col="dark",cex=1.6); 
#legend("topleft",paste("A) Infections  peak=",round(INF.quant[2,ipeak])," (",round(INF.quant[1,ipeak]),",",round(INF.quant[3,ipeak]),") on ",as.Date(DateDay[ipeak],format="%d/%m/%Y"),sep=""),bty="n",cex=1.35);
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
legend("topleft","b",bty="n",col="dark",cex=1.6);   
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
legend("topleft","c",bty="n",col="dark",cex=1.36);		
 mtext("Deaths", side=2, line=2.5, cex=0.9,las=0, col="black")
axis(2L,cex.axis=1.2); 


##### ----------------------- Recovery -----------------------------------------------------------
plot(as.Date(DateDay),RecoveryD.quant[2,],col="black",xlab="Dates",xaxt = "n",yaxt = "n",ylab="Daily no. of recoveries",lty=1,type="b",ylim=c(0,max(Recovery,RecoveryD.quant[3,]))); 
 polygon(c(as.Date(DateDay),rev(as.Date(DateDay))),c(RecoveryD.quant[1,],rev(RecoveryD.quant[3,])),col="grey",border=NA); 

abline(v=DateDay[round(t_start_delta+SampleT_brk_delta[i_MLE])],col="red",lty=1);		
abline(v=DateDay[round(t_start_omicron+SampleT_brk_omicron[i_MLE])],col="purple",lty=2);	

points(as.Date(DateDay[1:max(which(Recovery>0))]),Recovery[1:max(which(Recovery>0))],,cex=.8,pch=17,col="blue")
   lines(as.Date(DateDay),RecoveryD.quant[2,],lty="solid",col="black",lwd=2.5)
legend("topleft","d",bty="n",col="dark",cex=1.6);		
axis.Date(1, at=seq(min(as.Date(DateDay)), max(as.Date(DateDay)), by="1 mon"), format="%d/%m/%Y",cex.axis=1.2);
axis(2L,cex.axis=1.2); 
 mtext("Dates",side=1,cex.lab=1.25,outer=TRUE,line=2.2);  mtext("The recovered", side=2, line=2.5, cex=0.9,las=0, col="black")

dev.off();



