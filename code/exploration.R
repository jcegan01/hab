#Author: Kaitlin Volk; with edits by Jeff Cegan
#Date: May 19, 2021
#Project HAB with Jodi Ryder

library(dplyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(grid)
library(tidyr)

###----Data exploration----
##====Load and format data====
load("data/reservoir_wq.rda") #the text file of the results tab from the Louisville DASLER

names(dat) <- c("loc_id", "sample_num", "lab_id","lab_snum", "storet_num", "units",   
                "value", "text_value", "detect","pquant", "pre_recov", "post_recov",
                "dilution", "lab_anlz", "dcl_flag", "date_entered", "revision", "val_qual", #I exported the text file without headers
                "val_id", "export_flag", "act_value","test_methd", "prep_mthd")             #so I'm adding them in manually here

target.param <- c("temp_c","DO_mgl","DO_perc","chlA_tri","chlA_spect",
                  "P_total","N_total","flow","secchi_in","TSS",           #A list of the priority parameters to be used when 
                  "TDS", "ammonia","nitrate","phosphate", "C_total")      #filtering in the next dataframe

dat2 <- dat %>%   
  mutate(param = recode(storet_num, "00010"= "temp_c",  #Replace the storet codes of the priority params with descriptive text
                                        "00299" = "DO_mgl","00301" = "DO_perc",
                                        "32210" = "chlA_tri","32211" = "chlA_spect",
                                        "00665" = "P_total","00600" = "N_total",
                                        "74073" = "flow","00077" = "secchi_in",
                                        "00515" = "TSS","00530" = "TDS",
                                        "00610" = "ammonia","00620" = "nitrate",
                                        "00650" = "phosphate","00680" = "C_total"),
         sample_num=recode(sample_num, "210808162000000.0" = "201808162000000.0"),  #fix clearly wrong date
         yr = str_sub(sample_num,1,4),  #Year column
         mo=str_sub(sample_num,5,6),    #Month column
         dy = str_sub(sample_num,7,8),  #day column
         sample_date=make_date(yr, mo, dy), #full date column
         depth = as.numeric(str_sub(sample_num,-5,-1)), #depth of sample
         reservoir = str_sub(loc_id,1,3)) %>% #reservior of sample
  select(c(1,30,2,28, 25, 26,29,24,7,6)) %>% #extract just the location, date, depth, and paramater data
  filter(param %in% target.param, !is.na(value)) #extract rows with the priority parameters only and remove NA's


##====Explore ChlA data==== 
#ChlAis an important variable for identifying HABs and for the correlations in 1.6 of the SOW
#Two methods have been used: trichromatic and Spectrophotomatic. Is one better? Should they both be used?
dat.chla <- dat2 %>% filter(param == "chlA_tri"| param =="chlA_spect",!is.na(value)) #pull just the chla dat
ggplot(dat.chla, aes(x=sample_date, y=value, color=param)) + geom_line() +  #graph chla over time for the two methods
  ylim(0,200) + facet_wrap( ~ reservoir, nrow = 6)                          #facet by reservior to see differences

dat.chla.dif <- dat.chla %>%    
  pivot_wider(names_from=param, values_from=value, values_fn=mean) %>% 
  select(1,3,9,10) %>%  mutate(difference=chlA_tri-chlA_spect) %>%    #calculate difference in methods for samples with values for both methods
  filter(!is.na(difference)) %>% arrange(-difference)                   #2% of values are in disagreement

#Given the small percentage of differences, combine the two chlA columns into 1
dat3 <- dat2 %>% mutate(param = recode(param, "chlA_tri" = "chlA", "chlA_spect" = "chlA"))
dat3w <- dat3 %>% select(-10) %>% pivot_wider(names_from= param, values_from = value, values_fn=mean) 
rm(dat2) #remove previous version

  
##====Look for outliers/weirdness in the variables====
ggplot(dat3w) + geom_qq(aes(sample=temp_c)) #exclude below -2 and above 40
ggplot(dat3w) + geom_qq(aes(sample=DO_mgl)) #eclude below 0 and above 50
ggplot(dat3w) + geom_qq(aes(sample=flow))   #one outlier at over 13000 cfu (next highest 5388)
ggplot(dat3w) + geom_qq(aes(sample=DO_perc)) #excle above 300
ggplot(dat3w) + geom_qq(aes(sample=TSS)) #exclude above 2000
ggplot(dat3w) + geom_qq(aes(sample=C_total)) #exclude above 100
ggplot(dat3w) + geom_qq(aes(sample=P_total)) # gonna need to mess around with this b/c of units
ggplot(dat3w) + geom_qq(aes(sample=depth)) #exclude above 200 
ggplot(dat3w) + geom_qq(aes(sample=depth)) + facet_wrap(~reservoir, nrow=6) #depth by reservoir confirm that I should exclude depths over 200
#Tested and looks passable: chla, phoshate (only 7 values), nitrate, N_total, secchi_in,TDS, ammonia

#Dive into total phosphorous (P_total)
dat.p <- dat3 %>% filter(param=="P_total") %>%                               #Looking at P data only
  mutate(units=recode(units, "MG/L" = "mg/L", "mg/l" = "mg/L", "pH"="ug/L")) #simplify variations so you only have ug/L and mg/L

dat.p %>% filter(value<1) %>% ggplot(aes(x=sample_date, y=value, color=units)) + 
  geom_point() + facet_wrap(~reservoir, nrow=6)  #Look at P by reservior and units for values under 1
dat.p %>% filter(value<1000 & value >1) %>% ggplot(aes(x=sample_date, y=value, color=units)) + 
  geom_point() + facet_wrap(~reservoir, nrow=6) #Look at P by reservior and units for values over 1 but below 1000 
dat.p %>% filter(value>1000) %>% ggplot(aes(x=sample_date, y=value, color=units)) + 
  geom_point() + facet_wrap(~reservoir, nrow=6) #Look at P by reservior and for values over 1000 (the rounded upper limit in my brief lit search)
#If x is less than 1, convert to ug/l by multiplying value by 1000 
#if x > 1 but < 1000, recode mg/l to be uq/l
#if x > 1000, exclude regardless of unit

#Clean the data
fine <- c("chlA", "phosphate", "nitrate","N_total", "secchi_in", "TDS", "ammonia")   #parameters that won't be cleaned

dat4 <- dat3 %>% filter((param=="temp_c" & value <40 & value> -2)|      #exclude outlying data that
                          (param=="DO_mgl" & value>0 & value <50)|      #we determined above from visualiztion
                          (param=="flow" & value <6000)|                            
                          (param=="DO_perc" & value<300)|                           
                          (param=="TSS" & value < 2000)|                             
                          (param=="C_total" & value <100)|                          
                          (param=="P_total" & value <1000)|                       
                          (param%in%fine), 
                        depth <= 200) %>% 
  mutate(value=case_when(param=="P_total" & value <1 ~ value*1000, TRUE~value), #complete proper conversion/recoding for P data
         units=case_when(param=="P_total"~ "ug/L", TRUE~as.character(units)))
rm(dat3) #remove previous version

dat4w <- dat4%>% select(-10) %>% pivot_wider(names_from= param, values_from = value, values_fn=mean)
rm(dat3w) #remove previous version

#Visualize the cleaned data
ggplot(dat4w) + geom_qq(aes(sample=temp_c))
ggplot(dat4w) + geom_qq(aes(sample=DO_mgl)) 
ggplot(dat4w) + geom_qq(aes(sample=flow))   
ggplot(dat4w) + geom_qq(aes(sample=DO_perc)) 
ggplot(dat4w) + geom_qq(aes(sample=TSS)) 
ggplot(dat4w) + geom_qq(aes(sample=C_total)) 
ggplot(dat4w) + geom_qq(aes(sample=P_total))
ggplot(dat4w) + geom_qq(aes(sample=depth)) + facet_wrap(~reservoir, nrow=6)

#Curious about how many reserviors and sites per reservoir there are:
reservoir.summary <- dat4 %>% group_by(reservoir) %>% summarize(n.sites = n_distinct(loc_id))  #Get counts for the number of samples sites per reservoir


##==== 1.2 Time series plots for a single location====
#They're pretty messy, especially temp and DO if we don't control for depth.
dat4 %>% filter(loc_id=="TAR20001") %>%             #choose one sample site
  ggplot(aes(x=sample_date,y=value, color=depth)) + #graph values over time, with color displayign depth
  geom_point() + geom_smooth(method = "loess") +    #add loess curve to look at trend
  facet_wrap(~param,nrow=4,scales = "free_y")       #facet by paramter

#What if we do data for a single reservoir instead of a single site? -> Messier
dat4 %>% filter(reservoir=="TAR") %>%               #Choose single reservoir instead of single site
  ggplot(aes(x=sample_date,y=value, color=depth)) + #do the same as above 
  geom_point() + geom_smooth(method = "loess") + facet_wrap(~param,nrow=4,scales = "free_y")
#Surface only
test.graph <-dat4 %>% filter(reservoir=="TAR" & depth==0) %>% #use only data where the depth=0 (surface) 
  ggplot(aes(x=sample_date,y=value)) + geom_point() + geom_smooth(method = "loess") + 
  facet_wrap(~param,nrow=4,scales = "free_y")
test.graph 

#I'm thinking we could have a default (like all data color coded for depth or just surface) 
#and then the ability for the users to select/input a depth or depth range

#The SOW wants a label of how many data points are displayed on the graph. I make a table 
#of these values below but then the code I looked up to add it to a facet with different y axes
#doesn't work and I don't know why
site.nsamples <- dat4 %>% group_by(loc_id, param) %>% summarize(count=n()) %>% ungroup() %>% #Get counts for the number of observations per param per site
  pivot_wider(names_from = param, values_from=count) %>% mutate_all(~replace(., is.na(.), 0)) #give each param a column and replace na's with 0's

reservoir.nsamples<- dat4 %>%  group_by(reservoir, param) %>% summarize(count=n(),   #Get counts for the number of observations per param per reservoir
                                                                        ystar=(0.8 * max(value)))  #y placement for graph (doesn't seem to work though) 

test.graph + geom_text(data=reservoir.nsamples,   #counts from the reservoir.nsamples table
                       aes(y=ystar, x=ymd(2010-01-01), label=paste("N=", count)), #placement and content of label
                       colour="red", inherit.aes=FALSE, parse=TRUE)         #Doesn't work...

#Also adding # samples for each parameter to the reserveroir summary table so it's all in one place
reservoir.summary2 <- reservoir.nsamples %>% select(-4) %>% pivot_wider(names_from = param, values_from=count) %>% #remove graphing column and spread 
  mutate_all(~replace(., is.na(.), 0)) %>% right_join(reservoir.summary, by="reservoir") %>% #replace NAs with 0's
  relocate(n.sites, .after=reservoir) #put number of sites before # of param samples
rm(reservoir.summary) #remove previous version


##==== 1.3 Characterizing reservoirs based on last five years of summer months (may-sept)====
#             Oligo   Meso    Eutro
#TP (ug/L)    8       27      84   
#TN (mg/L)    .66     .75     1.9
#ChlA (ug/L)  1.7     4.7     14

trophic.param <- c("P_total", "N_total","chlA")   #Paramters that we'll use to determine trophic state
summer.months <- c("05", "06", "07", "08", "09")  #Summer months as specified in SOW (May-Sept)

trophicstates <- dat4 %>% 
  filter(param %in% trophic.param & mo %in% summer.months & yr > 2013) %>% #keep only trophic parameters from the summer for the last 5 years
  pivot_wider(names_from= param, values_from = value,values_fn=mean) %>%   #go wide
  group_by(reservoir) %>% 
  summarize(P_total = mean(P_total, na.rm=TRUE),                           #mean total phosphorous for each reservoir
            N_total = mean(N_total, na.rm=TRUE),                           #mean total nitrogen for each reservoir
            chlA = mean(chlA, na.rm=TRUE)) %>%                             #mean chlA for each reservoir
  mutate(troph_P = as.numeric(case_when(P_total < 27 ~ "1",         #use thresholds to determine trophic state based on P
                             P_total >= 27 & P_total < 84 ~ "2",     #1= Oligo, 2= Meso, 3=Eutro
                             P_total > 84 ~ "3",
                             TRUE ~ NA_character_)),
         troph_N = as.numeric(case_when(N_total < .75 ~ "1",         #use thresholds to determine trophic state based on N
                             N_total >= .75 & N_total < 1.9 ~ "2",
                             N_total >= 1.9 ~ "3",
                             TRUE ~ NA_character_)),
         troph_chla= as.numeric(case_when(chlA < 4.7 ~ "1",            #use thresholds to determine trophic state based on chlA
                                 chlA >= 4.7 & chlA < 14 ~ "2",
                                 chlA >= 14 ~ "3",
                                 TRUE ~ NA_character_))) %>% 
  rowwise() %>% mutate(troph_combo =(median(c(troph_P, troph_N, troph_chla), na.rm = T))) %>% #Find th median trophic state
  ungroup %>% mutate(trophic_state= case_when(troph_combo==1 ~"Oligotrophic",   #Recode median trophic state to words for final output
                                              troph_combo==1.5 ~ "Oligo-Mesotrophic",
                                              troph_combo==2 ~ "Mesotrophic",
                                              troph_combo==2.5 ~ "Meso-Eutrotrophic",
                                              troph_combo==3 ~ "Eutrotrophic"))

#Q: Do we want to take the mean or the median of the paramters to base the trophic states on?

#Q: What do we want to do when one parameters disagree on trophic state (i.e. it's eutro for P
# but meso for N and chlA)? We coud take the median or use hyphens (e.g. meso-eutro) 
# I did a combo of both above. 

#Add trophic state to the reservoir summary table with a left join since not all reservoirs have data for a trophic state
reservoir.summary3 <- left_join(reservoir.summary2, trophicstates, by="reservoir") %>% 
  select(-(17:23)) %>% relocate(trophic_state, .after = n.sites) 
rm(reservoir.summary2) #remove previous iteration; reminder, the values for the different parameters are counts of samples


##==== 1.4 Assign reservoirs to mixing regime====


##==== 1.5 Apply Mann-Kendall test to time series to find increasingtrends====


##==== 1.6 Cross correlation of time series (esp. chla, our bloom indicator )

