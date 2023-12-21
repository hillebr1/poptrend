library(tidyverse)
library(vegan)
library(reshape2)
library(lubridate)

data <- read_delim("data/PPKTcount_noHT_28092022.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE)
names(data)

# STEP 1: If your variable names do not conform to our standard (StationID, year, species, abu), rename and delete empty rows 
data<-data %>% rename(species = Species, abu = biovolume_l) 
# in our case StationID is okay, year needs to be calculated
data<-data[!is.na(data$abu),]
data<-data[!is.na(data$sampledate),]

# STEP 2: If sample date but not year is given, extract the year. Library lubridate allows to first change the dates to numeric values and then extract the year
data$date.numeric<-dmy(data$sampledate) 
# dmy is command if data are in format dd.mm.yyyy, ymd is for yyyy.mm.dd
data$year <- year(data$date.numeric) 
# extracts the year

# STEP 3: If your data can contain the same species more than once per sample (because different size classes or lifestages were distinguished) sum them
data<-data  %>%
  group_by(StationID, sampledate, year, species) %>%  
  summarise(abu= sum(abu, na.rm = T))
# obviously sampledate and STEP 4 can be ignored if only annual data exist

# STEP 4: Create annual means and the number of annual samples 
data<-data  %>%
  group_by(StationID, year, species) %>% 
  summarise(abu= mean(abu, na.rm = T),
            N.sample=length(abu))

# STEP 5: Further reductions
# if you need to exclude some StationID because they are not part of the WaddenSea or for other reasons
data<-data %>%
  filter(StationID %in% c("BOOMKDP", "DANTZGT", "DOOVBWT","GROOTGND","BOCHTVWTM",
                          "HUIBGOT", "MARSDND", "ROTTMPT3", "TERSLG10",
                          "Nney_W_2", "Bork_W_1", "WeMu_W_1", "JaBu_W_1"))

# In other data sets it might be needed to exclude certain species or years.

# Our example is already in long format, if your data is in wide format consider the "melt" command as in mentioned here, remove inital # to run
#data <- melt(data, 
#                 id.vars=c("StationID", "year"), # all the variables to keep but not split apart
#                 measure.vars=c(x:y), # number of the variables containing species, that shall be melted into one
#                 variable.name="species", # name for new variable that contains the names of the melded variables
#                 value.name="abu") # name of the new variables with the values

 
length.spec <-data %>% 
  group_by(StationID, species) %>% 
  summarise(N.spec=length(unique(year)))

length.stat <-data %>% 
  group_by(StationID) %>% 
  summarise(N.stat=length(unique(year)))
length<-merge(length.spec,length.stat,by="StationID")
length$threshold<-length$N.spec/length$N.stat
data<-merge(data,length,by=c("StationID","species"))

data<-data[data$N.spec>4,]
data<-data[data$threshold>.75,]

#get slope b of LN(abu) ~ year from linear regression

results.logabu = data %>% 
  group_by(StationID, species) %>% 
  do(tidy(lm(log(abu)~year, data = .)))%>% 
  filter(term == "year")
results.logabu<-results.logabu %>% rename(b.linear = estimate,
                                          SE.b.linear = std.error,
                                          pval.lin = p.value)

results.logabu<-results.logabu[,c("StationID","species","b.linear","SE.b.linear", "pval.lin")]
glance.logabu = data %>% 
  group_by(StationID, species) %>% 
  do(glance(lm(log(abu)~year, data = .))) 
glance.logabu<-glance.logabu %>% rename(AIC.linear = AIC,
                                        nobs.linear = nobs,
                                        r2.lin = adj.r.squared)

glance.logabu<-glance.logabu[,c("StationID","species","AIC.linear","nobs.linear","r2.lin")]
results.logabu<-merge(results.logabu,glance.logabu, by=c("StationID","species"))

results.year = data %>% 
  group_by(StationID, species) %>% 
  summarize(min.year=min(year),
            max.year=max(year))

results.logabu<-merge(results.logabu,results.year, by=c("StationID","species"))

#get slope b of LN(abu) ~ poly (year,2) from quadratic regression

results.logabu2.lin = data %>% 
  group_by(StationID,species) %>% 
  do(tidy(lm(log(abu)~poly(year,2), data = .)))%>% 
  filter(term == "poly(year, 2)1") 
results.logabu2.lin<-results.logabu2.lin %>% rename(b.poly = estimate,
                                                    SE.b.poly = std.error,
                                                    pval.b.poly = p.value)

results.logabu2.lin<-results.logabu2.lin[,c("StationID","species","b.poly","SE.b.poly", "pval.b.poly")]

#get slope c of LN(abu) ~ poly (year,2) from quadratic regression
results.logabu2.quad = data %>% 
  group_by(StationID,species) %>% 
  do(tidy(lm(log(abu)~poly(year,2), data = .)))%>% 
  filter(term == "poly(year, 2)2") 

results.logabu2.quad<-results.logabu2.quad %>% rename(c.poly = estimate,
                                                      SE.c.poly = std.error,
                                                      pval.c.poly = p.value)


results.logabu2.quad<-results.logabu2.quad[,c("StationID","species","c.poly","SE.c.poly", "pval.c.poly")]

glance.logabu2 = data %>% 
  group_by(StationID,species) %>% 
  do(glance(lm(log(abu)~poly(year,2), data = .)))
glance.logabu2<-glance.logabu2 %>% rename(AIC.poly = AIC,
                                          nobs.poly = nobs,
                                          r2.poly = adj.r.squared,
                                          pval.poly = p.value)

glance.logabu2<-glance.logabu2[,c("StationID","species","AIC.poly","nobs.poly","r2.poly","pval.poly")]

results.logabu2<-merge(results.logabu2.lin,results.logabu2.quad, by=c("StationID","species"))
results.logabu2<-merge(results.logabu2,glance.logabu2, by=c("StationID","species"))


#merge linear and quadratic regression
output<-merge(results.logabu,results.logabu2, by=c("StationID","species"))

 
#create the votes according to Figure 2

output$vote<-NA
# none of the regressions is significant
output$vote[output$pval.lin>=0.1 & output$pval.poly>=0.1]<-"no"
# only linear is significant
output$vote[output$pval.lin<0.1 & output$pval.poly>0.1 &
              output$b.linear<0]<-"negative linear"
output$vote[output$pval.lin<0.1 & output$pval.poly>0.1 &
              output$b.linear>0]<-"positive linear"

# only quadratic is significant
output$vote[output$pval.lin>0.1 & output$pval.poly<0.1 &
              output$b.poly<0 & output$c.poly<0]<-"negative accelerating"
output$vote[output$pval.lin>0.1 & output$pval.poly<0.1 &
              output$b.poly>0 & output$c.poly>0]<-"positive accelerating"
output$vote[output$pval.lin>0.1 & output$pval.poly<0.1 &
              output$b.poly<0 & output$c.poly>0]<-"negative to positive"
output$vote[output$pval.lin>0.1 & output$pval.poly<0.1 &
              output$b.poly>0 & output$c.poly<0]<-"positive to negative"

#both are significant, requiring involving AIC, first if AIC is lower for the polynomial
output$vote[output$pval.lin<0.1 & output$pval.poly<0.1 &
              output$AIC.linear>output$AIC.poly & output$b.poly<0 & output$c.poly<0]<-"negative accelerating"
output$vote[output$pval.lin<0.1 & output$pval.poly<0.1 &
              output$AIC.linear>output$AIC.poly & output$b.poly>0 & output$c.poly>0]<-"positive accelerating"
output$vote[output$pval.lin<0.1 & output$pval.poly<0.1 &
              output$AIC.linear>output$AIC.poly & output$b.poly<0 & output$c.poly>0]<-"negative to positive"
output$vote[output$pval.lin<0.1 & output$pval.poly<0.1 &
              output$AIC.linear>output$AIC.poly & output$b.poly>0 & output$c.poly<0]<-"positive to negative"
# same but when AIC is lower for linear regression
output$vote[output$pval.lin<0.1 & output$pval.poly<0.1 &
              output$AIC.linear<output$AIC.poly & output$b.linear<0]<-"negative linear"
output$vote[output$pval.lin<0.1 & output$pval.poly<0.1 &
              output$AIC.linear<output$AIC.poly & output$b.linear>0]<-"positive linear"
output$vote[output$pval.lin<0.1 & is.na(output$pval.poly) &
              output$b.linear<0]<-"negative linear"
output$vote[output$pval.lin<0.1 & is.na(output$pval.poly) &
              output$b.linear>0]<-"positive linear"

#sometimes, the polynomial or even the linear do not produce results, such that the p-value is NA, make these neutral trends as well
output$vote[is.na(output$vote) & output$pval.lin>0.1 & is.na(output$pval.poly)]<-"no"
output$vote[is.na(output$vote) & is.na(output$pval.lin)] <-"no"

 
#create a unique case identifier
output$UCI<-paste(output$StationID,output$species) 
data$UCI<-paste(data$StationID,data$species) 
UCI<-unique(output$UCI[output$vote=="negative to positive"|
                         output$vote=="positive to negative"])

# as there might be cases when there is no such case, i do the loop as in if-else clause
if(length(UCI)>0){
  MOS.logabu<-data.frame()
  for(i in 1:length(UCI)){
    temp<-data[data$UCI==UCI[i], ]
    if(dim(temp)[1]>2){
      MOS<-MOStest(temp$year,log(temp$abu))
      MOSout<-MOS$isBracketed
      MOS.logabu<-rbind(MOS.logabu,data.frame(UCI=temp$UCI[1],MOSout))
      rm(temp)
    }
  }
} else{
  MOS.logabu<-data.frame(matrix("", ncol = 1, nrow = dim(output)[1]))
  MOS.logabu$UCI<-output$UCI
  MOS.logabu$MOSout<-NA}


output<-merge(output,MOS.logabu,by="UCI",all=T)

# change votes if MOS Test was FALSE

output$vote[output$vote=="positive to negative" &
              output$MOSout=="FALSE"]<-"positive decelerating"
output$vote[output$vote=="negative to positive" &
              output$MOSout=="FALSE"]<-"negative decelerating"

summary(as.factor(output$vote))
