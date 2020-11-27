###########################################
#   Behavior of precip before landslides  #
#       Pauline Rivoire 23.11.2020        #
###########################################
library(lubridate)

excel_landslides <- readxl::read_xlsx(path = "../PrecipClustering_BottomUpAnalysis_LOCAL/Data/Landslides_Date.xlsx")

dates_landslides <- as.matrix(excel_landslides)

#Load the data, change the path if necessary
excel_precip<-readxl::read_xlsx(path = "Data/Precipitation_data_1950_2008_IB02.xlsx")

#Name columns properly
colnames(excel_precip) <- c("year","month","day","gridpoint_W","gridpoint_E")

#Convert data of the gridpoint we want to a vector
precip_GP_West <- as.numeric(as.matrix(excel_precip[-c(1,2,21277),"gridpoint_W"]))
date_matrix <- as.matrix(excel_precip[-c(1,2,21277),c("year","month","day")])

ref_date_ls <- numeric()
for (ls in 1:length(dates_landslides)) {
  ref_date_ls[ls] <- which(date_matrix[,1]==year(dates_landslides[ls]) &
                             date_matrix[,2]==month(dates_landslides[ls]) &
                             date_matrix[,3]==day(dates_landslides[ls]))
}#end for ls



##### Mean number events before landslide #####


mean_number_events_before_landslide <- matrix(data=NA, nrow = length((15:90)), ncol = 1) #time window length: 15:90
rownames(mean_number_events_before_landslide) <- as.character((15:90))
quantile_75 <- 10.5
for (t_w in 1:length((15:90))) {
  w_length <- (15:90)[t_w]
  
  precip_w_before <- matrix(data = NA, nrow = w_length, ncol = length(ref_date_ls))
  events_w_before <- matrix(data = NA, nrow = w_length, ncol = length(ref_date_ls))
  
  for (ls in 1:length(dates_landslides)) {
    precip_w_before[,ls] <- precip_GP_West[(ref_date_ls[ls]-w_length+1):ref_date_ls[ls]]
  }#end for ls
  
  for (ls in 1:length(dates_landslides)) {
    gt_1 <- as.numeric(precip_w_before[,ls]>quantile_75)
    gt_1_nocorr <- gt_1
    for (day in 2:w_length) {
      if(gt_1_nocorr[day-1] == 1 & gt_1_nocorr[day] == 1){gt_1_nocorr[day] <- 0}
    }#end for day
    events_w_before[,ls] <- gt_1_nocorr
  }#end for ls
  
  nb_events_w <- apply(X=events_w_before, MARGIN = 2, FUN = sum)
  
  mean_number_events_before_landslide[t_w] <- mean(nb_events_w)
}#end for t_w




#################### Precip time series : compute events ####################


#we exclude the two first rows and the last one because they don't contain data
All_years <- unique(date_matrix[,"year"])
years_fullwinter <- All_years[-c(1,length(All_years))] #Exclude first and last winter because not complete

nrow(date_matrix)==length(precip_GP_West)#check if same dimension

quantile_75 <- 10.5 #/!\ seas quantile, Nov to March (not October)

#extract seasonal precip
precip_extended_winter <- precip_GP_West
days_to_exclude <- which(date_matrix[,"month"] %in% 4:10) #exclude days from April to October included
precip_extended_winter[days_to_exclude] <- rep(0, length(days_to_exclude))
#set those excluded days to 0 (can also be set to NA if needed)

cbind(date_matrix[,"month"], precip_GP_West, precip_extended_winter)

compute_events <- function(precip_vector, threshold){
  n_days <- length(precip_vector)
  events_vector <- numeric(length = n_days)
  events_vector[1] <- as.numeric(precip_vector[1]>=quantile_75)
  
  for (day in 2:n_days) {
    if(events_vector[day-1]==0 & precip_vector[day]>=quantile_75){
      events_vector[day] <- 1
    }#end if
  }#end for day
  return(events_vector)
}#end compute_events()


Events_75 <- compute_events(precip_vector = precip_extended_winter, threshold = quantile_75)


precip_per_winter <- function(precip_events, days, months, years, year_to_extract){
  days_to_extract <- c(which(years == year_to_extract & months%in%11:12),
                       which(years == (year_to_extract+1) & months%in%1:3))
  out <- cbind(paste0(as.character(years[days_to_extract]),"-", as.character(months[days_to_extract]),
                      "-",as.character(days[days_to_extract])),
               precip_events[days_to_extract])
  colnames(out) <- c("date", "precip")
  return(out)
}#end precip_per_winter_30() 

list_precip_winter <- list()
count <- 0
for (year in years_fullwinter) {
  count <- count + 1
  precip_winter <- as.matrix(as.numeric(precip_per_winter(year_to_extract = year,
                                                             precip_events = Events_75,
                                                             days = date_matrix[,"day"], months = date_matrix[,"month"],
                                                             years = date_matrix[,"year"])[,"precip"]))
  colnames(precip_winter) <- paste("winter Nov.",as.character(year))
  list_precip_winter[[count]] <- precip_winter
}#end for year
rm(count)




########## simu of mean number of events: negbinom and subsamples ########## 
nsimu <- 1000




##### Mean number of events following negbinom in 30-60-90 in winter #####
load("Fit/Fit_negative_binomial.Rdata")

nyears_per_simu <- length(dates_landslides)
MC_mean_nb_ev_nbinom <- matrix(data = NA, nrow = nsimu, ncol = 3)
colnames(MC_mean_nb_ev_nbinom) <- c("30", "60", "90")
for (simu in 1:nsimu) {
  MC_mean_nb_ev_nbinom[simu,"30"] <- mean(rnbinom(n=nyears_per_simu, size = Fit_neg_bino$fit_w30$estimate["size"],
                                                  mu = Fit_neg_bino$fit_w30$estimate["mu"]))
  MC_mean_nb_ev_nbinom[simu,"60"] <- mean(rnbinom(n=nyears_per_simu, size = Fit_neg_bino$fit_w60$estimate["size"],
                                                  mu = Fit_neg_bino$fit_w60$estimate["mu"]))
  MC_mean_nb_ev_nbinom[simu,"90"] <- mean(rnbinom(n=nyears_per_simu, size = Fit_neg_bino$fit_w90$estimate["size"],
                                                  mu = Fit_neg_bino$fit_w90$estimate["mu"]))
  
}#end for simu
boxplot(MC_mean_nb_ev_nbinom)

hist(MC_mean_nb_ev_nbinom[,"30"], main = "Monte-Carlo simulation (1000 repet.)\n Mean nb events in 30d wind. following neg.binom")

text(x=3.5, y=150, "before landslide")
text(x=3.2, y=130, "Mean=")
text(x=3.6, y=130, mean_number_events_before_landslide["30",], col = "red")


hist(MC_mean_nb_ev_nbinom[,"60"], main = "Monte-Carlo simulation (1000 repet.)\n Mean nb events in 60d wind. following neg.binom")

text(x=5.5, y=250, "before landslide")
text(x=5.2, y=230, "Mean=")
text(x=5.6, y=230, round(mean_number_events_before_landslide["60",],1), col = "red")


hist(MC_mean_nb_ev_nbinom[,"90"], main = "Monte-Carlo simulation (1000 repet.)\n Mean nb events in 90d wind. following neg.binom")

text(x=8.5, y=250, "before landslide")
text(x=8.2, y=230, "Mean=")
text(x=9, y=230, round(mean_number_events_before_landslide["90",],1), col = "red")

##### Mean number of events per t_w in winter #####
#function to extract the number of events in a random window of fixed size for a fix year
nb_events_random_window <- function(list_precip_winter, year, years_fullwinter, size_window){ #list_precip_winter = 0s and 1s
  year_ind <- which(years_fullwinter==year)
  rdm_start <- sample(1:(length(list_precip_winter[[year_ind]])-size_window+1),size=1)
  return(sum(list_precip_winter[[year_ind]][rdm_start:(rdm_start+size_window-1)]))
}#end nb_events_random_window()

#function to obtain a sample of the number of events per window of a fixed size, chosen randomly within all the years
sample_nb_events_per_window <- function(list_precip_winter, years_fullwinter, size_window, sample_size){
  rdm_years <- as.matrix(sample(years_fullwinter, size = sample_size, replace = T))

  out <- apply(X=rdm_years, MARGIN = 1, FUN = nb_events_random_window,
               list_precip_winter = list_precip_winter, years_fullwinter =years_fullwinter, size_window = size_window)
  return(out)
}#end sample_nb_events_per_window()

quantile_75 <- 10.5
nyears_per_simu <- length(dates_landslides)
MC_mean_nb_ev_subsamp <- matrix(data = NA, nrow = nsimu, ncol = length((15:90)))
colnames(MC_mean_nb_ev_subsamp) <- as.character((15:90))

for (simu in 1:nsimu) {
  
  for (t_w in 1:length((15:90))) {
    w_length <- (15:90)[t_w]
    
    
    MC_mean_nb_ev_subsamp[simu,t_w] <- mean(sample_nb_events_per_window(list_precip_winter = list_precip_winter,
                                                                        years_fullwinter = years_fullwinter,
                                                                        size_window = w_length,
                                                                        sample_size = nyears_per_simu))
  }#end for t_w
  
} #end for simu
boxplot(MC_mean_nb_ev_subsamp[,c("30", "60", "90")])
low_quantile <- apply(X=MC_mean_nb_ev_subsamp, MARGIN = 2,FUN = quantile, probs=0.025)
up_quantile <- apply(X=MC_mean_nb_ev_subsamp, MARGIN = 2,FUN = quantile, probs=0.975)

plot(15:90,mean_number_events_before_landslide, pch="+", xlab="window length", ylim = c(0,12),
     ylab = "mean number of events", main = "Bootstrap number of events per time window\n(1000 repet., 22 winds chosen randomly in all winters)")
polygon(c(15:90,rev(15:90)),c(up_quantile,rev(low_quantile)),col="thistle",border=NA)
legend("bottomright", legend = c("before a landslide", "95%CI"), pch=c(3,15), col=c("black", "thistle"))





########## Ripley's K function: negbinom and subsamples ########## 

source("ripleyk.R")

Ripley_k_tw <- matrix(data = NA, nrow = length(15:90))
rownames(Ripley_k_tw)=as.character(15:90)
for (ref_tw in 1:length(15:90)) {
  tw <- (15:90)[ref_tw]
  Ripley_k_tw[ref_tw] <- RipleyK(Listvectors = list_precip_winter, size = length(list_precip_winter), xvect = tw)[1]
}#end for tw


plot(15:90, Ripley_k_tw, ylab = "", xlab="window length")
points(15:90, mean_number_events_before_landslide, pch="+")
legend("bottomright", legend = c("mean nb events before landslide", "Ripley's K November-March"), pch=c(1,3))



#mydata <- list()
#for(i in 1:57) mydata[[i]] <- rbinom(n=30*5,size=1,prob=(sum(unlist(list_precip_winter))/length(unlist(list_precip_winter))))
#Listvectors <- mydata
#size <- 57 #number of year

xvect <- (15:90)[(15:90)%%2==0]

y <- Ripley_k_tw

year_sizes <- unlist(lapply(X=list_precip_winter, FUN=length))

obj <- SimuRipleyK(nsimu=100,xvect=xvect,size=57,yearsize=year_sizes,
                   n.events.avg=sum(unlist(list_precip_winter)))

plot(x=xvect,y=y[(15:90)%%2==0],type="l",ylab="Ripley",lwd=3,ylim=c(2,max(y)),xlab="time window length (days)")
points(xvect, mean_number_events_before_landslide[(15:90)%%2==0], pch="+")
lines(x=xvect,y=obj$avg,col="red")
lines(x=xvect,y=obj$lower,col="grey",lty=2)
lines(x=xvect,y=obj$upper,col="grey",lty=2)
legend("bottomright", legend = c("mean nb events before landslide",
                                 "Ripley's K November-March",
                                 "mean Ripley's K binomial distr.",
                                 "90% CI Ripley's K binomial distr."),
       pch=c("+", NA, NA, NA), lty = c(NA,1,1,2), col=c("black","black", "red", "grey"), lwd=c(NA, 2, 1, 1))
