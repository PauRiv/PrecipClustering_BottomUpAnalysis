###########################################
#   Behavior of precip before landslides  #
#       Pauline Rivoire 23.11.2020        #
###########################################
library(lubridate)

excel_landslides <- readxl::read_xlsx(path = "../PrecipClustering_BottomUpAnalysis_LOCAL/Data/Landslides_Date.xlsx")

dates_landslides <- as.matrix(excel_landslides)

dates_landslides <- rbind("1958-12-19", dates_landslides)
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


mean_number_events_before_landslide <- matrix(data=NA, nrow = length((2:90)), ncol = 1) #time window length: 2:90
number_events_before_landslide <- matrix(data=NA, nrow = length((2:90)), ncol = 23)
rownames(mean_number_events_before_landslide) <- as.character((2:90))
quantile_75 <- 10.5
acc_precip_before_ls <- matrix(nrow = length(2:90), ncol = length(dates_landslides))
for (t_w in 1:length((2:90))) {
  w_length <- (2:90)[t_w]
  
  precip_w_before <- matrix(data = NA, nrow = w_length, ncol = length(ref_date_ls))
  events_w_before <- matrix(data = NA, nrow = w_length, ncol = length(ref_date_ls))
  
  for (ls in 1:length(dates_landslides)) {
    precip_w_before[,ls] <- precip_GP_West[(ref_date_ls[ls]-w_length+1):ref_date_ls[ls]]
    acc_precip_before_ls[t_w,ls] <- sum(precip_w_before[,ls])
  }#end for ls
  
  for (ls in 1:length(dates_landslides)) {
    gt_1 <- as.numeric(precip_w_before[,ls]>quantile_75)
    gt_1_nocorr <- gt_1
    for (day in 2:w_length) {
      if(gt_1_nocorr[day-1] == 1 & gt_1_nocorr[day] == 1){gt_1_nocorr[day] <- 0}
    }#end for day
    events_w_before[,ls] <- gt_1_nocorr
    number_events_before_landslide[t_w,ls] <- sum(gt_1_nocorr)
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
list_precip_winter_RAW <- list()
count <- 0
for (year in years_fullwinter) {
  count <- count + 1
  precip_winter <- as.matrix(as.numeric(precip_per_winter(year_to_extract = year,
                                                          precip_events = Events_75,
                                                          days = date_matrix[,"day"], months = date_matrix[,"month"],
                                                          years = date_matrix[,"year"])[,"precip"]))
  colnames(precip_winter) <- paste("winter Nov.",as.character(year))
  list_precip_winter[[count]] <- precip_winter
  list_precip_winter_RAW[[count]] <- as.matrix(as.numeric(precip_per_winter(year_to_extract = year,
                                                                            precip_events = precip_extended_winter,
                                                                            days = date_matrix[,"day"], months = date_matrix[,"month"],
                                                                            years = date_matrix[,"year"])[,"precip"]))
  colnames(list_precip_winter_RAW[[count]]) <- paste("winter Nov.",as.character(year))
}#end for year
rm(count)




########## simu of mean number of events: subsamples ##########
nsimu <- 1000

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
MC_mean_nb_ev_subsamp <- matrix(data = NA, nrow = nsimu, ncol = length((2:90)))
colnames(MC_mean_nb_ev_subsamp) <- as.character((2:90))
set.seed(2020)
for (simu in 1:nsimu) {

  for (t_w in 1:length((2:90))) {
    w_length <- (2:90)[t_w]


    MC_mean_nb_ev_subsamp[simu,t_w] <- mean(sample_nb_events_per_window(list_precip_winter = list_precip_winter,
                                                                        years_fullwinter = years_fullwinter,
                                                                        size_window = w_length,
                                                                        sample_size = nyears_per_simu))
  }#end for t_w

} #end for simu

nsimu2 <- 10000

MC_nb_ev_subsamp <- matrix(data = NA, nrow = nsimu2, ncol = length((2:90)))
colnames(MC_nb_ev_subsamp) <- as.character((2:90))

for (t_w in 1:length((2:90))) {
  w_length <- (2:90)[t_w]


  MC_nb_ev_subsamp[,t_w] <- sample_nb_events_per_window(list_precip_winter = list_precip_winter,
                                                        years_fullwinter = years_fullwinter,
                                                        size_window = w_length,
                                                        sample_size = nsimu2)
}#end for t_w

boxplot(MC_mean_nb_ev_subsamp[,c("30", "60", "90")])
low_quantile <- apply(X=MC_mean_nb_ev_subsamp, MARGIN = 2,FUN = quantile, probs=0.025)
up_quantile <- apply(X=MC_mean_nb_ev_subsamp, MARGIN = 2,FUN = quantile, probs=0.975)

plot(2:90,mean_number_events_before_landslide, pch="+", xlab="window length", ylim = c(0,12),
     ylab = "mean number of events", main = "Bootstrap number of events per time window\n(1000 repet., 23 winds chosen randomly in all winters)")
polygon(c(2:90,rev(2:90)),c(up_quantile,rev(low_quantile)),col="thistle",border=NA)
legend("bottomright", legend = c("before a landslide", "95%CI"), pch=c(3,15), col=c("black", "thistle"))

win_critical_duration <- c(10,10,1,30,15,75,5,1,1,15,15,30,40,60,75,40,90,60,60,4,10,40,1)
low_quantile2 <- apply(X=MC_nb_ev_subsamp, MARGIN = 2,FUN = quantile, probs=0.05)
up_quantile2 <- apply(X=MC_nb_ev_subsamp, MARGIN = 2,FUN = quantile, probs=0.95)

for (ls in 1:23) {
  png(paste0("../PrecipClustering_BottomUpAnalysis_LOCAL/nb_events_before_landslide_",dates_landslides[ls],".png"))
  plot(2:90,number_events_before_landslide[,ls], pch="+", xlab="window before event",
       ylim = c(0,max(number_events_before_landslide[,ls],12)),
       ylab = "Number of days precip>75th perc", main = paste("Lanslide", dates_landslides[ls]))
  polygon(c(2:90,rev(2:90)),c(up_quantile,rev(low_quantile)),col="thistle",border=NA)
  legend("bottomright", legend = c("before the landslide", "95%CI mean in winter","critical duration in data"),
         pch=c(3,15,NA), col=c("black", "thistle", "black"), lty=c(NA,NA,2))
  points(2:90,number_events_before_landslide[,ls], pch="+")
  abline(v=win_critical_duration[ls], lty=2)
  dev.off()
}#end for ls



# Compare to fit before precip ####

load(file = "Fit/Fit_precip_before_precip.Rdata")

quantile_0025 <- numeric()
quantile_005 <- numeric()
quantile_01 <- numeric()
quantile_025 <- numeric()
quantile_05 <- numeric()
quantile_075 <- numeric()
quantile_09 <- numeric()
quantile_095 <- numeric()
quantile_0975 <- numeric()
for (wind in 2:90) {
  if(!is.na(fit_before_precip$Fit_NegBinom[[wind-1]])){
    quantile_0025[wind-1] <- qnbinom(p=0.025, size = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["size"],
                                 mu = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["mu"])
    quantile_005[wind-1] <- qnbinom(p=0.05, size = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["size"],
                                   mu = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["mu"])
    quantile_01[wind-1] <- qnbinom(p=0.01, size = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["size"],
                                   mu = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["mu"])
    quantile_025[wind-1] <- qnbinom(p=0.25, size = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["size"],
                                   mu = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["mu"])
    quantile_05[wind-1] <- qnbinom(p=0.5, size = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["size"],
                                  mu = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["mu"])
    quantile_075[wind-1] <- qnbinom(p=0.75, size = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["size"],
                                   mu = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["mu"])
    quantile_09[wind-1] <- qnbinom(p=0.9, size = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["size"],
                                  mu = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["mu"])
    quantile_095[wind-1] <- qnbinom(p=0.95, size = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["size"],
                                 mu = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["mu"])
    quantile_0975[wind-1] <- qnbinom(p=0.975, size = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["size"],
                                 mu = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["mu"])
  } else {
    quantile_0025[wind-1] = quantile_005[wind-1] = quantile_025[wind-1] = quantile_075[wind-1] = NA
    quantile_09[wind-1] = quantile_095[wind-1] = quantile_0975[wind-1] = NA
  }
  
}#end for wind


quantile_0025[6] <- quantile_0025[7]
quantile_005[6] <- quantile_005[7]
quantile_01[6] <- quantile_01[7]
quantile_025[6] <- quantile_025[7]
quantile_05[6] <- quantile_05[7]
quantile_075[6] <- quantile_075[7]
quantile_09[6] <- quantile_09[7]
quantile_095[6] <- quantile_095[7]
quantile_0975[6] <- quantile_0975[7]


# #Smooth one
# quantile_0025 <- numeric()
# quantile_005 <- numeric()
# quantile_01 <- numeric()
# quantile_025 <- numeric()
# quantile_05 <- numeric()
# quantile_075 <- numeric()
# quantile_09 <- numeric()
# quantile_095 <- numeric()
# quantile_0975 <- numeric()
# for (wind in 2:90) {
#   if(!is.na(fit_before_precip$Fit_NegBinom[[wind-1]])){
#     
#     RUN <- rnbinom(n = 10000, size = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["size"],
#                    mu = fit_before_precip$Fit_NegBinom[[wind-1]]$estimate["mu"])
#     
#     quantile_0025[wind-1] <- quantile(x=RUN, probs=0.025)
#     quantile_005[wind-1] <- quantile(x=RUN, probs=0.05)
#     quantile_01[wind-1] <- quantile(x=RUN, probs=0.01)
#     quantile_025[wind-1] <- quantile(x=RUN, probs=0.25)
#     quantile_05[wind-1] <- quantile(x=RUN, probs=0.5)
#     quantile_075[wind-1] <- quantile(x=RUN, probs=0.75)
#     quantile_09[wind-1] <- quantile(x=RUN, probs=0.9)
#     quantile_095[wind-1] <- quantile(x=RUN, probs=0.95)
#     quantile_0975[wind-1] <- quantile(x=RUN, probs=0.975)
#   } else {
#     quantile_0025[wind-1] = quantile_005[wind-1] = quantile_025[wind-1] = quantile_075[wind-1] = NA
#     quantile_09[wind-1] = quantile_095[wind-1] = quantile_0975[wind-1] = NA
#   }
#   
# }#end for wind
# 
# 
# quantile_0025[6] <- mean(quantile_0025[5],quantile_0025[7])
# quantile_005[6] <-  mean(quantile_005[5],quantile_005[7])
# quantile_01[6] <-  mean(quantile_01[5],quantile_01[7])
# quantile_025[6] <-  mean(quantile_025[5],quantile_025[7])
# quantile_05[6] <-  mean(quantile_05[5],quantile_05[7])
# quantile_075[6] <-  mean(quantile_075[5],quantile_075[7])
# quantile_09[6] <-  mean(quantile_09[5],quantile_09[7])
# quantile_095[6] <-  mean(quantile_095[5],quantile_095[7])
# quantile_0975[6] <-  mean(quantile_0975[5],quantile_0975[7])
win_critical_duration <- c(10,10,1,30,15,75,5,1,1,15,15,30,40,60,75,40,90,60,60,4,10,40,1)
significant <- fit_before_precip$test_distrib[,"pval chi2 Negabinom"]>0.05
significant[is.na(significant)] <- FALSE

for (ls in 1:23) {
  # png(paste0("../PrecipClustering_BottomUpAnalysis_LOCAL/nb_events_before_landslide_",dates_landslides[ls],"_fitting.png"))
  plot(2:90,number_events_before_landslide[,ls], pch="+", xlab="window before event",
       # ylim = c(0,max(number_events_before_landslide[,ls],12)),
       ylim = c(0,18), col=rgb(0.5,0.5,0.5),
       ylab = "Number of days precip>75th perc", main = paste("Lanslide", dates_landslides[ls]))
  
  
  polygon(c(2:90,rev(2:90)),c(quantile_0975,rev(quantile_0025)),col=rgb(0.9,0.9,1),border=NA)
  polygon(c(2:90,rev(2:90)),c(quantile_095,rev(quantile_005)),col=rgb(0.85,0.85,1),border=NA)
  polygon(c(2:90,rev(2:90)),c(quantile_09,rev(quantile_01)),col=rgb(0.8,0.8,1),border=NA)
  polygon(c(2:90,rev(2:90)),c(quantile_075,rev(quantile_025)),col=rgb(0.75,0.75,1),border=NA)
  abline(v=win_critical_duration[ls], lty=2)
  legend("topleft", legend = c("before landslide", "CI before 1 precip event\n(nbinom fit)","critical duration"),
         pch=c(3,15,NA), col=c("black", rgb(0.85,0.85,1), "black"), lty=c(NA,NA,2))
  points(2:90,number_events_before_landslide[,ls], pch="+", col=rgb(0.5,0.5,0.5))
  points(names(significant[significant]),number_events_before_landslide[significant,ls], pch="+", lwd=3)
  text(x=90,y=quantile_05[89],srt=90,"50%", col="purple",cex=1.4)
  text(x=90,y=12.5,srt=90,"80%", col="purple",cex=1.2)
  text(x=90,y=15,srt=90,"90%", col="purple",cex=1.2)
  text(x=90,y=17.5,srt=90,"95%", col="purple",cex=1.2)
  # dev.off()
}#end for ls





####### number of events before a precip event in winter #######
for (size_wind in c(10,15,20,25,30,40,50,60,70,80,90)) {
  nb_events_before_precip <- list()
  accum_precip_before_event <- list()
  for (ref_winter in 1:length(list_precip_winter)) {
    nb_events_before_precip[[ref_winter]] <- numeric()
    accum_precip_before_event[[ref_winter]] <- numeric()

    ind <- which(list_precip_winter[[ref_winter]]==1)
    ind <- ind[ind>=size_wind]
    counting <- 0
    for (ID in ind) {
      counting <- counting + 1
      nb_events_before_precip[[ref_winter]][counting] <- sum(list_precip_winter[[ref_winter]][(ID-size_wind+1):ID])
      accum_precip_before_event[[ref_winter]][counting] <- sum(list_precip_winter_RAW[[ref_winter]][(ID-size_wind+1):ID])
    }#end for ID

  }#end for ref winter
  
  
  # png(paste0("../PrecipClustering_BottomUpAnalysis_LOCAL/nb_events_",size_wind,"days_before_event.png"))
  boxplot(unlist(nb_events_before_precip), number_events_before_landslide[size_wind-1,],
          names = c(paste0("before a precip event (",length(unlist(nb_events_before_precip)),")"),"before a landslide (23)"),
          main=paste("Number of precip events\n", size_wind, "days before precip event or lansdlide"))
  # dev.off()
  # png(paste0("../PrecipClustering_BottomUpAnalysis_LOCAL/accumulated_precip_",size_wind,"days_before_event.png"))
  boxplot(unlist(accum_precip_before_event), acc_precip_before_ls[size_wind-1,],
          names = c(paste0("before a precip event (",length(unlist(nb_events_before_precip)),")"),"before a landslide (23)"),
          main=paste("Accumulated precipitation [mm]\n", size_wind, "days before precip event or lansdlide"))
  # dev.off()

}#end for size_wind

NB_bef_precip <- matrix()
for (size_wind in 2:90) {
  nb_events_before_precip <- list()
  for (ref_winter in 1:length(list_precip_winter)) {
    nb_events_before_precip[[ref_winter]] <- numeric()
    
    ind <- which(list_precip_winter[[ref_winter]]==1)
    ind <- ind[ind>=size_wind]
    counting <- 0
    for (ID in ind) {
      counting <- counting + 1
      nb_events_before_precip[[ref_winter]][counting] <- sum(list_precip_winter[[ref_winter]][(ID-size_wind+1):ID])
    }#end for ID
    
  }#end for ref winter
  
  boxplot(unlist(nb_events_before_precip), number_events_before_landslide[size_wind-1,],
          names = c(paste0("before a precip event (",length(unlist(nb_events_before_precip)),")"),"before a landslide (23)"),
          main=paste("Number of precip events\n", size_wind, "days before precip event or lansdlide"))
  
  
}#end for size_wind


# ######## Fitting number of events before a precip event in winter  #######
# name_distrib <- matrix(data = NA, nrow = length(2:90), ncol = 5)
# row.names(name_distrib) <- c(2:90)
# colnames(name_distrib) <- c("best likelikihood", "pval chi2 Normal", "pval chi2 Poisson",
#                             "pval chi2 Logistic", "pval chi2 Negabinom")
# Fit_Poisson <- list()
# Fit_NegBinom <- list()
# Fit_Normal <- list()
# 
# chi_square_test <- function(empirical_freq, theoritic.proba){
#   theoritic.freq <- theoritic.proba*sum(empirical_freq)
# 
#   out <- matrix(data = NA, nrow = 1, ncol = 2)
#   colnames(out) <- c("stat", "pval")
#   stat_chi <- sum(((empirical_freq-theoritic.freq)^2)/theoritic.freq)
#   out[,"stat"] <- stat_chi
#   out[,"pval"] <- 1-pchisq(q = stat_chi, df = length(empirical_freq)-1)
#   return(out)
# }#end for chi_square_test
# 
# for (size_wind in 2:90) {
#   nb_events_before_precip <- list()
#   for (ref_winter in 1:length(list_precip_winter)) {
#     nb_events_before_precip[[ref_winter]] <- numeric()
# 
#     ind <- which(list_precip_winter[[ref_winter]]==1)
#     ind <- ind[ind>=size_wind]
#     counting <- 0
#     for (ID in ind) {
#       counting <- counting + 1
#       nb_events_before_precip[[ref_winter]][counting] <- sum(list_precip_winter[[ref_winter]][(ID-size_wind+1):ID])
#     }#end for ID
# 
#   }#end for ref winter
# 
#   nb_events_before_precip <- unlist(nb_events_before_precip)
# 
#   freq_empirical <- numeric()
# 
#   for (nb_ev_ind in 1:length(unique(nb_events_before_precip))) {
#     nb_ev <- sort(unique(nb_events_before_precip))[nb_ev_ind]
#     freq_empirical[nb_ev_ind] <- length(which(nb_events_before_precip==nb_ev))
#   }#end for nb_ev
# 
#   # nb_events_before_precip <-nb_events_before_precip[1:100]
# 
#   Fit_poiss <- MASS::fitdistr(nb_events_before_precip, densfun = "Poisson")
#   Fit_Poisson[[size_wind-1]] <- Fit_poiss
#   Fit_normal <- MASS::fitdistr(nb_events_before_precip, densfun = "Normal")
#   Fit_Normal[[size_wind-1]] <- Fit_normal
# 
#   if(length(unique(nb_events_before_precip))>1){
#     prob_poiss <- dpois(x=sort(unique(nb_events_before_precip)), lambda = Fit_poiss$estimate["lambda"])
# 
#     name_distrib[size_wind-1,"pval chi2 Poisson"] <- round(chi_square_test(empirical_freq = freq_empirical,
#                                                                            theoritic.proba = prob_poiss)[,"pval"],5)
# 
#     prob_norm <- dnorm(x=sort(unique(nb_events_before_precip)),
#                         mean = Fit_normal$estimate["mean"], sd = Fit_normal$estimate["sd"])
# 
#     name_distrib[size_wind-1,"pval chi2 Normal"] <- round(chi_square_test(empirical_freq = freq_empirical,
#                                                                           theoritic.proba = prob_norm)[,"pval"],5)
#   }
# 
# 
#   if(class(try(MASS::fitdistr(nb_events_before_precip, densfun = "negative binomial")))=="try-error" &
#      class(try(Fit_logi<- MASS::fitdistr(nb_events_before_precip, densfun = "logistic")))=="try-error" ){
#     Fit_NegBinom[[size_wind-1]] <- NA
#     max_loglik <- max(c(Fit_poiss$loglik,Fit_normal$loglik))
# 
#     if (max_loglik==Fit_poiss$loglik){name_distrib[size_wind-1,1]<-"Poisson"}
#     if (max_loglik==Fit_normal$loglik){name_distrib[size_wind-1,1]<-"Normal"}
#   } else {
#     if( class(try(Fit_logi<- MASS::fitdistr(nb_events_before_precip, densfun = "logistic")))=="try-error" ){
#       Fit_negabino <- MASS::fitdistr(nb_events_before_precip, densfun = "negative binomial")
#       Fit_NegBinom[[size_wind-1]] <- Fit_negabino
#       max_loglik <- max(c(Fit_poiss$loglik,Fit_normal$loglik, Fit_negabino$loglik))
# 
#       if (max_loglik==Fit_poiss$loglik){name_distrib[size_wind-1,1]<-"Poisson"}
#       if (max_loglik==Fit_normal$loglik){name_distrib[size_wind-1,1]<-"Normal"}
#       if (max_loglik==Fit_negabino$loglik){name_distrib[size_wind-1,1]<-"Negabinom"}
#     }
#     if(class(try(MASS::fitdistr(nb_events_before_precip, densfun = "negative binomial")))=="try-error"){
#       Fit_NegBinom[[size_wind-1]] <- NA
#       Fit_logi<- MASS::fitdistr(nb_events_before_precip, densfun = "logistic")
#       max_loglik <- max(c(Fit_poiss$loglik,Fit_normal$loglik, Fit_logi$loglik))
#       if (max_loglik==Fit_poiss$loglik){name_distrib[size_wind-1,1]<-"Poisson"}
#       if (max_loglik==Fit_normal$loglik){name_distrib[size_wind-1,1]<-"Normal"}
#       if (max_loglik==Fit_logi$loglik){name_distrib[size_wind-1,1]<-"Logistic"}
#     }
#   }
# 
#   if(class(try(MASS::fitdistr(nb_events_before_precip, densfun = "negative binomial")))!="try-error" &
#      class(try(Fit_logi<- MASS::fitdistr(nb_events_before_precip, densfun = "logistic")))!="try-error" ){
#     Fit_negabino <- MASS::fitdistr(nb_events_before_precip, densfun = "negative binomial")
#     Fit_NegBinom[[size_wind-1]] <- Fit_negabino
#     max_loglik <- max(c(Fit_poiss$loglik,Fit_normal$loglik, Fit_negabino$loglik))
#     Fit_logi<- MASS::fitdistr(nb_events_before_precip, densfun = "logistic")
#     max_loglik <- max(c(Fit_poiss$loglik,Fit_normal$loglik, Fit_logi$loglik, Fit_negabino$loglik))
#     if (max_loglik==Fit_poiss$loglik){name_distrib[size_wind-1,1]<-"Poisson"}
#     if (max_loglik==Fit_normal$loglik){name_distrib[size_wind-1,1]<-"Normal"}
#     if (max_loglik==Fit_logi$loglik){name_distrib[size_wind-1,1]<-"Logistic"}
#     if (max_loglik==Fit_negabino$loglik){name_distrib[size_wind-1,1]<-"Negabinom"}
#   }
# 
# 
#   if(class(try(MASS::fitdistr(nb_events_before_precip, densfun = "negative binomial")))!="try-error" &
#      length(unique(nb_events_before_precip))>1){
#   prob_nbin <- dnbinom(x=sort(unique(nb_events_before_precip)),
#                         size =  Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"])
# 
#   name_distrib[size_wind-1,"pval chi2 Negabinom"] <- round(chi_square_test(empirical_freq = freq_empirical,
#                                                                            theoritic.proba = prob_nbin)[,"pval"],5)
#   }
# 
#   if(class(try(Fit_logi<- MASS::fitdistr(nb_events_before_precip, densfun = "logistic")))!="try-error" &
#      length(unique(nb_events_before_precip))>1){
#     prob_logi <- dlogis(x=sort(unique(nb_events_before_precip)),
#                           location = Fit_logi$estimate["location"],
#                          scale = Fit_logi$estimate["scale"])
# 
#     name_distrib[size_wind-1,"pval chi2 Logistic"] <- round(chi_square_test(empirical_freq = freq_empirical,
#                                                                             theoritic.proba = prob_logi)[,"pval"],5)
#   }
# 
# }#end for size_wind
# 
# fit_before_precip <- list(Fit_NegBinom=Fit_NegBinom, Fit_Poisson=Fit_Poisson, Fit_Normal=Fit_Normal, test_distrib = name_distrib)
# 
# save(fit_before_precip, file = "Fit/Fit_precip_before_precip.Rdata")
# 
