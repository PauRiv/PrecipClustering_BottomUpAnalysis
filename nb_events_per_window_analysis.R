# Pauline Rivoire 04.11.2020
#Preliminary analysis of the precipitation time series at the grid point

###### Extract the precipitation data #####
#Load the data, change the path if necessary
excel_precip<-readxl::read_xlsx(path = "Data/Precipitation_data_1950_2008_IB02.xlsx")

#Name columns properly
colnames(excel_precip) <- c("year","month","day","gridpoint_W","gridpoint_E")

#Convert data of the gridpoint we want to a vector
precip_GP_West <- as.numeric(as.matrix(excel_precip[-c(1,2,21277),"gridpoint_W"]))
date_matrix <- as.matrix(excel_precip[-c(1,2,21277),c("year","month","day")])
#we exclude the two first rows and the last one because they don't contain data

nrow(date_matrix)==length(precip_GP_West)#check if same dimension

#extract seasonal precip
precip_extended_winter <- precip_GP_West
days_to_exclude <- which(date_matrix[,"month"] %in% 4:10) #exclude days from April to October included
precip_extended_winter[days_to_exclude] <- rep(0, length(days_to_exclude))
#set those excluded days to 0 (can also be set to NA if needed)

cbind(date_matrix[,"month"], precip_GP_West, precip_extended_winter)

quantile_75 <- quantile(x=precip_extended_winter[precip_extended_winter>=1], probs = 0.75) 


##### Compute the events above the 75th perc ####
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
}#end precip_per_winter() 

All_years <- unique(date_matrix[,"year"])
years_fullwinter <- All_years[-c(1,length(All_years))]

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




##### Extract number of events per window #####
nb_evts_30 <- numeric(length = 5*length(years_fullwinter))
nb_evts_60 <- numeric(length = 2*length(years_fullwinter))
nb_evts_90 <- numeric(length = length(years_fullwinter))

for (ind_yr in 1:length(years_fullwinter)) {
  ndays <- length(list_precip_winter[[ind_yr]])
  for (w_30 in 1:5) {#select 5 time windows of size 30 during the winter of this year
    nb_evts_30[(ind_yr-1)*5+w_30] <- sum(list_precip_winter[[ind_yr]][((w_30-1)*30+1):(w_30*30)])
  } #end for w_30
  
  for (w_60 in 1:2) {#select 2 time windows of size 60 during the winter of this year
    nb_evts_60[(ind_yr-1)*2+w_60] <- sum(list_precip_winter[[ind_yr]][((w_60-1)*60+1):(w_60*60)])
  } #end for w_60
  
  #select 1 time window of size 90 during the winter of this year
  nb_evts_90[ind_yr] <- sum(list_precip_winter[[ind_yr]][1:90])
}

plot(ecdf(nb_evts_30))
Fit_poiss <- MASS::fitdistr(nb_evts_30, densfun = "Poisson")
lines(ppois(0:30, lambda = Fit_poiss$estimate), col="chartreuse4")
#ks.test(x = nb_evts_30+rnorm(n=5*length(years_fullwinter),sd=0.00003), y = "ppois", lambda = Fit_poiss$estimate)

Fit_normal <- MASS::fitdistr(nb_evts_30, densfun = "normal")
lines(pnorm(0:30, mean = Fit_normal$estimate["mean"], sd = Fit_normal$estimate["sd"]), col="blue")
#ks.test(x = nb_evts_30+rnorm(n=5*length(years_fullwinter),sd=0.00003), y = "pnorm", mean = Fit_normal$estimate["mean"], sd = Fit_normal$estimate["sd"])

Fit_logi<- MASS::fitdistr(nb_evts_30, densfun = "logistic")
lines(plogis(0:30, location = Fit_logi$estimate["location"], scale = Fit_logi$estimate["scale"]), col="orange")
#ks.test(x = nb_evts_30+rnorm(n=5*length(years_fullwinter),sd=0.00003), y = "plogis", location = Fit_logi$estimate["location"], scale = Fit_logi$estimate["scale"])

Fit_negabino<- MASS::fitdistr(nb_evts_30, densfun = "negative binomial")
lines(pnbinom(0:30, size = Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"]), col="red")
KS <- ks.test(x = nb_evts_30+rnorm(n=5*length(years_fullwinter),sd=0.00003), y = "pnbinom", size = Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"])
text(x=8.5, y= 0.5, paste("pval ks-test nbinom", round(KS$p.value, 6)), col="red")

legend("bottomright",title="Log-likelihood",
       legend = c(paste(round(Fit_poiss$loglik,2), "    Poisson"),
                  paste(round(Fit_normal$loglik,2), "    Normal"),
                  paste(round(Fit_logi$loglik,2), "    Logistic"),
                  paste(round(Fit_negabino$loglik,2), "    Negative binomial")),
       text.col=c("chartreuse4", "blue", "orange", "red"))


plot(ecdf(nb_evts_60))
Fit_poiss <- MASS::fitdistr(nb_evts_60, densfun = "Poisson")
lines(ppois(0:30, lambda = Fit_poiss$estimate), col="chartreuse4")
#ks.test(x = nb_evts_60+rnorm(n=2*length(years_fullwinter),sd=0.00003), y = "ppois", lambda = Fit_poiss$estimate)

Fit_normal <- MASS::fitdistr(nb_evts_60, densfun = "normal")
lines(pnorm(0:30, mean = Fit_normal$estimate["mean"], sd = Fit_normal$estimate["sd"]), col="blue")
#ks.test(x = nb_evts_60+rnorm(n=2*length(years_fullwinter),sd=0.00003), y = "pnorm", mean = Fit_normal$estimate["mean"], sd = Fit_normal$estimate["sd"])

Fit_logi<- MASS::fitdistr(nb_evts_60, densfun = "logistic")
lines(plogis(0:30, location = Fit_logi$estimate["location"], scale = Fit_logi$estimate["scale"]), col="orange")
#ks.test(x = nb_evts_60+rnorm(n=2*length(years_fullwinter),sd=0.00003), y = "plogis", location = Fit_logi$estimate["location"], scale = Fit_logi$estimate["scale"])

Fit_negabino<- MASS::fitdistr(nb_evts_60, densfun = "negative binomial")
lines(pnbinom(0:30, size = Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"]), col="red")
KS <- ks.test(x = nb_evts_60+rnorm(n=2*length(years_fullwinter),sd=0.00003), y = "pnbinom", size = Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"])
text(x=11.5, y= 0.5, paste("pval ks-test nbinom", round(KS$p.value, 6)), col="red")

legend("bottomright",title="Log-likelihood",
       legend = c(paste(round(Fit_poiss$loglik,2), "    Poisson"),
                  paste(round(Fit_normal$loglik,2), "    Normal"),
                  paste(round(Fit_logi$loglik,2), "    Logistic"),
                  paste(round(Fit_negabino$loglik,2), "    Negative binomial")),
       text.col=c("chartreuse4", "blue", "orange", "red"))

plot(ecdf(nb_evts_90))
Fit_poiss <- MASS::fitdistr(nb_evts_90, densfun = "Poisson")
lines(ppois(0:30, lambda = Fit_poiss$estimate), col="chartreuse4")
#ks.test(x = nb_evts_90+rnorm(n=length(years_fullwinter),sd=0.00003), y = "ppois", lambda = Fit_poiss$estimate)

Fit_normal <- MASS::fitdistr(nb_evts_90, densfun = "normal")
lines(pnorm(0:30, mean = Fit_normal$estimate["mean"], sd = Fit_normal$estimate["sd"]), col="blue")
#ks.test(x = nb_evts_90+rnorm(n=length(years_fullwinter),sd=0.00003), y = "pnorm", mean = Fit_normal$estimate["mean"], sd = Fit_normal$estimate["sd"])

Fit_logi<- MASS::fitdistr(nb_evts_90, densfun = "logistic")
lines(plogis(0:30, location = Fit_logi$estimate["location"], scale = Fit_logi$estimate["scale"]), col="orange")
#ks.test(x = nb_evts_90+rnorm(n=length(years_fullwinter),sd=0.00003), y = "plogis", location = Fit_logi$estimate["location"], scale = Fit_logi$estimate["scale"])

Fit_negabino<- MASS::fitdistr(nb_evts_90, densfun = "negative binomial")
lines(pnbinom(0:30, size = Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"]), col="red")
KS <- ks.test(x = nb_evts_90+rnorm(n=length(years_fullwinter),sd=0.00003), y = "pnbinom", size = Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"])
text(x=16.4, y= 0.5, paste("pval ks-test nbinom", round(KS$p.value, 6)), col="red")
legend("bottomright",title="Log-likelihood",
       legend = c(paste(round(Fit_poiss$loglik,2), "    Poisson"),
                  paste(round(Fit_normal$loglik,2), "    Normal"),
                  paste(round(Fit_logi$loglik,2), "    Logistic"),
                  paste(round(Fit_negabino$loglik,2), "    Negative binomial")),
       text.col=c("chartreuse4", "blue", "orange", "red"))








# ##### Select randomly the windows #####
# 
# #function to extract the number of events in a random window of fixed size for a fix year
# nb_events_random_window <- function(list_precip_winter, year, years_fullwinter, size_window){
#   year_ind <- which(years_fullwinter==year)
#   rdm_start <- sample(1:(length(list_precip_winter[[year_ind]])-size_window+1),size=1)
#   return(sum(list_precip_winter[[year_ind]][rdm_start:(rdm_start+size_window-1)]))
# }#end nb_events_random_window()
# 
# #function to obtain a sample of the number of events per window of a fixed size, chosen randomly within all the years
# sample_nb_events_per_window <- function(list_precip_winter, years_fullwinter, size_window, sample_size){
#   rdm_years <- as.matrix(sample(years_fullwinter, size = sample_size, replace = T))
# 
#   out <- apply(X=rdm_years, MARGIN = 1, FUN = nb_events_random_window,
#                list_precip_winter = list_precip_winter, years_fullwinter =years_fullwinter, size_window = size_window)
#   return(out)
# }#end sample_nb_events_per_window()
# 
# samp_size <- 100
# 
# nb_events_30 <- sample_nb_events_per_window(list_precip_winter = list_precip_winter,
#                                             years_fullwinter = years_fullwinter, size_window = 30, sample_size = samp_size)
# nb_events_60 <- sample_nb_events_per_window(list_precip_winter = list_precip_winter,
#                                             years_fullwinter = years_fullwinter, size_window = 60, sample_size = samp_size)
# nb_events_90 <- sample_nb_events_per_window(list_precip_winter = list_precip_winter,
#                                             years_fullwinter = years_fullwinter, size_window = 90, sample_size = samp_size)
# plot(ecdf(nb_events_30))
# Fit_poiss <- MASS::fitdistr(nb_events_30, densfun = "Poisson")
# lines(ppois(0:30, lambda = Fit_poiss$estimate), col="chartreuse4")
# 
# Fit_normal <- MASS::fitdistr(nb_events_30, densfun = "normal")
# lines(pnorm(0:30, mean = Fit_normal$estimate["mean"], sd = Fit_normal$estimate["sd"]), col="blue")
# 
# Fit_logi<- MASS::fitdistr(nb_events_30, densfun = "logistic")
# lines(plogis(0:30, location = Fit_logi$estimate["location"], scale = Fit_logi$estimate["scale"]), col="orange")
# 
# Fit_negabino<- MASS::fitdistr(nb_events_30, densfun = "negative binomial")
# lines(pnbinom(0:30, size = Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"]), col="red")
# KS <- ks.test(x = nb_events_30+rnorm(n=length(samp_size),sd=0.00003), y = "pnbinom", size = Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"])
# text(x=6.4, y= 0.5, paste("pval ks-test nbinom", round(KS$p.value, 6)), col="red")
# 
# 
# legend("bottomright",title="Log-likelihood",
#        legend = c(paste(round(Fit_poiss$loglik,2), "    Poisson"),
#                   paste(round(Fit_normal$loglik,2), "    Normal"),
#                   paste(round(Fit_logi$loglik,2), "    Logistic"),
#                   paste(round(Fit_negabino$loglik,2), "    Negative binomial")),
#        text.col=c("chartreuse4", "blue", "orange", "red"))
# 
# 
# plot(ecdf(nb_events_60))
# Fit_poiss <- MASS::fitdistr(nb_events_60, densfun = "Poisson")
# lines(ppois(0:30, lambda = Fit_poiss$estimate), col="chartreuse4")
# Fit_normal <- MASS::fitdistr(nb_events_60, densfun = "normal")
# lines(pnorm(0:30, mean = Fit_normal$estimate["mean"], sd = Fit_normal$estimate["sd"]), col="blue")
# 
# Fit_logi<- MASS::fitdistr(nb_events_60, densfun = "logistic")
# lines(plogis(0:30, location = Fit_logi$estimate["location"], scale = Fit_logi$estimate["scale"]), col="orange")
# 
# Fit_negabino<- MASS::fitdistr(nb_events_60, densfun = "negative binomial")
# lines(pnbinom(0:30, size = Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"]), col="red")
# KS <- ks.test(x = nb_events_60+rnorm(n=length(samp_size),sd=0.00003), y = "pnbinom", size = Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"])
# text(x=12.4, y= 0.5, paste("pval ks-test nbinom", round(KS$p.value, 6)), col="red")
# 
# legend("bottomright",title="Log-likelihood",
#        legend = c(paste(round(Fit_poiss$loglik,2), "    Poisson"),
#                   paste(round(Fit_normal$loglik,2), "    Normal"),
#                   paste(round(Fit_logi$loglik,2), "    Logistic"),
#                   paste(round(Fit_negabino$loglik,2), "    Negative binomial")),
#        text.col=c("chartreuse4", "blue", "orange", "red"))
# 
# plot(ecdf(nb_events_90))
# Fit_poiss <- MASS::fitdistr(nb_events_90, densfun = "Poisson")
# lines(ppois(0:30, lambda = Fit_poiss$estimate), col="chartreuse4")
# Fit_normal <- MASS::fitdistr(nb_events_90, densfun = "normal")
# lines(pnorm(0:30, mean = Fit_normal$estimate["mean"], sd = Fit_normal$estimate["sd"]), col="blue")
# 
# Fit_logi<- MASS::fitdistr(nb_events_90, densfun = "logistic")
# lines(plogis(0:30, location = Fit_logi$estimate["location"], scale = Fit_logi$estimate["scale"]), col="orange")
# 
# Fit_negabino<- MASS::fitdistr(nb_events_90, densfun = "negative binomial")
# lines(pnbinom(0:30, size = Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"]), col="red")
# KS <- ks.test(x = nb_events_90+rnorm(n=length(samp_size),sd=0.00003), y = "pnbinom", size = Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"])
# text(x=15.4, y= 0.5, paste("pval ks-test nbinom", round(KS$p.value, 6)), col="red")
# 
# legend("bottomright",title="Log-likelihood",
#        legend = c(paste(round(Fit_poiss$loglik,2), "    Poisson"),
#                   paste(round(Fit_normal$loglik,2), "    Normal"),
#                   paste(round(Fit_logi$loglik,2), "    Logistic"),
#                   paste(round(Fit_negabino$loglik,2), "    Negative binomial")),
#        text.col=c("chartreuse4", "blue", "orange", "red"))
