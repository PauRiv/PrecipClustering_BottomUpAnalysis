#####################################################
#  Test of autocorrelation in precipitation events  #
#                                                   #
#            Pauline Rivoire 18.01.2021             #
#####################################################


library(lubridate)


#Load the data, change the path if necessary
excel_precip<-readxl::read_xlsx(path = "Data/Precipitation_data_1950_2008_IB02.xlsx")

#Name columns properly
colnames(excel_precip) <- c("year","month","day","gridpoint_W","gridpoint_E")

#Convert data of the gridpoint we want to a vector
precip_GP_West <- as.numeric(as.matrix(excel_precip[-c(1,2,21277),"gridpoint_W"]))
date_matrix <- as.matrix(excel_precip[-c(1,2,21277),c("year","month","day")])

# DATA <- data.frame(precip=precip_GP_West,precip_extended_winter = precip_extended_winter2,  Time=date_matrix)
# save(DATA, file = "/Users/admin/Desktop/Data_precip_portugal.Rdata")

# Precip time series : compute events -------------------------------------

#we exclude the two first rows and the last one because they don't contain data
All_years <- unique(date_matrix[,"year"])
years_fullwinter <- All_years[-c(1,length(All_years))] #Exclude first and last winter because not complete

nrow(date_matrix)==length(precip_GP_West)#check if same dimension

quantile_75 <- 10.5 #/!\ seas quantile, Nov to March (not October)

#extract seasonal precip
precip_extended_winter <- precip_GP_West
days_to_exclude <- which(date_matrix[,"month"] %in% 4:10) #exclude days from April to October included
precip_extended_winter[days_to_exclude] <- rep(0, length(days_to_exclude))
precip_extended_winter2 <- precip_GP_West
precip_extended_winter2[days_to_exclude] <- rep(NA, length(days_to_exclude))
winter_precip_event <- as.numeric(precip_extended_winter2>quantile_75)
#set those excluded days to 0 (can also be set to NA if needed)

cbind(date_matrix[,"month"], precip_GP_West, precip_extended_winter)

# compute_events <- function(precip_vector, threshold){
#   n_days <- length(precip_vector)
#   events_vector <- numeric(length = n_days)
#   events_vector[1] <- as.numeric(precip_vector[1]>=quantile_75)
#   
#   for (day in 2:n_days) {
#     if(events_vector[day-1]==0 & precip_vector[day]>=quantile_75){
#       events_vector[day] <- 1
#     }#end if
#   }#end for day
#   return(events_vector)
# }#end compute_events()
# 
# Events_75 <- compute_events(precip_vector = precip_extended_winter, threshold = quantile_75)


precip_per_winter <- function(precip_events, days, months, years, year_to_extract){
  days_to_extract <- c(which(years == year_to_extract & months%in%11:12),
                       which(years == (year_to_extract+1) & months%in%1:3))
  out <- cbind(paste0(as.character(years[days_to_extract]),"-", as.character(months[days_to_extract]),
                      "-",as.character(days[days_to_extract])),
               precip_events[days_to_extract])
  colnames(out) <- c("date", "precip")
  return(out)
}#end precip_per_winter_30() 

# list_precip_winter <- list()
list_precip_winter2 <- list()
list_precip_winter_RAW <- list()
count <- 0
for (year in years_fullwinter) {
  count <- count + 1
  # precip_winter <- as.matrix(as.numeric(precip_per_winter(year_to_extract = year,
  #                                                         precip_events = Events_75,
  #                                                         days = date_matrix[,"day"], months = date_matrix[,"month"],
  #                                                         years = date_matrix[,"year"])[,"precip"]))
  # colnames(precip_winter) <- paste("winter Nov.",as.character(year))
  # list_precip_winter[[count]] <- precip_winter
  
  list_precip_winter2[[count]] <- as.matrix(as.numeric(precip_per_winter(year_to_extract = year,
                                                                         precip_events = winter_precip_event,
                                                                         days = date_matrix[,"day"], months = date_matrix[,"month"],
                                                                         years = date_matrix[,"year"])[,"precip"]))
  
  list_precip_winter_RAW[[count]] <- as.matrix(as.numeric(precip_per_winter(year_to_extract = year,
                                                                            precip_events = precip_extended_winter,
                                                                            days = date_matrix[,"day"], months = date_matrix[,"month"],
                                                                            years = date_matrix[,"year"])[,"precip"]))
  colnames(list_precip_winter_RAW[[count]]) <- paste("winter Nov.",as.character(year))
}#end for year
rm(count)


# for (Y in 1:length(list_precip_winter)) {
#   acf(list_precip_winter2[[Y]],
#       plot = F,
#       lag.max = 5)
# }




# Autocorrelation test ----------------------------------------------------

acf(unlist(list_precip_winter2), lag.max = 10, main="All winters in one time series")

for (Y in 1:length(list_precip_winter2)) {
  
  acf(list_precip_winter2[[Y]], lag.max = 10, main=paste("winter",years_fullwinter[Y]))
}#end for Y




# Modify the extraml index function  --------

# to not take into account the long dependancies

require(extRemes)

dVec <- seq.Date(as.Date('1950-1-1'),as.Date('2008-3-30'),by='day')
djf_threshold <- 10.5
x = precip_extended_winter
# sum(x>=djf_threshold) # 587 values

extremalindex_mod = function (x, threshold, max_lag){
  # x = data vector
  # threshold = threshold to define POT events
  # max_lag = discrad inter-exceedance times above that time lag
  xout <- x
  data.name <- deparse(substitute(x))
  dname <- as.character(substitute(x))
  if (length(dname) > 1){
    dname <- c(dname[2], dname[length(dname)]) 
  } else {dname <- data.name}
  u <- threshold
  # x <- na.action(x)
  n <- length(x)
  eid <- x > u
  # N.u <- sum(eid)
  if (sum(eid) == 0) {
    res <- numeric(0)
    attr(res, "theta.msg") <- "extremalindex: No values of x exceed the threshold."
    class(res) <- "extremalindex"
    return(res)
  }
  # Intervals method
  T.u <- diff((1:n)[eid])
  T.u <- T.u[T.u<max_lag]
  N.u <- length(T.u)+1
  if (!any(T.u > 2)){
    hold <- colSums(cbind(T.u, T.u^2))
    theta.type <- "theta.hat"
  }else{
    hold <- colSums(cbind(T.u - 1, (T.u - 1) * (T.u - 2)))
    theta.type <- "theta.tilde"
  }
  theta <- 2 * (hold[1])^2/((N.u - 1) * hold[2])
  if (theta > 1) {
    theta <- 1
    theta.msg <- "Extremal index estimate > 1.  Re-set to 1."
  }else{theta.msg <- "Valid extremal index estimated"}
  K <- ifelse(round(theta * N.u, digits = 0) != theta * 
                N.u, ceiling(theta * N.u), theta * N.u)
  o <- order(T.u, na.last = TRUE, decreasing = TRUE)
  oT.u <- T.u[o]
  r <- oT.u[K]
  T.u2 <- c(diff(oT.u), 0)
  ind0 <- T.u2 == 0
  ind1 <- !ind0
  if ((K > 1) && (T.u2[K - 1] == 0)) {
    ind2.0 <- (1:(N.u - 1))[ind0]
    ind2.1 <- (1:(N.u - 1))[ind1]
    ind3 <- (ind2.0 < K) & (ind2.0 > max(ind2.1[ind2.1 < 
                                                  K]))
    K <- min(ind2.0[ind3])
    r <- oT.u[K]
  }
  res <- c(theta, K, r)
  attr(res, "theta.type") <- theta.type
  names(res) <- c("extremal.index", "number.of.clusters", "run.length")
  attr(res, "data") <- xout
  attr(res, "data.name") <- dname
  attr(res, "data.call") <- data.name
  attr(res, "call") <- match.call()
  attr(res, "na.action") <- na.action
  attr(res, "threshold") <- u
  # class(res) <- "extremalindex"
  return(res)
}

ext_index <- as.numeric(extremalindex_mod(x,djf_threshold,max_lag=60)[1])
avg_cluster_size <- 1/ext_index

# declustering
y = decluster(x,10.5,r=2)
inds = which(x>=10.5)
w = attr(y,"clusters") # for each threshold exceedance, w indicates the cluster it belongs to
declustered_precip <-as.numeric(y > 10.5)
cbind(x,y)[670:1170,]


Decl_precip <- cbind(date_matrix, declustered_precip)
# save(Decl_precip, file = "Data/declustered_precipitation.Rdata")






List_precip_events_winter01 <- list()
count <- 0
for (year in years_fullwinter) {
  count <- count + 1
  precip_winter <- as.matrix(as.numeric(precip_per_winter(year_to_extract = year,
                                                          precip_events = Decl_precip[,"declustered_precip"],
                                                          days = date_matrix[,"day"],
                                                          months = date_matrix[,"month"],
                                                          years = date_matrix[,"year"])[,"precip"]))
  colnames(precip_winter) <- paste("winter Nov.",as.character(year))
  List_precip_events_winter01[[count]] <- precip_winter
}#end for year
rm(count)



######## Fitting number of events before a precip event in winter  #######
name_distrib <- matrix(data = NA, nrow = length(2:90), ncol = 5)
row.names(name_distrib) <- c(2:90)
colnames(name_distrib) <- c("best likelikihood", "pval chi2 Normal", "pval chi2 Poisson",
                            "pval chi2 Logistic", "pval chi2 Negabinom")
Fit_Poisson <- list()
Fit_NegBinom <- list()
Fit_Normal <- list()

chi_square_test <- function(empirical_freq, theoritic.proba){
  theoritic.freq <- theoritic.proba*sum(empirical_freq)
  
  out <- matrix(data = NA, nrow = 1, ncol = 2)
  colnames(out) <- c("stat", "pval")
  stat_chi <- sum(((empirical_freq-theoritic.freq)^2)/theoritic.freq)
  out[,"stat"] <- stat_chi
  out[,"pval"] <- 1-pchisq(q = stat_chi, df = length(empirical_freq)-1)
  return(out)
}#end for chi_square_test

for (size_wind in 2:90) {
  nb_events_before_precip <- list()
  for (ref_winter in 1:length(List_precip_events_winter01)) {
    nb_events_before_precip[[ref_winter]] <- numeric()
    
    ind <- which(List_precip_events_winter01[[ref_winter]]==1)
    ind <- ind[ind>=size_wind]
    counting <- 0
    for (ID in ind) {
      counting <- counting + 1
      nb_events_before_precip[[ref_winter]][counting] <- sum(List_precip_events_winter01[[ref_winter]][(ID-size_wind+1):ID])
    }#end for ID
    
  }#end for ref winter
  
  nb_events_before_precip <- unlist(nb_events_before_precip)
  
  freq_empirical <- numeric()
  
  for (nb_ev_ind in 1:length(unique(nb_events_before_precip))) {
    nb_ev <- sort(unique(nb_events_before_precip))[nb_ev_ind]
    freq_empirical[nb_ev_ind] <- length(which(nb_events_before_precip==nb_ev))
  }#end for nb_ev
  
  # nb_events_before_precip <-nb_events_before_precip[1:100]
  
  Fit_poiss <- MASS::fitdistr(nb_events_before_precip, densfun = "Poisson")
  Fit_Poisson[[size_wind-1]] <- Fit_poiss
  Fit_normal <- MASS::fitdistr(nb_events_before_precip, densfun = "Normal")
  Fit_Normal[[size_wind-1]] <- Fit_normal
  
  if(length(unique(nb_events_before_precip))>1){
    prob_poiss <- dpois(x=sort(unique(nb_events_before_precip)), lambda = Fit_poiss$estimate["lambda"])
    
    name_distrib[size_wind-1,"pval chi2 Poisson"] <- round(chi_square_test(empirical_freq = freq_empirical,
                                                                           theoritic.proba = prob_poiss)[,"pval"],5)
    
    prob_norm <- dnorm(x=sort(unique(nb_events_before_precip)),
                       mean = Fit_normal$estimate["mean"], sd = Fit_normal$estimate["sd"])
    
    name_distrib[size_wind-1,"pval chi2 Normal"] <- round(chi_square_test(empirical_freq = freq_empirical,
                                                                          theoritic.proba = prob_norm)[,"pval"],5)
  }
  
  
  if(class(try(MASS::fitdistr(nb_events_before_precip, densfun = "negative binomial")))=="try-error" &
     class(try(Fit_logi<- MASS::fitdistr(nb_events_before_precip, densfun = "logistic")))=="try-error" ){
    Fit_NegBinom[[size_wind-1]] <- NA
    max_loglik <- max(c(Fit_poiss$loglik,Fit_normal$loglik))
    
    if (max_loglik==Fit_poiss$loglik){name_distrib[size_wind-1,1]<-"Poisson"}
    if (max_loglik==Fit_normal$loglik){name_distrib[size_wind-1,1]<-"Normal"}
  } else {
    if( class(try(Fit_logi<- MASS::fitdistr(nb_events_before_precip, densfun = "logistic")))=="try-error" ){
      Fit_negabino <- MASS::fitdistr(nb_events_before_precip, densfun = "negative binomial")
      Fit_NegBinom[[size_wind-1]] <- Fit_negabino
      max_loglik <- max(c(Fit_poiss$loglik,Fit_normal$loglik, Fit_negabino$loglik))
      
      if (max_loglik==Fit_poiss$loglik){name_distrib[size_wind-1,1]<-"Poisson"}
      if (max_loglik==Fit_normal$loglik){name_distrib[size_wind-1,1]<-"Normal"}
      if (max_loglik==Fit_negabino$loglik){name_distrib[size_wind-1,1]<-"Negabinom"}
    }
    if(class(try(MASS::fitdistr(nb_events_before_precip, densfun = "negative binomial")))=="try-error"){
      Fit_NegBinom[[size_wind-1]] <- NA
      Fit_logi<- MASS::fitdistr(nb_events_before_precip, densfun = "logistic")
      max_loglik <- max(c(Fit_poiss$loglik,Fit_normal$loglik, Fit_logi$loglik))
      if (max_loglik==Fit_poiss$loglik){name_distrib[size_wind-1,1]<-"Poisson"}
      if (max_loglik==Fit_normal$loglik){name_distrib[size_wind-1,1]<-"Normal"}
      if (max_loglik==Fit_logi$loglik){name_distrib[size_wind-1,1]<-"Logistic"}
    }
  }
  
  if(class(try(MASS::fitdistr(nb_events_before_precip, densfun = "negative binomial")))!="try-error" &
     class(try(Fit_logi<- MASS::fitdistr(nb_events_before_precip, densfun = "logistic")))!="try-error" ){
    Fit_negabino <- MASS::fitdistr(nb_events_before_precip, densfun = "negative binomial")
    Fit_NegBinom[[size_wind-1]] <- Fit_negabino
    max_loglik <- max(c(Fit_poiss$loglik,Fit_normal$loglik, Fit_negabino$loglik))
    Fit_logi<- MASS::fitdistr(nb_events_before_precip, densfun = "logistic")
    max_loglik <- max(c(Fit_poiss$loglik,Fit_normal$loglik, Fit_logi$loglik, Fit_negabino$loglik))
    if (max_loglik==Fit_poiss$loglik){name_distrib[size_wind-1,1]<-"Poisson"}
    if (max_loglik==Fit_normal$loglik){name_distrib[size_wind-1,1]<-"Normal"}
    if (max_loglik==Fit_logi$loglik){name_distrib[size_wind-1,1]<-"Logistic"}
    if (max_loglik==Fit_negabino$loglik){name_distrib[size_wind-1,1]<-"Negabinom"}
  }
  
  
  if(class(try(MASS::fitdistr(nb_events_before_precip, densfun = "negative binomial")))!="try-error" &
     length(unique(nb_events_before_precip))>1){
    prob_nbin <- dnbinom(x=sort(unique(nb_events_before_precip)),
                         size =  Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"])
    
    name_distrib[size_wind-1,"pval chi2 Negabinom"] <- round(chi_square_test(empirical_freq = freq_empirical,
                                                                             theoritic.proba = prob_nbin)[,"pval"],5)
  }
  
  if(class(try(Fit_logi<- MASS::fitdistr(nb_events_before_precip, densfun = "logistic")))!="try-error" &
     length(unique(nb_events_before_precip))>1){
    prob_logi <- dlogis(x=sort(unique(nb_events_before_precip)),
                        location = Fit_logi$estimate["location"],
                        scale = Fit_logi$estimate["scale"])
    
    name_distrib[size_wind-1,"pval chi2 Logistic"] <- round(chi_square_test(empirical_freq = freq_empirical,
                                                                            theoritic.proba = prob_logi)[,"pval"],5)
  }
  
}#end for size_wind


fit_before_precip <- list(Fit_NegBinom=Fit_NegBinom, Fit_Poisson=Fit_Poisson, Fit_Normal=Fit_Normal, test_distrib = name_distrib)

for (s in 2:90) {
  if(!is.na(Fit_NegBinom[[s-1]])){
    print(paste("prob=",as.numeric(((Fit_NegBinom[[s-1]]$estimate["size"])/(Fit_NegBinom[[s-1]]$estimate["mu"]+Fit_NegBinom[[s-1]]$estimate["size"])))))
  }
  
}

