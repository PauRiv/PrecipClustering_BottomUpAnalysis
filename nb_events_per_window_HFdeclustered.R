########################################
#   Number of events per time window   #
#     High-frequency declustered data  #
#     Pauline Rivoire   21.01.2021     #
########################################

rm(list=ls())

load("Data/declustered_precipitation.Rdata")
Decl_precip[,"declustered_precip"]

excel_precip<-readxl::read_xlsx(path = "Data/Precipitation_data_1950_2008_IB02.xlsx")
colnames(excel_precip) <- c("year","month","day","gridpoint_W","gridpoint_E")
precip_GP_West <- as.numeric(as.matrix(excel_precip[-c(1,2,21277),"gridpoint_W"]))
date_matrix <- as.matrix(excel_precip[-c(1,2,21277),c("year","month","day")])

#extract seasonal precip
precip_extended_winter <- precip_GP_West
days_to_exclude <- which(date_matrix[,"month"] %in% 4:10) #exclude days from April to October included
precip_extended_winter[days_to_exclude] <- rep(0, length(days_to_exclude))

excel_landslides <- readxl::read_xlsx(path = "../PrecipClustering_BottomUpAnalysis_LOCAL/Data/Landslides_Date.xlsx")
dates_landslides <- as.matrix(excel_landslides)
dates_landslides <- rbind("1958-12-19", dates_landslides)
ref_date_ls <- numeric()
for (ls in 1:length(dates_landslides)) {
  ref_date_ls[ls] <- which(date_matrix[,1]==year(dates_landslides[ls]) &
                             date_matrix[,2]==month(dates_landslides[ls]) &
                             date_matrix[,3]==day(dates_landslides[ls]))
}#end for ls

All_years <- unique(date_matrix[,"year"])
years_fullwinter <- All_years[-c(1,length(All_years))] #Exclude first and last winter because not complete

precip_per_winter <- function(precip_events, days, months, years, year_to_extract){
  days_to_extract <- c(which(years == year_to_extract & months%in%11:12),
                       which(years == (year_to_extract+1) & months%in%1:3))
  out <- cbind(paste0(as.character(years[days_to_extract]),"-", as.character(months[days_to_extract]),
                      "-",as.character(days[days_to_extract])),
               precip_events[days_to_extract])
  colnames(out) <- c("date", "precip")
  return(out)
}#end precip_per_winter_30() 


list_precip_events_winter <- list()
list_precip_winter_RAW <- list()
count <- 0
for (year in years_fullwinter) {
  count <- count + 1
  precip_winter <- as.matrix(as.numeric(precip_per_winter(year_to_extract = year,
                                                                      precip_events = Decl_precip[,"declustered_precip"],
                                                                      days = date_matrix[,"day"],
                                                                      months = date_matrix[,"month"],
                                                                      years = date_matrix[,"year"])[,"precip"]))
  colnames(precip_winter) <- paste("winter Nov.",as.character(year))
  list_precip_events_winter[[count]] <- precip_winter
  list_precip_winter_RAW[[count]] <- as.matrix(as.numeric(precip_per_winter(year_to_extract = year,
                                                                            precip_events = precip_extended_winter,
                                                                            days = date_matrix[,"day"], months = date_matrix[,"month"],
                                                                            years = date_matrix[,"year"])[,"precip"]))
  colnames(list_precip_winter_RAW[[count]]) <- paste("winter Nov.",as.character(year))
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
  for (ref_winter in 1:length(list_precip_events_winter)) {
    nb_events_before_precip[[ref_winter]] <- numeric()

    ind <- which(list_precip_events_winter[[ref_winter]]==1)
    ind <- ind[ind>=size_wind]
    counting <- 0
    for (ID in ind) {
      counting <- counting + 1
      nb_events_before_precip[[ref_winter]][counting] <- sum(list_precip_events_winter[[ref_winter]][(ID-size_wind+1):ID])
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
  print(as.numeric(((Fit_NegBinom[[s-1]]$estimate["size"])/(Fit_NegBinom[[s-1]]$estimate["mu"]+Fit_NegBinom[[s-1]]$estimate["size"]))))
}
