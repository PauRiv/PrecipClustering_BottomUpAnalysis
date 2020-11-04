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


#Plot the obtain logical vector for winter 1993/1994
min_days_to_plot <- min(which((date_matrix[,"month"] %in% 11:12 & date_matrix[,"year"]==1993) |
                        (date_matrix[,"month"] %in% 1:3 & date_matrix[,"year"]==1994)))
max_days_to_plot <- max(which((date_matrix[,"month"] %in% 11:12 & date_matrix[,"year"]==1993) |
                                (date_matrix[,"month"] %in% 1:3 & date_matrix[,"year"]==1994)))
days_to_plot <- min_days_to_plot:max_days_to_plot


df.bar <-barplot(precip_extended_winter[days_to_plot],
                 col="blue", main = "winter 1993/1994")
abline(h=quantile_75, lty=3, lwd=1)
points(x = df.bar, y=Events_75[days_to_plot]*quantile_75, pch="+", col="red", lwd=4)



##### Study the number of events per winter #####

#we denote by winter of the year A the extended winter starting in Nov of year A en ending in March of year A+1
nb_events_per_winter <- function(precip_events, months, years, year_to_extract){
  days_to_extract <- c(which(years == year_to_extract & months%in%11:12),
                       which(years == (year_to_extract+1) & months%in%1:3))
  return(sum(precip_events[days_to_extract]))
}#end number_of_events_per_winter() 

nb_events_per_winter(precip_events = Events_75,
                     months = date_matrix[,"month"], years = date_matrix[,"year"],
                     year_to_extract = 1993)


All_years <- unique(date_matrix[,"year"])
years_fullwinter <- All_years[-c(1,length(All_years))]

NB_events_winters <- apply(X=as.matrix(years_fullwinter), MARGIN = 1,
                           FUN = nb_events_per_winter,
                           precip_events = Events_75, months = date_matrix[,"month"], years = date_matrix[,"year"])

winter_events <- cbind(years_fullwinter, NB_events_winters)
colnames(winter_events) <- c("year","nb_events")


plot(ecdf(winter_events[,"nb_events"]), xlab="Number of events", main="ecdf of the number of events per winter")

Fit_poiss <- MASS::fitdistr(winter_events[,"nb_events"], densfun = "Poisson")
lines(ppois(0:30, lambda = Fit_poiss$estimate), col="chartreuse4")

Fit_normal <- MASS::fitdistr(winter_events[,"nb_events"], densfun = "normal")
lines(pnorm(0:30, mean = Fit_normal$estimate["mean"], sd = Fit_normal$estimate["sd"]), col="blue")

Fit_logi<- MASS::fitdistr(winter_events[,"nb_events"], densfun = "logistic")
lines(plogis(0:30, location = Fit_logi$estimate["location"], scale = Fit_logi$estimate["scale"]), col="orange")

Fit_negabino<- MASS::fitdistr(winter_events[,"nb_events"], densfun = "negative binomial")
lines(pnbinom(0:30, size = Fit_negabino$estimate["size"], mu = Fit_negabino$estimate["mu"]), col="red")

legend("bottomright",title="Log-likelihood",
       legend = c(paste(round(Fit_poiss$loglik,2), "    Poisson"),
                  paste(round(Fit_normal$loglik,2), "    Normal"),
                  paste(round(Fit_logi$loglik,2), "    Logistic"),
                  paste(round(Fit_negabino$loglik,2), "    Negative binomial")),
       text.col=c("chartreuse4", "blue", "orange", "red"))
