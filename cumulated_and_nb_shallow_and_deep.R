###########################################
#  Number of events and cumulated precip  #
#    before shallow or deep landslides    #
#       Pauline Rivoire 07.01.2021        #
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
date_vector <- apply(X=date_matrix, FUN=paste, collapse="-",MARGIN = 1)

ref_date_ls <- numeric()
for (ls in 1:length(dates_landslides)) {
  ref_date_ls[ls] <- which(date_matrix[,1]==year(dates_landslides[ls]) &
                             date_matrix[,2]==month(dates_landslides[ls]) &
                             date_matrix[,3]==day(dates_landslides[ls]))
}#end for ls




# Mean number events before landslide -------------------------------------



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
list_date <- list()
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
  list_date[[count]] <- as.Date(precip_per_winter(year_to_extract = year,
                                         precip_events = date_vector,
                                         days = date_matrix[,"day"], months = date_matrix[,"month"],
                                         years = date_matrix[,"year"])[,"precip"])
}#end for year
rm(count)



# Compare to fit before precip --------------------------------------------


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

win_critical_duration <- c(10,10,1,30,15,75,5,1,1,15,15,30,40,60,75,40,90,60,60,4,10,40,1)
significant <- fit_before_precip$test_distrib[,"pval chi2 Negabinom"]>0.05
significant[is.na(significant)] <- FALSE

landslide_type <- c("S","S","S","D","S","D","S","S","S","S","S","D","D",
                    "D","D","D","D","D","D","S","S","D","S")

for (ls in 1:23) {
  # png(paste0("../PrecipClustering_BottomUpAnalysis_LOCAL/nb_events_before_landslide_",dates_landslides[ls],"_fitting.png"))
  plot(2:90,number_events_before_landslide[,ls], pch="+", xlab="window before event",
       # ylim = c(0,max(number_events_before_landslide[,ls],12)),
       ylim = c(0,18), col=rgb(0.5,0.5,0.5),
       ylab = "Number of days precip>75th perc",
       main = paste0("Lanslide ", dates_landslides[ls],", type=",landslide_type[ls]))
  
  
  polygon(c(2:90,rev(2:90)),c(quantile_0975,rev(quantile_0025)),col=rgb(0.75,0.75,1),border=NA)
  polygon(c(2:90,rev(2:90)),c(quantile_095,rev(quantile_005)),col=rgb(0.8,0.8,1),border=NA)
  polygon(c(2:90,rev(2:90)),c(quantile_09,rev(quantile_01)),col=rgb(0.85,0.85,1),border=NA)
  polygon(c(2:90,rev(2:90)),c(quantile_075,rev(quantile_025)),col=rgb(0.9,0.9,1),border=NA)
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






# Cumulated precip before deep and shallow landslides ------------------------

ind_deep <- which(year(dates_landslides)==1996)
ind_shallow <- which(year(dates_landslides)==2006 & month(dates_landslides)<11)

landslide_type <- c("S","S","S","D","S","D","S","S","S","S","S","D","D",
                    "D","D","D","D","D","D","S","S","D","S")

plot(((-88):0),rev(acc_precip_before_ls[,1]), type="l", col="green", lwd=2,
     ylim = c(0,800), xlab = "Days before landslide", ylab = "cumulated precip (backward, mm/day)")

for (ind in 2:23) {
  if(landslide_type[ind]=="S"){COL <- "green"}else{COL <- "chartreuse4"}
  lines(((-88):0),rev(acc_precip_before_ls[,ind]), type="l", col=COL, lwd=2)
}
legend("topright", legend = c("Deep", "Shallow"), title="Landlside type", col=c("chartreuse4", "green"),
       lty=1, lwd=2)


plot(((-88):0),rev(acc_precip_before_ls[,1]), type="l", col="green", lwd=2,
     ylim = c(0,800), xlab = "Days before landslide", ylab = "cumulated precip (backward, mm/day)")

for (ind in 2:23) {
  if(landslide_type[ind]=="S"){COL <- "green"}else{COL <- "chartreuse4"}
  lines(((-88):0),rev(acc_precip_before_ls[,ind]), type="l", col=COL, lwd=2)
}

for (ind in ind_deep) {
  lines(((-88):0),rev(acc_precip_before_ls[,ind]), type="l", col="darkgreen", lwd=3)
}
for (ind in ind_shallow) {
  lines(((-88):0),rev(acc_precip_before_ls[,ind]), type="l", col="darkolivegreen1", lwd=3)
}

legend("topright", legend = c("Deep", "Shallow", "1996 (deep)", "2006 (shallow)"), title="Landlside type",
       col=c("chartreuse4", "green","darkgreen", "darkolivegreen1"), lty=1, lwd=c(2,2,3,3))




# Winters without landslides ----------------------------------------------

winter_landslide_Oct <- c(1958,1958,1967,1968,1977,1978,1981,1983,1986,1989,1989,1989,1989,1995,
                          1995,1995,1995,2000,2000,2005,2006,2006,2007)
length(unique(winter_landslide_Oct))

list_precip_winter_RAW[[1]]

cumul_precip_winter <- as.matrix(unlist(lapply(X=list_precip_winter_RAW, FUN=sum)))
rownames(cumul_precip_winter) <- years_fullwinter
sort(cumul_precip_winter, decreasing = T)[1:14]

unique(winter_landslide_Oct[which(landslide_type=="D")])

Is_landslide <- years_fullwinter[which(cumul_precip_winter %in% sort(cumul_precip_winter, decreasing = T)[1:7])] %in% unique(winter_landslide_Oct)

Years_no_landsl_lot_precip <- years_fullwinter[which(cumul_precip_winter %in% sort(cumul_precip_winter, decreasing = T)[1:7])][!Is_landslide]


#Info in a table
sorted_winters <- years_fullwinter[sort(cumul_precip_winter, decreasing = T, index.return = T)$ix]
sorted_precip <- sort(cumul_precip_winter, decreasing = T)
Is_landslide <- years_fullwinter[sort(cumul_precip_winter, decreasing = T, index.return = T)$ix] %in% winter_landslide_Oct
landslide_T <- numeric(length = length(years_fullwinter))
for (Y in 1:length(years_fullwinter)) {
  if(Is_landslide[Y]){
    if(length(landslide_type[sorted_winters[Y]==winter_landslide_Oct])>1){ 
      if(("D" %in% landslide_type[sorted_winters[Y]==winter_landslide_Oct])){
        landslide_T[Y] <- "D"
      } else {landslide_T[Y] <- "S"}
    } else {
      landslide_T[Y] <- landslide_type[sorted_winters[Y]==winter_landslide_Oct]
    }
  } 
}#end for Y

Table_sorted_precip <- cbind(1:length(years_fullwinter),sorted_winters, sorted_precip, Is_landslide, landslide_T)

colnames(Table_sorted_precip) <- c("Rank", "Nov year", "cumulated precip (mm, Nov-Mar)",
                                   "Landslide occurrence", "Largest landslide type")

write.csv(Table_sorted_precip,
          file = "/Users/admin/Documents/PrecipClustering_BottomUpAnalysis_LOCAL/winters_sorted_cumul_precip.csv")

# 1955

plot(list_precip_winter_RAW[[which(years_fullwinter==Years_no_landsl_lot_precip[1])]],
     main=paste("Precipitation winter", Years_no_landsl_lot_precip[1]), ylab="precip (mm/day)")
abline(h=10.5, col="blue", lty=2)

plot(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[1])]],
     main=paste("Precipitation winter", Years_no_landsl_lot_precip[1]), ylab="day with precip > 75th perc",
     pch="+")
rect(xleft = 0, ybottom = 0.99, xright = 78.5, ytop = 1.01, border="blue")

plot(2:78,cumsum(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[1])]][78:1])[-1], pch="+", xlab="window before event",
     # ylim = c(0,max(number_events_before_landslide[,ls],12)),
     ylim = c(0,18), col=rgb(0.5,0.5,0.5),
     ylab = "Number of days precip>75th perc", main = paste("1955 1st part"))

polygon(c(2:78,rev(2:78)),c(quantile_0975[1:77],rev(quantile_0025[1:77])),col=rgb(0.75,0.75,1),border=NA)
polygon(c(2:78,rev(2:78)),c(quantile_095[1:77],rev(quantile_005[1:77])),col=rgb(0.8,0.8,1),border=NA)
polygon(c(2:78,rev(2:78)),c(quantile_09[1:77],rev(quantile_01[1:77])),col=rgb(0.85,0.85,1),border=NA)
polygon(c(2:78,rev(2:78)),c(quantile_075[1:77],rev(quantile_025[1:77])),col=rgb(0.9,0.9,1),border=NA)

legend("topleft", legend = c("before last event in rect.", "CI before 1 precip event\n(nbinom fit)"),
       pch=c(3,15,NA), col=c("black", rgb(0.85,0.85,1), "black"), lty=c(NA,NA,2))
points(2:78,cumsum(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[1])]][78:1])[-1], pch="+", col="black")
text(x=79,y=quantile_05[89],srt=90,"50%", col="purple",cex=1.4)
text(x=79,y=12.5,srt=90,"80%", col="purple",cex=1.2)
text(x=79,y=15,srt=90,"90%", col="purple",cex=1.2)
text(x=79,y=17.5,srt=90,"95%", col="purple",cex=1.2)

plot(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[1])]],
     main=paste("Precipitation winter", Years_no_landsl_lot_precip[1]), ylab="day with precip > 75th perc",
     pch="+")
rect(xleft = 60.5, ybottom = 0.99, xright = 149.5, ytop = 1.01, border="blue")

plot(2:90,cumsum(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[1])]][149:60])[-1], pch="+", xlab="window before event",
     # ylim = c(0,max(number_events_before_landslide[,ls],12)),
     ylim = c(0,18), col=rgb(0.5,0.5,0.5),
     ylab = "Number of days precip>75th perc", main = paste("1955 2nd part"))

polygon(c(2:90,rev(2:90)),c(quantile_0975,rev(quantile_0025)),col=rgb(0.75,0.75,1),border=NA)
polygon(c(2:90,rev(2:90)),c(quantile_095,rev(quantile_005)),col=rgb(0.8,0.8,1),border=NA)
polygon(c(2:90,rev(2:90)),c(quantile_09,rev(quantile_01)),col=rgb(0.85,0.85,1),border=NA)
polygon(c(2:90,rev(2:90)),c(quantile_075,rev(quantile_025)),col=rgb(0.9,0.9,1),border=NA)

legend("topleft", legend = c("before last event in rect.", "CI before 1 precip event\n(nbinom fit)"),
       pch=c(3,15,NA), col=c("black", rgb(0.85,0.85,1), "black"), lty=c(NA,NA,2))
points(2:90,cumsum(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[1])]][149:60])[-1], pch="+", col="black")
text(x=90,y=quantile_05[89],srt=90,"50%", col="purple",cex=1.4)
text(x=90,y=12.5,srt=90,"80%", col="purple",cex=1.2)
text(x=90,y=15,srt=90,"90%", col="purple",cex=1.2)
text(x=90,y=17.5,srt=90,"95%", col="purple",cex=1.2)


# 1963

plot(list_precip_winter_RAW[[which(years_fullwinter==Years_no_landsl_lot_precip[2])]],
     main=paste("Precipitation winter", Years_no_landsl_lot_precip[2]), ylab="precip (mm/day)")
abline(h=10.5, col="blue", lty=2)

# 1st part
plot(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[2])]],
     main=paste("Precipitation winter", Years_no_landsl_lot_precip[2]), ylab="day with precip > 75th perc",
     pch="+")
rect(xleft = 0, ybottom = 0.99, xright = 76.5, ytop = 1.01, border="blue")

plot(2:76,cumsum(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[2])]][76:1])[-1], pch="+", xlab="window before event",
     # ylim = c(0,max(number_events_before_landslide[,ls],12)),
     ylim = c(0,18), col=rgb(0.5,0.5,0.5),
     ylab = "Number of days precip>75th perc", main = paste("1963 1st part"))

polygon(c(2:76,rev(2:76)),c(quantile_0975[1:75],rev(quantile_0025[1:75])),col=rgb(0.75,0.75,1),border=NA)
polygon(c(2:76,rev(2:76)),c(quantile_095[1:75],rev(quantile_005[1:75])),col=rgb(0.8,0.8,1),border=NA)
polygon(c(2:76,rev(2:76)),c(quantile_09[1:75],rev(quantile_01[1:75])),col=rgb(0.85,0.85,1),border=NA)
polygon(c(2:76,rev(2:76)),c(quantile_075[1:75],rev(quantile_025[1:75])),col=rgb(0.9,0.9,1),border=NA)

legend("topleft", legend = c("before last event in rect.", "CI before 1 precip event\n(nbinom fit)"),
       pch=c(3,15,NA), col=c("black", rgb(0.85,0.85,1), "black"), lty=c(NA,NA,2))
points(2:76,cumsum(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[2])]][76:1])[-1], pch="+", col="black")
text(x=77,y=quantile_05[89],srt=90,"50%", col="purple",cex=1.4)
text(x=77,y=12,srt=90,"80%", col="purple",cex=1.2)
text(x=77,y=14,srt=90,"90%", col="purple",cex=1.2)
text(x=77,y=16,srt=90,"95%", col="purple",cex=1.2)

#2nd part

plot(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[2])]],
     main=paste("Precipitation winter", Years_no_landsl_lot_precip[2]), ylab="day with precip > 75th perc",
     pch="+")
rect(xleft = 28.5, ybottom = 0.99, xright = 117.5, ytop = 1.01, border="blue")

plot(2:90,cumsum(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[2])]][117:28])[-1], pch="+", xlab="window before event",
     # ylim = c(0,max(number_events_before_landslide[,ls],12)),
     ylim = c(0,18), col=rgb(0.5,0.5,0.5),
     ylab = "Number of days precip>75th perc", main = paste("1963 2nd part"))

polygon(c(2:90,rev(2:90)),c(quantile_0975,rev(quantile_0025)),col=rgb(0.75,0.75,1),border=NA)
polygon(c(2:90,rev(2:90)),c(quantile_095,rev(quantile_005)),col=rgb(0.8,0.8,1),border=NA)
polygon(c(2:90,rev(2:90)),c(quantile_09,rev(quantile_01)),col=rgb(0.85,0.85,1),border=NA)
polygon(c(2:90,rev(2:90)),c(quantile_075,rev(quantile_025)),col=rgb(0.9,0.9,1),border=NA)

legend("topleft", legend = c("before last event in rect.", "CI before 1 precip event\n(nbinom fit)"),
       pch=c(3,15,NA), col=c("black", rgb(0.85,0.85,1), "black"), lty=c(NA,NA,2))
points(2:90,cumsum(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[2])]][117:28])[-1], pch="+", col="black")
text(x=90,y=quantile_05[89],srt=90,"50%", col="purple",cex=1.4)
text(x=90,y=12.5,srt=90,"80%", col="purple",cex=1.2)
text(x=90,y=15,srt=90,"90%", col="purple",cex=1.2)
text(x=90,y=17.5,srt=90,"95%", col="purple",cex=1.2)

#3rd part

plot(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[2])]],
     main=paste("Precipitation winter", Years_no_landsl_lot_precip[2]), ylab="day with precip > 75th perc",
     pch="+")
rect(xleft = 42.5, ybottom = 0.99, xright = 131.5, ytop = 1.01, border="blue")


plot(2:90,cumsum(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[2])]][131:42])[-1], pch="+", xlab="window before event",
     # ylim = c(0,max(number_events_before_landslide[,ls],12)),
     ylim = c(0,18), col=rgb(0.5,0.5,0.5),
     ylab = "Number of days precip>75th perc", main = paste("1963 3rd part"))

polygon(c(2:90,rev(2:90)),c(quantile_0975,rev(quantile_0025)),col=rgb(0.75,0.75,1),border=NA)
polygon(c(2:90,rev(2:90)),c(quantile_095,rev(quantile_005)),col=rgb(0.8,0.8,1),border=NA)
polygon(c(2:90,rev(2:90)),c(quantile_09,rev(quantile_01)),col=rgb(0.85,0.85,1),border=NA)
polygon(c(2:90,rev(2:90)),c(quantile_075,rev(quantile_025)),col=rgb(0.9,0.9,1),border=NA)

legend("topleft", legend = c("before last event in rect.", "CI before 1 precip event\n(nbinom fit)"),
       pch=c(3,15,NA), col=c("black", rgb(0.85,0.85,1), "black"), lty=c(NA,NA,2))
points(2:90,cumsum(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[2])]][131:42])[-1], pch="+", col="black")
text(x=90,y=quantile_05[89],srt=90,"50%", col="purple",cex=1.4)
text(x=90,y=12.5,srt=90,"80%", col="purple",cex=1.2)
text(x=90,y=15,srt=90,"90%", col="purple",cex=1.2)
text(x=90,y=17.5,srt=90,"95%", col="purple",cex=1.2)


#4th part
plot(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[2])]],
     main=paste("Precipitation winter", Years_no_landsl_lot_precip[2]), ylab="day with precip > 75th perc",
     pch="+")
rect(xleft = 62.5, ybottom = 0.99, xright = 151.5, ytop = 1.01, border="blue")

plot(2:90,cumsum(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[2])]][151:62])[-1], pch="+", xlab="window before event",
     # ylim = c(0,max(number_events_before_landslide[,ls],12)),
     ylim = c(0,18), col=rgb(0.5,0.5,0.5),
     ylab = "Number of days precip>75th perc", main = paste("1963 4th part"))

polygon(c(2:90,rev(2:90)),c(quantile_0975,rev(quantile_0025)),col=rgb(0.75,0.75,1),border=NA)
polygon(c(2:90,rev(2:90)),c(quantile_095,rev(quantile_005)),col=rgb(0.8,0.8,1),border=NA)
polygon(c(2:90,rev(2:90)),c(quantile_09,rev(quantile_01)),col=rgb(0.85,0.85,1),border=NA)
polygon(c(2:90,rev(2:90)),c(quantile_075,rev(quantile_025)),col=rgb(0.9,0.9,1),border=NA)

legend("topleft", legend = c("before last event in rect.", "CI before 1 precip event\n(nbinom fit)"),
       pch=c(3,15,NA), col=c("black", rgb(0.85,0.85,1), "black"), lty=c(NA,NA,2))
points(2:90,cumsum(list_precip_winter[[which(years_fullwinter==Years_no_landsl_lot_precip[2])]][151:62])[-1], pch="+", col="black")
text(x=90,y=quantile_05[89],srt=90,"50%", col="purple",cex=1.4)
text(x=90,y=12.5,srt=90,"80%", col="purple",cex=1.2)
text(x=90,y=15,srt=90,"90%", col="purple",cex=1.2)
text(x=90,y=17.5,srt=90,"95%", col="purple",cex=1.2)


plot(list_precip_winter_RAW[[which(years_fullwinter==1995)]],
     main=paste("Precipitation winter 1995"), ylab="precip (mm/day)")

plot(list_precip_winter[[which(years_fullwinter==1995)]],
     main=paste("Precipitation winter 1995"), ylab="day with precip > 75th perc", pch="+")




# Cumul precip before landslide -------------------------------------------


precip_w_before_landslide <- matrix(nrow = 91, ncol = length(dates_landslides))

for (ls in 1:length(dates_landslides)) {
  precip_w_before_landslide[1,ls] <- precip_GP_West[ref_date_ls[ls]]
  for (t_w in 2:91) {
    precip_w_before_landslide[t_w,ls] <- precip_GP_West[(ref_date_ls[ls]-t_w+1)]
  }
}#end for ls


Id_ls_sorted <- sort(apply(X=precip_w_before_landslide, FUN = sum, MARGIN = 2), index.return=T,decreasing = T)$ix
date_ls_sorted <- dates_landslides[Id_ls_sorted]
cumulated_precip_90days_prior_ls_sorted <- apply(X=precip_w_before_landslide, FUN = sum, MARGIN = 2)[Id_ls_sorted]
ls_type_sorted <- landslide_type[Id_ls_sorted]

Table_sorted_landslides <- cbind(1:length(Id_ls_sorted),date_ls_sorted,
                                 cumulated_precip_90days_prior_ls_sorted, ls_type_sorted)

colnames(Table_sorted_landslides) <- c("Rank", "Date", "cumul. precip 91 days prior (mm)",
                                       "Landslide type")

write.csv(Table_sorted_landslides,
          file = "/Users/admin/Documents/PrecipClustering_BottomUpAnalysis_LOCAL/landslides_sorted_cumul_precip.csv")






# Clustering test where no landslide --------------------------------------
detect_clustering <- as.Date(character())
count <- 0
count2 <- 0
for (id_y in 1:length(years_fullwinter)) {
  for (D in 1:length(list_precip_winter[[id_y]])) {
    if(list_precip_winter[[id_y]][D]==1){
      wind_length <- min(as.numeric(list_date[[id_y]][D]-list_date[[id_y]][1]), 90)
      
      test_in_ls <- 0
      for (DAY in (D-wind_length):D) {
        test_in_ls <- sum(dates_landslides == list_date[[id_y]][DAY]) + test_in_ls
      }#end for DAY
      
      if(test_in_ls==0 & wind_length>1){
        count <- count + 1
        
        if(sum(cumsum(list_precip_winter[[id_y]][(D-wind_length):(D-1)])[-1]>=quantile_095[1:(wind_length-1)] &
               cumsum(list_precip_winter[[id_y]][(D-wind_length):(D-1)])[-1]>quantile_09[1:(wind_length-1)])){
          count2 <- count2 + 1
          detect_clustering[count2] <- list_date[[id_y]][D]
          
          png(paste0("../PrecipClustering_BottomUpAnalysis_LOCAL/No_landslide_but_clusteredprecip/nb_events_before_event_",list_date[[id_y]][D],"_clustering_nolandslide.png"))
          
          plot(2:wind_length,cumsum(list_precip_winter[[id_y]][(D-wind_length):(D-1)])[-1], pch="+", xlab="window before event",
               # ylim = c(0,max(number_events_before_landslide[,ls],12)),
               ylim = c(0,18), col=rgb(0.5,0.5,0.5),
               ylab = "Number of days precip>75th perc", main = paste("Before", list_date[[id_y]][D]))

          polygon(c(2:wind_length,rev(2:wind_length)),c(quantile_0975[1:(wind_length-1)],rev(quantile_0025[1:(wind_length-1)])),
                  col=rgb(0.75,0.75,1),border=NA)
          polygon(c(2:wind_length,rev(2:wind_length)),c(quantile_095[1:(wind_length-1)],rev(quantile_005[1:(wind_length-1)])),
                  col=rgb(0.8,0.8,1),border=NA)
          polygon(c(2:wind_length,rev(2:wind_length)),c(quantile_09[1:(wind_length-1)],rev(quantile_01[1:(wind_length-1)])),
                  col=rgb(0.85,0.85,1),border=NA)
          polygon(c(2:wind_length,rev(2:wind_length)),c(quantile_075[1:(wind_length-1)],rev(quantile_025[1:(wind_length-1)])),
                  col=rgb(0.9,0.9,1),border=NA)

          legend("topleft", legend = c("before last event in rect.", "CI before 1 precip event\n(nbinom fit)"),
                 pch=c(3,15,NA), col=c("black", rgb(0.85,0.85,1), "black"), lty=c(NA,NA,2))
          points(2:wind_length,cumsum(list_precip_winter[[id_y]][(D-wind_length):(D-1)])[-1], pch="+", col="black")
          
          dev.off()
        }
      

      }
      
    }#end if event
  }#end for D
}#end for id_y
