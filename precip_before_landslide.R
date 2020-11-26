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



precip_90days_before <- matrix(data = NA, nrow = 90, ncol = length(ref_date_ls))
precip_60days_before <- matrix(data = NA, nrow = 60, ncol = length(ref_date_ls))
precip_30days_before <- matrix(data = NA, nrow = 30, ncol = length(ref_date_ls))

for (ls in 1:length(dates_landslides)) {
  precip_90days_before[,ls] <- precip_GP_West[(ref_date_ls[ls]-89):ref_date_ls[ls]]
  precip_60days_before[,ls] <- precip_GP_West[(ref_date_ls[ls]-59):ref_date_ls[ls]]
  precip_30days_before[,ls] <- precip_GP_West[(ref_date_ls[ls]-29):ref_date_ls[ls]]
}#end for ls

quantile_75 <- 10.5

events_90_days_before <- matrix(data = NA, nrow = 90, ncol = length(ref_date_ls))
events_60_days_before <- matrix(data = NA, nrow = 60, ncol = length(ref_date_ls))
events_30_days_before <- matrix(data = NA, nrow = 30, ncol = length(ref_date_ls))

for (ls in 1:length(dates_landslides)) {
  gt_1 <- as.numeric(precip_90days_before[,ls]>quantile_75)
  gt_1_nocorr <- gt_1
  for (day in 2:90) {
    if(gt_1_nocorr[day-1] == 1 & gt_1_nocorr[day] == 1){gt_1_nocorr[day] <- 0}
  }#end for day
  events_90_days_before[,ls] <- gt_1_nocorr
  
  gt_1 <- as.numeric(precip_60days_before[,ls]>quantile_75)
  gt_1_nocorr <- gt_1
  for (day in 2:60) {
    if(gt_1_nocorr[day-1] == 1 & gt_1_nocorr[day] == 1){gt_1_nocorr[day] <- 0}
  }#end for day
  events_60_days_before[,ls] <- gt_1_nocorr
  
  gt_1 <- as.numeric(precip_30days_before[,ls]>quantile_75)
  gt_1_nocorr <- gt_1
  for (day in 2:30) {
    if(gt_1_nocorr[day-1] == 1 & gt_1_nocorr[day] == 1){gt_1_nocorr[day] <- 0}
  }#end for day
  events_30_days_before[,ls] <- gt_1_nocorr
}#end for ls

dim(events_90_days_before)
nb_events_90 <- apply(X=events_90_days_before, MARGIN = 2, FUN = sum)
nb_events_60 <- apply(X=events_60_days_before, MARGIN = 2, FUN = sum)
nb_events_30 <- apply(X=events_30_days_before, MARGIN = 2, FUN = sum)


cumul_precip_90 <- apply(X=precip_90days_before, MARGIN = 2, FUN = sum)
cumul_precip_60 <- apply(X=precip_60days_before, MARGIN = 2, FUN = sum)
cumul_precip_30 <- apply(X=precip_30days_before, MARGIN = 2, FUN = sum)


plot(nb_events_90,cumul_precip_90)
plot(nb_events_60,cumul_precip_60)
plot(nb_events_30,cumul_precip_30)
 
load("Fit/Fit_negative_binomial.Rdata")

plot(ecdf(nb_events_30))
lines(0:12,pnbinom(0:12, size = Fit_neg_bino$fit_w30$estimate["size"], mu = Fit_neg_bino$fit_w30$estimate["mu"]), col="red")
legend("bottomright",
       legend = c("30-wind. before landslide", "30-wind. anytime in winter"),
       text.col=c("black", "red"))

plot(ecdf(nb_events_60))
lines(0:20,pnbinom(0:20, size = Fit_neg_bino$fit_w60$estimate["size"], mu = Fit_neg_bino$fit_w60$estimate["mu"]), col="red")
legend("bottomright",
       legend = c("60-wind. before landslide", "60-wind. anytime in winter"),
       text.col=c("black", "red"))

plot(ecdf(nb_events_90))
lines(0:22,pnbinom(0:22, size = Fit_neg_bino$fit_w60$estimate["size"], mu = Fit_neg_bino$fit_w60$estimate["mu"]), col="red")
legend("bottomright",
       legend = c("90-wind. before landslide", "90-wind. anytime in winter"),
       text.col=c("black", "red"))
