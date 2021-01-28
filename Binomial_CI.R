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



# Precip time series : compute events -------------------------------------

#we exclude the two first rows and the last one because they don't contain data
All_years <- unique(date_matrix[,"year"])
years_fullwinter <- All_years[-c(1,length(All_years))] #Exclude first and last winter because not complete

nrow(date_matrix)==length(precip_GP_West)#check if same dimension

frequency(precip_GP_West[precip_GP_West>1])



TIES <- rownames(table(precip_GP_West))[which(table(precip_GP_West)>1)]

precip_wo_ties <- precip_GP_West

# set.seed(2021)
# for (tie in 17:length(TIES)) {
#   where_to_change <- which(precip_GP_West==TIES[tie])
#   precip_wo_ties[where_to_change] <- precip_GP_West[where_to_change] + (runif(n = length(where_to_change))-0.5)*0.25
# }#end for tie

# cbind(precip_wo_ties, precip_GP_West)


quantile(precip_GP_West[precip_GP_West>1], probs = 0.95)


quant_precip_75 <- 10.5
quant_precip_90 <- 18
quant_precip_95 <- 24

Q_chosen <- quant_precip_95
Qname <- "95th"


#extract seasonal precip
precip_extended_winter <- precip_GP_West
days_to_exclude <- which(date_matrix[,"month"] %in% 4:10) #exclude days from April to October included
precip_extended_winter[days_to_exclude] <- rep(0, length(days_to_exclude))
precip_extended_winter2 <- precip_GP_West
precip_extended_winter2[days_to_exclude] <- rep(NA, length(days_to_exclude))
winter_precip_event <- as.numeric(precip_extended_winter2>Q_chosen)
#set those excluded days to 0 (can also be set to NA if needed)

cbind(date_matrix[,"month"], precip_GP_West, precip_extended_winter, winter_precip_event)


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
list_precip_winter <- list()
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
  
  list_precip_winter[[count]] <- as.matrix(as.numeric(precip_per_winter(year_to_extract = year,
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

acf(unlist(list_precip_winter), lag.max = 10, main="All winters in one time series")

for (Y in 1:length(list_precip_winter)) {
  
  acf(list_precip_winter[[Y]], lag.max = 10, main=paste("winter",years_fullwinter[Y]))
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

ext_index <- as.numeric(extremalindex_mod(x,djf_threshold,max_lag=90)[1])
avg_cluster_size <- 1/ext_index

# declustering
y = decluster(precip_extended_winter,Q_chosen,r=2)
inds = which(x>=10.5)
w = attr(y,"clusters") # for each threshold exceedance, w indicates the cluster it belongs to
declustered_precip <-as.numeric(y > 10.5)
# cbind(x,y)[670:1170,]
# 
# # Another method of declustering:
# library(evd)
# CL <- clusters(data = precip_extended_winter, u=Q_chosen,r=2, cmax = T)
# EVENTs_evd <- numeric(length = length(precip_GP_West))
# EVENTs_evd[as.numeric(names(CL))] <- rep(1, length(CL))
# 
# cbind(precip_extended_winter2, declustered_precip, EVENTs_evd)[670:1170,]
# sum(declustered_precip[which(EVENTs_evd==1)]==0)
# sum(EVENTs_evd[which(declustered_precip==1)]==0)

Decl_precip <- cbind(date_matrix, declustered_precip)
# save(Decl_precip, file = "Data/declustered_precipitation.Rdata")

# load("Data/declustered_precipitation.Rdata")

List_precip_events_winter01 <- list()
count <- 0
for (year in years_fullwinter) {
  count <- count + 1
  precip_winter <- as.matrix(as.numeric(precip_per_winter(year_to_extract = year,
                                                          precip_events = declustered_precip,
                                                          days = date_matrix[,"day"],
                                                          months = date_matrix[,"month"],
                                                          years = date_matrix[,"year"])[,"precip"]))
  colnames(precip_winter) <- paste("winter Nov.",as.character(year))
  List_precip_events_winter01[[count]] <- precip_winter
}#end for year
rm(count)





# Model a binomial --------------------------------------------------------

P_bin <- sum(unlist(List_precip_events_winter01))/length(unlist(List_precip_events_winter01))


rbinom(p = P_bin, n = 100, size = 30)
qbinom(prob = P_bin, size = 2:90, p=0.99)

quantile_80 <- qbinom(prob = P_bin, size = 2:90, p=0.80)
quantile_90 <- qbinom(prob = P_bin, size = 2:90, p=0.90)
quantile_95 <- qbinom(prob = P_bin, size = 2:90, p=0.95)
quantile_975 <- qbinom(prob = P_bin, size = 2:90, p=0.975)
zeros <-rep(0,length(2:90))

# Precip events before landslide ------------------------------------------

excel_landslides <- readxl::read_xlsx(path = "../PrecipClustering_BottomUpAnalysis_LOCAL/Data/Landslides_Date.xlsx")
dates_landslides <- as.matrix(excel_landslides)
dates_landslides <- rbind("1958-12-19", dates_landslides)

ref_date_ls <- numeric()
for (ls in 1:length(dates_landslides)) {
  ref_date_ls[ls] <- which(date_matrix[,1]==year(dates_landslides[ls]) &
                             date_matrix[,2]==month(dates_landslides[ls]) &
                             date_matrix[,3]==day(dates_landslides[ls]))
}#end for ls

number_events_before_landslide <- matrix(data=NA, nrow = length((2:90)),
                                         ncol = length(dates_landslides))
rownames(number_events_before_landslide) <- as.character((2:90))
acc_precip_before_ls <- matrix(nrow = length(2:90), ncol = length(dates_landslides))

for (t_w in 1:length((2:90))) {
  w_length <- (2:90)[t_w]
  
  precip_w_before <- matrix(data = NA, nrow = w_length, ncol = length(ref_date_ls))
  events_w_before <- matrix(data = NA, nrow = w_length, ncol = length(ref_date_ls))
  
  for (ls in 1:length(dates_landslides)) {
    precip_w_before[,ls] <- precip_GP_West[(ref_date_ls[ls]-w_length+1):ref_date_ls[ls]]
    acc_precip_before_ls[t_w,ls] <- sum(precip_w_before[,ls])
    precip_nocorr = decluster(precip_w_before[,ls],Q_chosen,r=2)
    events_w_before[,ls] <- as.numeric(precip_nocorr > 10.5)
    number_events_before_landslide[t_w,ls] <- sum(events_w_before[,ls])
  }#end for ls
}#end for t_w


win_critical_duration <- c(10,10,1,30,15,75,5,1,1,15,15,30,40,60,75,40,90,60,60,4,10,40,1)

landslide_type <- c("S","S","S","D","S","D","S","S","S","S","S","D","D",
                    "D","D","D","D","D","D","S","S","D","S")
max_nb_ev <- max(number_events_before_landslide)
for (ls in 1:23) {
  # png(paste0("../PrecipClustering_BottomUpAnalysis_LOCAL/nb_events_before_landslide_",dates_landslides[ls],"_fitting.png"))
  plot(2:90,number_events_before_landslide[,ls], pch="+", xlab="window before event",
       # ylim = c(0,max(number_events_before_landslide[,ls],12)),
       ylim = c(0,max_nb_ev), col=rgb(0.5,0.5,0.5),
       ylab = paste0("Number of days precip>",Qname,"perc"),
       main = paste0("Lanslide ", dates_landslides[ls],", type=",landslide_type[ls]))
  
  
  polygon(c(2:90,rev(2:90)),c(quantile_975,rev(zeros)),col=rgb(0.75,0.75,1),border=NA)
  polygon(c(2:90,rev(2:90)),c(quantile_95,rev(zeros)),col=rgb(0.8,0.8,1),border=NA)
  polygon(c(2:90,rev(2:90)),c(quantile_90,rev(zeros)),col=rgb(0.85,0.85,1),border=NA)
  polygon(c(2:90,rev(2:90)),c(quantile_80,rev(zeros)),col=rgb(0.9,0.9,1),border=NA)
  abline(v=win_critical_duration[ls], lty=2)
  legend("topleft", legend = c("before landslide", "CI before 1 precip event\n(binom)","critical duration"),
         pch=c(3,15,NA), col=c("black", rgb(0.85,0.85,1), "black"), lty=c(NA,NA,2))
  points(2:90,number_events_before_landslide[,ls], pch="+", col=rgb(0.5,0.5,0.5))
  # points(names(significant[significant]),number_events_before_landslide[significant,ls], pch="+", lwd=3)
  # text(x=90,y=quantile_05[89],srt=90,"50%", col="purple",cex=1.4)
  # text(x=90,y=12.5,srt=90,"80%", col="purple",cex=1.2)
  # text(x=90,y=15,srt=90,"90%", col="purple",cex=1.2)
  # text(x=90,y=17.5,srt=90,"95%", col="purple",cex=1.2)
  # dev.off()
}#end for ls

