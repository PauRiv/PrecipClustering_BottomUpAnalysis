## The Ripley's K function adapted for one dimension
## Clément Chevalier - 2014-02-01
## Pauline Rivoire - 27-11-2020

RipleyK <- function( Listvectors, size , xvect  ){
  
  # Listvectors : a list containing vectors of zeros and ones. Each data in a vector represents one day.
  # 0 => no extreme; 1 => extreme event
  # size: the total number of vectors
  # xvect: all the length (in days) of the considered window.
  
  # The function calculates the average number of events in the interval [ti-x, ti+x]
  # where ti the time of a randomly extreme event (i.e. at time ti, we have a "1")
  
  # First step, get informations of extreme events
  yearExtreme <- NULL
  dayExtreme <- NULL
  length.year <- rep(0,times=size)
  
  for(i in 1:size){
    year.i <- Listvectors[[i]]
    n.i <- length(year.i)
    length.year[i] <- n.i
    index <- which(x=as.logical(year.i),arr.ind=TRUE)
    if(length(index)!=0){
      yearExtreme <- c(yearExtreme,rep(i,times=length(index)))
      dayExtreme <- c(dayExtreme,index)
    }
  }
  n.events <- length(yearExtreme)
  
  #second step, we calculate the kfunction
  
  bigres <- rep(0,times=length(xvect))
  
  for(i in 1:n.events){
    year.myevent <- yearExtreme[i]
    day.myevent <- dayExtreme[i]
    index <- which(yearExtreme==year.myevent)
    if( length(index)!= 1 ){
      #there are other events this year
      dist_in_days <- sort(abs(dayExtreme[index] - day.myevent))[-1]
      
      tmp <- findInterval(x=dist_in_days,vec=xvect)
      res.i <- rep(0,times=length(xvect))
      for(k in 1:length(tmp)) res.i[tmp[k]+1] <- res.i[tmp[k]+1]+1
      res.i <- cumsum(res.i)
      bigres <- bigres+res.i
    }
  }
  bigres <- bigres/n.events
  return(bigres)
}

SimuRipleyK <- function(nsimu, xvect, size , yearsize, n.events.avg ){
  total.size <- sum(yearsize)
  proba <- n.events.avg / total.size
  
  result1 <- matrix(c(0),nrow=nsimu,ncol=length(xvect))
  
  for(i in 1:nsimu){
#    print(paste("simu",i,"out of",nsimu))
    Listvectors <- list()
    for(j in 1:size) Listvectors[[j]] <- rbinom(n=yearsize[j],size=1,prob=proba)
    for (tw in 1:length(xvect)) { #Modif 3 lines
      TW <- xvect[tw]
      result1[i,tw] <- RipleyK(Listvectors=Listvectors,size=size,xvect=TW)[1]
    }#end for tw
    
  }
  
  result.avg <- apply(X=result1,MARGIN=2,FUN=mean)
  result.upper <- apply(X=result1,MARGIN=2,FUN=quantile,probs=0.95) # probs=0.975
  result.lower <- apply(X=result1,MARGIN=2,FUN=quantile,probs=0.05) # probs=0.025
  
  return(list(avg=result.avg,upper=result.upper,lower=result.lower,allsimu=result1))
  
}


### RIPLEY'S K EXAMPLE OF USE ###
#################################

# mydata <- list()
# for(i in 1:40) mydata[[i]] <- rbinom(n=91,size=1,prob=0.04)
# Listvectors <- mydata
# size <- 40 #number of year
# xvect <- c(1:91) #One season
# 
# y <- RipleyK(Listvectors=Listvectors,size=size,xvect=xvect)
# 
# obj <- SimuRipleyK(nsimu=1000,xvect=xvect,size=40,yearsize=rep(91,times=40),n.events.avg=140)
# 
# plot(x=xvect,y=y,type="l",ylab="Ripley",lwd=3,ylim=c(0,5),xlab="days")
# lines(x=xvect,y=obj$avg,col="red")
# lines(x=xvect,y=obj$lower,col="grey",lty=2)
# lines(x=xvect,y=obj$upper,col="grey",lty=2)

# mydata <- list()
# mydata[[1]] <- rbinom(n=18262,size=1,prob=0.01)
# Listvectors <- mydata
# size <- 1
# xvect <- c(1:18262)
# 
# y <- RipleyK(Listvectors=Listvectors,size=size,xvect=xvect)
# 
# obj <- SimuRipleyK(nsimu=10,xvect=xvect,size=1,yearsize=rep(18262,times=1),n.events.avg=length(which(Listvectors[[1]]>0)))
# 
# plot(x=xvect,y=y,type="l",ylab="Ripley",lwd=3,xlab="days") #ylim=c(0,5)
# lines(x=xvect,y=obj$avg,col="red")
# lines(x=xvect,y=obj$lower,col="grey",lty=2)
# lines(x=xvect,y=obj$upper,col="grey",lty=2)
