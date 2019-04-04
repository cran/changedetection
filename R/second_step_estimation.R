recursiveDetection <- function(dataset, s, e, tao, responseIndex,settings){

  # get the center of max energy distance
  maxCenter <- getMaximumCenter(dataset,s,e,tao,responseIndex,settings)

  s <- s + (maxCenter-1)*tao
  e <- s + 2*tao
  #drill down
  updatedTao <- floor(tao*settings$gamma)
  if (updatedTao > 10*(length(dataset[1,])-settings$q)){
    #if there is still enough variables to estimate models well
    #rule: 10 observations per variable (Frank Harrell's book, Regression Modeling Strategies)
    recursiveDetection(dataset,s,e,updatedTao,responseIndex,settings)
  } else {
    #exact detection
    result <- getPreciseChangeLocation(dataset,s,tao,responseIndex,settings)
  }
}

getMaximumCenter<-function(dataset,s,e,tao,responseIndex,settings){
  # data length
  T <- e-s
  # number of periods to observe
  settings$numberOfPeriods <- round(T/tao)
  cat("Number of periods = ",settings$numberOfPeriods,"\n")
  stats = rep(0.0,settings$numberOfPeriods-1)
  for (i in 1:(settings$numberOfPeriods-1))
  {
    # get updated left hand variable indexes
    left = dataset[((i-1)*tao+1+s):(i*tao+s),]

    #update variable indexes
    indexes <- getVariableIndexes(left,settings,responseIndex)
    variables <- c(1:settings$q,indexes+settings$q)

    #cat(c("variable indexes: ", paste(variables,sep=" "),"\n"))

    left <- dataset[((i-1)*tao+1+s):(i*tao+s),variables]
    right <- dataset[(i*tao+1+s):((i+1)*tao+s),variables]

    stats[i] <- calculateEnergyDistance(left, right,settings, responseIndex)

    #cat(c("Periods [",(i-1)*tao+1+s,i*tao+s,"] and [",i*tao+1+s,(i+1)*tao+s,"], StatValue=",stats[i],"\n"))
  }

  maximumIndex <- which.max(stats)
  return (maximumIndex)

}

#point - location of a point indicating small area of a change
getPreciseChangeLocation <-function(dataset,s,tao,responseIndex,settings)
{
  errors = rep(0,2*tao)
  leftBoundary <- s

  cat(c("Left boundary=",leftBoundary,"\n"))

  for (i in leftBoundary:(leftBoundary+2*tao))
  {
    # get updated left hand variable indexes
    if (i-tao<0) start=1
    else start=i-tao


    left = dataset[start:i,]

    #update variable indexes
    indexes <- getVariableIndexes(left,settings,responseIndex)
    variables <- c(1:settings$q,indexes+settings$q)

    #cat(c("variable indexes: ", paste(variables,sep=" "),"\n"))

    left = dataset[start:i,variables]
    right = dataset[(i+1):(3*tao+leftBoundary),variables]


    errors[i-leftBoundary+1] = calculateEnergyDistance(left, right,settings, responseIndex)


    #cat(c("Periods [",i-tao,i,"] and [",i+1,3*tao+leftBoundary,"], energy=",errors[i-leftBoundary+1],"\n"))

  }

  #find the minimum energy distance
  result <- which.max(errors)
  cat(c("Found optimal index within the period ",result,"\n"))
  result = result + s -1
  return (result)
}
