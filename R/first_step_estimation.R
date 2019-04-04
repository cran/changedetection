# First step estimation. Get the suspected change areas centers
#
# @param dataset A matrix containing both responses and predictors.
# @param settings algorithm settings.
#
# @return The set of approximate change centers indexes \code{flags}.

getInitialFlags <- function(dataset, settings){

  tao <-settings$tao

  #estimate valuable variables for the whole dtaset per each period separately
  settings$variables <- estimateValuableVariables(dataset,settings)
  flags = rep(FALSE,settings$numberOfPeriods-1)
  stats = rep(0.0,settings$numberOfPeriods-1)

  # assess the statistics flags for all neighbour pairs of periods
  for (i in 1:(settings$numberOfPeriods-1))
  {
    left <- dataset[((i-1)*tao+1):(i*tao),settings$variables[[i]]]
    right <- dataset[(i*tao+1):((i+1)*tao),settings$variables[[i]]]

    flags[i] <- calculateNPTestResult(left, right, settings)
    stats[i] <- calculateEnergyDistance(left, right, settings)


    #cat(c("Periods [",(i-1)*tao+1,i*tao,"] and [",i*tao+1,(i+1)*tao,"], Flag=",flags[i],"; StatValue=",stats[i],"\n"))
  }

  # exclude excess flags
  for (i in 2:(settings$numberOfPeriods-1)){
    # if a test showed a change in two consequential points
    if(flags[i]&&flags[i-1]){
      # if the current change is stronger than previous
      if (stats[i]> stats[i-1]){
        flags[i-1] <- FALSE
      } else { #otherwise
        flags[i] <- FALSE
      }
    }
  }

  return (flags)
}

# First step estimation. Get the suspected change areas centers
#
# @param dataset A matrix containing both responses and predictors.
# @param changes a set of approximate change points.
# @param settings algorithm settings.
#
# @return The set of response indexes having biggest changes.

getResponseIndexes <- function(dataset,changes,settings)
{
  result <- vector()
  tao <- settings$tao

  for (i in 1:(length(changes)))
  {
    # get updated left hand variable indexes
    left = dataset[(changes[i]*tao - tao):(changes[i]*tao-1),]

    #update valuable variables
    indexes<-getVariableIndexes(left,settings)
    variables <- c(1:settings$q,indexes+settings$q)

    # cut relevant data pieces
    left = dataset[(changes[i]*tao - tao):(changes[i]*tao-1),variables]
    right = dataset[(changes[i]*tao):(changes[i]*tao+tao),variables]
    stats = rep(0.0,settings$q)

    for (j in 1:settings$q)
    {
      stats[j] = calculateEnergyDistance(left,right,settings,j)
    }

    result <- c(result,which.max(stats))
  }

  return (result)
}

#estimate valuable variables for the whole dtaset per each period separately
estimateValuableVariables <- function(dataset, settings)
{
  cat("Variable selection step \n")
  finalIndexes <- NULL
  tao <- settings$tao

  variables <- vector("list", length = settings$numberOfPeriods)

  for (i in 1:settings$numberOfPeriods)
  {
    cat(c("Period ",i,"[",((i-1)*tao+1),",",(i*tao),"]:"))

    left = dataset[((i-1)*tao+1):(i*tao),]

    indexes<-getVariableIndexes(left,settings)
    cat(c(" variables: ", paste(indexes,sep=" "),"\n"))

    # add response indexes to the final set of variables
    variables[[i]] <- c(1:settings$q,indexes+settings$q)

  }

  return (variables)

}
