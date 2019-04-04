
# Estimate multivariate linear model of data.
#
# @param data multivariate dataset (matrix).
# @param settings algorithm settings.
# @return set of linear models.

trainMultivariateModel <- function(data, settings)
{
  # resulting array of models
  models <- matrix(0.0, nrow = settings$q, ncol = length(data[1,])-settings$q+1)
  #all x from (q+1) column to the end of the matrix
  X <- data[,(settings$q+1):ncol(data)]

  for (i in 1:settings$q)
  {
    # i-th column is i-th response
    Y <- data[,i]
    model <- l1fit(X, Y, intercept = TRUE)
    models[i,] <- model$coefficients
  }

  return (models)
}


# Estimate univariate linear model of data.
#
# @param data multivariate dataset (matrix).
# @param i the index of a response to be analyzed.
# @param settings algorithm settings.
# @return a linear model
trainUnivariateModel <- function(data,i,settings)
{
  if (missing(i)){i=1}
  Y <- data[,i] #i-th column is y
  X <- data[,(settings$q+1):ncol(data)]

  model <- l1fit(X, Y, intercept = TRUE)

  return (model)
}

# Calculate a vector of multivariate model residuals.
#
# @param data multivariate dataset (matrix).
# @param models matrix of linear models coefficients.
# @param settings algorithm settings.
# @return resulting vector of residuals
getResidualsMultivariate <-function(data, models, settings)
{
  matrixOfResiduals <- getMatrixOfResiduals(data, models, settings)
  result <- vector()
  for (i in 1:length(matrixOfResiduals[1,]))
  {
    result[i] <- sum(abs(matrixOfResiduals[,i])^2)^(1/2)#sum(matrixOfResiduals[,i])
  }
  return (result)
}

# Calculate a vector of univariate model residuals.
#
# @param data multivariate dataset (matrix).
# @param model vactor of linear models coefficients.
# @param settings algorithm settings.
# @param i the index of a response to be analyzed.
# @return resulting vector of residuals
getResidualsUnivariate <-function(data,model,i,settings)
{
  if (missing(i)){i=1}
  # get independent variables
  X <- data[,(settings$q+1):ncol(data)]

  # multiply weights by variables
  predictions <- as.matrix(X) %*% model$coefficients[2:(length(data[1,])-settings$q+1)]  + model$coefficients[1]

  # calculate absolute residuals
  result <- abs(data[,i]-predictions)

  return(result)

}


# Calculate matrix of multivariate model residuals.
#
# @param data multivariate dataset (matrix).
# @param models matrix of linear models coefficients
# @param settings algorithm settings.
# @return resulting matrix of residuals
getMatrixOfResiduals <-function(data, models, settings)
{
  # prepare the matrix variable
  result <- matrix(0.0, nrow = settings$q, ncol = length(data[,1]))
  # all x from (q+1) column to the end of the matrix
  X <- data[,(settings$q+1):ncol(data)]

  for (i in 1:settings$q)
  {
    #multiply weights by variables
    predictions <- as.matrix(X) %*% models[i,2:length(models[1,])]  + models[i,1]

    #calculate absolute residuals
    result[i,] <- abs(data[,i]-predictions)
  }

  return (result)
}

# Estimate variable indexes for the given dataset.
#
# @param data multivariate dataset (matrix).
# @param settings algorithm settings.
# @param index index of y if univariate estimation needed (default value is -1)
# @return resulting list of variables
getVariableIndexes <- function(data,settings,index = -1)
{

  # extract x matrix
  XX <- data[,(settings$q+1):ncol(data)]
  X <- as.matrix(XX)
  # max number of variables to select
  l <- settings$l

  if (length(X[1,])!=l)
  {
    allIndexes <- 0
    # define the responses to analyze
      rangeOfResponses <- c(1:settings$q)
    if (index!=-1)
      rangeOfResponses <- c(index)

    # iterate over the range
    for (i in rangeOfResponses)
    {
      # extract single response (i.e. i-th column (i-th y))
      Y <- data[,i]
      # run Lasso
      glmnet1 <- cv.glmnet(X,Y)
      ##ic.glmnet(X, Y, crit = "bic", alpha = 1, nlambda = 1000)

      # estimated coefficients
      c<-coef(glmnet1,s='lambda.min',exact=TRUE)

      # get rid of bias
      c<-c[-1]

      # select non-zero contributed variables
      index <- which(c!=0, arr.ind = TRUE)
      values <- c[index]

      # sort absolute contributions od all vars
      orderedArray <- sort(abs(c), decreasing=T)

      # get rid of excess variables (by the setup count l)
      if (length(index)>l && orderedArray[l]!=0)
      {
        boundary <- orderedArray[l]
        index <- which(abs(c) >= boundary, arr.ind=TRUE)
        values <- c[abs(c)>= boundary]
      }

      # conctenate with results for previous responses
      allIndexes <- c(allIndexes,index)
    }

    # get rid of first zero
    allIndexes<-allIndexes[-1]

    # extract unique indexes out of all estimated
    uniqueIndexes <- unique(allIndexes)

    result <- sort(uniqueIndexes)
  }
  else
    result <- c(1:length(X[1,]))


  return (result)
}





