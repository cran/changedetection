#' Sample dataset
#'
#' This is artificial data following changing multivariate liner structure.
#'
#' @format A data frame with 600 rows and 8 columns:
#' \itemize{
#'   \item{\eqn{y_1, y_2, y_3}}{ response variables}
#'   \item{\eqn{x_1, x_2, x_3, x_4, x_5}}{ predictors}
#' }
#'
#'
#' @details
#' The underlying model is divided into 4 parts:
#' \tabular{llllllll}{
#' 1) for rows 1-59:\cr
#'       \tab \eqn{y_1 =  x_1 +  x_2  + x_3}\cr
#'       \tab \eqn{y_2 =  x_1 + 3x_2  + x_3}\cr
#'       \tab \eqn{y_3 = 3x_1 + 3x_2  + x_3}\cr
#' 2) for rows 60-299:\cr
#'       \tab \eqn{y_1 =  x_1 + 3x_2  + x_3}\cr
#'       \tab \eqn{y_2 =  x_1 + 3x_2  + x_3}\cr
#'       \tab \eqn{y_3 = 3x_1 + 3x_2  + x_3}\cr
#' 3) for rows 300-479:\cr
#'       \tab \eqn{y_1 =  x_1 + 3x_2  + x_3}\cr
#'       \tab \eqn{y_2 = 3x_1 + 3x_2  + x_3}\cr
#'       \tab \eqn{y_3 = 3x_1 + 3x_2  + x_3}\cr
#' 4) for rows 480-600:\cr
#'       \tab \eqn{y_1 =  x_1 + 3x_2  + x_3}\cr
#'       \tab \eqn{y_2 = 3x_1 + 3x_2  + x_3}\cr
#'       \tab \eqn{y_3 = 5x_1 + 3x_2  + x_3}\cr
#'       }
#'
#' By construction, the \code{sampledataset} has three structural changes happened at points 60, 300 and 480, which makes it useful to be applied to various functions in the current package.
#'
#' @name sampledataset

"sampledataset"