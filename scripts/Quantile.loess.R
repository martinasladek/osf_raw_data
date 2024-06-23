# The code was published on: http://www.r-statistics.com/2010/04/quantile-lowess-combining-a-moving-quantile-window-with-lowess-r-function
# Based on the algorithm published here: http://www.e-publications.org/ims/submission/index.php/AOAS/user/submissionFile/4295?confirm=37ca4b72
# http://en.wikipedia.org/wiki/Moving_average

# install.packages("zoo")
library(zoo) # rollmean
#?rollapply

Quantile.loess	 <- function(Y, X = NULL, 
                            number.of.splits = NULL,
                            window.size = 20,
                            percent.of.overlap.between.two.windows = NULL,
                            the.distance.between.each.window = NULL,
                            the.quant = .95,
                            window.alignment = c("center"), 
                            window.function = function(x) {quantile(x, the.quant)},
                            ...)
{
  # input: Y and X, and smothing parameters
  # output: new y and x
  
  # Extra parameter "..." goes to the loess	
  
  # window.size ==  the number of observation in the window (not the window length!)
  
  # "number.of.splits" will override "window.size"
  # let's compute the window.size:	
  if(!is.null(number.of.splits)) {window.size <- ceiling(length(Y)/number.of.splits)}
  
  # If the.distance.between.each.window is not specified, let's make the distances fully distinct
  if(is.null(the.distance.between.each.window)) {the.distance.between.each.window <- window.size}
  
  # If percent.of.overlap.between.windows is not null, it will override the.distance.between.each.window 
  if(!is.null(percent.of.overlap.between.two.windows)) 
  {
    the.distance.between.each.window <- window.size * (1-percent.of.overlap.between.two.windows)
  }
  
  
  
  # loading zoo
  if(!require(zoo)) 	
  {
    print("zoo is not installed - please install it.")
    install.packages("zoo")
  }
  
  
  
  if(is.null(X)) {X <- index(Y)} # if we don't have any X, then Y must be ordered, in which case, we can use the indexes of Y as X.
  
  # creating our new X and Y
  zoo.Y <- zoo(x = Y, order.by = X)
  #zoo.X <- attributes(zoo.Y)$index
  
  new.Y <- rollapply(zoo.Y, width = window.size, 
                     FUN = window.function,
                     by = the.distance.between.each.window,
                     align = window.alignment)
  new.X <- attributes(new.Y)$index
  
  new.Y.loess <- loess(new.Y~new.X, family = "sym", ...)$fitted 
  
  return(list(y = new.Y, x = new.X, y.loess = new.Y.loess))
}









moving.quantile.demo <- function()
{
  # example
  
  data(airquality)
  attach(airquality)
  
  no.na <- !(is.na(Ozone) | is.na(Temp))
  Ozone.2 <- Ozone[no.na]
  Temp.2 <- jitter(Temp[no.na])
  
  plot(Ozone ~ Temp, main = "Predicting the 95% Ozone level accourding Temperature")
  
  # fitting the Quantile regression
  require(quantreg)
  abline(rq(Ozone ~ Temp, tau = .95), col = "red")
  
  # fitting the Quantile LOWESS
  QL <- Quantile.loess(Y = Ozone.2, X = Temp.2, 
                       the.quant = .95,
                       window.size = 5,
                       window.alignment = c("center"))
  points(QL$y.loess ~ QL$x, type = "l", col = "green")
  
  legend("topleft",legend = c("95% Quantile regression", "95% Quantile LOWESS"), fill = c("red","green"))
  
}

# moving.quantile.demo()

