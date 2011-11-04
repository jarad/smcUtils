ess.weights <-
function(weights) {
  return(1/(sum(weights^2)))
  
  # Could compute this as 
  # return(length(weights)/(1+cov.weights(weights)))
  # as in Liu and Chen 1995.
}

