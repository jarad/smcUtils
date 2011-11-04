cov.weights <-
function(weights) {
  return(var(weights)/mean(weights)^2)
  
  # Could compute this as
  # return(mean((length(weights)*weights-1)^2))
  # as in Liu and Chen 1995.
}

