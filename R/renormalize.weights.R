renormalize.weights <-
function(weights,log=FALSE) {
  # If provided log-weights, exponentiate.
  if (log) weights = exp(weights-max(weights)) # max(weights) for numerical stability
  
  return(weights/sum(weights))
}

