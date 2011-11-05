renormalize.weights <-
function(weights, log=FALSE, useC=FALSE) {
  if (useC) return(.C("renormalize_wrap", 
                      weights=as.double(weights),
                      as.integer(length(weights)),
                      as.integer(log))$weights)

  # If provided log-weights, exponentiate.
  if (log) weights = exp(weights-max(weights)) # max(weights) for numerical stability
  
  return(weights/sum(weights))
}

