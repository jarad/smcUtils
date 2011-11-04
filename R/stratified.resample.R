stratified.resample <-
function(weights,num.samples=length(weights)) {
  lower.bound = (1:num.samples-1)/num.samples
  uniforms    = runif(num.samples,lower.bound,lower.bound+1/num.samples)

  return(inverse.cdf.weights(weights,uniforms))
}

