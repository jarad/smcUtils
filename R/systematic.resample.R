systematic.resample <-
function(weights,num.samples=length(weights)) {
  uniforms                = runif(1,0,1/num.samples)
  uniforms[2:num.samples] = uniforms[1]+(1:(num.samples-1))/num.samples

  return(inverse.cdf.weights(weights,uniforms))
}

