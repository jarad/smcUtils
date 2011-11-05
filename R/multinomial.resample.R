multinomial.resample <-
function(weights, num.samples=length(weights)) {
  return(sample(length(weights), num.samples, replace=TRUE, prob=weights))
}

