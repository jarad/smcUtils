multinomial.resample <-
function(weights,num.samples=length(weights)) {
  return(sample(1:length(weights),num.samples,replace=TRUE,prob=weights))
}

