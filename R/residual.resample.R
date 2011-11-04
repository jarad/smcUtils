residual.resample <-
function(weights,num.samples=length(weights),resample.function=multinomial.resample) {
  expected.num.samples = num.samples*weights
  deterministic.reps   = floor(expected.num.samples)
  sum.reps             = sum(deterministic.reps)
  residual             = num.samples-sum.reps

  deterministic.ids = rep2id(deterministic.reps)
  random.ids        = resample.function((expected.num.samples-deterministic.reps)/residual,residual)

  return(c(deterministic.ids,random.ids))
}

