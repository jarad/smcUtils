branching.resample <-
function(weights,num.samples=length(weights)) {
  expected.num.samples = weights*num.samples
  deterministic.reps   = floor(expected.num.samples)
  random.reps          = rbinom(num.samples,1,expected.num.samples-deterministic.reps)
  
  return(rep2id(deterministic.reps+random.reps))
}

