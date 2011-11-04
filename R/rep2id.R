rep2id <-
function(rep) {
  ids           = integer(sum(rep))
  current.index = integer(1)
  
  for (i in 1:length(rep)) {
    if (rep[i] != 0) {
      ids[current.index+1:rep[i]] = i
      current.index = current.index+rep[i]
    }
  }
  return(ids)
}

