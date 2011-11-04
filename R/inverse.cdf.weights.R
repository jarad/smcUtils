inverse.cdf.weights <-
function(weights,uniforms=runif(length(weights))) {
  if (is.unsorted(uniforms)) uniforms = sort(uniforms)

  num.samples   = length(uniforms)
  ids           = integer(num.samples)
  weights.cusum = cumsum(weights)

  index = as.integer(1)
  one   = as.integer(1)
  for (i in 1:num.samples) {
    found = FALSE
    while (!found) {
      if (uniforms[i] > weights.cusum[index]) {
        index = index+one
      } else {
        found = TRUE
      }
    }
    ids[i] = index
  }
  return(ids)
}

