ent.weights <-
function(weights) {
  return(-sum(weights*log2(weights+.Machine$double.eps)))

  # Could take the maximum of this number and 0 to avoid negative results.
  # Others define this with log (ln) rather than log2.
}

