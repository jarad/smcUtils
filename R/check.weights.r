
check.weights = function(weights,log=F,normalized=T) 
{
    stopifnot(length(weights)>0)

    if (all(log,normalized))
        warning("logged and normalized weights are unusual")

    if (!normalized) 
    {
        weights = renormalize.weights(weights, log, engine="R") # returns unlogged
        log=F
    }

    if (log) weights=exp(weights)

    stopifnot(all(weights>=0))
}


