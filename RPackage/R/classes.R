setClass("withoutReplacementResult", slots = list(estimate = "mpfr", call = "call", start = "POSIXct", end = "POSIXct", n = "integer", seed = "integer", rowSums = "integer", columnSums = "integer", keepTables = "logical", tables = "list"))
setClass("withoutReplacementMergingResult", contains = "withoutReplacementResult", slots = list(mergeFrequency = "integer"))
setClass("conditionalPoissonResult", slots = list(varianceEstimate = "mpfr", estimate = "mpfr", call = "call", start = "POSIXct", end = "POSIXct", n = "integer", seed = "integer", rowSums = "integer", columnSums = "integer"))
