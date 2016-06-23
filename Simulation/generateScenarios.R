methods <- c("WOR", "CP")
scenariosExact12x12 <- expand.grid(method = methods, scenario = "exact12x12", sampleSize = 10000, replications = 100, stringsAsFactors = FALSE)
scenarios50x50_1 <- expand.grid(method = methods ,scenario = "50x50_1", sampleSize = 1000, replications = 100, stringsAsFactors=FALSE)

scenarios <- rbind(scenariosExact12x12, scenarios50x50_1)
scenarios$file <- paste0(scenarios$method, "-", scenarios$scenario, "-", scenarios$sampleSize, ".RData")
