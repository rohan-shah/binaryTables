scenariosExact12x12 <- expand.grid(method = c("WOR", "WOR-Merge1", "WOR-Merge12"), scenario = "exact12x12", sampleSize = 10000, replications = 400, stringsAsFactors = FALSE)

methods <- c("WOR", "CP", "WOR-Merge1", "WOR-Merge50")
scenarios50x50_1 <- expand.grid(method = methods, scenario = "50x50_1", sampleSize = 1000, replications = 400, stringsAsFactors=FALSE)
scenarios50x50_2 <- expand.grid(method = methods, scenario = "50x50_2", sampleSize = 1000, replications = 400, stringsAsFactors=FALSE)

methods <- c("WOR", "CP", "WOR-Merge1", "WOR-Merge200")
scenariosPathological200 <- expand.grid(method = methods, scenario = "pathological_200", sampleSize = 100, replications = 400, stringsAsFactors=FALSE)

scenarios <- rbind(scenariosExact12x12, scenarios50x50_1, scenarios50x50_2, scenariosPathological200)
scenarios$file <- paste0(scenarios$method, "-", scenarios$scenario, "-", scenarios$sampleSize, ".RData")
