scenariosExact12x12 <- expand.grid(method = c("WOR", "WOR-Merge1", "WOR-Merge12"), scenario = "exact12x12", sampleSize = 10000, replications = 400, stringsAsFactors = FALSE)

methods <- c("WOR", "CP", "WOR-Merge1", "WOR-Merge50")
scenarios50x50_1 <- expand.grid(method = methods, scenario = "50x50_1", sampleSize = c(1000, 10000), replications = 400, stringsAsFactors=FALSE)

methods <- c("WOR", "CP", "WOR-Merge1", "WOR-Merge200")
scenariosPathological200 <- expand.grid(method = methods, scenario = "pathological_200", sampleSize = c(100, 200), replications = 400, stringsAsFactors=FALSE)

methods <- c("WOR", "CP", "WOR-Merge1", "WOR-Merge100")
scenariosPathological100 <- expand.grid(method = methods, scenario = "pathological_100", sampleSize = 100, replications = 400, stringsAsFactors=FALSE)

methods <- c("WOR", "CP", "WOR-Merge1", "WOR-Merge300")
scenariosPathological300 <- expand.grid(method = methods, scenario = "pathological_300", sampleSize = c(100, 300), replications = 400, stringsAsFactors=FALSE)

scenarios <- rbind(scenariosExact12x12, scenarios50x50_1, scenariosPathological300, scenariosPathological200, scenariosPathological100)
scenarios$file <- paste0(scenarios$method, "-", scenarios$scenario, "-", scenarios$sampleSize, ".RData")
