# Unset java options
Sys.unsetenv("_JAVA_OPTIONS")

# Define Java options to suppress warnings
options(java.parameters = "-Djava.util.logging.ConsoleHandler.level=SEVERE")

# Load required libraries
library(nlrx)
library(future)  # Optimized for SLURM clusters
library(ggplot2)
library(dplyr)
library(docstring)

# Define paths
netlogopath <- "NetLogo-6.2.2"
modelpath <- file.path(
  netlogopath, "app", "models", "Code Examples",
  "Extensions Examples", "gis", "swine_model_asf_v8_fencing_10042025.nlogo"
)
modelname <- basename(modelpath)
outpath <- "out/"

# Check if paths exist
check_paths <- function(netlogopath, modelpath, outpath) {
  #' @title Check paths
  #' @description Check NetLogo installation, model file, and output folder
  if (!dir.exists(netlogopath)) stop(paste("? NetLogo path not found:", netlogopath))
  message(paste("? NetLogo path found:", netlogopath))

  if (!file.exists(modelpath)) stop(paste("? Model file not found:", modelpath))
  message(paste("? Model file found:", modelpath))

  if (!dir.exists(outpath)) {
    message(paste("?? Output directory does not exist. Creating:", outpath))
    dir.create(outpath, recursive = TRUE)
  } else {
    message(paste("? Output directory exists:", outpath))
  }
}
check_paths(netlogopath, modelpath, outpath)

# Detect SLURM environment variables
num_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(num_cpus) || num_cpus <= 0) num_cpus <- 1  # Default to 1 CPU if not set

slurm_mem <- as.numeric(Sys.getenv("SLURM_MEM_PER_CPU"))  
if (is.na(slurm_mem)) {
  slurm_mem <- as.numeric(Sys.getenv("SLURM_MEM_PER_NODE")) / num_cpus  
}
jvm_memory <- min(4096, slurm_mem)  # Cap at 4GB per CPU

# Initialize nlrx object
nl <- nl(
  nlversion = "6.2.2",
  nlpath = netlogopath,
  modelpath = modelpath,
  jvmmem = jvm_memory
)

# Print nlrx object to verify
print(nl)

# Attach an experiment
nl@experiment <- experiment(
  expname = "asf-model-nl",
  outpath = outpath,
  repetition = 1,
  tickmetrics = "true",
  idsetup = "setup-all-in",
  idgo = "go",
  runtime = 53,
  metrics = c(
    "avg-distance-traveled-by-infection", 
    "count patches with [pcolor = pink]", 
    "count farms with [infected-farm? = true]", 
    "r0-farm", "r0-patch"
  ),
  variables = list(),
  constants = list(
    "num-patches-to-infect" = 1,
    "termination-t" = 53,
    "cull?" = "false",
    "cull-infected-farm?" = "false",
    "stop-at-t?" = "true",
    "infect-random-wildboar-area-per-time-step?" = "false",
    "cull-n-random-wildboar-area?" = "false",
    "restrict-farm-connection?" = "false",
    "stop-when-I-plateaus?" = "false",
    "remove-wb-carcass?" = "false",
	"build-fence?" = "false",
    "alpha-patch-patch" = 1.611032,
    "beta-patch-patch" = 1.198028,
	"fencing-speed" = 0,
	"proba-of-fence-failing" = 0
  )
)

# Attach a simulation design
nseeds <- 100
nl@simdesign <- simdesign_simple(nl=nl, nseeds=nseeds)

# Validate nlrx object
eval_variables_constants(nl)
print(nl)

# Parallelization in local machine using future
library(future)

# Set up parallelization using the number of CPUs allocated by SLURM
plan(cluster, workers = num_cpus)

# Run simulations in parallel on SLURM
results <- run_nl_all(nl)

# Save results
setsim(nl, "simoutput") <- results
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
filename <- paste0("ASFModelSim_DeMeinweg_BASE", modelname, "_seeds", nseeds, "_", timestamp, ".rds")
saveRDS(nl, file.path(nl@experiment@outpath, filename))

message("? Simulation complete! Results saved to ", filename)