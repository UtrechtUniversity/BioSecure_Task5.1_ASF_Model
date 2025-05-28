#---
#  title: "Parameter calibration of the NetLogo ASF model"
#  output: html_notebook
#---
  
#Using Latin Hypercube Sampling for Approximate bayesian computation
#Step 1: Create a nl object:
  
library(nlrx)
library(ggplot2)
library(dplyr)

# Define paths
netlogopath <- "NetLogo-6.2.2"
modelpath <- file.path(
  netlogopath, "app", "models", "Code Examples",
  "Extensions Examples", "gis", "swine_model_asf_v6_anon_WB_ONLY.nlogo"
)
modelname <- basename(modelpath)
outpath <- "out/"

# Check if paths exist
check_paths <- function(netlogopath, modelpath, outpath) {
  #' @title Check paths
  #' @description This function check the paths of the netlogo installation,
  #' the netlogo model file, and the output folder
  #' @param netlogopath Path of the NetLogo installation
  #' @param modelpath Path of the NetLogo model file
  #' @param outpath Path of the output directory for nlrx
  #' @return Prompts whether the paths are found or not.

  if (dir.exists(netlogopath)) {
    message(paste("NetLogo path found:", netlogopath))
  } else {
    stop(paste("NetLogo path not found:", netlogopath))
  }

  # Check if model file exists
  if (file.exists(modelpath)) {
    message(paste("Model file found:", modelpath))
  } else {
    stop(paste("Model file not found:", modelpath))
  }

  # Check if output directory exists, create if it doesn't
  if (dir.exists(outpath)) {
    message(paste("Output directory exists:", outpath))
  } else {
    warning(paste("Output directory does not exist. Creating:", outpath))
    dir.create(outpath, recursive = TRUE)
  }
}

check_paths(netlogopath, modelpath, outpath)

# Define Java options to suppress warnings
options(java.parameters = "-Djava.util.logging.ConsoleHandler.level=SEVERE")

# Initialize nlrx object
nl <- nl(nlversion = "6.2.2",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024
)

# Print nlrx object to verify
print(nl)

# Step 2: Attach an experiment
# We want to run a Latin Hypercube sampling, thus we need to define proper variable ranges. We also need to define our output metrics. These metrics are also used for the rejection sampling later on.

nl@experiment <- experiment(expname="asf-model-nl",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="false",
                            idsetup="setup-all-in",
                            idgo="go",
                            runtime=53,
                            metrics=c("avg-distance-traveled-by-infection"),
                            variables = list(#"num-patches-to-infect" = list(min=0, max=20, step=1, qfun="qunif"),
                                             #"alpha-farm-farm" = list(min=0, max=20, step=1, qfun="qunif"),
                                             #"alpha-farm-patch" = list(min=0, max=20, step=1, qfun="qunif"),
                                             #"alpha-patch-farm" = list(min=0, max=20, step=1, qfun="qunif"),
                                             "alpha-patch-patch" = list(min=0, max=3, step=0.1, qfun="qunif"),
                                             #"beta-farm-farm" = list(min=0, max=20, step=1, qfun="qunif"),
                                             #"beta-farm-patch" = list(min=0, max=20, step=1, qfun="qunif"),
                                             #"beta-patch-farm" = list(min=0, max=20, step=1, qfun="qunif"),
                                             "beta-patch-patch" = list(min=0, max=3, step=0.1, qfun="qunif")),
                            constants = list("num-patches-to-infect" = 1,
                                             "termination-t" = 53,
                                             "cull?" = "false",
                                             "cull-infected-farm?" = "false",
                                             "stop-at-t?" = "true",
                                             "infect-random-wildboar-area-per-time-step?" = "false",
                                             "cull-n-random-wildboar-area?" = "false",
                                             "restrict-farm-connection?" = "false",
                                             "stop-when-I-plateaus?" = "false"
                            ),
)

# Step 3: Attach a simulation design
# While the experiment defines the variables and specifications of the model, the simulation design creates a parameter input table based on these model specifications and the chosen simulation design method. nlrx provides a bunch of different simulation designs, such as full-factorial, latin-hypercube, sobol, morris and eFast (see “Simdesign Examples” vignette for more information on simdesigns). All simdesign helper functions need a properly defined nl object with a valid experiment design. Each simdesign helper also allows to define a number of random seeds that are randomly generated and can be used to execute repeated simulations of the same parameter matrix with different random-seeds (see “Advanced configuration” vignette for more information on random-seed and repetition management). A simulation design is attached to a nl object by using one of the simdesign helper functions:
  
samples <- 10000
nseeds <- 1
precision <- 3
nl@simdesign <- simdesign_lhs(nl=nl,
                              samples=samples,
                              nseeds=nseeds,
                              precision=precision)

# Step 4: Run simulations
# All information that is needed to run the simulations is now stored within the nl object. The run_nl_one() function allows to run one specific simulation from the siminput parameter table. The run_nl_all() function runs a loop over all simseeds and rows of the parameter input table siminput. The loops are constructed in a way that allows easy parallelisation, either locally or on remote HPC machines (see “Advanced configuration” vignette for more information on parallelisation). Before running your simulations you might want to check your current nl object setup. eval_variables_constants(nl) evaluates if the defined variables and constants are correctly defined and are consistent with the attached model. print(nl) prints a complete summary of the provided nl object including checkmarks that might help to indicate potential problems.

# Evaluate nl object:
eval_variables_constants(nl)
print(nl)

# Parallelization in local machine using future
library(future)

# Get the number of CPUs allocated by SLURM
num_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(num_cpus) || num_cpus <= 0) {
  num_cpus <- 1  # Default to 1 if SLURM_CPUS_PER_TASK is not valid
}

# Set up parallelization using the number of CPUs allocated by SLURM
plan(cluster, workers = num_cpus)

# Unset java options
Sys.unsetenv("_JAVA_OPTIONS")

# Execute simulations
results <- run_nl_all(nl = nl)

# Step 5: Investigate output
# We first attach the output results to our nl object and store a copy of the nl object on disk.
setsim(nl, "simoutput") <- results
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
filename <- paste0("ABClhs_", modelname,"_samples", samples,"_seeds", nseeds, "_", timestamp, ".rds") #init1 means 1 initial infection
saveRDS(nl, file.path(nl@experiment@outpath, filename))