# Central path configuration. Every R script sources this file:
#
#   source(file.path(Sys.getenv("ADSEX_ROOT", "."), "config", "paths.R"))
#
# Run scripts from the repository root, or set ADSEX_ROOT to it.
# Data and results can live outside the repo (e.g. shared storage on the
# server) by setting ADSEX_DATA / ADSEX_RESULTS; see config/README.md.

PROJECT_ROOT <- normalizePath(Sys.getenv("ADSEX_ROOT", "."), mustWork = TRUE)
if (!file.exists(file.path(PROJECT_ROOT, "config", "paths.R"))) {
  stop("Project root not found. Run from the repository root or set ADSEX_ROOT.",
       call. = FALSE)
}

DATA_DIR    <- normalizePath(Sys.getenv("ADSEX_DATA",    file.path(PROJECT_ROOT, "data")),    mustWork = FALSE)
RESULTS_DIR <- normalizePath(Sys.getenv("ADSEX_RESULTS", file.path(PROJECT_ROOT, "results")), mustWork = FALSE)

# Results folder of an analysis module, e.g. module_results("07_master_regulators")
module_results <- function(module) file.path(RESULTS_DIR, module)
